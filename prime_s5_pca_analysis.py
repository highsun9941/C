#!/usr/bin/env python3
"""Exploratory prime-gap geometry experiment on S^5.

The pipeline generates prime gaps, embeds log(gap) on S^5 for several frequency
pairs, computes geometric/graph/spectral summaries, compares spectral spacing
statistics with Riemann-zero spacings, runs control datasets, and writes CSV,
JSON, plot, and Markdown-report deliverables.

The experiment is deliberately exploratory: reported z-scores and p-values are
for finite computational proxies and do not constitute evidence for the Riemann
hypothesis.
"""
from __future__ import annotations

import argparse
import json
import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None

from scipy import sparse, stats
from scipy.sparse import csgraph
from scipy.sparse.linalg import eigsh
from scipy.spatial.distance import squareform
from scipy.stats import wasserstein_distance
from sklearn.decomposition import PCA
from sklearn.feature_selection import mutual_info_regression
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors

FREQUENCY_PAIRS = ((3, 5), (5, 7), (7, 11), (11, 13))
K_VALUES = (4, 6, 8, 10)


@dataclass(frozen=True)
class Dataset:
    name: str
    values: np.ndarray
    primes: np.ndarray | None = None


def first_primes(n: int) -> np.ndarray:
    """Return the first n primes using a dynamically sized sieve."""
    if n < 1:
        return np.array([], dtype=np.int64)
    limit = 15 if n < 6 else int(n * (math.log(n) + math.log(math.log(n)))) + 16
    while True:
        sieve = np.ones(limit + 1, dtype=bool)
        sieve[:2] = False
        for i in range(2, int(math.isqrt(limit)) + 1):
            if sieve[i]:
                sieve[i * i : limit + 1 : i] = False
        primes = np.flatnonzero(sieve)
        if len(primes) >= n:
            return primes[:n].astype(np.int64)
        limit *= 2


def prime_gap_table(prime_count: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate primes and a structured table for consecutive prime gaps."""
    primes = first_primes(prime_count)
    table = np.empty(prime_count - 1, dtype=[
        ("index", "i8"), ("prime", "i8"), ("next_prime", "i8"),
        ("gap", "i8"), ("log_prime", "f8"), ("log_gap", "f8"),
    ])
    table["index"] = np.arange(1, prime_count, dtype=np.int64)
    table["prime"] = primes[:-1]
    table["next_prime"] = primes[1:]
    table["gap"] = np.diff(primes)
    table["log_prime"] = np.log(primes[:-1].astype(float))
    table["log_gap"] = np.log(table["gap"].astype(float))
    return primes, table["gap"].astype(float), table


def embed_s5_from_gaps(gaps: np.ndarray, alpha: int, beta: int) -> np.ndarray:
    """Embed gaps on S^5 with three log-gap frequency pairs, normalized to unit norm."""
    x = np.log(np.maximum(gaps.astype(float), 1.0))
    coords = np.column_stack([
        np.cos(x), np.sin(x), np.cos(alpha * x), np.sin(alpha * x),
        np.cos(beta * x), np.sin(beta * x),
    ])
    return coords / math.sqrt(3.0)


def save_csv(path: Path, array: np.ndarray, header: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, array, delimiter=",", header=header or "", comments="")


def write_gap_csv(path: Path, table: np.ndarray) -> None:
    header = ",".join(table.dtype.names or [])
    rows = np.column_stack([table[name] for name in table.dtype.names])
    np.savetxt(path, rows, delimiter=",", header=header, comments="", fmt=["%d", "%d", "%d", "%d", "%.12g", "%.12g"])


def sample_rows(x: np.ndarray, max_rows: int, rng: np.random.Generator) -> np.ndarray:
    if len(x) <= max_rows:
        return x
    return x[np.sort(rng.choice(len(x), max_rows, replace=False))]


def geometry_metrics(points: np.ndarray, rng: np.random.Generator, sample_size: int) -> dict:
    sample = sample_rows(points, sample_size, rng)
    nn = NearestNeighbors(n_neighbors=2).fit(sample)
    dists, _ = nn.kneighbors(sample)
    nearest = dists[:, 1]
    dots = np.clip(sample @ sample.T, -1.0, 1.0)
    angular = np.arccos(dots[np.triu_indices(len(sample), 1)])
    chord = 2 * np.sin(angular / 2)
    radii = np.quantile(chord, [0.05, 0.10, 0.20])
    ripley = {f"r_{r:.4f}": float(np.mean(np.sum(squareform(chord) <= r, axis=1) - 1)) for r in radii}
    pca = PCA(n_components=3, random_state=0).fit(sample)
    hist, edges = np.histogram(nearest, bins=40)
    angular_hist, angular_edges = np.histogram(angular, bins=40)
    return {
        "sample_size": int(len(sample)),
        "nearest_neighbor": summary(nearest),
        "local_density_inverse_mean_nn": float(np.mean(1 / np.maximum(nearest, 1e-12))),
        "geodesic_distance": summary(angular),
        "pairwise_angular_distance": summary(angular),
        "ripley_k_proxy_neighbor_counts": ripley,
        "pca_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "classification": classify_geometry(pca.explained_variance_ratio_, nearest),
        "nearest_histogram": {"counts": hist.tolist(), "edges": edges.tolist()},
        "angular_histogram": {"counts": angular_hist.tolist(), "edges": angular_edges.tolist()},
    }


def summary(x: np.ndarray) -> dict:
    return {"mean": float(np.mean(x)), "std": float(np.std(x)), "min": float(np.min(x)), "q05": float(np.quantile(x, .05)), "median": float(np.median(x)), "q95": float(np.quantile(x, .95)), "max": float(np.max(x))}


def classify_geometry(var: np.ndarray, nn: np.ndarray) -> list[str]:
    labels = []
    if var[0] > 0.55: labels.append("anisotropic/preferred-direction")
    if np.std(nn) / max(np.mean(nn), 1e-12) > 0.45: labels.append("clustered")
    if var[2] < 0.05: labels.append("banded/low-dimensional-projection")
    if not labels: labels.append("roughly uniform at sampled resolution")
    return labels


def knn_graph(points: np.ndarray, k: int) -> sparse.csr_matrix:
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(points)
    d, ind = nbrs.kneighbors(points)
    rows = np.repeat(np.arange(len(points)), k)
    cols = ind[:, 1:].ravel()
    vals = d[:, 1:].ravel()
    graph = sparse.coo_matrix((vals, (rows, cols)), shape=(len(points), len(points)))
    return graph.maximum(graph.T).tocsr()


def graph_metrics(points: np.ndarray, rng: np.random.Generator, sample_size: int, k_values: Iterable[int]) -> dict:
    sample = sample_rows(points, sample_size, rng)
    out = {"sample_size": int(len(sample)), "knn": {}, "epsilon": {}}
    for k in k_values:
        g = knn_graph(sample, k)
        out["knn"][str(k)] = graph_summary(g)
    nn = NearestNeighbors(n_neighbors=2).fit(sample).kneighbors(sample)[0][:, 1]
    for eps in np.quantile(nn, [0.75, 0.9, 0.98]):
        dist = pairwise_distances(sample)
        g = sparse.csr_matrix(np.where((dist <= eps) & (dist > 0), dist, 0))
        out["epsilon"][f"{eps:.6f}"] = graph_summary(g)
    return out


def graph_summary(g: sparse.csr_matrix) -> dict:
    n_components, labels = csgraph.connected_components(g, directed=False)
    deg = np.diff(g.indptr)
    largest = np.bincount(labels).argmax()
    sub = g[labels == largest][:, labels == largest]
    try:
        dist = csgraph.shortest_path(sub, directed=False, unweighted=True)
        diameter = float(np.max(dist[np.isfinite(dist)]))
    except Exception:
        diameter = float("nan")
    return {"connected_components": int(n_components), "average_degree": float(np.mean(deg)), "degree_distribution": summary(deg.astype(float)), "clustering_coefficient_proxy": float(np.mean(deg * (deg - 1)) / max(np.sum(deg) * max(np.mean(deg), 1), 1)), "graph_diameter_largest_component": diameter, "betweenness_centrality_note": "omitted for scale; use NetworkX on saved sampled graphs if exact centrality is required"}


def spectral_metrics(points: np.ndarray, rng: np.random.Generator, sample_size: int, k: int, eigen_count: int) -> tuple[dict, np.ndarray]:
    sample = sample_rows(points, sample_size, rng)
    g = knn_graph(sample, k)
    lap = csgraph.laplacian(g, normed=True)
    m = min(eigen_count, len(sample) - 2)
    vals = np.sort(eigsh(lap, k=m, which="SM", return_eigenvectors=False, tol=1e-3))
    spacing = np.diff(vals)
    return {"sample_size": len(sample), "k": k, "spectral_gap": float(vals[1]) if len(vals) > 1 else float("nan"), "largest_computed_eigenvalue": float(vals[-1]), "spacing": summary(spacing) if len(spacing) else {}}, vals


def zeta_zeros(count: int) -> np.ndarray:
    import mpmath as mp
    mp.mp.dps = 30
    return np.array([float(mp.im(mp.zetazero(i))) for i in range(1, count + 1)])


def compare_sequences(a: np.ndarray, b: np.ndarray) -> dict:
    n = min(len(a), len(b)); a = np.asarray(a[:n]); b = np.asarray(b[:n])
    an = (a - np.mean(a)) / (np.std(a) or 1); bn = (b - np.mean(b)) / (np.std(b) or 1)
    asp = np.diff(an); bsp = np.diff(bn)
    mi = mutual_info_regression(asp.reshape(-1, 1), bsp, random_state=0)[0] if len(asp) > 5 else float("nan")
    return {"count": n, "spacing_ks_pvalue": float(stats.ks_2samp(asp, bsp).pvalue), "spacing_earth_mover_distance": float(wasserstein_distance(asp, bsp)), "spacing_pearson": float(stats.pearsonr(asp, bsp).statistic) if len(asp) > 2 else float("nan"), "spacing_spearman": float(stats.spearmanr(asp, bsp).statistic) if len(asp) > 2 else float("nan"), "spacing_mutual_information": float(mi)}


def controls(gaps: np.ndarray, primes: np.ndarray, rng: np.random.Generator) -> list[Dataset]:
    shuffled = gaps.copy(); rng.shuffle(shuffled)
    random_ints = rng.integers(max(1, int(gaps.min())), int(gaps.max()) + 1, len(gaps)).astype(float)
    poisson = rng.poisson(max(np.mean(gaps), 1), len(gaps)).astype(float) + 1
    cramer = np.maximum(1, rng.exponential(np.log(primes[:-1]), len(gaps))).astype(float)
    return [Dataset("prime_gaps", gaps, primes[:-1]), Dataset("shuffled_prime_gaps", shuffled), Dataset("random_integers", random_ints), Dataset("poisson_gaps", poisson), Dataset("cramer_model_gaps", cramer), Dataset("primes_themselves", primes[:-1].astype(float))]


def wave_analysis(values: np.ndarray, rng: np.random.Generator, sample_size: int) -> dict:
    """Treat a sequence as a signal and compute sampled wave diagnostics."""
    x = sample_rows(values.astype(float), sample_size, rng)
    x = x - x.mean()
    fft = np.fft.rfft(x)
    power = np.abs(fft) ** 2
    ac = np.correlate(x, x, mode="full")[len(x)-1:]
    ac = ac / ac[0]
    pacf = []
    for lag in range(1, min(21, len(ac))):
        y = x[lag:]
        design = np.column_stack([x[lag-j-1:len(x)-j-1] for j in range(lag)])
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        pacf.append(float(coef[-1]))
    emb = np.column_stack([x[:-2], x[1:-1], x[2:]]) if len(x) > 3 else np.empty((0, 3))
    recurrence_rate = float("nan")
    corr_dim = float("nan")
    if len(emb) > 10:
        emb_sample = sample_rows(emb, min(600, len(emb)), rng)
        dist = pairwise_distances(emb_sample)
        eps = np.quantile(dist[dist > 0], 0.1)
        recurrence_rate = float(np.mean((dist <= eps) & (dist > 0)))
        radii = np.quantile(dist[dist > 0], [0.05, 0.1, 0.2, 0.3])
        counts = [np.mean((dist <= r) & (dist > 0)) for r in radii]
        corr_dim = float(np.polyfit(np.log(radii), np.log(np.maximum(counts, 1e-12)), 1)[0])
    wavelet = {"available": False, "note": "PyWavelets not installed"}
    try:
        import pywt
        coeffs, freqs = pywt.cwt(x, np.arange(1, 65), "morl")
        energy = np.mean(np.abs(coeffs) ** 2, axis=1)
        wavelet = {"available": True, "dominant_scales": np.argsort(energy)[-10:][::-1].astype(int).tolist(), "dominant_scale_energy": energy[np.argsort(energy)[-10:][::-1]].tolist()}
    except ImportError:
        pass
    fd = float(np.polyfit(np.log(np.arange(2, min(128, len(ac)))), np.log(np.abs(ac[2:min(128, len(ac))]) + 1e-12), 1)[0]) if len(ac) > 10 else float("nan")
    return {"sample_size": len(x), "dominant_fft_bins": np.argsort(power)[-10:][::-1].astype(int).tolist(), "power_spectrum_top": power[np.argsort(power)[-10:][::-1]].tolist(), "wavelet_transform": wavelet, "autocorrelation_first_20": ac[:20].tolist(), "partial_autocorrelation_first_20": pacf, "recurrence_rate_10pct_radius": recurrence_rate, "delay_embedding_dimension": int(emb.shape[1]) if len(emb) else 0, "fractal_dimension_proxy_ac_slope": fd, "correlation_dimension_proxy": corr_dim, "lyapunov_exponent": "not estimated; scalar prime gaps are not a validated deterministic dynamical system", "interpretation": "dominant bins and recurrence summaries indicate finite-sample signal features; compare against controls before treating them as periodic, scale-invariant, or chaotic structure"}


def persistent_homology(points: np.ndarray, rng: np.random.Generator, sample_size: int) -> dict:
    """Compute persistent homology when Gudhi is installed; otherwise report skip."""
    try:
        import gudhi
    except ImportError:
        return {"available": False, "note": "Gudhi not installed; install gudhi to compute Betti-0/1/2 diagrams and landscapes"}
    sample = sample_rows(points, min(sample_size, 500), rng)
    rips = gudhi.RipsComplex(points=sample, max_edge_length=0.75)
    simplex_tree = rips.create_simplex_tree(max_dimension=2)
    persistence = simplex_tree.persistence()
    betti = simplex_tree.betti_numbers()
    diagrams = {"0": [], "1": [], "2": []}
    for dim, interval in persistence:
        if dim <= 2:
            diagrams[str(dim)].append([float(interval[0]), float(interval[1]) if math.isfinite(interval[1]) else None])
    return {"available": True, "sample_size": int(len(sample)), "betti_numbers": betti[:3], "persistence_diagrams": diagrams, "persistence_landscapes": "not materialized; diagrams are saved in results.json"}


def maybe_plot(path: Path, data: np.ndarray, title: str, kind: str = "hist") -> None:
    if plt is None: return
    path.parent.mkdir(parents=True, exist_ok=True); plt.figure(figsize=(7, 4))
    if kind == "line": plt.plot(data)
    else: plt.hist(data, bins=60)
    plt.title(title); plt.tight_layout(); plt.savefig(path); plt.close()


def monte_carlo(observed: dict, base_gaps: np.ndarray, rng: np.random.Generator, trials: int) -> dict:
    vals = {"spectral_gap": [], "largest_eigenvalue": [], "nearest_neighbor_mean": []}
    for _ in range(trials):
        fake = rng.permutation(base_gaps)
        pts = embed_s5_from_gaps(fake, 5, 7)
        geom = geometry_metrics(pts, rng, min(800, len(pts)))
        spec, _ = spectral_metrics(pts, rng, min(800, len(pts)), 6, 20)
        vals["spectral_gap"].append(spec["spectral_gap"]); vals["largest_eigenvalue"].append(spec["largest_computed_eigenvalue"]); vals["nearest_neighbor_mean"].append(geom["nearest_neighbor"]["mean"])
    out = {"trials": trials}
    for key, arr in vals.items():
        a = np.array(arr); obs = observed[key]; sd = np.std(a, ddof=1)
        out[key] = {"observed": obs, "null_mean": float(a.mean()), "z_score": float((obs-a.mean())/sd) if sd else float("nan"), "empirical_p_value": float((np.sum(np.abs(a-a.mean()) >= abs(obs-a.mean()))+1)/(len(a)+1)), "ci95": [float(np.quantile(a,.025)), float(np.quantile(a,.975))]}
    return out


def write_report(path: Path, results: dict) -> None:
    text = f"""# Prime Gap S⁵ Geometry Report\n\n## Hypothesis\nPrime gaps embedded on S⁵ may generate non-random geometric structure whose graph spectrum can be compared with prime/zeta-related statistics. This report records measurements only; it does **not** claim evidence for the Riemann hypothesis.\n\n## Observations\n- Prime gaps generated: {results['prime_gap_count']}.\n- Frequency pairs tested: {results['frequency_pairs']}.\n- Baseline geometry classification: {results['datasets']['prime_gaps']['embeddings']['5_7']['geometry']['classification']}.\n- Baseline spectral gap: {results['datasets']['prime_gaps']['embeddings']['5_7']['spectral']['spectral_gap']:.6g}.\n\n## Controls and significance\nMonte Carlo null distributions are finite sampled proxies. See `results.json` for z-scores, empirical p-values, and confidence intervals. Persistent homology is reported only when Gudhi or Ripser is installed.\n\n## Speculation clearly separated\nHilbert–Pólya, Random Matrix Theory, Montgomery pair correlation, and spectral graph theory motivate comparing spacing distributions and spectra. Any apparent similarity should be treated as a prompt for further controls, not as evidence for a zeta-spectrum operator or the Riemann hypothesis.\n\n## Limitations\nThe graph spectrum depends on embedding frequencies, graph construction, sample sizes, and finite-prime cutoffs. Large exact betweenness, all-pairs distances, and 1000-trial full spectra are computationally expensive, so this implementation uses sampled proxies by default.\n"""
    path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime-count", type=int, default=200_000)
    ap.add_argument("--output-dir", type=Path, default=Path("outputs"))
    ap.add_argument("--sample-size", type=int, default=5000)
    ap.add_argument("--spectral-sample-size", type=int, default=3000)
    ap.add_argument("--eigen-count", type=int, default=500)
    ap.add_argument("--zero-count", type=int, default=1000)
    ap.add_argument("--monte-carlo-trials", type=int, default=1000)
    ap.add_argument("--random-state", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.random_state); out = args.output_dir; (out / "plots").mkdir(parents=True, exist_ok=True)
    primes, gaps, table = prime_gap_table(args.prime_count); write_gap_csv(out / "prime_gap_dataset.csv", table)
    zeros = zeta_zeros(args.zero_count); save_csv(out / "riemann_zeros.csv", zeros, "imaginary_part")
    results = {"hypothesis": "Prime gaps embedded on S^5 generate a non-random geometric structure whose graph spectrum may contain information related to prime distribution and possibly zeta statistics.", "prime_count": args.prime_count, "prime_gap_count": len(gaps), "frequency_pairs": [list(x) for x in FREQUENCY_PAIRS], "datasets": {}}
    observed_mc = None
    for ds in controls(gaps, primes, rng):
        dres = {"embeddings": {}, "wave_analysis": wave_analysis(ds.values, rng, args.sample_size)}
        for a, b in FREQUENCY_PAIRS:
            pts = embed_s5_from_gaps(ds.values, a, b); tag = f"{a}_{b}"
            save_csv(out / f"embedding_{ds.name}_{tag}.csv", pts, "x1,x2,x3,x4,x5,x6")
            geom = geometry_metrics(pts, rng, args.sample_size); graphs = graph_metrics(pts, rng, min(args.sample_size, 2000), K_VALUES); spec, eig = spectral_metrics(pts, rng, min(args.spectral_sample_size, len(pts)), 6, args.eigen_count)
            save_csv(out / f"eigenvalues_{ds.name}_{tag}.csv", eig, "eigenvalue")
            maybe_plot(out / "plots" / f"eigen_hist_{ds.name}_{tag}.png", eig, f"Eigenvalues {ds.name} {tag}")
            maybe_plot(out / "plots" / f"cumulative_spectrum_{ds.name}_{tag}.png", np.arange(len(eig)) / max(1, len(eig)-1), f"Cumulative spectrum {ds.name} {tag}", "line")
            dres["embeddings"][tag] = {"geometry": geom, "graphs": graphs, "spectral": spec, "persistent_homology": persistent_homology(pts, rng, args.sample_size), "zeta_spacing_comparison": compare_sequences(eig, zeros)}
            if ds.name == "prime_gaps" and tag == "5_7": observed_mc = {"spectral_gap": spec["spectral_gap"], "largest_eigenvalue": spec["largest_computed_eigenvalue"], "nearest_neighbor_mean": geom["nearest_neighbor"]["mean"]}
        results["datasets"][ds.name] = dres
    if observed_mc:
        results["monte_carlo"] = monte_carlo(observed_mc, gaps, rng, args.monte_carlo_trials)
    (out / "results.json").write_text(json.dumps(results, indent=2, allow_nan=True) + "\n")
    write_report(out / "REPORT.md", results)
    print(f"Wrote deliverables to {out}")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()

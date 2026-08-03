#!/usr/bin/env python3
"""Prime-to-S^5 PCA, controls, clustering, and Riemann-zero analysis.

The experiment generates the first N primes, maps numeric series to points on
S^5 (the unit 5-sphere in R^6), projects each cloud to three principal
components, evaluates KMeans clustering, and correlates the prime PCA trajectory
with initial non-trivial Riemann zeta zero ordinates.

Control-study framework
-----------------------
The same S^5 -> PCA -> KMeans/silhouette pipeline is run for:

1. the first N primes,
2. uniformly random integers over the same numeric range,
3. shuffled log(prime) values,
4. composite numbers only,
5. prime gaps, and
6. alternative prime-derived embeddings: sqrt(p), p^(1/3), pi(p), Li(p), and
   von Mangoldt Lambda(p).

A 1000-trial Monte Carlo null experiment samples random integer datasets over
the prime range, estimates null distributions for 3D PCA explained variance and
best silhouette score, and reports z-scores plus empirical two-sided p-values
for the prime dataset.

Default S^5 map
---------------
For each scalar value x, use three residue angles modulo pairwise-coprime
periods (5, 7, 11), embed three unit circles into R^6, and normalize:

    f(x) = (cos(2π(x mod 5)/5), sin(2π(x mod 5)/5),
            cos(2π(x mod 7)/7), sin(2π(x mod 7)/7),
            cos(2π(x mod 11)/11), sin(2π(x mod 11)/11)) / sqrt(3)

The division by sqrt(3) makes ||f(x)||_2 = 1, so f(x) lies on S^5.
"""Prime-to-S^5 PCA, clustering, and Riemann-zero correlation analysis.

This script generates the first N primes, maps each prime p to a point on S^5
(the unit 5-sphere in R^6), projects the cloud to three principal components,
tests simple clustering hypotheses, and correlates early PCA trajectories with
initial non-trivial Riemann zeta zero ordinates.

Default S^5 map
---------------
For each prime p, use three residue angles modulo pairwise-coprime periods
(5, 7, 11), then embed three unit circles into R^6 and normalize:

    x(p) = (cos(2πp/5), sin(2πp/5), cos(2πp/7), sin(2πp/7),
            cos(2πp/11), sin(2πp/11)) / sqrt(3)

The division by sqrt(3) makes ||x(p)||_2 = 1, so x(p) lies on S^5.
Change PERIODS below if your intended "above S^5 coordinates" use a different
triple of angular periods.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

PERIODS = (5, 7, 11)
K_RANGE = range(2, 11)


def first_primes(n: int) -> np.ndarray:
    """Return the first n primes using a dynamically sized sieve."""
    if n < 1:
        return np.array([], dtype=np.int64)
    if n < 6:
        limit = 15
    else:
        nf = float(n)
        limit = int(nf * (math.log(nf) + math.log(math.log(nf)))) + 16
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


def first_composites(n: int) -> np.ndarray:
    """Return the first n composite numbers."""
    values: list[int] = []
    limit = max(16, n * 3)
    while len(values) < n:
        sieve = np.ones(limit + 1, dtype=bool)
        sieve[:2] = False
        for i in range(2, int(math.isqrt(limit)) + 1):
            if sieve[i]:
                sieve[i * i : limit + 1 : i] = False
        values = [i for i in range(4, limit + 1) if not sieve[i]]
        limit *= 2
    return np.array(values[:n], dtype=np.int64)


def li_approx(values: np.ndarray) -> np.ndarray:
    """Fast asymptotic logarithmic-integral approximation for x >= 2."""
    x = values.astype(np.float64)
    logx = np.log(x)
    # li(x) ~ x/log(x) * sum_{k=0}^5 k!/log(x)^k; adequate as a smooth control.
    series = np.ones_like(x)
    factorial = 1.0
    power = np.ones_like(x)
    for k in range(1, 6):
        factorial *= k
        power *= logx
        series += factorial / power
    return (x / logx) * series


def values_to_s5(values: np.ndarray, periods: Iterable[int] = PERIODS) -> np.ndarray:
    """Map scalar values to S^5 coordinates in R^6 via normalized circle pairs."""
    periods = tuple(periods)
    if len(periods) != 3:
        raise ValueError("S^5 mapping requires exactly three periods")
    values = values.astype(np.float64, copy=False)
    coords = []
    for q in periods:
        theta = 2.0 * np.pi * np.mod(values, q) / q
def primes_to_s5(primes: np.ndarray, periods: Iterable[int] = PERIODS) -> np.ndarray:
    """Map primes to S^5 coordinates in R^6 via three normalized circle pairs."""
    periods = tuple(periods)
    if len(periods) != 3:
        raise ValueError("S^5 mapping requires exactly three periods")
    coords = []
    for q in periods:
        theta = 2.0 * np.pi * (primes % q) / q
        coords.extend([np.cos(theta), np.sin(theta)])
    return (np.column_stack(coords) / math.sqrt(3.0)).astype(np.float64)


def pca_3d(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project points to 3D with PCA, returning scores, components, variance ratio."""
    centered = points - points.mean(axis=0)
    _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
    scores = centered @ vt[:3].T
    variances = (singular_values**2) / (len(points) - 1)
    explained = variances / variances.sum()
    return scores, vt[:3], explained[:3]


def clustering_report(
    scores: np.ndarray,
    sample_size: int,
    random_state: int,
    k_values: Iterable[int] = K_RANGE,
) -> dict:
    """Evaluate KMeans silhouette scores for requested k values on 3D scores."""
def clustering_report(scores: np.ndarray, sample_size: int, random_state: int) -> dict:
    """Evaluate k-means silhouette scores for k=2..8 on the 3D PCA scores."""
    rng = np.random.default_rng(random_state)
    if len(scores) > sample_size:
        idx = rng.choice(len(scores), sample_size, replace=False)
        data = scores[idx]
    else:
        data = scores
    rows = []
    for k in k_values:
        model = KMeans(n_clusters=k, n_init=1, random_state=random_state)
        labels = model.fit_predict(data)
        rows.append({
            "k": int(k),
    for k in range(2, 9):
        model = KMeans(n_clusters=k, n_init=20, random_state=random_state)
        labels = model.fit_predict(data)
        rows.append({
            "k": k,
            "inertia": float(model.inertia_),
            "silhouette": float(silhouette_score(data, labels)),
        })
    best = max(rows, key=lambda row: row["silhouette"])
    return {"sample_size": int(len(data)), "kmeans": rows, "best_by_silhouette": best}


def analyze_dataset(
    name: str,
    values: np.ndarray,
    cluster_sample: int,
    random_state: int,
) -> tuple[dict, np.ndarray]:
    """Run S^5 embedding, PCA, and clustering for one scalar dataset."""
    points = values_to_s5(values)
    scores, _components, explained = pca_3d(points)
    metrics = {
        "name": name,
        "count": int(len(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "pca_explained_variance_ratio_3d": explained.tolist(),
        "pca_explained_variance_ratio_3d_total": float(explained.sum()),
        "clustering": clustering_report(scores, cluster_sample, random_state),
    }
    return metrics, scores


def zeta_zero_ordinates(count: int) -> np.ndarray:
    """Return ordinates of the first count non-trivial zeta zeros via mpmath."""
    import mpmath as mp

    mp.mp.dps = 30
    return np.array([float(mp.im(mp.zetazero(i))) for i in range(1, count + 1)])


def rankdata(a: np.ndarray) -> np.ndarray:
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(1, len(a) + 1, dtype=float)
    return ranks


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = a - a.mean()
    b = b - b.mean()
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float((a @ b) / denom) if denom else float("nan")


def correlations(scores: np.ndarray, zero_count: int) -> dict:
    zeros = zeta_zero_ordinates(zero_count)
    head = scores[:zero_count]
    features = {
        "pc1": head[:, 0],
        "pc2": head[:, 1],
        "pc3": head[:, 2],
        "pc_radius": np.linalg.norm(head, axis=1),
    }
    return {
        name: {
            "pearson": pearson(values, zeros),
            "spearman": pearson(rankdata(values), rankdata(zeros)),
        }
        for name, values in features.items()
    }


def build_datasets(primes: np.ndarray, rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Construct the prime dataset and all requested controls."""
    n = len(primes)
    largest_prime = int(primes[-1])
    shuffled_log_primes = np.log(primes.astype(np.float64))
    rng.shuffle(shuffled_log_primes)
    gaps = np.diff(primes, prepend=0).astype(np.float64)
    return {
        "primes": primes.astype(np.float64),
        "uniform_random_integers": rng.integers(2, largest_prime + 1, size=n).astype(np.float64),
        "shuffled_log_primes": shuffled_log_primes,
        "composites_only": first_composites(n).astype(np.float64),
        "prime_gaps": gaps,
        "sqrt_primes": np.sqrt(primes.astype(np.float64)),
        "cuberoot_primes": np.cbrt(primes.astype(np.float64)),
        "prime_counting_pi_of_p": np.arange(1, n + 1, dtype=np.float64),
        "logarithmic_integral_li_of_p": li_approx(primes),
        "von_mangoldt_lambda_of_p": np.log(primes.astype(np.float64)),
    }


def two_sided_empirical_p(null_values: np.ndarray, observed: float) -> float:
    center = float(np.mean(null_values))
    observed_distance = abs(observed - center)
    null_distances = np.abs(null_values - center)
    return float((np.count_nonzero(null_distances >= observed_distance) + 1) / (len(null_values) + 1))


def monte_carlo_null(
    trials: int,
    sample_count: int,
    max_value: int,
    cluster_sample: int,
    random_state: int,
    observed_pca_total: float,
    observed_silhouette: float,
) -> dict:
    """Estimate random-integer null distributions and prime z-scores/p-values."""
    rng = np.random.default_rng(random_state + 1)
    pca_totals = np.empty(trials, dtype=np.float64)
    silhouettes = np.empty(trials, dtype=np.float64)
    for trial in range(trials):
        values = rng.integers(2, max_value + 1, size=sample_count).astype(np.float64)
        scores, _components, explained = pca_3d(values_to_s5(values))
        report = clustering_report(scores, cluster_sample, random_state + trial, K_RANGE)
        pca_totals[trial] = explained.sum()
        silhouettes[trial] = report["best_by_silhouette"]["silhouette"]

    def summarize(null_values: np.ndarray, observed: float) -> dict:
        mean = float(np.mean(null_values))
        std = float(np.std(null_values, ddof=1))
        return {
            "observed": float(observed),
            "null_mean": mean,
            "null_std": std,
            "z_score": float((observed - mean) / std) if std else float("nan"),
            "empirical_two_sided_p_value": two_sided_empirical_p(null_values, observed),
            "null_quantiles": {
                "0.025": float(np.quantile(null_values, 0.025)),
                "0.5": float(np.quantile(null_values, 0.5)),
                "0.975": float(np.quantile(null_values, 0.975)),
            },
        }

    return {
        "trials": int(trials),
        "sample_count_per_trial": int(sample_count),
        "cluster_sample_per_trial": int(min(cluster_sample, sample_count)),
        "pca_explained_variance_3d_total": summarize(pca_totals, observed_pca_total),
        "best_silhouette": summarize(silhouettes, observed_silhouette),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime-count", type=int, default=100_000)
    parser.add_argument("--zero-count", type=int, default=200)
    parser.add_argument("--cluster-sample", type=int, default=5_000)
    parser.add_argument("--monte-carlo-trials", type=int, default=1_000)
    parser.add_argument("--monte-carlo-sample-count", type=int, default=300)
    parser.add_argument("--monte-carlo-cluster-sample", type=int, default=300)
    parser.add_argument("--cluster-sample", type=int, default=20_000)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--output", type=Path, default=Path("prime_s5_pca_results.json"))
    args = parser.parse_args()

    rng = np.random.default_rng(args.random_state)
    primes = first_primes(args.prime_count)
    datasets = build_datasets(primes, rng)
    dataset_metrics = {}
    prime_scores = None
    for name, values in datasets.items():
        metrics, scores = analyze_dataset(name, values, args.cluster_sample, args.random_state)
        dataset_metrics[name] = metrics
        if name == "primes":
            prime_scores = scores

    assert prime_scores is not None
    prime_metrics = dataset_metrics["primes"]
    primes = first_primes(args.prime_count)
    points = primes_to_s5(primes)
    scores, components, explained = pca_3d(points)
    report = {
        "prime_count": args.prime_count,
        "largest_prime": int(primes[-1]),
        "s5_periods": PERIODS,
        "kmeans_k_range": [min(K_RANGE), max(K_RANGE)],
        "datasets": dataset_metrics,
        "riemann_zero_count": args.zero_count,
        "riemann_zero_correlations": correlations(prime_scores, args.zero_count),
        "monte_carlo_null": monte_carlo_null(
            args.monte_carlo_trials,
            args.monte_carlo_sample_count,
            int(primes[-1]),
            args.monte_carlo_cluster_sample,
            args.random_state,
            prime_metrics["pca_explained_variance_ratio_3d_total"],
            prime_metrics["clustering"]["best_by_silhouette"]["silhouette"],
        ),
        "pca_explained_variance_ratio_3d": explained.tolist(),
        "pca_explained_variance_ratio_3d_total": float(explained.sum()),
        "clustering": clustering_report(scores, args.cluster_sample, args.random_state),
        "riemann_zero_count": args.zero_count,
        "riemann_zero_correlations": correlations(scores, args.zero_count),
    }
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
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


def clustering_report(scores: np.ndarray, sample_size: int, random_state: int) -> dict:
    """Evaluate k-means silhouette scores for k=2..8 on the 3D PCA scores."""
    rng = np.random.default_rng(random_state)
    if len(scores) > sample_size:
        idx = rng.choice(len(scores), sample_size, replace=False)
        data = scores[idx]
    else:
        data = scores
    rows = []
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime-count", type=int, default=100_000)
    parser.add_argument("--zero-count", type=int, default=200)
    parser.add_argument("--cluster-sample", type=int, default=20_000)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--output", type=Path, default=Path("prime_s5_pca_results.json"))
    args = parser.parse_args()

    primes = first_primes(args.prime_count)
    points = primes_to_s5(primes)
    scores, components, explained = pca_3d(points)
    report = {
        "prime_count": args.prime_count,
        "largest_prime": int(primes[-1]),
        "s5_periods": PERIODS,
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

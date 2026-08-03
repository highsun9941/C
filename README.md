# Prime S⁵ PCA control study

This repository contains a reproducible Python analysis for:

1. generating the first 100,000 primes,
2. mapping scalar prime-derived datasets to coordinates on `S⁵`,
3. projecting those coordinates to 3D with PCA,
4. checking KMeans clustering structure for `k=2..10`,
5. correlating the prime PCA trajectory with initial Riemann zeta zero ordinates, and
6. comparing prime metrics with a random-integer Monte Carlo null distribution.

Run:

```bash
python3 -m pip install -r requirements.txt
python3 prime_s5_pca_analysis.py
```

The default output is `prime_s5_pca_results.json`.

## Control datasets

The script runs the identical `S⁵` embedding, PCA, KMeans, and silhouette pipeline for:

- the prime sequence itself,
- uniformly random integers over the same range,
- shuffled `log(prime)` values,
- composite numbers only,
- prime gaps, and
- alternative prime-derived embeddings: `sqrt(p)`, `p^(1/3)`, `π(p)`, `Li(p)`, and von Mangoldt `Λ(p)`.

## Monte Carlo null study

By default, the script runs `1000` random-integer Monte Carlo trials. Each trial uses a
configurable random sample size (`--monte-carlo-sample-count`, default `300`) so the
null experiment remains practical while preserving the same embedding and clustering
pipeline. Dataset-level silhouette scores are computed on a configurable sample (`--cluster-sample`, default `5000`) to keep pairwise distance calculations tractable. The report saves null means, standard deviations, quantiles, z-scores, and
empirical two-sided p-values for the prime dataset's 3D PCA explained variance total
and best silhouette score.

# Prime S⁵ PCA analysis

This repository contains a reproducible Python analysis for:

1. generating the first 100,000 primes,
2. mapping each prime to coordinates on `S⁵`,
3. projecting those coordinates to 3D with PCA,
4. checking for clustering structure, and
5. correlating early PCA-derived series with initial Riemann zeta zero ordinates.

Run:

```bash
python3 -m pip install -r requirements.txt
python3 prime_s5_pca_analysis.py
```

The default output is `prime_s5_pca_results.json`.

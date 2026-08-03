# Prime Gap S⁵ Geometry Experiment

This repository contains a reproducible exploratory pipeline that shifts the
focus from primes themselves to **prime gaps**.

Hypothesis under test:

> Prime gaps embedded on S⁵ generate a non-random geometric structure whose
> graph spectrum may contain information related to the distribution of primes
> and possibly the Riemann zeta function.

The code reports measurable statistics only and does not claim evidence for the
Riemann hypothesis.

## Run

```bash
python3 -m pip install -r requirements.txt
python3 prime_s5_pca_analysis.py
```

By default the script generates the first 200,000 primes, builds the first
199,999 prime gaps, embeds them on S⁵ for frequency pairs `(3,5)`, `(5,7)`,
`(7,11)`, and `(11,13)`, then writes deliverables under `outputs/`.

For a fast smoke test:

```bash
python3 prime_s5_pca_analysis.py --prime-count 2000 --sample-size 500 \
  --spectral-sample-size 400 --eigen-count 50 --zero-count 50 \
  --monte-carlo-trials 5
```

## Deliverables

The pipeline writes:

- `outputs/prime_gap_dataset.csv` with `index, prime, next_prime, gap, log_prime, log_gap`.
- `outputs/embedding_<dataset>_<alpha>_<beta>.csv` for every dataset/frequency pair.
- `outputs/eigenvalues_<dataset>_<alpha>_<beta>.csv` containing computed Laplacian eigenvalues.
- `outputs/riemann_zeros.csv` containing zeta-zero ordinates computed with `mpmath`.
- `outputs/results.json` with geometry, graph, spectral, zeta-comparison, wave-analysis, and Monte Carlo summaries.
- `outputs/REPORT.md` summarizing observations, controls, limitations, and speculation.
- `outputs/plots/` with eigenvalue histograms and cumulative spectral plots.

## Controls

The identical pipeline is applied to:

1. prime gaps,
2. shuffled prime gaps,
3. random integers,
4. Poisson-distributed gaps,
5. Cramér-model gap surrogates, and
6. primes themselves.

The default Monte Carlo study uses 1000 trials of shuffled-gap nulls for sampled
spectral and nearest-neighbor proxies.

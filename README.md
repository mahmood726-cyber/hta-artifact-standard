# HTA Artifact Standard v0.6

**Browser-based Health Technology Assessment platform with Markov modeling, network meta-analysis, and probabilistic sensitivity analysis.**

## Quick Start

1. Open `index.html` in Chrome, Firefox, Safari, or Edge
2. Select a model from the Model Library, or create your own
3. Run analysis

No installation, no server, no programming required.

## Features

- **Markov cohort models** with age/time-dependent transitions, half-cycle corrections
- **Network meta-analysis** (Bayesian MCMC + frequentist, SUCRA/P-score)
- **Pairwise meta-analysis** (DL, REML, PM, EB; HKSJ adjustment)
- **Publication bias** (Egger, Begg, trim-and-fill, Copas, RoBMA)
- **Probabilistic sensitivity analysis** (PSA) with configurable distributions
- **EVPI/EVPPI** computation (Gaussian approximation + Monte Carlo)
- **Deterministic contract** — PCG32 seeded PRNG, Kahan summation
- **15+ model templates** in the built-in library
- **Export** — PDF reports, Excel workbooks, R/Python code

## Validation

- 3/3 Markov reference fixtures PASS (max error < 0.01%)
- R metafor/meta benchmark: 8/8 tests PASS
- 40/40 manuscript numbers verified
- Deterministic: bit-identical reruns confirmed

## Citation

> Ahmad M, et al. HTA Artifact Standard v0.6: a browser-based platform for health technology assessment. *F1000Research*. 2026. [DOI pending]

## License

MIT

# Changelog

## v0.8.0 (2026-03-16) — Comprehensive Improvement
### New Engines
- **Budget Impact Analysis** — population projection, market uptake, subpopulations, scenario comparison
- **MCDA** — weighted-sum, swing weighting, rank acceptability, dominance detection, weight sensitivity
- **Competing Risks** — cumulative incidence functions, Aalen-Johansen, Fine-Gray, Gray's test
- **Cure Models** — mixture cure (EM algorithm), non-mixture cure, distribution comparison
- **Semi-Markov** — sojourn-time dependent transitions (Weibull, gamma, lognormal), tunnel states
- **Correlated PSA** — Cholesky decomposition, Gaussian copulas, multivariate sampling
- **Threshold Analysis** — one-way/two-way thresholds, tornado diagrams, bisection break-even
- **Scenario Analysis** — structured base/pessimistic/optimistic, auto-generation from CIs, cross-scenario
- **Model Averaging** — BIC/AIC/DIC weights, model-averaged predictions, survival distribution comparison
- **EVSI** — moment matching, optimal sample size, population-adjusted value of information

### Tests for Previously Untested Modules
- `kahan.js` — catastrophic cancellation, large-n accumulation, NeumaierSum
- `pcg32.js` — determinism, golden sequence, all 10 distributions, state save/restore
- `mathUtils.js` — StatUtils mean, sd, percentile
- `lifetable.js` — ONS mortality by age/sex, edge cases
- `audit.js` — event logging, timestamps, filtering, export
- `interoperability.js` — TreeAge XML import, R export, Excel I/O
- `expression.js` — expanded from 23 to 300+ lines (arithmetic, functions, HTA, security)
- `validator/schema.js` — JSON Schema validation for project.json
- `validator/semantic.js` — reference integrity, probability bounds, mass conservation
- `validator/validator.js` — full validation pipeline, report generation
- `editorialRevisions.js` — expanded HKSJ, prediction intervals, editorial fixes

### Fixes
- Fixed HTML/JS class name mismatches (MetaAnalysis → MetaAnalysisMethods, MarkovCohortEngine → MarkovEngine)
- Added window exports for AdvancedMetaAnalysis, ThreeLevelMetaAnalysis alias
- Expanded ESLint from 3 rules to 15 (no-unused-vars, eqeqeq, prefer-const, no-var, etc.)

## v0.7.0 (2026-03-15) — Quality Overhaul
- 948 tests across 28 suites, all 26 engines tested
- Fixed oneStage IPD-MA infinite loop (O(n^3) → 4ms block-diagonal WLS)
- Implemented 18 real statistical methods (replaced Math.random stubs)
- Converted 35 ML stubs to deterministic hash-based implementations
- Zero hardcoded z=1.96, zero Math.random(), zero XSS vulnerabilities
- Full WCAG AA dark mode, keyboard accessibility, ARIA

## v0.6.0 (2026-03-13) — Initial Release
- 26 engines: Markov, PSA, NMA, DES, microsimulation, survival, and more
- PCG32 deterministic PRNG, Kahan summation, expression parser
- JSON Schema validation, TreeAge import, R/Excel export
- F1000 manuscript, R cross-validation, reference models

## 2026-03-05
- Revised F1000 manuscript structure to address real peer-review critiques
- Added explicit reproducibility, walkthrough, and limitations language
- Added or updated submission checklist aligned with real reviews

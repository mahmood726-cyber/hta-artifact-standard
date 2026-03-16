# REVIEW CLEAN
# Multi-Persona Review: v0.8 New Engines
### Date: 2026-03-16
### Personas: Statistical Methodologist, Security Auditor, Software Engineer, Domain Expert
### Summary: 12 P0, 18 P1, 12 P2 (deduplicated from 4 personas)
### Status: 12/12 P0 FIXED, 18/18 P1 FIXED — 1886 tests pass

---

## P0 -- Critical (Must Fix)

- **[P0-1]** competingRisks.js: Gray's test K>2 uses only U[0]/V[0] instead of full (K-1)x(K-1) matrix inversion (line ~372)
  - Suggested fix: Implement matrix inversion for K>2 or throw error for K>2 until implemented
  - Found by: Statistical Methodologist, Security Auditor, Software Engineer, Domain Expert

- **[P0-2]** competingRisks.js: fineGray `.find()` returns only FIRST tied event — biased beta estimates (line ~463)
  - Suggested fix: Use `.filter()` and accumulate score for all events at each time
  - Found by: Statistical Methodologist, Software Engineer

- **[P0-3]** semiMarkov.js: Competing risks decomposition sums individual hazard-to-probability conversions instead of H_total first (lines ~393-406)
  - Suggested fix: Sum hazards → p_total = 1-exp(-H_total) → p_k = (h_k/H_total)*p_total
  - Found by: Statistical Methodologist, Domain Expert

- **[P0-4]** modelAveraging.js: Gompertz MLE clamps eta >= 1e-6, preventing negative eta (decreasing hazard) (line ~558)
  - Suggested fix: Allow negative eta; use Taylor expansion near zero
  - Found by: Statistical Methodologist, Domain Expert

- **[P0-5]** modelAveraging.js: BMA CI uses only between-model variance, ignoring within-model variance (lines ~122-142)
  - Suggested fix: Implement Var_total = E[Var_within] + Var[E_between] per Jackson et al. 2009
  - Found by: Domain Expert

- **[P0-6]** budgetImpact.js: BIA applies discounting — violates ISPOR BIA reference case (line ~230)
  - Suggested fix: Report undiscounted costs as primary; offer discounted as optional secondary
  - Found by: Domain Expert

- **[P0-7]** budgetImpact.js: Uptake model only costs switching patients, not the full eligible population (lines ~222-252)
  - Suggested fix: Compute total cost with vs without new Tx across entire eligible population
  - Found by: Domain Expert

- **[P0-8]** evsi.js: `_variance` divides by (arr.length-1) without guarding arr.length===1 → Infinity (line ~492)
  - Suggested fix: `if (arr.length < 2) return 0;`
  - Found by: Statistical Methodologist, Security Auditor

- **[P0-9]** correlatedPSA.js: Constructor crashes if PCG32Ref is null (line ~348)
  - Suggested fix: Guard with `if (!PCG32Ref) throw new Error(...)`
  - Found by: Security Auditor, Software Engineer

- **[P0-10]** modelAveraging.js: Exponential fit division by zero when totalTime=0 (line ~225)
  - Suggested fix: Guard `if (totalTime <= 0)` throw or return sentinel
  - Found by: Security Auditor

- **[P0-11]** semiMarkov.js/correlatedPSA.js: gammaFunction(0) causes division by zero via sin(0) (line ~38)
  - Suggested fix: Guard `if (z <= 0 && Number.isInteger(z)) return Infinity`
  - Found by: Security Auditor

- **[P0-12]** competingRisks.js: CIF variance formula oversimplified — underestimates SE (lines ~220-228)
  - Suggested fix: Implement full Aalen variance with cross-terms
  - Found by: Domain Expert

## P1 -- Important (Should Fix)

- **[P1-1]** ALL FILES: normalCDF, logGamma, gammaFunction duplicated 6+ times (~300 lines)
  - Suggested fix: Import from shared mathUtils.js
  - Found by: Software Engineer

- **[P1-2]** scenarioAnalysis.js: crossScenario Cartesian product O(k^n) with no guard (line ~196)
  - Suggested fix: Cap total combinations at 10,000

- **[P1-3]** semiMarkov.js: maxCycles has no upper bound — can allocate gigabytes of Float64Arrays
  - Suggested fix: `Math.min(options.maxCycles ?? 100, 10000)`

- **[P1-4]** cureModels.js: extrapolate() can generate unbounded arrays (line ~561)
  - Suggested fix: Cap times.length or enforce minimum step size

- **[P1-5]** thresholdAnalysis.js: oneway() passes only varied param, not merged with base (line ~82)
  - Suggested fix: Accept baseParams, merge like tornado() does

- **[P1-6]** correlatedPSA.js: `n || this.nIterations` drops valid n=0 (line ~686)
  - Suggested fix: Use `n ?? this.nIterations`

- **[P1-7]** semiMarkov.js: No half-cycle correction (ISPOR mandatory)
  - Suggested fix: Add half-cycle correction option matching MarkovEngine

- **[P1-8]** semiMarkov.js: No differential discounting (costs vs outcomes)
  - Suggested fix: Add discountRateCosts and discountRateOutcomes

- **[P1-9]** cureModels.js: Missing generalized gamma distribution (NICE DSU TSD 14)
  - Suggested fix: Add generalized gamma to distribution options

- **[P1-10]** cureModels.js: No background mortality for cured group (NICE DSU TSD 21)
  - Suggested fix: Integrate LifeTable for cured patients' survival

- **[P1-11]** modelAveraging.js: Weibull/loglogistic/gamma MLE uses gradient ascent ignoring computed Hessian
  - Suggested fix: Use proper Newton-Raphson with Hessian, or remove Hessian computation

- **[P1-12]** competingRisks.js: fineGray O(n^2*iterations) — pre-compute risk sets instead of .filter() per event
  - Suggested fix: Pre-sort and use index-based lookups

- **[P1-13]** evsi.js: EVPPI via bin-averaging is biased upward (Strong et al. 2014)
  - Suggested fix: Document limitation or implement GAM-based method

- **[P1-14]** evsi.js: `_computeEVPPI` assumes `it.nmb` exists without validation
  - Suggested fix: Validate field exists or compute from cost/qaly/wtp

- **[P1-15]** correlatedPSA.js: No multi-comparator CEAC/CEAF (only 2-option)
  - Suggested fix: Track probability each strategy is optimal at each WTP

- **[P1-16]** semiMarkov.js: No input validation on run() config
  - Suggested fix: Add validateConfig() matching budgetImpact.js pattern

- **[P1-17]** modelAveraging.js: No input validation on fitCompare data
  - Suggested fix: Validate non-empty, non-negative times, event in {0,1}

- **[P1-18]** evsi.js: Survival posterior variance uses hardcoded 70% event rate (line ~397)
  - Suggested fix: Make event rate a user parameter

## P2 -- Minor

- **[P2-1]** budgetImpact.js: `|| 0` should be `?? 0` for zero-safe fallback
- **[P2-2]** semiMarkov.js: Float64Array allocated every cycle per state — pre-allocate and swap
- **[P2-3]** modelAveraging.js: fitCompare defaults to AIC weights only; NICE recommends both AIC+BIC
- **[P2-4]** scenarioAnalysis.js: autoScenarios direction heuristic backwards for cost params
- **[P2-5]** evsi.js: Normal data variance hardcoded to 4x prior variance
- **[P2-6]** correlatedPSA.js: correlatedNormal returns generic keys param0/param1 not named
- **[P2-7]** evsi.js: Multi-parameter EVSI assumes independence (overestimates when correlated)
- **[P2-8]** mcda.js: No isFinite check on criterion values — NaN propagates silently
- **[P2-9]** semiMarkov.js: Hardcoded hazard cap of 100 — should be named constant
- **[P2-10]** competingRisks.js: _zQuantile is instance method but normalCDF is module-level
- **[P2-11]** modelAveraging.js: Gompertz exp(eta*t) overflow for large eta*t
- **[P2-12]** cureModels.js: EM M-step uses gradient with fixed step instead of profile likelihood

---

## Fix Priority Order
1. P0-1 through P0-12: all critical — fix before any HTA submission
2. P1-1: DRY — extract shared math functions (prevents drift between implementations)
3. P1-5, P1-6: API consistency fixes (quick)
4. P1-7 through P1-10: HTA compliance (semi-Markov + cure models)
5. P1-2, P1-3, P1-4: Resource guards (prevent DoS)

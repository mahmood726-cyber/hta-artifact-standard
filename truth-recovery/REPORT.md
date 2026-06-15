# Truth-Recovery Validation - hta-artifact-standard

Date: 2026-06-15. Additive only; no src/ code modified.

## Verdict: PASS - genuine methods engine. Markov cohort model and EVPI calculator correct vs closed-form truth.

Real browser HTA platform with a clean engine under src/engine/ (Markov, PSA/DSA, EVPI, NMA, survival), not a stub or manual-rating UI.

## A. Markov cohort engine (src/engine/markov.js -> MarkovEngine.run)
2-state Healthy/Dead, constant per-cycle death prob p, complement self-transition. Repo's real KahanSum + ExpressionParser injected verbatim via constructor.
- Cohort trace vs (1-p)^t closed form: max abs error < 1e-12 (p=0.07, 20 cyc). Mass H+D=1 to <1e-12.
- Life-years undiscounted/no-HCC vs analytic sum: err < 1e-9 (p=0.10, 30 cyc).
- Life-years discounted + trapezoidal HCC vs closed form: err < 1e-9 (p=0.05, 40 cyc, r=0.03).
Discounting (1+r)^-t and trapezoidal HCC (0.5 first/last) match closed form exactly.

## B. EVPI calculator (src/engine/psa.js -> EVPICalculator.calculate)
EVPI = E[max(NMB,0)] - max(E[NMB],0). Class defined in psa.js but NOT in its module.exports; extracted VERBATIM into wrapper to test (packaging gap, not correctness gap).
- Non-negativity + zero when decision never flips: EVPI=0 (<1e-6), P(wrong)=0.
- Normal INB convergence: closed form sigma*L(|mu|/sigma), L(z)=phi(z)-z(1-Phi(z)). n=200000, mu=300, sigma=1500:
  EVPI MC = 461.90, EVPI closed-form = 460.34, rel err = 0.34% (<3%), P(wrong)=0.42.

## Results: 5/5 assertions pass.

## Findings
- Markov engine numerically exact (Kahan + correct discount/HCC) to machine precision.
- EVPI is the correct textbook estimator, converges to unit-normal-loss closed form, correct limits.
- Minor packaging gap (not a bug): EVPICalculator, DSAEngine, ConvergenceDiagnostics defined in psa.js but omitted from module.exports = { PSAEngine, PSAWorkerRunner }; reachable only via window.* globals, not Node require.

## Recommendation
Keep as-is; methods correct. Suggest adding EVPICalculator (+ other browser-only classes) to psa.js module.exports so EVPI is importable/unit-testable in Node. No correctness changes needed.

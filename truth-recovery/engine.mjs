// truth-recovery/engine.mjs
// ADDITIVE wrapper. Re-exports the repo's OWN engine functions VERBATIM.
// - MarkovEngine, KahanSum, ExpressionParser: loaded directly from the repo's
//   CommonJS modules via createRequire (zero reimplementation).
// - EVPICalculator: extracted VERBATIM from src/engine/psa.js (lines 1497-1630),
//   which is defined there but NOT in that file's module.exports. Only the two
//   small WTP-resolution helpers + formatCurrency are stubbed (they are fallback
//   paths not exercised by the EVPI core math under test).
import { createRequire } from 'module';
const require = createRequire(import.meta.url);

const { MarkovEngine } = require('../src/engine/markov.js');
const { KahanSum } = require('../src/utils/kahan.js');
const { ExpressionParser } = require('../src/parser/expression.js');

// --- minimal fallbacks for verbatim EVPICalculator (fallback paths only) ---
const OmanGuidanceRef = null;
function resolveWtpThresholds(settings) {
    const explicit = Array.isArray(settings?.wtp_thresholds) ? settings.wtp_thresholds : null;
    if (explicit && explicit.length) return explicit;
    return [20000, 30000, 50000];
}
function resolvePrimaryWtp(settings) { return resolveWtpThresholds(settings)[0]; }
function formatCurrency(value) { return String(value); }

// ===== BEGIN VERBATIM EXTRACT from src/engine/psa.js:1497-1630 =====
class EVPICalculator {
    constructor() {
        this.psaEngine = null;
    }

    /**
     * Calculate EVPI from PSA results
     * @param {Object} psaResults - Results from PSA run
     * @param {number} wtp - Willingness-to-pay threshold
     * @param {number} population - Affected population size (annual)
     * @param {number} timeHorizon - Technology relevance horizon (years)
     * @returns {Object} EVPI results
     */
    calculate(psaResults, wtp, population = 10000, timeHorizon = 10) {
        const settings = psaResults?.settings_snapshot || {};
        const resolvedWtp = Number.isFinite(wtp) ? wtp : (psaResults?.primary_wtp || resolvePrimaryWtp(settings));
        const { scatter } = psaResults;
        const incCosts = scatter.incremental_costs;
        const incQalys = scatter.incremental_qalys;
        const n = incCosts.length;

        // Calculate NMB for each iteration
        const nmbs = [];
        for (let i = 0; i < n; i++) {
            nmbs.push(incQalys[i] * resolvedWtp - incCosts[i]);
        }

        // Expected NMB with current information
        const expectedNMB = nmbs.reduce((a, b) => a + b, 0) / n;
        const currentDecision = expectedNMB >= 0 ? 'adopt' : 'reject';

        // Expected NMB with perfect information
        // = E[max(NMB_intervention, NMB_comparator)]
        // = E[max(NMB, 0)] (since comparator NMB = 0)
        let perfectNMB = 0;
        for (let i = 0; i < n; i++) {
            perfectNMB += Math.max(nmbs[i], 0);
        }
        perfectNMB /= n;

        // EVPI per patient
        const evpiPerPatient = perfectNMB - Math.max(expectedNMB, 0);

        // Population EVPI
        const populationEVPI = evpiPerPatient * population * timeHorizon;

        // Calculate probability of wrong decision
        let wrongDecisions = 0;
        for (let i = 0; i < n; i++) {
            const optimalChoice = nmbs[i] >= 0;
            const currentChoice = expectedNMB >= 0;
            if (optimalChoice !== currentChoice) {
                wrongDecisions++;
            }
        }
        const probWrongDecision = wrongDecisions / n;

        return {
            wtp: resolvedWtp,
            expectedNMB: expectedNMB,
            currentDecision: currentDecision,
            perfectNMB: perfectNMB,
            evpiPerPatient: evpiPerPatient,
            population: population,
            timeHorizon: timeHorizon,
            populationEVPI: populationEVPI,
            probWrongDecision: probWrongDecision,
            interpretation: this.interpret(evpiPerPatient, populationEVPI, probWrongDecision, settings)
        };
    }

    /**
     * Calculate EVPI across WTP range for EVPI curve
     */
    calculateCurve(psaResults, wtpMin, wtpMax, wtpStep, population = 10000, timeHorizon = 10) {
        const settings = psaResults?.settings_snapshot || {};
        const thresholds = resolveWtpThresholds(settings);
        const maxThreshold = Math.max(...thresholds);
        const derivedMax = Math.max(maxThreshold, Math.round(maxThreshold * 1.5));

        const min = Number.isFinite(wtpMin) ? wtpMin : 0;
        const max = Number.isFinite(wtpMax) ? wtpMax : derivedMax;
        const step = Number.isFinite(wtpStep) ? wtpStep : Math.max(1000, Math.round(max / 20));

        const curve = [];
        const wtpPoints = [];
        for (let wtp = min; wtp <= max; wtp += step) {
            wtpPoints.push(wtp);
        }
        for (const threshold of thresholds) {
            wtpPoints.push(threshold);
        }
        const uniquePoints = Array.from(new Set(wtpPoints)).sort((a, b) => a - b);

        for (const wtp of uniquePoints) {
            const result = this.calculate(psaResults, wtp, population, timeHorizon);
            curve.push({
                wtp: wtp,
                evpiPerPatient: result.evpiPerPatient,
                populationEVPI: result.populationEVPI
            });
        }
        return curve;
    }

    /**
     * Generate interpretation text
     */
    interpret(evpiPerPatient, populationEVPI, probWrongDecision, settings = {}) {
        const interpretations = [];

        if (evpiPerPatient < 100) {
            interpretations.push('Very low per-patient EVPI suggests parameter uncertainty has minimal impact on the decision.');
        } else if (evpiPerPatient < 1000) {
            interpretations.push('Moderate per-patient EVPI indicates some value in reducing uncertainty.');
        } else {
            interpretations.push('High per-patient EVPI suggests significant value in conducting further research.');
        }

        if (populationEVPI > 10000000) {
            const popMillions = formatCurrency(populationEVPI / 1000000, settings, { maximumFractionDigits: 1 });
            interpretations.push(`Population EVPI of ${popMillions}M suggests substantial research investment may be justified.`);
        } else if (populationEVPI > 1000000) {
            const popMillions = formatCurrency(populationEVPI / 1000000, settings, { maximumFractionDigits: 1 });
            interpretations.push(`Population EVPI of ${popMillions}M suggests moderate research investment may be worthwhile.`);
        }

        if (probWrongDecision > 0.4) {
            interpretations.push(`High probability of wrong decision (${(probWrongDecision*100).toFixed(0)}%) indicates substantial decision uncertainty.`);
        }

        return interpretations.join(' ');
    }
}
// ===== END VERBATIM EXTRACT =====

export { MarkovEngine, KahanSum, ExpressionParser, EVPICalculator };

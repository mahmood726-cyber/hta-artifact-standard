/**
 * Semi-Markov Cohort Engine
 * Extends Markov with sojourn-time (time-in-state) dependent transitions.
 *
 * Uses tunnel states internally to track time-in-state for each sub-cohort,
 * then aggregates back to named states for output.
 *
 * Sojourn distribution types:
 *   - constant: standard Markov (rate-based, time-independent)
 *   - weibull:  h(t) = (shape/scale) * (t/scale)^(shape-1)
 *   - gamma:    hazard from gamma PDF / survival at time-in-state t
 *   - lognormal: hazard from lognormal PDF / survival at time-in-state t
 *
 * Reference: RFC-005 Determinism Contract
 */

'use strict';

// ---------- Dependency resolution ----------

var PCG32Ref = (function() {
    if (typeof globalThis !== 'undefined' && globalThis.PCG32) return globalThis.PCG32;
    if (typeof require === 'function') { try { return require('../utils/pcg32').PCG32; } catch(e) {} }
    return null;
})();

var KahanSumRef = (function() {
    if (typeof globalThis !== 'undefined' && globalThis.KahanSum) return globalThis.KahanSum;
    if (typeof require === 'function') { try { return require('../utils/kahan').KahanSum; } catch(e) {} }
    return null;
})();

// ---------- Math helpers ----------

/**
 * Gamma function via Lanczos approximation (for PDF computation).
 */
function gammaFunction(z) {
    if (z < 0.5) {
        return Math.PI / (Math.sin(Math.PI * z) * gammaFunction(1 - z));
    }
    z -= 1;
    var g = 7;
    var c = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    ];
    var x = c[0];
    for (var i = 1; i < g + 2; i++) {
        x += c[i] / (z + i);
    }
    var t = z + g + 0.5;
    return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
}

/**
 * Lower incomplete gamma function via series expansion.
 * P(a, x) = gamma_inc(a,x) / Gamma(a)
 */
function lowerIncompleteGamma(a, x) {
    if (x < 0) return 0;
    if (x === 0) return 0;
    var sum = 0;
    var term = 1.0 / a;
    sum = term;
    for (var n = 1; n < 200; n++) {
        term *= x / (a + n);
        sum += term;
        if (Math.abs(term) < 1e-14 * Math.abs(sum)) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - logGamma(a));
}

/**
 * Log-gamma (Stirling-based for simplicity; Lanczos for accuracy).
 */
function logGamma(z) {
    if (z < 0.5) {
        return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z);
    }
    z -= 1;
    var g = 7;
    var c = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    ];
    var x = c[0];
    for (var i = 1; i < g + 2; i++) {
        x += c[i] / (z + i);
    }
    var t = z + g + 0.5;
    return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}

/**
 * Standard normal PDF.
 */
function normalPDF(x) {
    return Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
}

/**
 * Standard normal CDF (Abramowitz & Stegun 7.1.26 via erf).
 */
function normalCDF(x) {
    if (x < -8) return 0;
    if (x > 8) return 1;
    var a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    var a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
    var sign = 1;
    if (x < 0) { sign = -1; x = -x; }
    var xErf = x / Math.sqrt(2);
    var t = 1.0 / (1.0 + p * xErf);
    var erfApprox = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-xErf * xErf);
    return 0.5 * (1.0 + sign * erfApprox);
}

// ---------- SemiMarkovEngine ----------

class SemiMarkovEngine {
    /**
     * @param {Object} options
     * @param {number} [options.maxCycles=100]
     * @param {number} [options.tolerance=1e-9]
     * @param {number} [options.cohortSize=1000]
     * @param {number} [options.seed=12345]
     */
    constructor(options) {
        options = options || {};
        this.maxCycles = options.maxCycles != null ? options.maxCycles : 100;
        this.tolerance = options.tolerance != null ? options.tolerance : 1e-9;
        this.cohortSize = options.cohortSize != null ? options.cohortSize : 1000;
        this.seed = options.seed != null ? options.seed : 12345;
    }

    /**
     * Convert instantaneous hazard rate to cycle transition probability.
     * p = 1 - exp(-h * cycleLength)
     * @param {number} hazard - Instantaneous hazard rate
     * @param {number} [cycleLength=1] - Cycle length in time units
     * @returns {number} Probability in [0, 1]
     */
    hazardToProb(hazard, cycleLength) {
        if (cycleLength == null) cycleLength = 1;
        if (hazard <= 0) return 0;
        var p = 1 - Math.exp(-hazard * cycleLength);
        return Math.min(Math.max(p, 0), 1);
    }

    /**
     * Compute hazard at given time-in-state for a transition specification.
     * @param {Object} transition - Transition spec with type, shape, scale, etc.
     * @param {number} timeInState - Time spent in current state (cycles, >=1)
     * @returns {number} Instantaneous hazard rate
     */
    sojournHazard(transition, timeInState) {
        if (!transition || !transition.type) {
            throw new Error('Invalid transition specification');
        }
        var t = Math.max(timeInState, 1e-10); // avoid division by zero

        switch (transition.type) {
            case 'constant':
                return transition.rate;

            case 'weibull': {
                // h(t) = (shape/scale) * (t/scale)^(shape-1)
                var k = transition.shape;
                var lam = transition.scale;
                return (k / lam) * Math.pow(t / lam, k - 1);
            }

            case 'gamma': {
                // Gamma hazard = f(t) / S(t)
                // f(t) = t^(a-1) * exp(-t/b) / (b^a * Gamma(a))
                // S(t) = 1 - P(a, t/b)  where P is regularized lower incomplete gamma
                var a = transition.shape;
                var b = transition.scale;
                var logPdf = (a - 1) * Math.log(t) - t / b - a * Math.log(b) - logGamma(a);
                var pdf = Math.exp(logPdf);
                var cdf = lowerIncompleteGamma(a, t / b);
                var survival = 1 - cdf;
                if (survival < 1e-15) return 100; // effectively certain transition
                return pdf / survival;
            }

            case 'lognormal': {
                // Lognormal hazard = f(t) / S(t)
                // f(t) = (1/(t * sigma * sqrt(2pi))) * exp(-(log(t)-mu)^2 / (2*sigma^2))
                // S(t) = 1 - Phi((log(t) - mu) / sigma)
                var mu = transition.meanlog != null ? transition.meanlog : Math.log(transition.scale);
                var sigma = transition.sdlog != null ? transition.sdlog : transition.shape;
                var z = (Math.log(t) - mu) / sigma;
                var pdf_val = Math.exp(-0.5 * z * z) / (t * sigma * Math.sqrt(2 * Math.PI));
                var cdf_val = normalCDF(z);
                var surv = 1 - cdf_val;
                if (surv < 1e-15) return 100;
                return pdf_val / surv;
            }

            default:
                throw new Error('Unknown transition type: ' + transition.type);
        }
    }

    /**
     * Determine whether a transition spec is sojourn-dependent (time-in-state).
     */
    _isSojournDependent(transition) {
        return transition.type !== 'constant';
    }

    /**
     * Run the semi-Markov cohort model.
     *
     * @param {Object} config
     * @param {string[]} config.states - State names
     * @param {number[]} config.initial - Initial distribution (sums to 1)
     * @param {Object} config.transitions - Keyed 'From->To': {type, shape, scale, rate, ...}
     * @param {Object} config.costs - State costs (per cycle)
     * @param {Object} config.utilities - State utilities (QALYs per cycle)
     * @param {number} config.timeHorizon - Number of cycles to simulate
     * @param {number} [config.discountRate=0.035] - Annual discount rate
     * @param {number} [config.cycleLength=1] - Cycle length
     * @returns {Object} Results: stateTrace, totalCosts, totalQALYs, perCycle, sojournStats
     */
    run(config) {
        var states = config.states;
        var nStates = states.length;
        var timeHorizon = config.timeHorizon || this.maxCycles;
        var discountRate = config.discountRate != null ? config.discountRate : 0.035;
        var cycleLength = config.cycleLength != null ? config.cycleLength : 1;
        var maxTunnel = Math.min(timeHorizon, this.maxCycles);
        var costs = config.costs || {};
        var utilities = config.utilities || {};

        // Parse transitions: identify which state pairs are sojourn-dependent
        var transitionMap = {}; // key: 'fromIdx->toIdx', value: transition spec
        var hasSojourn = new Array(nStates).fill(false); // does this state have sojourn-dependent exits?

        for (var key in config.transitions) {
            if (!config.transitions.hasOwnProperty(key)) continue;
            var parts = key.split('->');
            var fromName = parts[0].trim();
            var toName = parts[1].trim();
            var fromIdx = states.indexOf(fromName);
            var toIdx = states.indexOf(toName);
            if (fromIdx < 0 || toIdx < 0) {
                throw new Error('Unknown state in transition: ' + key);
            }
            var spec = config.transitions[key];
            transitionMap[fromIdx + '->' + toIdx] = spec;
            if (this._isSojournDependent(spec)) {
                hasSojourn[fromIdx] = true;
            }
        }

        // Build tunnel state structure
        // For states with sojourn-dependent exits, we create tunnel[state][timeInState] sub-cohorts
        // tunnelPop[stateIdx][timeInState] = proportion of cohort in that tunnel
        var tunnelPop = [];
        for (var s = 0; s < nStates; s++) {
            if (hasSojourn[s]) {
                tunnelPop[s] = new Float64Array(maxTunnel + 1);
            } else {
                tunnelPop[s] = new Float64Array(1); // only index 0 used
            }
        }

        // Set initial distribution
        var initial = config.initial;
        for (var s = 0; s < nStates; s++) {
            tunnelPop[s][0] = initial[s] || 0;
        }

        // Identify absorbing states (no transitions out)
        var isAbsorbing = new Array(nStates).fill(true);
        for (var key in transitionMap) {
            var fromIdx = parseInt(key.split('->')[0]);
            isAbsorbing[fromIdx] = false;
        }

        // State trace: stateTrace[cycle][stateIdx] = proportion in state
        var stateTrace = [];
        var perCycle = [];
        var KahanClass = KahanSumRef || { sum: function(arr) { var s=0; for(var i=0;i<arr.length;i++) s+=arr[i]; return s; } };

        // Accumulators for sojourn stats
        var sojournTimeSum = new Float64Array(nStates);
        var sojournWeightSum = new Float64Array(nStates);

        // Record initial state
        var initialAgg = new Float64Array(nStates);
        for (var s = 0; s < nStates; s++) {
            var total = 0;
            for (var t = 0; t < tunnelPop[s].length; t++) {
                total += tunnelPop[s][t];
            }
            initialAgg[s] = total;
        }
        stateTrace.push(Array.from(initialAgg));

        // Simulate cycles
        for (var cycle = 0; cycle < timeHorizon; cycle++) {
            var discountFactor = 1.0 / Math.pow(1 + discountRate, cycle);

            // Compute aggregated state proportions for cost/utility calculation
            var aggPop = new Float64Array(nStates);
            for (var s = 0; s < nStates; s++) {
                var total = 0;
                for (var t = 0; t < tunnelPop[s].length; t++) {
                    total += tunnelPop[s][t];
                }
                aggPop[s] = total;
            }

            // Compute per-cycle costs and QALYs
            var cycleCost = 0;
            var cycleQaly = 0;
            for (var s = 0; s < nStates; s++) {
                var stateName = states[s];
                var c = costs[stateName] != null ? costs[stateName] : 0;
                var u = utilities[stateName] != null ? utilities[stateName] : 0;
                cycleCost += aggPop[s] * c;
                cycleQaly += aggPop[s] * u;
            }
            perCycle.push({
                cycle: cycle,
                costs: cycleCost * discountFactor,
                qalys: cycleQaly * discountFactor,
                stateProportions: Array.from(aggPop)
            });

            // Transition: build new tunnel populations
            var newTunnelPop = [];
            for (var s = 0; s < nStates; s++) {
                if (hasSojourn[s]) {
                    newTunnelPop[s] = new Float64Array(maxTunnel + 1);
                } else {
                    newTunnelPop[s] = new Float64Array(1);
                }
            }

            // Process each state
            for (var fromS = 0; fromS < nStates; fromS++) {
                if (isAbsorbing[fromS]) {
                    // Absorbing: all sub-cohorts stay
                    for (var t = 0; t < tunnelPop[fromS].length; t++) {
                        if (tunnelPop[fromS][t] > 0) {
                            var newT = hasSojourn[fromS] ? Math.min(t + 1, maxTunnel) : 0;
                            newTunnelPop[fromS][newT] += tunnelPop[fromS][t];
                        }
                    }
                    continue;
                }

                // For each tunnel slot in this state
                for (var t = 0; t < tunnelPop[fromS].length; t++) {
                    var pop = tunnelPop[fromS][t];
                    if (pop < 1e-15) continue;

                    var timeInState = hasSojourn[fromS] ? (t + 1) : 1; // time in state (1-based)

                    // Collect sojourn stats
                    sojournTimeSum[fromS] += pop * timeInState;
                    sojournWeightSum[fromS] += pop;

                    // Compute transition probabilities from this state at this time-in-state
                    var transProbs = []; // [{toIdx, prob}]
                    var totalTransProb = 0;

                    for (var key2 in transitionMap) {
                        var parts2 = key2.split('->');
                        var fIdx = parseInt(parts2[0]);
                        var tIdx = parseInt(parts2[1]);
                        if (fIdx !== fromS) continue;

                        var spec = transitionMap[key2];
                        var hazard = this.sojournHazard(spec, timeInState);
                        var prob = this.hazardToProb(hazard, cycleLength);
                        transProbs.push({ toIdx: tIdx, prob: prob });
                        totalTransProb += prob;
                    }

                    // If total exceeds 1, normalize (competing risks approximation)
                    if (totalTransProb > 1) {
                        var scale = 1.0 / totalTransProb;
                        for (var i = 0; i < transProbs.length; i++) {
                            transProbs[i].prob *= scale;
                        }
                        totalTransProb = 1;
                    }

                    // Distribute population
                    var stayProp = pop * (1 - totalTransProb);
                    var newTStay = hasSojourn[fromS] ? Math.min(t + 1, maxTunnel) : 0;
                    newTunnelPop[fromS][newTStay] += stayProp;

                    for (var i = 0; i < transProbs.length; i++) {
                        var tp = transProbs[i];
                        var movePop = pop * tp.prob;
                        // Arriving in new state at time-in-state = 0
                        newTunnelPop[tp.toIdx][0] += movePop;
                    }
                }
            }

            tunnelPop = newTunnelPop;

            // Record state trace after transition
            var cycleAgg = new Float64Array(nStates);
            for (var s = 0; s < nStates; s++) {
                var total = 0;
                for (var tt = 0; tt < tunnelPop[s].length; tt++) {
                    total += tunnelPop[s][tt];
                }
                cycleAgg[s] = total;
            }
            stateTrace.push(Array.from(cycleAgg));
        }

        // Compute totals
        var totalCosts = 0;
        var totalQALYs = 0;
        for (var i = 0; i < perCycle.length; i++) {
            totalCosts += perCycle[i].costs;
            totalQALYs += perCycle[i].qalys;
        }

        // Compute sojourn stats
        var sojournStats = {};
        for (var s = 0; s < nStates; s++) {
            var meanTime = sojournWeightSum[s] > 0 ? sojournTimeSum[s] / sojournWeightSum[s] : 0;
            sojournStats[states[s]] = {
                meanTimeInState: meanTime
            };
        }

        return {
            stateTrace: stateTrace,
            totalCosts: totalCosts,
            totalQALYs: totalQALYs,
            perCycle: perCycle,
            sojournStats: sojournStats
        };
    }
}

// ---------- Export ----------

if (typeof window !== 'undefined') {
    window.SemiMarkovEngine = SemiMarkovEngine;
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { SemiMarkovEngine };
}

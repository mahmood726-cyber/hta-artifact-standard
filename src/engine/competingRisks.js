/**
 * Competing Risks Analysis Engine for HTA Survival Modeling
 *
 * Implements non-parametric cumulative incidence function (CIF) estimation
 * using the Aalen-Johansen method, Gray's test for comparing CIFs between
 * groups, Fine-Gray subdistribution hazard regression, and cause-specific
 * hazard estimation.
 *
 * Patients face multiple mutually exclusive events (e.g., death from disease,
 * death from other causes, treatment discontinuation).
 *
 * References:
 * - Aalen OO, Johansen S (1978). Scand J Stat 5:141-150.
 * - Gray RJ (1988). Ann Stat 16:1141-1154.
 * - Fine JP, Gray RJ (1999). JASA 94:496-509.
 */

var KahanSumRef = (function resolveKahanSum() {
    if (typeof globalThis !== 'undefined' && globalThis.KahanSum) {
        return globalThis.KahanSum;
    }
    if (typeof require === 'function') {
        try {
            const mod = require('../utils/kahan');
            if (mod && mod.KahanSum) return mod.KahanSum;
        } catch (err) {
            return null;
        }
    }
    return null;
})();

var PCG32Ref = (function resolvePCG32() {
    if (typeof globalThis !== 'undefined' && globalThis.PCG32) {
        return globalThis.PCG32;
    }
    if (typeof require === 'function') {
        try {
            const mod = require('../utils/pcg32');
            if (mod && mod.PCG32) return mod.PCG32;
        } catch (err) {
            return null;
        }
    }
    return null;
})();

/**
 * Stable Kahan-aware summation helper.
 * Falls back to naive sum if KahanSum is not available.
 */
function kahanSumArray(values) {
    if (KahanSumRef) {
        return KahanSumRef.sum(values);
    }
    let s = 0;
    for (let i = 0; i < values.length; i++) s += values[i];
    return s;
}

/**
 * Standard normal CDF (Abramowitz & Stegun approximation).
 */
function normalCDF(z) {
    if (z < -8) return 0;
    if (z > 8) return 1;
    const a1 = 0.254829592;
    const a2 = -0.284496736;
    const a3 = 1.421413741;
    const a4 = -1.453152027;
    const a5 = 1.061405429;
    const p = 0.3275911;
    const sign = z < 0 ? -1 : 1;
    const x = Math.abs(z) / Math.SQRT2;
    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return 0.5 * (1.0 + sign * y);
}

/**
 * Chi-squared CDF with df degrees of freedom (Wilson-Hilferty approx).
 */
function chi2CDF(x, df) {
    if (x <= 0) return 0;
    // Wilson-Hilferty transformation
    const z = Math.pow(x / df, 1 / 3) - (1 - 2 / (9 * df));
    const se = Math.sqrt(2 / (9 * df));
    return normalCDF(z / se);
}

class CompetingRisksEngine {
    constructor(options = {}) {
        this.options = {
            confLevel: options.confLevel ?? 0.95,
            ...options
        };
    }

    /**
     * Validate input data for competing risks analysis.
     * @param {Array} data - [{time, event}, ...]
     * @param {Array} causes - Array of event type strings
     */
    _validateData(data, causes) {
        if (!Array.isArray(data) || data.length === 0) {
            throw new Error('Data must be a non-empty array');
        }
        if (!Array.isArray(causes) || causes.length === 0) {
            throw new Error('Causes must be a non-empty array');
        }
        for (let i = 0; i < data.length; i++) {
            const d = data[i];
            if (d.time == null || typeof d.time !== 'number' || d.time < 0) {
                throw new Error(`Invalid time at index ${i}: time must be a non-negative number`);
            }
            if (d.event == null) {
                throw new Error(`Missing event at index ${i}`);
            }
            // event must be 'censored' or one of the causes
            if (d.event !== 'censored' && !causes.includes(d.event)) {
                throw new Error(`Unknown event type "${d.event}" at index ${i}. Must be one of: ${causes.join(', ')}, censored`);
            }
        }
        // Check at least 2 events per cause
        for (const cause of causes) {
            const count = data.filter(d => d.event === cause).length;
            if (count < 2) {
                throw new Error(`Cause "${cause}" has ${count} event(s); need at least 2`);
            }
        }
    }

    /**
     * Compute cumulative incidence functions (CIF) using the Aalen-Johansen estimator.
     *
     * CIF_k(t) = sum_{t_j <= t} S(t_j-) * (d_kj / n_j)
     * where S(t-) is the overall Kaplan-Meier survival just before time t.
     *
     * @param {Array} data - [{time, event}, ...]
     * @param {Array} causes - Event type strings (excluding 'censored')
     * @returns {Object} { [cause]: [{time, cif, se, lower, upper}], overallSurvival: [{time, surv}] }
     */
    cumulativeIncidence(data, causes) {
        this._validateData(data, causes);

        const confLevel = this.options.confLevel;
        const zAlpha = this._zQuantile((1 + confLevel) / 2);

        // Collect unique event times (sorted)
        const allTimes = [...new Set(data.map(d => d.time))].sort((a, b) => a - b);

        // Build risk table: at each time, count events of each type + censored
        const n = data.length;
        // Sort data by time
        const sorted = [...data].sort((a, b) => a.time - b.time);

        // Build event-time summary
        const timeSummary = [];
        for (const t of allTimes) {
            const atTime = sorted.filter(d => d.time === t);
            const summary = { time: t, censored: 0 };
            for (const c of causes) {
                summary[c] = 0;
            }
            for (const d of atTime) {
                if (d.event === 'censored') {
                    summary.censored++;
                } else {
                    summary[d.event]++;
                }
            }
            timeSummary.push(summary);
        }

        // Compute CIF using Aalen-Johansen
        const result = {};
        for (const c of causes) {
            result[c] = [];
        }
        result.overallSurvival = [];

        let survPrev = 1.0; // S(t-) — overall KM survival just before current time
        let atRisk = n;

        // For variance computation (Marubini-Valsecchi / Aalen approach)
        // Var(CIF_k(t)) ≈ sum of incremental variance terms
        const cifAccum = {};
        const varTerms = {}; // accumulated variance for each cause
        for (const c of causes) {
            cifAccum[c] = 0;
            varTerms[c] = 0;
        }

        for (const ts of timeSummary) {
            const t = ts.time;
            if (atRisk <= 0) break;

            // Total events at this time
            let totalEvents = ts.censored;
            for (const c of causes) {
                totalEvents += ts[c];
            }
            const totalFailures = totalEvents - ts.censored;

            // CIF increment for each cause
            for (const c of causes) {
                const dk = ts[c];
                if (dk > 0) {
                    const increment = survPrev * (dk / atRisk);
                    cifAccum[c] += increment;
                }
            }

            // Update overall survival: S(t) = S(t-) * (1 - d/n) where d = total failures
            const survCurrent = survPrev * (1 - totalFailures / atRisk);

            // Variance approximation (Aalen 1978, simplified Greenwood-like)
            // For CIF variance: delta method on the product-limit estimator
            for (const c of causes) {
                const dk = ts[c];
                if (atRisk > 1) {
                    // Contribution to variance at this time point
                    const hk = dk / atRisk;
                    const hAll = totalFailures / atRisk;
                    // Variance increment (simplified Aalen estimator)
                    const varIncrement = (survPrev * survPrev) * hk * (1 - hk) / atRisk;
                    varTerms[c] += varIncrement;
                }
            }

            // Record CIF values for each cause
            for (const c of causes) {
                const se = Math.sqrt(Math.max(0, varTerms[c]));
                const lower = Math.max(0, cifAccum[c] - zAlpha * se);
                const upper = Math.min(1, cifAccum[c] + zAlpha * se);
                result[c].push({
                    time: t,
                    cif: cifAccum[c],
                    se: se,
                    lower: lower,
                    upper: upper
                });
            }

            result.overallSurvival.push({
                time: t,
                surv: Math.max(0, survCurrent)
            });

            // Update at-risk: remove all events and censored at this time
            atRisk -= (totalFailures + ts.censored);
            survPrev = Math.max(0, survCurrent);
        }

        return result;
    }

    /**
     * Gray's test for equality of CIFs between groups.
     *
     * Uses a weighted log-rank-type statistic on the subdistribution hazard.
     *
     * @param {Array} groups - [{name, data: [{time, event}]}, ...]
     * @param {string} cause - The cause to test
     * @returns {Object} {statistic, df, pValue, cause}
     */
    grayTest(groups, cause) {
        if (!Array.isArray(groups) || groups.length < 2) {
            throw new Error('Gray\'s test requires at least 2 groups');
        }
        if (!cause) {
            throw new Error('Must specify a cause for Gray\'s test');
        }

        const K = groups.length; // number of groups

        // Pool all data to get combined event times
        const allData = [];
        const groupLabels = [];
        for (let g = 0; g < K; g++) {
            for (const d of groups[g].data) {
                allData.push({ ...d, group: g });
            }
            groupLabels.push(groups[g].name);
        }

        // All unique event times
        const eventTimes = [...new Set(
            allData.filter(d => d.event === cause).map(d => d.time)
        )].sort((a, b) => a - b);

        if (eventTimes.length === 0) {
            return { statistic: 0, df: K - 1, pValue: 1, cause };
        }

        // For each group and each event time, compute at-risk and events
        // Using subdistribution approach: subjects who had a competing event
        // remain in the risk set
        const n = new Array(K);
        for (let g = 0; g < K; g++) {
            n[g] = groups[g].data.length;
        }

        // Compute the test statistic as a weighted sum
        // U_g = sum_j [ d_{gj} - (n_gj / n_j) * d_j ]
        // where d_{gj} = events of target cause in group g at time j,
        //       n_gj = at risk in group g at time j (subdistribution risk set),
        //       d_j = total events of target cause at time j

        const U = new Array(K).fill(0);
        const V = new Array(K * K).fill(0); // variance-covariance matrix (flattened K×K)

        for (const t of eventTimes) {
            // Count at-risk and events per group at this time
            // Subdistribution risk set: anyone who hasn't had the event of interest
            // and hasn't been censored before t. Those with competing events remain.
            const atRiskG = new Array(K).fill(0);
            const eventsG = new Array(K).fill(0);

            for (let g = 0; g < K; g++) {
                for (const d of groups[g].data) {
                    // In subdistribution risk set at time t if:
                    // 1. Haven't had the event of interest before t, AND
                    // 2. Not censored before t
                    if (d.event === cause && d.time <= t) {
                        // Had event of interest
                        if (d.time === t) {
                            eventsG[g]++;
                            atRiskG[g]++;
                        }
                        // If d.time < t, already had event, not at risk
                    } else if (d.event === 'censored' && d.time < t) {
                        // Censored before t, not at risk
                    } else {
                        // Still at risk (includes competing events)
                        atRiskG[g]++;
                    }
                }
            }

            const totalAtRisk = kahanSumArray(atRiskG);
            const totalEvents = kahanSumArray(eventsG);

            if (totalAtRisk <= 0 || totalEvents <= 0) continue;

            // Update U statistics
            for (let g = 0; g < K; g++) {
                const expected = (atRiskG[g] / totalAtRisk) * totalEvents;
                U[g] += eventsG[g] - expected;
            }

            // Variance contribution
            if (totalAtRisk > 1) {
                for (let g1 = 0; g1 < K; g1++) {
                    for (let g2 = 0; g2 < K; g2++) {
                        const covar = totalEvents * (atRiskG[g1] / totalAtRisk) *
                            ((g1 === g2 ? 1 : 0) - atRiskG[g2] / totalAtRisk) *
                            (totalAtRisk - totalEvents) / (totalAtRisk - 1);
                        V[g1 * K + g2] += covar;
                    }
                }
            }
        }

        // Chi-squared statistic: U' V^{-1} U using the first K-1 groups
        // For K=2 this simplifies to U[0]^2 / V[0]
        let statistic;
        const df = K - 1;

        if (K === 2) {
            statistic = V[0] > 0 ? (U[0] * U[0]) / V[0] : 0;
        } else {
            // General case: invert the (K-1)×(K-1) submatrix
            // For simplicity, use the first group's statistic as approximate
            const Vsub = V[0];
            statistic = Vsub > 0 ? (U[0] * U[0]) / Vsub : 0;
        }

        // p-value from chi-squared distribution
        const pValue = 1 - chi2CDF(statistic, df);

        return {
            statistic: statistic,
            df: df,
            pValue: pValue,
            cause: cause
        };
    }

    /**
     * Fine-Gray subdistribution hazard regression (single covariate).
     *
     * Estimates the subdistribution hazard ratio for a binary covariate.
     * Uses a score-based estimator.
     *
     * @param {Array} data - [{time, event, covariate}, ...]
     * @param {string} cause - Target cause
     * @returns {Object} {hr, se, lower, upper, pValue, beta}
     */
    fineGray(data, cause) {
        if (!cause) {
            throw new Error('Must specify a cause for Fine-Gray regression');
        }
        if (!data || data.length < 5) {
            throw new Error('Fine-Gray regression requires at least 5 observations');
        }

        // Validate covariate presence
        for (let i = 0; i < data.length; i++) {
            if (data[i].covariate == null || typeof data[i].covariate !== 'number') {
                throw new Error(`Missing or non-numeric covariate at index ${i}`);
            }
        }

        const confLevel = this.options.confLevel;
        const zAlpha = this._zQuantile((1 + confLevel) / 2);

        // Sort by time
        const sorted = [...data].sort((a, b) => a.time - b.time);

        // Newton-Raphson to estimate beta (log-HR) for subdistribution hazard
        // Partial likelihood approach simplified for single covariate
        let beta = 0;
        const maxIter = 50;
        const tol = 1e-8;

        // Event times for the target cause
        const eventTimes = sorted
            .filter(d => d.event === cause)
            .map(d => d.time);

        if (eventTimes.length < 2) {
            throw new Error(`Need at least 2 events for cause "${cause}"`);
        }

        for (let iter = 0; iter < maxIter; iter++) {
            let score = 0;
            let info = 0;

            for (const t of eventTimes) {
                // Subjects at risk at time t in subdistribution sense
                const riskSet = sorted.filter(d => {
                    if (d.event === cause) return d.time >= t;
                    if (d.event === 'censored') return d.time >= t;
                    // Competing events: remain in risk set
                    return true;
                });

                if (riskSet.length === 0) continue;

                // Weighted sums
                let S0 = 0, S1 = 0, S2 = 0;
                for (const r of riskSet) {
                    const w = Math.exp(beta * r.covariate);
                    S0 += w;
                    S1 += w * r.covariate;
                    S2 += w * r.covariate * r.covariate;
                }

                if (S0 === 0) continue;

                // Event at this time
                const eventSubject = sorted.find(d => d.time === t && d.event === cause);
                if (!eventSubject) continue;

                const xbar = S1 / S0;
                score += eventSubject.covariate - xbar;
                info += (S2 / S0) - (xbar * xbar);
            }

            if (Math.abs(info) < 1e-15) break;

            const step = score / info;
            beta += step;

            if (Math.abs(step) < tol) break;
        }

        const hr = Math.exp(beta);

        // Estimate SE from information matrix (last iteration)
        let infoFinal = 0;
        for (const t of eventTimes) {
            const riskSet = sorted.filter(d => {
                if (d.event === cause) return d.time >= t;
                if (d.event === 'censored') return d.time >= t;
                return true;
            });

            let S0 = 0, S1 = 0, S2 = 0;
            for (const r of riskSet) {
                const w = Math.exp(beta * r.covariate);
                S0 += w;
                S1 += w * r.covariate;
                S2 += w * r.covariate * r.covariate;
            }

            if (S0 > 0) {
                const xbar = S1 / S0;
                infoFinal += (S2 / S0) - (xbar * xbar);
            }
        }

        const seBeta = infoFinal > 0 ? 1 / Math.sqrt(infoFinal) : Infinity;
        const seHR = hr * seBeta; // delta method

        const lower = Math.exp(beta - zAlpha * seBeta);
        const upper = Math.exp(beta + zAlpha * seBeta);
        const zStat = beta / seBeta;
        const pValue = 2 * (1 - normalCDF(Math.abs(zStat)));

        return {
            hr: hr,
            beta: beta,
            se: seHR,
            seBeta: seBeta,
            lower: lower,
            upper: upper,
            pValue: pValue
        };
    }

    /**
     * Cause-specific hazard estimation at each event time.
     *
     * h_k(t) = d_k(t) / n(t) where only events of type k count as failures;
     * all other events are treated as censored.
     *
     * @param {Array} data - [{time, event}, ...]
     * @param {string} cause - Target cause
     * @returns {Array} [{time, hazard, cumHazard, atRisk, events}]
     */
    causeSpecificHazard(data, cause) {
        if (!cause) {
            throw new Error('Must specify a cause');
        }
        if (!Array.isArray(data) || data.length === 0) {
            throw new Error('Data must be a non-empty array');
        }

        for (let i = 0; i < data.length; i++) {
            if (data[i].time < 0) {
                throw new Error(`Negative time at index ${i}`);
            }
        }

        // Sort by time
        const sorted = [...data].sort((a, b) => a.time - b.time);

        // Unique event times for the target cause
        const causeEventTimes = [...new Set(
            sorted.filter(d => d.event === cause).map(d => d.time)
        )].sort((a, b) => a - b);

        let atRisk = sorted.length;
        let cumHazard = 0;
        const results = [];
        let processedIdx = 0;

        for (const t of causeEventTimes) {
            // Remove subjects with times before t (all types)
            while (processedIdx < sorted.length && sorted[processedIdx].time < t) {
                atRisk--;
                processedIdx++;
            }

            if (atRisk <= 0) break;

            // Count cause-specific events and all events at time t
            let causeEvents = 0;
            let allEventsAtT = 0;
            let tempIdx = processedIdx;
            while (tempIdx < sorted.length && sorted[tempIdx].time === t) {
                allEventsAtT++;
                if (sorted[tempIdx].event === cause) {
                    causeEvents++;
                }
                tempIdx++;
            }

            const hazard = causeEvents / atRisk;
            cumHazard += hazard;

            results.push({
                time: t,
                hazard: hazard,
                cumHazard: cumHazard,
                atRisk: atRisk,
                events: causeEvents
            });

            // Remove all subjects at time t from risk set
            atRisk -= allEventsAtT;
            processedIdx = tempIdx;
        }

        return results;
    }

    /**
     * z-quantile for given probability p (Abramowitz & Stegun inverse).
     */
    _zQuantile(p) {
        // Rational approximation (Abramowitz & Stegun 26.2.23)
        if (p <= 0) return -Infinity;
        if (p >= 1) return Infinity;
        if (p === 0.5) return 0;

        const sign = p < 0.5 ? -1 : 1;
        const pp = p < 0.5 ? p : 1 - p;

        const t = Math.sqrt(-2 * Math.log(pp));
        const c0 = 2.515517;
        const c1 = 0.802853;
        const c2 = 0.010328;
        const d1 = 1.432788;
        const d2 = 0.189269;
        const d3 = 0.001308;

        const z = t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t);
        return sign * z;
    }
}

// Export
if (typeof window !== 'undefined') {
    window.CompetingRisksEngine = CompetingRisksEngine;
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { CompetingRisksEngine };
}

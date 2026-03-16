/**
 * Model Averaging Engine
 * Bayesian Model Averaging for structural uncertainty in HTA models.
 *
 * Features:
 * - BIC/AIC/DIC-based posterior model weights
 * - Weighted model averaging with uncertainty
 * - Parametric survival distribution fitting (MLE via Newton-Raphson)
 * - Model-averaged survival predictions
 */

var OmanGuidanceRef = (function resolveOmanGuidance() {
    if (typeof globalThis !== 'undefined' && globalThis.OmanHTAGuidance) {
        return globalThis.OmanHTAGuidance;
    }
    if (typeof require === 'function') {
        try {
            return require('../utils/omanGuidance');
        } catch (err) {
            return null;
        }
    }
    return null;
})();

// ============ NAMED CONSTANTS ============
const NR_MAX_ITER = 200;
const NR_TOLERANCE = 1e-8;
const NR_STEP_SHRINK = 0.5;
const NORMAL_Z_975 = 1.959964;

class ModelAveragingEngine {
    constructor(options = {}) {
        this.options = options;
    }

    /**
     * Compute information-criterion weights from an array of IC values.
     * Generic helper used by bicWeights, aicWeights, dicWeights.
     *
     * @param {Object[]} models - [{name, ic}, ...] where ic = BIC/AIC/DIC
     * @param {string} icField - field name for the IC value
     * @returns {Object[]} [{name, ic, deltaIC, weight}, ...]
     */
    _icWeights(models, icField) {
        if (!models || models.length === 0) return [];

        const icValues = models.map(m => m[icField]);
        const minIC = Math.min(...icValues);

        const deltas = icValues.map(v => v - minIC);
        const rawWeights = deltas.map(d => Math.exp(-0.5 * d));
        const sumWeights = rawWeights.reduce((a, b) => a + b, 0);

        return models.map((m, i) => ({
            name: m.name,
            [icField]: m[icField],
            [`delta${icField.toUpperCase()}`]: deltas[i],
            weight: sumWeights > 0 ? rawWeights[i] / sumWeights : 1 / models.length
        }));
    }

    /**
     * Compute posterior model weights from BIC values.
     * Weight_i = exp(-0.5 * deltaBIC_i) / sum(exp(-0.5 * deltaBIC_j))
     *
     * @param {Object[]} models - [{name, bic}, ...]
     * @returns {Object[]} [{name, bic, deltaBIC, weight}, ...]
     */
    bicWeights(models) {
        return this._icWeights(models, 'bic');
    }

    /**
     * Compute posterior model weights from AIC values.
     *
     * @param {Object[]} models - [{name, aic}, ...]
     * @returns {Object[]} [{name, aic, deltaAIC, weight}, ...]
     */
    aicWeights(models) {
        return this._icWeights(models, 'aic');
    }

    /**
     * Compute posterior model weights from DIC values.
     *
     * @param {Object[]} models - [{name, dic}, ...]
     * @returns {Object[]} [{name, dic, deltaDIC, weight}, ...]
     */
    dicWeights(models) {
        return this._icWeights(models, 'dic');
    }

    /**
     * Model-averaged predictions: weighted combination of model outputs.
     *
     * @param {Object[]} models - [{name, predictions: [y1, y2, ...], weight}, ...]
     * @returns {Object} {averaged, uncertainty, credibleInterval}
     */
    modelAverage(models) {
        if (!models || models.length === 0) {
            return { averaged: [], uncertainty: [], credibleInterval: [] };
        }

        const nPred = models[0].predictions.length;
        const averaged = new Array(nPred).fill(0);
        const uncertainty = new Array(nPred).fill(0);

        // Weighted mean
        let totalWeight = 0;
        for (const m of models) {
            totalWeight += m.weight;
        }

        for (const m of models) {
            const w = totalWeight > 0 ? m.weight / totalWeight : 1 / models.length;
            for (let j = 0; j < nPred; j++) {
                averaged[j] += w * m.predictions[j];
            }
        }

        // Between-model variance (weighted)
        for (const m of models) {
            const w = totalWeight > 0 ? m.weight / totalWeight : 1 / models.length;
            for (let j = 0; j < nPred; j++) {
                const diff = m.predictions[j] - averaged[j];
                uncertainty[j] += w * diff * diff;
            }
        }

        // Convert variance to SD
        for (let j = 0; j < nPred; j++) {
            uncertainty[j] = Math.sqrt(uncertainty[j]);
        }

        // 95% credible interval using normal approximation
        const credibleInterval = averaged.map((avg, j) => ({
            lower: avg - NORMAL_Z_975 * uncertainty[j],
            upper: avg + NORMAL_Z_975 * uncertainty[j]
        }));

        return { averaged, uncertainty, credibleInterval };
    }

    // ========================================================================
    // Parametric survival distribution fitting (MLE)
    // ========================================================================

    /**
     * Fit parametric distributions to survival data and return comparison table.
     *
     * @param {Object[]} data - [{time, event}, ...] where event = 1 (failed) or 0 (censored)
     * @param {string[]} distributions - e.g. ['weibull', 'lognormal', 'loglogistic', 'exponential']
     * @returns {Object[]} [{name, params, logLik, aic, bic, weight}, ...]
     */
    fitCompare(data, distributions) {
        const results = [];
        const n = data.length;

        for (const dist of distributions) {
            try {
                const fit = this._fitDistribution(data, dist);
                const nParams = Object.keys(fit.params).length;
                const aic = -2 * fit.logLik + 2 * nParams;
                const bic = -2 * fit.logLik + nParams * Math.log(n);

                results.push({
                    name: dist,
                    params: fit.params,
                    logLik: fit.logLik,
                    aic,
                    bic,
                    weight: 0  // will be computed below
                });
            } catch (e) {
                // Skip failed fits
                results.push({
                    name: dist,
                    params: null,
                    logLik: -Infinity,
                    aic: Infinity,
                    bic: Infinity,
                    weight: 0,
                    error: e.message
                });
            }
        }

        // Compute AIC weights
        const validResults = results.filter(r => isFinite(r.aic));
        if (validResults.length > 0) {
            const minAIC = Math.min(...validResults.map(r => r.aic));
            const rawWeights = validResults.map(r => Math.exp(-0.5 * (r.aic - minAIC)));
            const sumW = rawWeights.reduce((a, b) => a + b, 0);
            for (let i = 0; i < validResults.length; i++) {
                validResults[i].weight = sumW > 0 ? rawWeights[i] / sumW : 1 / validResults.length;
            }
        }

        return results;
    }

    /**
     * Fit a single distribution to survival data via MLE.
     */
    _fitDistribution(data, dist) {
        switch (dist) {
            case 'exponential': return this._fitExponential(data);
            case 'weibull': return this._fitWeibull(data);
            case 'lognormal': return this._fitLognormal(data);
            case 'loglogistic': return this._fitLoglogistic(data);
            case 'gamma': return this._fitGamma(data);
            case 'gompertz': return this._fitGompertz(data);
            default: throw new Error(`Unknown distribution: ${dist}`);
        }
    }

    /**
     * Exponential distribution MLE. S(t) = exp(-lambda * t)
     * Closed-form: lambda = sum(event) / sum(time)
     */
    _fitExponential(data) {
        const totalEvents = data.reduce((s, d) => s + d.event, 0);
        const totalTime = data.reduce((s, d) => s + d.time, 0);
        const lambda = totalEvents / totalTime;

        // Log-likelihood: sum[ event*log(lambda) - lambda*t ]
        let logLik = 0;
        for (const d of data) {
            logLik += d.event * Math.log(lambda) - lambda * d.time;
        }

        return { params: { lambda }, logLik };
    }

    /**
     * Weibull distribution MLE. S(t) = exp(-(t/scale)^shape)
     * h(t) = (shape/scale) * (t/scale)^(shape-1)
     * Newton-Raphson on (log_shape, log_scale) for positivity.
     */
    _fitWeibull(data) {
        // Initial estimates
        const totalEvents = data.reduce((s, d) => s + d.event, 0);
        const totalTime = data.reduce((s, d) => s + d.time, 0);
        const meanTime = totalTime / data.length;
        let shape = 1.0;
        let scale = meanTime;

        // Newton-Raphson on log-parameterization
        let logShape = Math.log(shape);
        let logScale = Math.log(scale);

        for (let iter = 0; iter < NR_MAX_ITER; iter++) {
            shape = Math.exp(logShape);
            scale = Math.exp(logScale);

            let dll_dshape = 0, dll_dscale = 0;
            let d2ll_dshape2 = 0, d2ll_dscale2 = 0, d2ll_dshapescale = 0;

            for (const d of data) {
                const t = d.time;
                const e = d.event;
                if (t <= 0) continue;

                const logT = Math.log(t);
                const z = t / scale;
                const zk = Math.pow(z, shape);
                const logZ = logT - logScale;

                dll_dshape += e * (1 / shape + logZ) - zk * logZ;
                dll_dscale += (-e * shape / scale) + shape * zk / scale;

                d2ll_dshape2 += -e / (shape * shape) - zk * logZ * logZ;
                d2ll_dscale2 += e * shape / (scale * scale) - shape * (shape + 1) * zk / (scale * scale);
                d2ll_dshapescale += -e / scale + zk * (1 + shape * logZ) / scale;
            }

            // Chain rule for log-parameterization
            const g1 = dll_dshape * shape;
            const g2 = dll_dscale * scale;

            const norm = Math.abs(g1) + Math.abs(g2);
            if (norm < NR_TOLERANCE) break;

            // Simple gradient step (avoid full Hessian inversion complexity)
            const stepSize = 0.01;
            logShape += stepSize * g1 / (Math.abs(g1) + 1);
            logScale += stepSize * g2 / (Math.abs(g2) + 1);

            // Clamp for numerical stability
            logShape = Math.max(-5, Math.min(5, logShape));
            logScale = Math.max(-10, Math.min(20, logScale));
        }

        shape = Math.exp(logShape);
        scale = Math.exp(logScale);

        // Compute log-likelihood
        let logLik = 0;
        for (const d of data) {
            if (d.time <= 0) continue;
            const z = d.time / scale;
            const zk = Math.pow(z, shape);
            logLik += d.event * (Math.log(shape) - Math.log(scale) + (shape - 1) * Math.log(z)) - zk;
        }

        return { params: { shape, scale }, logLik };
    }

    /**
     * Lognormal distribution MLE. log(T) ~ N(mu, sigma^2)
     * Closed-form for uncensored; iterative with censoring.
     */
    _fitLognormal(data) {
        // For simplicity, use closed-form on uncensored observations
        // and EM-style adjustment for censored
        const logTimes = data.filter(d => d.time > 0).map(d => ({
            logT: Math.log(d.time),
            event: d.event
        }));

        const events = logTimes.filter(d => d.event === 1);
        if (events.length === 0) {
            return { params: { mu: 0, sigma: 1 }, logLik: -Infinity };
        }

        // Initial: MLE from uncensored only
        let mu = events.reduce((s, d) => s + d.logT, 0) / events.length;
        let sigma2 = events.reduce((s, d) => s + (d.logT - mu) ** 2, 0) / events.length;
        let sigma = Math.sqrt(Math.max(sigma2, 0.01));

        // Refine with censored data using gradient ascent
        for (let iter = 0; iter < NR_MAX_ITER; iter++) {
            let g_mu = 0, g_sigma = 0;

            for (const d of logTimes) {
                const z = (d.logT - mu) / sigma;
                if (d.event === 1) {
                    g_mu += z / sigma;
                    g_sigma += (z * z - 1) / sigma;
                } else {
                    // Censored: contribution from survival function
                    const phi = this._normalPDF(z);
                    const Phi = this._normalCDF(z);
                    const survProb = 1 - Phi;
                    if (survProb > 1e-15) {
                        const ratio = phi / survProb;
                        g_mu += -ratio / sigma;
                        g_sigma += -ratio * z / sigma;
                    }
                }
            }

            const norm = Math.abs(g_mu) + Math.abs(g_sigma);
            if (norm < NR_TOLERANCE) break;

            const step = 0.01;
            mu += step * g_mu;
            sigma += step * g_sigma;
            sigma = Math.max(sigma, 0.01);
        }

        // Log-likelihood
        let logLik = 0;
        for (const d of logTimes) {
            const z = (d.logT - mu) / sigma;
            if (d.event === 1) {
                logLik += -0.5 * Math.log(2 * Math.PI) - Math.log(sigma) - d.logT - 0.5 * z * z;
            } else {
                const survProb = 1 - this._normalCDF(z);
                logLik += Math.log(Math.max(survProb, 1e-300));
            }
        }

        return { params: { mu, sigma }, logLik };
    }

    /**
     * Log-logistic distribution MLE. S(t) = 1 / (1 + (t/alpha)^beta)
     */
    _fitLoglogistic(data) {
        const events = data.filter(d => d.event === 1 && d.time > 0);
        if (events.length === 0) {
            return { params: { alpha: 1, beta: 1 }, logLik: -Infinity };
        }

        // Initial estimates from log-logistic relationship
        const logTimes = events.map(d => Math.log(d.time));
        const meanLogT = logTimes.reduce((a, b) => a + b, 0) / logTimes.length;
        let alpha = Math.exp(meanLogT);
        let beta = 1.0;

        let logAlpha = Math.log(alpha);
        let logBeta = Math.log(beta);

        for (let iter = 0; iter < NR_MAX_ITER; iter++) {
            alpha = Math.exp(logAlpha);
            beta = Math.exp(logBeta);

            let g_alpha = 0, g_beta = 0;

            for (const d of data) {
                if (d.time <= 0) continue;
                const z = Math.pow(d.time / alpha, beta);
                const logZ = beta * (Math.log(d.time) - logAlpha);

                if (d.event === 1) {
                    // f(t) = (beta/alpha)(t/alpha)^(beta-1) / (1+(t/alpha)^beta)^2
                    g_alpha += -beta / alpha + 2 * beta * z / (alpha * (1 + z));
                    g_beta += 1 / beta + Math.log(d.time / alpha) - 2 * z * Math.log(d.time / alpha) / (1 + z);
                } else {
                    // S(t) = 1/(1+z)
                    g_alpha += beta * z / (alpha * (1 + z));
                    g_beta += -z * Math.log(d.time / alpha) / (1 + z);
                }
            }

            const norm = Math.abs(g_alpha * alpha) + Math.abs(g_beta * beta);
            if (norm < NR_TOLERANCE) break;

            const step = 0.01;
            logAlpha += step * g_alpha * alpha / (Math.abs(g_alpha * alpha) + 1);
            logBeta += step * g_beta * beta / (Math.abs(g_beta * beta) + 1);

            logAlpha = Math.max(-10, Math.min(20, logAlpha));
            logBeta = Math.max(-5, Math.min(5, logBeta));
        }

        alpha = Math.exp(logAlpha);
        beta = Math.exp(logBeta);

        // Log-likelihood
        let logLik = 0;
        for (const d of data) {
            if (d.time <= 0) continue;
            const z = Math.pow(d.time / alpha, beta);
            if (d.event === 1) {
                logLik += Math.log(beta) - Math.log(alpha) + (beta - 1) * Math.log(d.time / alpha)
                          - 2 * Math.log(1 + z);
            } else {
                logLik += -Math.log(1 + z);
            }
        }

        return { params: { alpha, beta }, logLik };
    }

    /**
     * Gamma distribution MLE. f(t) = (lambda^k / Gamma(k)) * t^(k-1) * exp(-lambda*t)
     * Uses method of moments for initial values + gradient ascent.
     */
    _fitGamma(data) {
        const events = data.filter(d => d.event === 1 && d.time > 0);
        if (events.length === 0) {
            return { params: { shape: 1, rate: 1 }, logLik: -Infinity };
        }

        const times = events.map(d => d.time);
        const mean = times.reduce((a, b) => a + b, 0) / times.length;
        const variance = times.reduce((a, b) => a + (b - mean) ** 2, 0) / times.length;

        let shape = variance > 0 ? (mean * mean) / variance : 1;
        let rate = variance > 0 ? mean / variance : 1;

        // Gradient ascent on (log_shape, log_rate)
        let logShape = Math.log(Math.max(shape, 0.01));
        let logRate = Math.log(Math.max(rate, 0.01));

        for (let iter = 0; iter < NR_MAX_ITER; iter++) {
            shape = Math.exp(logShape);
            rate = Math.exp(logRate);

            let g_shape = 0, g_rate = 0;

            for (const d of data) {
                if (d.time <= 0) continue;
                if (d.event === 1) {
                    g_shape += Math.log(rate) + Math.log(d.time) - this._digamma(shape);
                    g_rate += shape / rate - d.time;
                } else {
                    // Censored: use survival function gradient (approximate)
                    const gammaVal = this._upperIncompleteGammaRatio(shape, rate * d.time);
                    if (gammaVal > 1e-15) {
                        // Simplified gradient for censored observations
                        g_rate += -d.time * Math.exp(-rate * d.time) * Math.pow(rate * d.time, shape - 1)
                                  / (this._gammaFn(shape) * gammaVal);
                    }
                }
            }

            const norm = Math.abs(g_shape * shape) + Math.abs(g_rate * rate);
            if (norm < NR_TOLERANCE) break;

            const step = 0.005;
            logShape += step * g_shape * shape / (Math.abs(g_shape * shape) + 1);
            logRate += step * g_rate * rate / (Math.abs(g_rate * rate) + 1);

            logShape = Math.max(-5, Math.min(10, logShape));
            logRate = Math.max(-10, Math.min(10, logRate));
        }

        shape = Math.exp(logShape);
        rate = Math.exp(logRate);

        // Log-likelihood
        let logLik = 0;
        for (const d of data) {
            if (d.time <= 0) continue;
            if (d.event === 1) {
                logLik += shape * Math.log(rate) - this._logGamma(shape)
                          + (shape - 1) * Math.log(d.time) - rate * d.time;
            } else {
                const gammaRatio = this._upperIncompleteGammaRatio(shape, rate * d.time);
                logLik += Math.log(Math.max(gammaRatio, 1e-300));
            }
        }

        return { params: { shape, rate }, logLik };
    }

    /**
     * Gompertz distribution MLE. h(t) = b * exp(eta * t)
     * S(t) = exp(-(b/eta)*(exp(eta*t) - 1))
     */
    _fitGompertz(data) {
        const events = data.filter(d => d.event === 1 && d.time > 0);
        if (events.length === 0) {
            return { params: { b: 0.01, eta: 0.01 }, logLik: -Infinity };
        }

        let b = 0.01;
        let eta = 0.01;

        for (let iter = 0; iter < NR_MAX_ITER; iter++) {
            let g_b = 0, g_eta = 0;

            for (const d of data) {
                if (d.time <= 0) continue;
                const expEtaT = Math.exp(eta * d.time);
                const cumHaz = (b / eta) * (expEtaT - 1);

                if (d.event === 1) {
                    g_b += 1 / b - (expEtaT - 1) / eta;
                    g_eta += d.time - b * (d.time * expEtaT * eta - expEtaT + 1) / (eta * eta);
                } else {
                    g_b += -(expEtaT - 1) / eta;
                    g_eta += -b * (d.time * expEtaT * eta - expEtaT + 1) / (eta * eta);
                }
            }

            const norm = Math.abs(g_b) + Math.abs(g_eta);
            if (norm < NR_TOLERANCE) break;

            const step = 0.0005;
            b += step * g_b;
            eta += step * g_eta;

            b = Math.max(1e-6, b);
            eta = Math.max(1e-6, eta);
        }

        // Log-likelihood
        let logLik = 0;
        for (const d of data) {
            if (d.time <= 0) continue;
            const expEtaT = Math.exp(eta * d.time);
            const cumHaz = (b / eta) * (expEtaT - 1);

            if (d.event === 1) {
                logLik += Math.log(b) + eta * d.time - cumHaz;
            } else {
                logLik += -cumHaz;
            }
        }

        return { params: { b, eta }, logLik };
    }

    /**
     * Model-averaged survival prediction.
     *
     * @param {Object[]} fittedModels - [{name, params, weight, distribution}, ...]
     * @param {number[]} times - time points for prediction
     * @returns {Object} {times, survival, ci}
     */
    survivalPrediction(fittedModels, times) {
        const nT = times.length;
        const survCurves = [];

        for (const model of fittedModels) {
            const curve = times.map(t => this._survivalFunction(model.distribution ?? model.name, model.params, t));
            survCurves.push({ curve, weight: model.weight });
        }

        // Weighted average
        const survival = new Array(nT).fill(0);
        const variance = new Array(nT).fill(0);
        let totalWeight = survCurves.reduce((s, c) => s + c.weight, 0);

        for (const sc of survCurves) {
            const w = totalWeight > 0 ? sc.weight / totalWeight : 1 / survCurves.length;
            for (let j = 0; j < nT; j++) {
                survival[j] += w * sc.curve[j];
            }
        }

        // Between-model variance
        for (const sc of survCurves) {
            const w = totalWeight > 0 ? sc.weight / totalWeight : 1 / survCurves.length;
            for (let j = 0; j < nT; j++) {
                variance[j] += w * (sc.curve[j] - survival[j]) ** 2;
            }
        }

        const ci = survival.map((s, j) => {
            const sd = Math.sqrt(variance[j]);
            return {
                lower: Math.max(0, s - NORMAL_Z_975 * sd),
                upper: Math.min(1, s + NORMAL_Z_975 * sd)
            };
        });

        return { times, survival, ci };
    }

    /**
     * Evaluate survival function S(t) for a given distribution.
     */
    _survivalFunction(dist, params, t) {
        if (t <= 0) return 1.0;

        switch (dist) {
            case 'exponential':
                return Math.exp(-params.lambda * t);
            case 'weibull':
                return Math.exp(-Math.pow(t / params.scale, params.shape));
            case 'lognormal': {
                const z = (Math.log(t) - params.mu) / params.sigma;
                return 1 - this._normalCDF(z);
            }
            case 'loglogistic': {
                const z = Math.pow(t / params.alpha, params.beta);
                return 1 / (1 + z);
            }
            case 'gamma': {
                return this._upperIncompleteGammaRatio(params.shape, params.rate * t);
            }
            case 'gompertz': {
                return Math.exp(-(params.b / params.eta) * (Math.exp(params.eta * t) - 1));
            }
            default:
                return 0;
        }
    }

    // ========================================================================
    // Mathematical helper functions
    // ========================================================================

    /** Standard normal PDF */
    _normalPDF(z) {
        return Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);
    }

    /** Standard normal CDF (Abramowitz & Stegun approximation) */
    _normalCDF(z) {
        if (z < -8) return 0;
        if (z > 8) return 1;

        const a1 = 0.254829592;
        const a2 = -0.284496736;
        const a3 = 1.421413741;
        const a4 = -1.453152027;
        const a5 = 1.061405429;
        const p = 0.3275911;

        const sign = z < 0 ? -1 : 1;
        const x = Math.abs(z) / Math.sqrt(2);
        const t = 1.0 / (1.0 + p * x);
        const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

        return 0.5 * (1.0 + sign * y);
    }

    /** Log-gamma function (Stirling approximation for large values, Lanczos for small) */
    _logGamma(x) {
        if (x <= 0) return Infinity;

        // Lanczos approximation
        const g = 7;
        const c = [
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

        let sum = c[0];
        for (let i = 1; i < g + 2; i++) {
            sum += c[i] / (x + i - 1);
        }

        const t = x + g - 0.5;
        return 0.5 * Math.log(2 * Math.PI) + (x - 0.5) * Math.log(t) - t + Math.log(sum);
    }

    /** Gamma function */
    _gammaFn(x) {
        return Math.exp(this._logGamma(x));
    }

    /** Digamma function (psi) — series approximation */
    _digamma(x) {
        if (x <= 0) return -Infinity;

        let result = 0;
        // Shift x to be large enough for asymptotic expansion
        while (x < 6) {
            result -= 1 / x;
            x += 1;
        }

        // Asymptotic expansion
        result += Math.log(x) - 1 / (2 * x);
        const x2 = x * x;
        result -= 1 / (12 * x2);
        result += 1 / (120 * x2 * x2);
        result -= 1 / (252 * x2 * x2 * x2);

        return result;
    }

    /**
     * Upper incomplete gamma ratio: Q(a, x) = 1 - P(a, x) = Gamma(a,x)/Gamma(a)
     * Uses series or continued fraction depending on values.
     */
    _upperIncompleteGammaRatio(a, x) {
        if (x < 0) return 1;
        if (x === 0) return 1;
        if (a <= 0) return 0;

        // Use series for x < a+1, continued fraction otherwise
        if (x < a + 1) {
            return 1 - this._lowerIncompleteGammaSeries(a, x);
        } else {
            return this._upperIncompleteGammaCF(a, x);
        }
    }

    /** Lower incomplete gamma ratio via series expansion */
    _lowerIncompleteGammaSeries(a, x) {
        const logGammaA = this._logGamma(a);
        let sum = 1.0 / a;
        let term = 1.0 / a;

        for (let n = 1; n < 200; n++) {
            term *= x / (a + n);
            sum += term;
            if (Math.abs(term) < 1e-15 * Math.abs(sum)) break;
        }

        return sum * Math.exp(-x + a * Math.log(x) - logGammaA);
    }

    /** Upper incomplete gamma ratio via continued fraction (Lentz's method) */
    _upperIncompleteGammaCF(a, x) {
        const logGammaA = this._logGamma(a);
        let f = 1e-30;
        let c = 1e-30;
        let d = 1 / (x + 1 - a);
        let h = d;

        for (let n = 1; n < 200; n++) {
            const an = n * (a - n);
            const bn = x + 2 * n + 1 - a;
            d = bn + an * d;
            if (Math.abs(d) < 1e-30) d = 1e-30;
            c = bn + an / c;
            if (Math.abs(c) < 1e-30) c = 1e-30;
            d = 1 / d;
            const delta = d * c;
            h *= delta;
            if (Math.abs(delta - 1) < 1e-15) break;
        }

        return Math.exp(-x + a * Math.log(x) - logGammaA) * h;
    }
}

// Export
if (typeof window !== 'undefined') {
    window.ModelAveragingEngine = ModelAveragingEngine;
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { ModelAveragingEngine };
}

// truth-recovery/dgp-known-truth.mjs
// Standalone, seeded, KNOWN-TRUTH data-generating processes + closed-form oracles.
// No dependency on the repo under test -- these are independent analytic references.

// ---- Seeded PRNG (mulberry32) for reproducible PSA samples ----
export function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
        a |= 0; a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// Box-Muller standard normal from a uniform PRNG.
export function makeNormal(rng) {
    let spare = null;
    return function () {
        if (spare !== null) { const v = spare; spare = null; return v; }
        let u = 0, v = 0;
        while (u === 0) u = rng();
        while (v === 0) v = rng();
        const mag = Math.sqrt(-2 * Math.log(u));
        spare = mag * Math.sin(2 * Math.PI * v);
        return mag * Math.cos(2 * Math.PI * v);
    };
}

// ---- Markov oracle: 2-state Healthy/Dead with constant per-cycle death prob p ----
// Occupancy in Healthy at integer cycle t is exactly (1-p)^t (geometric decay).
export function markovHealthyTrace(p, cycles) {
    const trace = [];
    for (let t = 0; t <= cycles; t++) trace.push(Math.pow(1 - p, t));
    return trace;
}

// Discounted, half-cycle-corrected life-years for the 2-state model.
// LY = sum_{t=0..N} w_t * (1-p)^t * (1+r)^{-t}, with trapezoidal HCC giving
// w_0 = w_N = 0.5, else 1. cycle_length = 1.
export function markovLifeYears(p, cycles, r, hcc) {
    let ly = 0;
    for (let t = 0; t <= cycles; t++) {
        let w = 1;
        if (hcc === 'trapezoidal') { if (t === 0 || t === cycles) w = 0.5; }
        else if (hcc === 'start') { if (t === 0) w = 0.5; }
        else if (hcc === 'end') { if (t === cycles) w = 0.5; }
        ly += w * Math.pow(1 - p, t) * Math.pow(1 + r, -t);
    }
    return ly;
}

// ---- EVPI oracle: closed form for a Normal incremental net benefit ----
// If INB ~ Normal(mu, sigma^2), EVPI = sigma * L(|mu|/sigma) where the unit-normal
// loss L(z) = phi(z) - z*(1-Phi(z)).  (Standard partial-EVPI / value-of-info result.)
function stdNormalPdf(z) { return Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI); }
function stdNormalCdf(z) {
    // Abramowitz & Stegun 7.1.26 erf approximation.
    const t = 1 / (1 + 0.3275911 * Math.abs(z) / Math.SQRT2);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-z * z / 2);
    const cdf = 0.5 * (1 + (z >= 0 ? y : -y));
    return cdf;
}
export function evpiNormalClosedForm(mu, sigma) {
    if (sigma <= 0) return 0;
    const z = Math.abs(mu) / sigma;
    const L = stdNormalPdf(z) - z * (1 - stdNormalCdf(z));
    return sigma * L;
}

// Generate seeded PSA scatter whose INB = wtp*dQ - dC is Normal(muINB, sigmaINB^2).
// We set dQ ~ Normal(mq, sq^2), dC = 0 (deterministic) so INB = wtp*dQ,
// hence muINB = wtp*mq, sigmaINB = wtp*sq.
export function makePsaResultsNormalINB({ n, wtp, mq, sq, seed }) {
    const rng = mulberry32(seed);
    const norm = makeNormal(rng);
    const incremental_qalys = [];
    const incremental_costs = [];
    for (let i = 0; i < n; i++) {
        incremental_qalys.push(mq + sq * norm());
        incremental_costs.push(0);
    }
    return {
        settings_snapshot: { wtp_thresholds: [wtp] },
        primary_wtp: wtp,
        scatter: { incremental_qalys, incremental_costs }
    };
}

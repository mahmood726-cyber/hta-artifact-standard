// truth-recovery/harness.test.mjs
// Known-truth validation of the repo's OWN engine functions.
// Run: node --test truth-recovery/harness.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';

import { MarkovEngine, KahanSum, ExpressionParser, EVPICalculator } from './engine.mjs';
import {
    markovHealthyTrace, markovLifeYears,
    evpiNormalClosedForm, makePsaResultsNormalINB
} from './dgp-known-truth.mjs';

function build2State(p, { horizon = 20, r = 0, hcc = 'none' } = {}) {
    return {
        settings: {
            time_horizon: horizon, cycle_length: 1,
            discount_rate_costs: r, discount_rate_qalys: r,
            half_cycle_correction: hcc, starting_age: 50
        },
        parameters: {},
        states: {
            H: { label: 'Healthy', initial_probability: 1, cost: 0, utility: 1 },
            D: { label: 'Dead', type: 'absorbing', initial_probability: 0, cost: 0, utility: 0 }
        },
        transitions: {
            t1: { from: 'H', to: 'D', probability: p },
            t2: { from: 'H', to: 'H', probability: 'complement' },
            t3: { from: 'D', to: 'D', probability: 1 }
        }
    };
}

// 1) Cohort trace matches the closed-form geometric decay (1-p)^t EXACTLY.
test('Markov cohort trace == (1-p)^t closed form', () => {
    const eng = new MarkovEngine({ KahanSum, ExpressionParser });
    const p = 0.07, cycles = 20;
    const res = eng.run(build2State(p, { horizon: cycles }));
    const oracle = markovHealthyTrace(p, cycles);
    let maxErr = 0;
    for (let t = 0; t <= cycles; t++) {
        maxErr = Math.max(maxErr, Math.abs(res.trace.states.H[t] - oracle[t]));
    }
    assert.ok(maxErr < 1e-12, `max trace error ${maxErr} should be < 1e-12`);
    // mass conservation: H + D == 1 every cycle
    for (let t = 0; t <= cycles; t++) {
        assert.ok(Math.abs(res.trace.states.H[t] + res.trace.states.D[t] - 1) < 1e-12);
    }
});

// 2) Undiscounted, no-HCC life-years == analytic sum of occupancy.
test('Markov life-years == analytic sum (no discount, no HCC)', () => {
    const eng = new MarkovEngine({ KahanSum, ExpressionParser });
    const p = 0.1, cycles = 30;
    const res = eng.run(build2State(p, { horizon: cycles, r: 0, hcc: 'none' }));
    const oracle = markovLifeYears(p, cycles, 0, 'none');
    assert.ok(Math.abs(res.life_years - oracle) < 1e-9,
        `LY ${res.life_years} vs oracle ${oracle}`);
});

// 3) Discounted + trapezoidal HCC life-years match the discounted closed form.
test('Markov life-years == discounted + half-cycle-corrected closed form', () => {
    const eng = new MarkovEngine({ KahanSum, ExpressionParser });
    const p = 0.05, cycles = 40, r = 0.03;
    const res = eng.run(build2State(p, { horizon: cycles, r, hcc: 'trapezoidal' }));
    const oracle = markovLifeYears(p, cycles, r, 'trapezoidal');
    assert.ok(Math.abs(res.life_years - oracle) < 1e-9,
        `LY ${res.life_years} vs oracle ${oracle}`);
});

// 4) EVPI == 0 when the decision NEVER changes across PSA (all INB same sign).
test('EVPI == 0 when decision never changes (dominant adopt)', () => {
    const calc = new EVPICalculator();
    // mq large + small sd => INB = wtp*dQ always > 0 with seed -> all adopt.
    const psa = makePsaResultsNormalINB({ n: 5000, wtp: 30000, mq: 0.20, sq: 0.01, seed: 12345 });
    const r = calc.calculate(psa, 30000, 1, 1);
    // EVPI is non-negative by construction and effectively zero here.
    assert.ok(r.evpiPerPatient >= 0, 'EVPI must be >= 0');
    assert.ok(r.evpiPerPatient < 1e-6, `EVPI ${r.evpiPerPatient} should be ~0 when no decision flips`);
    assert.equal(r.probWrongDecision, 0);
});

// 5) EVPI converges to the unit-normal-loss closed form for a Normal INB.
test('EVPI converges to closed-form normal value-of-information', () => {
    const calc = new EVPICalculator();
    const wtp = 30000, mq = 0.01, sq = 0.05, n = 200000;
    const psa = makePsaResultsNormalINB({ n, wtp, mq, sq, seed: 2026 });
    const r = calc.calculate(psa, wtp, 1, 1);
    // INB ~ Normal(mu = wtp*mq, sigma = wtp*sq)
    const mu = wtp * mq, sigma = wtp * sq;
    const oracle = evpiNormalClosedForm(mu, sigma);
    const relErr = Math.abs(r.evpiPerPatient - oracle) / oracle;
    assert.ok(r.evpiPerPatient >= 0, 'EVPI must be >= 0');
    assert.ok(relErr < 0.03, `EVPI MC ${r.evpiPerPatient.toFixed(2)} vs closed form ${oracle.toFixed(2)} relErr ${(relErr*100).toFixed(2)}% should be < 3%`);
});

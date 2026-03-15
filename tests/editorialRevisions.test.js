/**
 * Jest-based tests for editorial revisions modules.
 * Replaces previous custom test runner to avoid process.exit side effects.
 */

'use strict';

const {
    HKSJMetaAnalysis,
    EVPPICalculator,
    SurvivalModelSelection,
    NetworkMetaAnalysis,
    PublicationBiasTests,
    NumericalValidation
} = require('../src/engine/editorialRevisions');

const { OptimizedAlgorithms, FastMath } = require('../src/engine/performanceWrapper');

describe('Editorial Revisions', () => {
    test('HKSJ returns structurally valid outputs and CI width relation', () => {
        const effects = [0.5, 0.3, 0.7, 0.4, 0.6];
        const variances = [0.04, 0.09, 0.06, 0.08, 0.05];

        const hksj = new HKSJMetaAnalysis({ method: 'REML' });
        const result = hksj.analyze(effects, variances);

        expect(result.effect).toBeGreaterThan(0);
        expect(result.se).toBeGreaterThan(0);
        expect(result.tau2).toBeGreaterThanOrEqual(0);
        expect(result.adjustment).toBe('HKSJ');

        const hksjWidth = result.ci[1] - result.ci[0];
        const stdWidth = result.ciStandard[1] - result.ciStandard[0];
        expect(hksjWidth).toBeGreaterThanOrEqual(stdWidth - 1e-6);
    });

    test('EVPPI computes results for at least one parameter', async () => {
        const modelFn = (params) => ({
            cost: 1000 + params.costMultiplier * 250,
            effect: 0.5 + params.effectModifier * 0.2,
            strategy: 'intervention'
        });

        const allParams = {
            costMultiplier: { distribution: 'normal', mean: 1, se: 0.2 },
            effectModifier: { distribution: 'normal', mean: 0, se: 0.1 }
        };

        const evppi = new EVPPICalculator({ outerSamples: 50, innerSamples: 50 });
        const result = await evppi.calculate(modelFn, allParams, ['costMultiplier'], 50000);

        expect(result.evpi).toBeDefined();
        expect(result.evppiResults.costMultiplier).toBeDefined();
        expect(result.evppiResults.costMultiplier.evppi).toBeGreaterThanOrEqual(-1e-6);
    });

    test('publication bias methods return bounded p-values', () => {
        const effects = [0.5, 0.4, 0.6, 0.3, 0.7, 0.55, 0.45];
        const ses = [0.1, 0.15, 0.12, 0.18, 0.08, 0.14, 0.11];
        const variances = ses.map(se => se * se);

        const bias = new PublicationBiasTests();
        const egger = bias.eggerTest(effects, ses);
        const begg = bias.beggTest(effects, ses);
        const peters = bias.petersTest(effects, variances);

        expect(egger.pValue).toBeGreaterThanOrEqual(0);
        expect(egger.pValue).toBeLessThanOrEqual(1);
        expect(begg.pValue).toBeGreaterThanOrEqual(0);
        expect(begg.pValue).toBeLessThanOrEqual(1);
        expect(peters.pValue).toBeGreaterThanOrEqual(0);
        expect(peters.pValue).toBeLessThanOrEqual(1);
    });

    test('network meta-analysis flags disconnected networks', () => {
        const studies = [
            { study: 'S1', treat1: 'A', treat2: 'B', effect: 0.3, se: 0.1 },
            { study: 'S2', treat1: 'C', treat2: 'D', effect: 0.2, se: 0.15 }
        ];

        const nma = new NetworkMetaAnalysis();
        const result = nma.analyze(studies);
        expect(result.error).toBe('Network is not connected');
    });

    test('survival model selection returns recommendation and averaged predictions', () => {
        const times = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15];
        const events = [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1];

        const selector = new SurvivalModelSelection();
        const result = selector.compare(times, events);

        expect(result.models.length).toBeGreaterThan(0);
        expect(result.recommended).toBeDefined();
        expect(result.modelAveraging.length).toBeGreaterThan(0);
        expect(result.modelAveraging[0].survival).toBeGreaterThanOrEqual(0);
        expect(result.modelAveraging[0].survival).toBeLessThanOrEqual(1);
    });
});

describe('Performance wrapper numerical correctness', () => {
    test('optimized REML includes prediction interval and I² CI', () => {
        const effects = [0.5, 0.3, 0.7, 0.4, 0.6];
        const variances = [0.04, 0.09, 0.06, 0.08, 0.05];

        const dl = OptimizedAlgorithms.derSimonianLaird(effects, variances, { hksj: true });
        const reml = OptimizedAlgorithms.reml(effects, variances);

        expect(dl.I2CI).toBeDefined();
        expect(dl.adjustment).toBe('HKSJ');
        expect(reml.predictionInterval).toBeDefined();
        expect(reml.predictionInterval.lower).toBeDefined();
        expect(reml.predictionInterval.upper).toBeDefined();
    });

    test('FastMath distribution helpers are within tolerance', () => {
        expect(FastMath.tQuantile(0.975, 10)).toBeCloseTo(2.228, 2);
        expect(FastMath.tQuantile(0.975, 30)).toBeCloseTo(2.042, 2);
        expect(FastMath.chiSquaredCDF(3.84, 1)).toBeCloseTo(0.95, 2);
    });

    test('NumericalValidation generates a report with passing tests', () => {
        const validator = new NumericalValidation();
        const report = validator.runAllValidations();

        expect(report.summary.total).toBeGreaterThan(0);
        expect(report.summary.passed).toBeGreaterThan(0);
        expect(report.summary.passRate).toBeDefined();
    });
});

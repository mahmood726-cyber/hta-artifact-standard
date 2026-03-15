/**
 * Regression tests for src/parser/expression.js
 */

'use strict';

const { ExpressionParser } = require('../../src/parser/expression');

describe('ExpressionParser lazy evaluation semantics', () => {
    test('if() evaluates only the selected branch', () => {
        expect(ExpressionParser.evaluate('if(1, 42, 1/0)', {})).toBe(42);
        expect(ExpressionParser.evaluate('if(0, 1/0, 7)', {})).toBe(7);
    });

    test('if() still throws when selected branch is invalid', () => {
        expect(() => ExpressionParser.evaluate('if(1, 1/0, 5)', {})).toThrow('Division by zero');
    });

    test('logical operators short-circuit', () => {
        expect(ExpressionParser.evaluate('0 and (1 / 0)', {})).toBe(0);
        expect(ExpressionParser.evaluate('1 or (1 / 0)', {})).toBe(1);
    });
});

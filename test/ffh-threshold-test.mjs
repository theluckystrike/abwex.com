import { createRequire } from 'node:module';
import assert from 'node:assert/strict';
import fs from 'node:fs';

const require = createRequire(import.meta.url);
const ffh = require('../assets/js/ffh-core.js');
function loadJson(file) {
  try { return JSON.parse(fs.readFileSync(file)); }
  catch (error) { throw new Error('Cannot parse FFH fixtures', { cause: error }); }
}
const fixtures = loadJson(new URL('../tools/fisher-freeman-halton-calculator/ffh-fixtures.json', import.meta.url)).fixtures;

for (const fixture of fixtures) {
  const result=ffh.calculate(fixture.table);
  const parts=fixture.p_value_fraction.split('/');
  assert.equal(result.pValueNumerator,parts[0]);
  assert.equal(result.pValueDenominator,parts[1]||'1');
}

const balanced = ffh.calculate([[2,2,2],[2,2,2]], { alphaText: '0.999999999999999' });
assert.equal(balanced.pValue, 1);
assert.equal(balanced.pValueNumerator, balanced.pValueDenominator);
assert.equal(balanced.exactDecision, 'insufficient');

const sparse = ffh.calculate([[1,4,0],[3,0,2]], { alphaText: '0.07936507936507935' });
assert.equal(sparse.pValueNumerator, '5');
assert.equal(sparse.pValueDenominator, '63');
assert.equal(sparse.exactDecision, 'insufficient');

assert.equal(ffh.compareRationalToDecimal('1','20','0.05'), 0);
assert.equal(ffh.compareRationalToDecimal('1','20','.05'), 0);
assert.equal(ffh.compareRationalToDecimal('1','20','5.0001e-2'), -1);
assert.equal(ffh.compareRationalToDecimal('1','20','0.049999'), 1);
assert.equal(ffh.compareRationalToDecimal('1','20','invalid'), null);

const limited = ffh.calculate([[8,2,1],[1,5,4]], { maxTables: 1, alphaText: '0.05' });
assert.equal(limited.complete, false);
assert.equal(limited.pValue, undefined);
assert.equal(limited.exactDecision, undefined);

process.stdout.write('ffh threshold tests passed\n');

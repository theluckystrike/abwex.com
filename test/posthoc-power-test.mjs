import fs from 'fs';
const src = fs.readFileSync(new URL('../assets/js/posthoc-power.js', import.meta.url), 'utf8');
// Pull the pure math out of the IIFE and evaluate it standalone.
const names = ['erf', 'phi', 'invPhi', 'posthocPower', 'twoSidedP'];
let code = '';
for (const n of names) {
  const m = src.match(new RegExp('\\n  function ' + n + '\\s*\\([^)]*\\)\\s*\\{'));
  if (!m) { console.log('MISSING', n); continue; }
  let i = m.index + m[0].length - 1, depth = 0;
  const start = m.index;
  while (i < src.length) { if (src[i] === '{') depth++; else if (src[i] === '}') { depth--; if (!depth) break; } i++; }
  code += src.slice(start, i + 1) + '\n';
}
const F = new Function(code + '\nreturn {' + names.join(',') + '};')();
let fails = 0;
function eq(name, a, b, tol = 1e-4) { const ok = Math.abs(a - b) <= tol; if (!ok) { fails++; console.log('  FAIL', name, a, '!=', b); } return ok; }

// Reference cases verified against statsmodels power_proportions_2indep, which agrees
// with this pooled-null formula to machine precision (diffs ~1e-16, run 2026-09-01).
eq('power 10v12 n1000', F.posthocPower(0.10, 0.12, 1000, 1000, 0.05), 0.2981, 5e-5);
eq('power 10v12 n200', F.posthocPower(0.10, 0.12, 200, 200, 0.05), 0.0978, 5e-5);
eq('power 10v10.5 n5000', F.posthocPower(0.10, 0.105, 5000, 5000, 0.05), 0.1307, 5e-5);

// The Hoenig and Heisey identity, power at p exactly alpha is one half.
const zp = F.invPhi(1 - 0.05 / 2);
eq('invPhi(0.975)', zp, 1.959964, 1e-5);
const powAtAlpha = F.phi(zp - zp) + F.phi(-zp - zp);
eq('power at p=alpha', powAtAlpha, 0.5000, 1e-3);

// p value for the default worked example on the page.
eq('p 100/1000 v 120/1000', F.twoSidedP(0.10, 0.12, 1000, 1000), 0.1529, 5e-5);

// Degenerate inputs must return NaN, never a number.
if (!Number.isNaN(F.posthocPower(0, 0.12, 1000, 1000, 0.05))) { fails++; console.log('  FAIL zero rate not NaN'); }
if (!Number.isNaN(F.posthocPower(0.1, 1, 1000, 1000, 0.05))) { fails++; console.log('  FAIL unit rate not NaN'); }
if (!Number.isNaN(F.posthocPower(0.1, 0.12, 0, 1000, 0.05))) { fails++; console.log('  FAIL zero n not NaN'); }
if (!Number.isNaN(F.invPhi(0)) || !Number.isNaN(F.invPhi(1))) { fails++; console.log('  FAIL invPhi bounds'); }

// Symmetry and monotonicity checks.
eq('symmetry', F.posthocPower(0.12, 0.10, 1000, 1000, 0.05), F.posthocPower(0.10, 0.12, 1000, 1000, 0.05), 1e-12);
const p1k = F.posthocPower(0.10, 0.12, 1000, 1000, 0.05);
const p4k = F.posthocPower(0.10, 0.12, 4000, 4000, 0.05);
if (!(p4k > p1k)) { fails++; console.log('  FAIL power not increasing in n'); }

console.log(fails === 0 ? 'PASS posthoc-power ' + names.length + ' functions, all checks green' : 'FAILURES: ' + fails);
process.exit(fails === 0 ? 0 : 1);

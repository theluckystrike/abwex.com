/* ABWex Cochran-Armitage trend test calculator. MIT licensed, github.com/theluckystrike/abwex.com.
 * Runs entirely in the browser. Reads the reader's inputs and the two injected dataset objects
 * (#ca-size-data, #ca-fixture-data); it never sends anything to a server.
 * Math is in pure functions (asymptotic, exact, nearestCell) so it can be tested outside a browser.
 * Written to NASA Power of 10 habits. No recursion, every loop has a fixed or input bounded ceiling,
 * const and let only, no continue, and every exported function checks its arguments. Results are
 * rendered with the site's own components (.card, .input-group, .stat-grid, .result-box).
 * Loaded as an ES module: its scope is private without a wrapper function. */
const MIN_GROUPS = 3;
const MAX_GROUPS = 8;
const MAX_VISITORS_PER_GROUP = 1000000000;
// The exact test's log gamma weights lose about 1e-15 x ln(n!) in absolute terms; at a million
// visitors that is near 1e-8 relative, at a billion near 1e-5, too coarse for six printed digits.
const MAX_EXACT_VISITORS = 1000000;
const STATE_BOUND = 4000000;      // (total conversions + 1) x (score sum range + 1)
const WORK_BOUND = 200000000;     // multiply-adds allowed for the exact dynamic programme
const ALPHA = 0.05;
const MAX_STUDY_CELLS = 1000;     // ceiling on size study rows read from the page
const MAX_FIXTURES = 50;          // ceiling on worked examples read from the page

// Complementary error function, Numerical Recipes erfcc, fractional error below 1.2e-7 everywhere,
// so small tail probabilities keep their relative accuracy instead of being lost to 1 - CDF.
function erfc(x) {
  if (typeof x !== 'number') return NaN;
  if (Number.isNaN(x)) return NaN;
  const z = Math.abs(x);
  const t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
    t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 +
    t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}
function twoSidedP(z) {
  if (typeof z !== 'number') return NaN;
  if (Number.isNaN(z)) return NaN;
  return Math.min(1, erfc(Math.abs(z) / Math.SQRT2));
}

// Lanczos log gamma, g = 7, nine coefficients. Only ever called on whole counts plus one, so
// z >= 1 and the reflection formula (which would recurse) is not needed.
const LANCZOS = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
  -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
function lgamma(z) {
  const zm = z - 1;
  let x = LANCZOS[0];
  for (let i = 1; i < 9; i++) x += LANCZOS[i] / (zm + i);
  const t = zm + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (zm + 0.5) * Math.log(t) - t + Math.log(x);
}
function lnChoose(n, k) { return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1); }
function binomPmf(n, k, q) { return Math.exp(lnChoose(n, k) + k * Math.log(q) + (n - k) * Math.log(1 - q)); }

function validGroup(g) {
  if (!g || typeof g !== 'object') return false;
  if (!Number.isFinite(g.score)) return false;
  if (!Number.isInteger(g.n) || g.n < 1 || g.n > MAX_VISITORS_PER_GROUP) return false;
  return Number.isInteger(g.x) && g.x >= 0 && g.x <= g.n;
}
// Returns '' for a usable list of groups, otherwise the reason it cannot be tested.
function checkGroups(groups) {
  if (!Array.isArray(groups)) return 'Groups must be a list.';
  if (groups.length < MIN_GROUPS || groups.length > MAX_GROUPS) return 'Enter between ' + MIN_GROUPS + ' and ' + MAX_GROUPS + ' ordered groups.';
  for (let i = 0; i < groups.length; i++) {
    if (!validGroup(groups[i])) return 'Group ' + (i + 1) + ' needs a numeric score, whole visitors and whole conversions no larger than visitors.';
  }
  return '';
}

// Sums over groups with every score shifted by `shift`. The trend statistic, both variances and
// Z are unchanged by a common shift, and summing shifted scores avoids the cancellation in
// St2 - St^2 / N that made scores like 1e8, 1e8 + 1, 1e8 + 2 print a wrong p value.
function totals(groups, shift) {
  const t = { N: 0, R: 0, St: 0, St2: 0, Sx: 0 }, c = shift || 0;
  for (let i = 0; i < groups.length; i++) {
    const g = groups[i], s = g.score - c;
    t.N += g.n; t.R += g.x; t.St += s * g.n; t.St2 += s * s * g.n; t.Sx += s * g.x;
  }
  return t;
}
function minScore(groups) {
  let m = Infinity;
  for (let i = 0; i < groups.length; i++) m = Math.min(m, groups[i].score);
  return m;
}

// Asymptotic Cochran-Armitage statistic in both variance forms.
function asymptotic(groups) {
  const bad = checkGroups(groups);
  if (bad) return { defined: false, reason: bad };
  const t = totals(groups, minScore(groups));
  if (t.R === 0 || t.R === t.N) return { defined: false, reason: 'Total conversions are ' + (t.R === 0 ? 'zero' : 'equal to total visitors') + ', so the variance is zero and no trend test is defined.' };
  const pbar = t.R / t.N;
  const spread = t.St2 - t.St * t.St / t.N;
  if (!Number.isFinite(spread)) return { defined: false, reason: 'The scores are too large to compute with in the browser.' };
  if (!(spread > 0)) return { defined: false, reason: 'Every group has the same score, so there is no ordering to test.' };
  const T = t.Sx - pbar * t.St;
  const varN = pbar * (1 - pbar) * spread;
  const varC = varN * t.N / (t.N - 1);
  const zN = T / Math.sqrt(varN);
  const zC = T / Math.sqrt(varC);
  return { defined: true, N: t.N, R: t.R, pbar: pbar, T: T, varN: varN, varC: varC, zN: zN, zC: zC,
    chi2N: zN * zN, chi2C: zC * zC, pN: twoSidedP(zN), pC: twoSidedP(zC) };
}

// Size of the dynamic programme, or the reason it will not run.
function exactSetup(groups, t) {
  let minScore = Infinity, maxScore = -Infinity;
  for (let i = 0; i < groups.length; i++) {
    if (!Number.isInteger(groups[i].score)) return { reason: 'The exact test needs whole number scores.' };
    if (groups[i].n > MAX_EXACT_VISITORS) return { reason: 'The exact test here takes at most ' + MAX_EXACT_VISITORS + ' visitors per group.' };
    minScore = Math.min(minScore, groups[i].score); maxScore = Math.max(maxScore, groups[i].score);
  }
  let sMax = 0;
  for (let i = 0; i < groups.length; i++) sMax += (groups[i].score - minScore) * groups[i].n;
  sMax = Math.min(sMax, (maxScore - minScore) * t.R);
  const width = sMax + 1;
  const states = (t.R + 1) * width;
  let work = 0;
  for (let i = 0; i < groups.length; i++) work += states * (Math.min(groups[i].n, t.R) + 1);
  if (states > STATE_BOUND || work > WORK_BOUND) return { reason: 'This table is larger than the browser enumeration bound, so the exact p value was not computed. No partial sum is shown.' };
  return { reason: '', minScore: minScore, sMax: sMax, width: width, states: states };
}

// Adds one occupied state (r conversions so far, score sum ss, weight v) times the group's pmf.
// The upper limit on xx replaces a break on score sums past sMax; scores are whole and s >= 0.
function addShifted(c, r, ss, v) {
  const lim = Math.min(c.xMax, c.R - r);
  const top = c.s > 0 ? Math.min(lim, Math.floor((c.sMax - ss) / c.s)) : lim;
  for (let xx = 0; xx <= top; xx++) {
    const s2 = ss + c.s * xx;
    c.nxt[(r + xx) * c.width + s2] += v * c.pmf[xx];
  }
}
function convolveGroup(c) {
  for (let r = 0; r <= c.rHi; r++) {
    for (let ss = 0; ss <= c.sHi; ss++) {
      const v = c.cur[r * c.width + ss];
      if (v !== 0) addShifted(c, r, ss, v);
    }
  }
}

// Exact conditional test given total conversions, by dynamic programming over groups.
// Each group contributes Binomial(n_i, q) weights with q = R / N; conditional on the total R the
// distribution of the score sum does not depend on q, so this is the permutation distribution.
function runDp(groups, setup, R, q) {
  let cur = new Float64Array(setup.states);
  let nxt = new Float64Array(setup.states);
  cur[0] = 1;
  let rHi = 0, sHi = 0;
  for (let i = 0; i < groups.length; i++) {
    const n = groups[i].n, s = groups[i].score - setup.minScore, xMax = Math.min(n, R);
    const pmf = new Float64Array(xMax + 1);
    for (let x = 0; x <= xMax; x++) pmf[x] = binomPmf(n, x, q);
    nxt.fill(0);
    const rTop = Math.min(R, rHi + xMax), sTop = Math.min(setup.sMax, sHi + s * xMax);
    convolveGroup({ cur: cur, nxt: nxt, pmf: pmf, s: s, xMax: xMax, R: R, rHi: rHi, sHi: sHi, sMax: setup.sMax, width: setup.width });
    const tmp = cur; cur = nxt; nxt = tmp;
    rHi = rTop; sHi = sTop;
  }
  return cur;
}

function rowMass(cur, rowStart, width) {
  let mass = 0, support = 0;
  for (let a = 0; a < width; a++) {
    const w = cur[rowStart + a];
    if (w > 0) { mass += w; support += 1; }
  }
  return { mass: mass, support: support };
}
function tailSums(cur, rowStart, width, mean, dObs) {
  const tol = 1e-9 * Math.max(1, dObs);
  let pGe = 0, pEq = 0;
  for (let a = 0; a < width; a++) {
    const wa = cur[rowStart + a];
    const d = Math.abs(a - mean);
    if (wa !== 0 && d >= dObs - tol) pGe += wa;
    if (wa !== 0 && Math.abs(d - dObs) <= tol) pEq += wa;
  }
  return { pGe: pGe, pEq: pEq };
}
function scoreSums(groups, minScore) {
  let visitors = 0, observed = 0;
  for (let i = 0; i < groups.length; i++) {
    visitors += (groups[i].score - minScore) * groups[i].n;
    observed += (groups[i].score - minScore) * groups[i].x;
  }
  return { visitors: visitors, observed: observed };
}

function exact(groups) {
  const bad = checkGroups(groups);
  if (bad) return { computed: false, reason: bad };
  const t = totals(groups);
  if (t.R === 0 || t.R === t.N) return { computed: false, reason: 'No test is defined when total conversions are zero or equal total visitors.' };
  const setup = exactSetup(groups, t);
  if (setup.reason) return { computed: false, reason: setup.reason };
  const q = t.R / t.N;
  const cur = runDp(groups, setup, t.R, q);
  const rowStart = t.R * setup.width;
  const m = rowMass(cur, rowStart, setup.width);
  const expectedMass = binomPmf(t.N, t.R, q);
  const sums = scoreSums(groups, setup.minScore);
  const mean = t.R * sums.visitors / t.N;
  const tail = tailSums(cur, rowStart, setup.width, mean, Math.abs(sums.observed - mean));
  const p = tail.pGe / m.mass, eq = tail.pEq / m.mass;
  return { computed: true, p: Math.min(1, p), mid: Math.min(1, p - 0.5 * eq), support: m.support,
    massRatio: m.mass / expectedMass };
}

function within(v, list) { return v >= Math.min.apply(null, list) && v <= Math.max.apply(null, list); }
// True when the design sits inside the grid the exact study covers.
function onGrid(g, k, nPerGroup, rate) {
  if (!g || !Array.isArray(g.k) || !Array.isArray(g.n_per_group) || !Array.isArray(g.p)) return false;
  return g.k.indexOf(k) >= 0 && within(nPerGroup, g.n_per_group) && within(rate, g.p);
}
// Number of tables with these group sizes and this total, the count the worked examples print.
// Whole counts summed in doubles with a sliding window, exact while every partial count stays at or
// below 2^53. If any partial count passes that, a small final count could be inexact, so the result
// is NaN (not computed) unless the final is clearly above 2^53, which is then reported as Infinity.
function tableCount(groups) {
  if (checkGroups(groups)) return NaN;
  const R = totals(groups).R;
  if (!(R <= STATE_BOUND)) return NaN;
  let cur = new Float64Array(R + 1), overflow = false;
  cur[0] = 1;
  for (let i = 0; i < groups.length; i++) {
    const next = new Float64Array(R + 1), n = groups[i].n;
    let window = 0;
    for (let r = 0; r <= R; r++) {
      if (r - n - 1 >= 0) window -= cur[r - n - 1];
      window += cur[r];
      next[r] = window;
      if (window > Number.MAX_SAFE_INTEGER) overflow = true;
    }
    cur = next;
  }
  if (!overflow) return cur[R];
  return cur[R] > 2 * Number.MAX_SAFE_INTEGER ? Infinity : NaN;
}

// Nearest cell of the exact size study; null when the design sits outside the grid.
function nearestCell(study, k, nPerGroup, rate) {
  if (!study || !Array.isArray(study.cells)) return null;
  if (!Number.isInteger(k) || !(nPerGroup > 0) || !(rate > 0)) return null;
  if (!onGrid(study.grid, k, nPerGroup, rate)) return null;
  let best = null, bestD = Infinity;
  const count = Math.min(study.cells.length, MAX_STUDY_CELLS);
  for (let i = 0; i < count; i++) {
    const c = study.cells[i];
    const d = c.k === k ? Math.abs(Math.log(c.n_per_group / nPerGroup)) + Math.abs(Math.log(c.p / rate)) : Infinity;
    if (d < bestD) { bestD = d; best = c; }
  }
  return best;
}

const api = { asymptotic: asymptotic, exact: exact, nearestCell: nearestCell, tableCount: tableCount, erfc: erfc, twoSidedP: twoSidedP,
  MIN_GROUPS: MIN_GROUPS, MAX_GROUPS: MAX_GROUPS, ALPHA: ALPHA };
globalThis.ABWexTrend = api;

// ---------------- page wiring ----------------
function el(tag, attrs, text) {
  const e = document.createElement(tag);
  if (attrs) Object.keys(attrs).forEach(function (k) { e.setAttribute(k, attrs[k]); });
  if (text !== undefined) e.textContent = text;
  return e;
}
// A p value of exactly 0 is an underflow, not a result. The two sided normal tail reaches 0 in
// doubles only below about 1e-323, and an exact tail only when every tail table weight underflows,
// so "below 1e-300" is true in every case that prints it.
function fmtP(p) {
  if (p < 1e-300) return 'below 1e-300';
  if (p < 0.0001) return p.toExponential(3);
  return p.toFixed(6);
}
function fmt(v, d) { return Number(v).toFixed(d); }
function readJson(id) {
  const node = document.getElementById(id);
  if (!node) return null;
  try { return JSON.parse(node.textContent); } catch { return null; }  // malformed data renders as "did not load"
}

// Inputs are found by data-field, not by class, so no class on the page lacks a CSS rule.
const FIELDS = [
  { key: 'label', label: 'Group label', type: 'text' },
  { key: 'score', label: 'Score', type: 'text' },
  { key: 'n', label: 'Visitors', type: 'number', step: '1', min: '1' },
  { key: 'x', label: 'Conversions', type: 'number', step: '1', min: '0' }
];
let rowSeq = 0;
// A new row continues the spacing of the last two scores (0, 2, 4 then 6), else takes its index.
function nextScore(container, idx) {
  const s = Array.prototype.map.call(container.querySelectorAll('[data-field="score"]'), function (i) { return num(i.value.trim()); });
  if (s.length < 2) return idx;
  const a = s[s.length - 2], b = s[s.length - 1];
  return Number.isFinite(a) && Number.isFinite(b) ? b + (b - a) : idx;
}
function defaultValue(key, idx, values, container) {
  if (values && values[key] !== undefined) return values[key];
  if (key === 'label') return 'Group ' + (idx + 1);
  return key === 'score' ? nextScore(container, idx) : '';
}
function addRow(container, values) {
  rowSeq += 1;
  const idx = container.querySelectorAll('.ca-row').length;
  const row = el('div', { 'class': 'input-group ca-row', role: 'group', 'aria-label': 'Group ' + (idx + 1) });
  FIELDS.forEach(function (f) {
    const id = 'ca-' + f.key + '-' + rowSeq;
    const field = el('div', { 'class': 'input-field' });
    const attrs = { id: id, 'data-field': f.key, type: f.type };
    if (f.step) { attrs.step = f.step; attrs.inputmode = 'decimal'; }
    if (f.min) attrs.min = f.min;
    const input = el('input', attrs);
    input.value = String(defaultValue(f.key, idx, values, container));
    field.appendChild(el('label', { 'for': id }, f.label));
    field.appendChild(input);
    row.appendChild(field);
  });
  container.appendChild(row);
}

function fieldValue(row, key) {
  const n = row.querySelector('[data-field="' + key + '"]');
  return n ? n.value.trim() : '';
}
// An empty field is missing, not zero. Number('') is 0, which used to turn a blank conversions
// field into a real looking result.
// Plain decimal numbers only: Number() would also take '0x10' as 16 or 'Infinity'.
const DECIMAL = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;
function num(raw) { return DECIMAL.test(raw) ? Number(raw) : NaN; }
function validVisitors(n) { return Number.isInteger(n) && n >= 1 && n <= MAX_VISITORS_PER_GROUP; }
function rowErrors(g) {
  const errors = [];
  if (!Number.isFinite(g.score)) errors.push(g.label + ': the score must be a number, with a point for decimals.');
  if (!validVisitors(g.n)) errors.push(g.label + ': visitors must be a whole number from 1 to ' + MAX_VISITORS_PER_GROUP + '.');
  if (!Number.isInteger(g.x) || g.x < 0) errors.push(g.label + ': conversions must be a whole number of zero or more.');
  else if (validVisitors(g.n) && g.x > g.n) errors.push(g.label + ': conversions cannot exceed visitors.');
  return errors;
}
function readGroups(container) {
  const rows = container.querySelectorAll('.ca-row');
  const groups = [];
  let errors = [];
  const count = Math.min(rows.length, MAX_GROUPS + 1);
  for (let i = 0; i < count; i++) {
    const g = { label: fieldValue(rows[i], 'label') || 'Group ' + (i + 1), score: num(fieldValue(rows[i], 'score')),
      n: num(fieldValue(rows[i], 'n')), x: num(fieldValue(rows[i], 'x')) };
    errors = errors.concat(rowErrors(g));
    groups.push(g);
  }
  if (rows.length < MIN_GROUPS || rows.length > MAX_GROUPS) errors.push('Enter between ' + MIN_GROUPS + ' and ' + MAX_GROUPS + ' ordered groups.');
  return { groups: groups, errors: errors };
}

// Header rows. A plain string is one column; {t, cs} spans columns over a second row and
// {t, rs: 2} spans both rows. Grouping keeps shared words out of every column header.
function headerRow(cells) {
  const tr = el('tr');
  cells.forEach(function (h) {
    const attrs = { scope: h.cs ? 'colgroup' : 'col' };
    if (h.cs) attrs.colspan = String(h.cs);
    if (h.rs) attrs.rowspan = String(h.rs);
    tr.appendChild(el('th', attrs, typeof h === 'string' ? h : h.t));
  });
  return tr;
}
// A table inside the site's .table-scroll wrapper, focusable so a keyboard can scroll it.
function table(headers, rows, label, opts) {
  const o = opts || {};
  const wrap = el('div', { 'class': 'table-scroll' + (o.cls ? ' ' + o.cls : ''), role: 'region', 'aria-label': label, tabindex: '0' });
  const t = el('table'), thead = el('thead');
  if (o.groups) thead.appendChild(headerRow(o.groups));
  thead.appendChild(headerRow(headers));
  t.appendChild(thead);
  const tb = el('tbody');
  rows.forEach(function (r) {
    const row = el('tr');
    r.forEach(function (c) {
      if (c instanceof Node) { const td = el('td'); td.appendChild(c); row.appendChild(td); } else row.appendChild(el('td', null, String(c)));
    });
    tb.appendChild(row);
  });
  t.appendChild(tb);
  wrap.appendChild(t);
  return wrap;
}
// The homepage's significance colors, same thresholds as sigClass in /assets/app.js.
function sigClass(p) {
  if (p < ALPHA) return 'sig-green';
  return p < 0.1 ? 'sig-yellow' : 'sig-red';
}
function statCard(value, label, sig) {
  const card = el('div', { 'class': 'stat-card' });
  const v = el('div', { 'class': 'stat-value' });
  if (sig) v.appendChild(el('span', { 'class': sig }, value)); else v.textContent = value;
  card.appendChild(v);
  card.appendChild(el('div', { 'class': 'stat-label' }, label));
  return card;
}
function pCard(p, label) {
  return typeof p === 'number' ? statCard(fmtP(p), label, sigClass(p)) : statCard('Not computed', label);
}
// Two grids, four p values then Z, chi-square and the verdict, as the homepage lays out its results.
function statGrids(a, ex, pRef) {
  const g1 = el('div', { 'class': 'stat-grid' }), g2 = el('div', { 'class': 'stat-grid' });
  g1.appendChild(pCard(a.pN, 'p, N form'));
  g1.appendChild(pCard(a.pC, 'p, conditional form'));
  g1.appendChild(pCard(ex.computed ? ex.p : null, 'Exact conditional p'));
  g1.appendChild(pCard(ex.computed ? ex.mid : null, 'Mid p'));
  g2.appendChild(statCard(fmt(a.zN, 6), 'Z, N form'));
  g2.appendChild(statCard(fmt(a.chi2N, 6), 'Chi-square, 1 df'));
  g2.appendChild(statCard(pRef < ALPHA ? 'p < ' + ALPHA : 'p \u2265 ' + ALPHA, ex.computed ? 'Verdict, exact p' : 'Verdict, conditional p', sigClass(pRef)));
  return [g1, g2];
}
function fmtCount(c) {
  if (c === Infinity) return 'above ' + Number.MAX_SAFE_INTEGER;
  return Number.isFinite(c) ? String(c) : 'Not computed';
}
function detailBox(a, ex, groups) {
  const box = el('div', { 'class': 'result-box' });
  const dl = el('dl', { 'class': 'ca-detail' });
  // A sentence as the value would take the whole value column and squeeze every label to one word
  // per line, so sentence values get the full row under their label.
  function item(k, v, wide) { dl.appendChild(el('dt', null, k)); dl.appendChild(el('dd', wide ? { 'class': 'ca-wide' } : null, v)); }
  item('Pooled conversion rate', fmt(a.pbar, 6) + ' (' + a.R + ' of ' + a.N + ')');
  item('Trend statistic T = sum of score x (conversions - visitors x pooled rate)', fmt(a.T, 6));
  item('Variance, N form = pooled x (1 - pooled) x (sum n t^2 - (sum n t)^2 / N)', fmt(a.varN, 6));
  item('Z and chi-square (1 df), N form, matches R prop.trend.test', fmt(a.zN, 6) + ' and ' + fmt(a.chi2N, 6));
  item('Two sided asymptotic p value, N form', fmtP(a.pN));
  item('Variance, conditional form = N form x N / (N - 1)', fmt(a.varC, 6));
  item('Z, conditional form, compared to statsmodels test_ordinal_association', fmt(a.zC, 6));
  item('Two sided asymptotic p value, conditional form', fmtP(a.pC));
  if (ex.computed) {
    item('Two sided exact conditional p value', fmtP(ex.p));
    item('Two sided mid p value', fmtP(ex.mid));
    item('Tables with these group sizes and total conversions', fmtCount(tableCount(groups)));
    item('Distinct score sums enumerated, and probability mass check (should be 1)', ex.support + ' and ' + fmt(ex.massRatio, 9));
  } else {
    item('Exact conditional p value', 'Not computed. ' + ex.reason, true);
  }
  box.appendChild(dl);
  return box;
}
function groupTable(groups) {
  // Rate, the only computed column, sits next to the group so it stays in view on a phone.
  return table(['Group', 'Rate', 'Score', 'Visitors', 'Conversions'], groups.map(function (g) {
    return [g.label, fmt(g.x / g.n, 4), g.score, g.n, g.x];
  }), 'Your groups');
}
function renderResult(out, groups, a, ex) {
  out.textContent = '';
  if (!a.defined) { out.appendChild(groupTable(groups)); out.appendChild(el('p', null, a.reason)); return; }
  const pRef = ex.computed ? ex.p : a.pC;
  const which = ex.computed ? 'exact conditional' : 'asymptotic conditional form';
  statGrids(a, ex, pRef).forEach(function (g) { out.appendChild(g); });
  const verdict = 'Test run, two sided Cochran-Armitage trend test, ' + which + ' p value ' + fmtP(pRef) +
    (pRef < ALPHA ? ', below' : ', not below') + ' alpha ' + ALPHA + '. This says whether the rate moves steadily with the score, not which arm is best, and it isn\'t the probability that any variant is better.';
  out.appendChild(el('p', { 'class': 'ca-verdict' }, verdict));
  // Screen readers hear this one sentence rather than every card and table in the result.
  const status = document.getElementById('ca-status');
  if (status) status.textContent = verdict;
  out.appendChild(detailBox(a, ex, groups));
  out.appendChild(groupTable(groups));
}

function renderNearest(node, study, groups, a) {
  node.textContent = '';
  if (!a.defined) return;
  const nAvg = a.N / groups.length;
  if (!study || !Array.isArray(study.cells)) { node.appendChild(el('p', null, 'The exact study data did not load, so no nearest cell is shown.')); return; }
  const cell = nearestCell(study, groups.length, nAvg, a.pbar);
  if (!cell) { node.appendChild(el('p', null, 'Your design, ' + groups.length + ' groups, about ' + Math.round(nAvg) + ' visitors per group and pooled rate ' + fmt(a.pbar, 4) + ', is outside the grid of the exact study, so no nearest cell is shown.')); return; }
  node.appendChild(el('p', null, 'The nearest cell of the exact study is ' + cell.k + ' groups, ' + cell.n_per_group +
    ' visitors per group and true rate ' + cell.p + '. It isn\'t your exact design. At nominal two sided alpha ' + study.alpha + ', each version of the test rejects a true null with the probability below.'));
  node.appendChild(table(['Version of the test', 'Probability'], [
    ['Asymptotic N form', cell.size_asym_n], ['Asymptotic conditional form', cell.size_asym_n1],
    ['Exact conditional test', cell.size_exact], ['Mid p test', cell.size_midp],
    ['No test defined at all, because no one or everyone converted', cell.p_undefined]
  ], 'Nearest cell of the exact study'));
}

// The group count picker is the homepage's .tab-bar: one .tab-btn per k, the active one filled.
function studyTabs(bar, onPick, ks) {
  ks.forEach(function (k, i) {
    const b = el('button', { type: 'button', 'class': 'tab-btn' + (i === 0 ? ' active' : ''), 'aria-pressed': i === 0 ? 'true' : 'false', 'data-k': String(k) }, k + ' groups');
    b.addEventListener('click', function () {
      bar.querySelectorAll('.tab-btn').forEach(function (o) { o.classList.toggle('active', o === b); o.setAttribute('aria-pressed', o === b ? 'true' : 'false'); });
      onPick(k);
    });
    bar.appendChild(b);
  });
}
function renderStudy(holder, study, k) {
  if (!study || !study.grid || !Array.isArray(study.cells)) { holder.textContent = 'The exact study data did not load.'; return; }
  const rows = study.cells.slice(0, MAX_STUDY_CELLS).filter(function (c) { return c.k === k; }).map(function (c) {
    return [c.n_per_group, c.p, c.size_asym_n, c.size_asym_n1, c.size_exact, c.size_midp, c.p_undefined];
  });
  holder.textContent = '';
  const groups = [{ t: 'Visitors per group', rs: 2 }, { t: 'True rate', rs: 2 }, { t: 'Rejection rate of a true null at alpha ' + study.alpha, cs: 4 }, { t: 'No test defined', rs: 2 }];
  holder.appendChild(table(['Asymptotic N form', 'Asymptotic conditional form', 'Exact conditional', 'Mid p'], rows, 'Exact false positive rates for ' + k + ' groups', { groups: groups }));
}

function renderFixtures(holder, fixtures, load) {
  if (!fixtures || !Array.isArray(fixtures.fixtures)) { holder.textContent = 'The worked examples did not load.'; return; }
  holder.textContent = '';
  const groups = [{ t: 'Load example', rs: 2 }, { t: 'Visitors', rs: 2 }, { t: 'Conversions', rs: 2 }, { t: 'Asymptotic p', cs: 3 }, { t: 'Exact conditional', cs: 3 }];
  holder.appendChild(table(['N form', 'conditional form', 'statsmodels', 'p', 'mid p', 'tables enumerated'],
    fixtures.fixtures.slice(0, MAX_FIXTURES).map(function (f, i) {
      const name = f.id.replace(/_/g, ' ');
      const b = el('button', { type: 'button', 'class': 'btn btn-secondary btn-sm', 'data-fixture': String(i), 'aria-label': 'Load ' + name }, name);
      b.addEventListener('click', function () { load(f); });
      return [b, f.n.join(', '), f.x.join(', '), f.p_asymptotic_n, f.p_asymptotic_n1, f.statsmodels_pvalue, f.p_exact_conditional, f.p_midp, f.tables_enumerated];
    }), 'Worked examples', { cls: 'ca-fixtures', groups: groups }));
  if (fixtures.check) holder.appendChild(el('p', { 'class': 'ca-caption' }, 'Cross check against ' + fixtures.check.reference + ', maximum absolute difference ' + fixtures.check.max_abs_diff + '. ' + (fixtures.check.note || '')));
}

// Fill empty visitor and conversion fields from the first published worked example, the way the
// homepage opens with example counts. Values come from the dataset, never typed into the page.
// The first worked example, or null when the record is missing or malformed.
function firstFixture(fixtures) {
  const f = fixtures && Array.isArray(fixtures.fixtures) ? fixtures.fixtures[0] : null;
  if (!f || !Array.isArray(f.n) || !Array.isArray(f.x)) return null;
  return f.n.length === f.x.length ? f : null;
}
function prefill(rows, fixtures) {
  const f = firstFixture(fixtures);
  const all = rows.querySelectorAll('.ca-row');
  if (!f || all.length !== f.n.length) return false;
  let filled = 0;
  for (let j = 0; j < all.length; j++) {
    const n = all[j].querySelector('[data-field="n"]'), x = all[j].querySelector('[data-field="x"]');
    if (n && x && n.value === '' && x.value === '') { n.value = String(f.n[j]); x.value = String(f.x[j]); filled += 1; }
  }
  return filled > 0;
}
// At the group limits Add or Remove is marked aria-disabled, so a click never silently does
// nothing. aria-disabled rather than disabled keeps keyboard focus on the button; a disabled
// button drops focus, and moving it elsewhere let a repeated Enter remove a row.
function syncLimits(rows) {
  const count = rows.querySelectorAll('.ca-row').length;
  const add = document.getElementById('ca-add'), remove = document.getElementById('ca-remove');
  if (!add || !remove) return count;
  add.setAttribute('aria-disabled', count >= MAX_GROUPS ? 'true' : 'false');
  remove.setAttribute('aria-disabled', count <= MIN_GROUPS ? 'true' : 'false');
  return count;
}
// After a calculation the reader asked for, bring the answer (or the error) into view when it is
// off screen. At 1440 by 900 the first result card started at y 845 of 900, and after a worked
// example's Load button it sat 3669px above the viewport.
const REVEAL_BOTTOM_GAP = 120;
// node is always an element: run() returns the error line or the first rendered result.
function reveal(node) {
  const top = node.getBoundingClientRect().top;
  if (top >= 0 && top <= window.innerHeight - REVEAL_BOTTOM_GAP) return false;
  const reduce = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  node.scrollIntoView({ behavior: reduce ? 'auto' : 'smooth', block: 'start' });
  return true;
}
function markStale(out) {
  if (!out.firstElementChild || out.querySelector('.ca-stale')) return false;
  out.insertBefore(el('p', { 'class': 'ca-note ca-stale' }, 'The inputs changed after this result. Press Calculate to update it.'), out.firstElementChild);
  return true;
}
function wireCalculator(ctx) {
  function run() {
    const r = readGroups(ctx.rows);
    ctx.err.textContent = r.errors.join(' ');
    const status = document.getElementById('ca-status');
    if (status) status.textContent = '';
    if (r.errors.length) { ctx.out.textContent = ''; if (ctx.nearest) ctx.nearest.textContent = ''; return ctx.err; }
    const a = asymptotic(r.groups);
    const ex = a.defined ? exact(r.groups) : { computed: false, reason: a.reason };
    renderResult(ctx.out, r.groups, a, ex);
    if (ctx.nearest) renderNearest(ctx.nearest, ctx.study, r.groups, a);
    return ctx.out.firstElementChild;
  }
  function load(f, moveFocus) {
    ctx.rows.textContent = '';
    for (let j = 0; j < f.n.length && j < MAX_GROUPS; j++) addRow(ctx.rows, { label: 'Group ' + (j + 1), score: f.scores ? f.scores[j] : j, n: f.n[j], x: f.x[j] });
    syncLimits(ctx.rows);
    reveal(run());
    // A worked example button sits far below the calculator; focus follows the answer up.
    if (moveFocus) ctx.out.focus({ preventScroll: true });
  }
  ctx.form.addEventListener('submit', function (e) { e.preventDefault(); reveal(run()); });
  // Editing a field after a calculation leaves an old answer on screen; say so above it.
  ctx.rows.addEventListener('input', function () { markStale(ctx.out); });
  // Chrome steps a focused number field on each wheel notch, so scrolling past it silently
  // changed a count (200 became 202 after two notches). Blurring it lets the wheel scroll the page.
  ctx.rows.addEventListener('wheel', function (e) {
    if (e.target && e.target.type === 'number' && document.activeElement === e.target) e.target.blur();
  }, { passive: true });
  const add = document.getElementById('ca-add'), remove = document.getElementById('ca-remove'), example = document.getElementById('ca-example');
  // Row changes clear the old error, which may name a group that no longer exists.
  function rowsChanged() { ctx.err.textContent = ''; markStale(ctx.out); syncLimits(ctx.rows); }
  if (add) add.addEventListener('click', function () { if (ctx.rows.querySelectorAll('.ca-row').length < MAX_GROUPS) { addRow(ctx.rows); rowsChanged(); } });
  if (remove) remove.addEventListener('click', function () { const all = ctx.rows.querySelectorAll('.ca-row'); if (all.length > MIN_GROUPS) { ctx.rows.removeChild(all[all.length - 1]); rowsChanged(); } });
  const hasFixtures = ctx.fixtures && Array.isArray(ctx.fixtures.fixtures) && ctx.fixtures.fixtures.length > 0;
  if (example && hasFixtures) example.addEventListener('click', function () { load(ctx.fixtures.fixtures[0], false); });
  return function (f) { load(f, true); };
}

function wireStudy(study) {
  const bar = document.getElementById('ca-ng'), holder = document.getElementById('ca-size-table');
  if (!bar || !holder) return false;
  const ks = study && study.grid && Array.isArray(study.grid.k) ? study.grid.k.slice(0, MAX_GROUPS) : [];
  studyTabs(bar, function (k) { renderStudy(holder, study, k); }, ks);
  renderStudy(holder, study, ks[0]);
  return true;
}
function init() {
  const form = document.getElementById('ca-form'), rows = document.getElementById('ca-rows');
  if (!form || !rows) return;
  const ctx = { form: form, rows: rows, out: document.getElementById('ca-out'), err: document.getElementById('ca-error'),
    nearest: document.getElementById('ca-nearest'), study: readJson('ca-size-data'), fixtures: readJson('ca-fixture-data') };
  if (!ctx.out || !ctx.err) return;
  // The page ships its first rows as static HTML so nothing shifts on load; new ids continue after them.
  rowSeq = rows.querySelectorAll('.ca-row').length;
  if (rowSeq === 0) for (let i = 0; i < MIN_GROUPS; i++) addRow(rows);
  const load = wireCalculator(ctx);
  // The controls ship disabled so that without the script, or before it arrives on a slow line,
  // Calculate cannot submit the form natively and reload the page over the reader's numbers.
  // They are enabled before anything that reads the datasets, so a bad record cannot leave them off.
  ['ca-add', 'ca-remove', 'ca-example', 'ca-calc'].forEach(function (id) { const b = document.getElementById(id); if (b) b.disabled = false; });
  syncLimits(rows);
  prefill(rows, ctx.fixtures);
  wireStudy(ctx.study);
  const fx = document.getElementById('ca-fixture-table');
  if (fx) renderFixtures(fx, ctx.fixtures, load);
}
// Module scripts run after parsing; outside a browser (the node tests) there is no page to wire.
if (typeof document !== 'undefined') {
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
}

export { asymptotic, exact, nearestCell, tableCount, erfc, twoSidedP, MIN_GROUPS, MAX_GROUPS, ALPHA };

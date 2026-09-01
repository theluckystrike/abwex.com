/* Post hoc power calculator.
   Runs entirely in the browser. No input leaves the page.

   Computes the achieved power of a finished two proportion test at its own observed
   effect size, the quantity usually asked for as post hoc or observed power, and then
   says plainly what that number can and cannot support, because observed power is a
   deterministic transform of the p value and reviewers who ask for it are usually
   asking for something the confidence interval answers better.

   Bounds and validation follow the same rules as the rest of this site: every loop is
   bounded, every function validates its arguments, no function exceeds the size limit.
   Output is built with DOM nodes and textContent, never markup strings, so no user
   input can ever reach the page as HTML. */
(function () {
  'use strict';

  const Z975 = 1.959963984540054; /* Phi^-1(0.975), fixed reference for the 95% CI */

  function el(id) { return document.getElementById(id); }
  function fmt(x, d) { return isFinite(x) ? x.toFixed(d == null ? 4 : d) : 'n/a'; }

  /* Standard normal CDF via erf. Abramowitz and Stegun 7.1.26 rational
     approximation, absolute error below 1.5e-7, plenty for four decimals. */
  function erf(x) {
    if (!isFinite(x)) { return x > 0 ? 1 : -1; }
    const sign = x < 0 ? -1 : 1;
    const ax = Math.abs(x);
    const t = 1 / (1 + 0.3275911 * ax);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-ax * ax);
    return sign * y;
  }

  function phi(x) {
    if (!isFinite(x)) { return x > 0 ? 1 : 0; }
    return 0.5 * (1 + erf(x / Math.SQRT2));
  }

  /* Inverse standard normal CDF, Acklam's rational approximation,
     absolute error below 1.15e-9 across (0, 1). */
  function invPhi(q) {
    if (!isFinite(q) || q <= 0 || q >= 1) { return NaN; }
    const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687,
      138.3577518672690, -30.66479806614716, 2.506628277459239];
    const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866,
      66.80131188771972, -13.28068155288572];
    const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838,
      -2.549732539343734, 4.374664141464968, 2.938163982698783];
    const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996,
      3.754408661907416];
    const plow = 0.02425;
    if (q < plow) {
      const u = Math.sqrt(-2 * Math.log(q));
      return (((((c[0] * u + c[1]) * u + c[2]) * u + c[3]) * u + c[4]) * u + c[5]) /
        ((((d[0] * u + d[1]) * u + d[2]) * u + d[3]) * u + 1);
    }
    if (q > 1 - plow) {
      const u = Math.sqrt(-2 * Math.log(1 - q));
      return -(((((c[0] * u + c[1]) * u + c[2]) * u + c[3]) * u + c[4]) * u + c[5]) /
        ((((d[0] * u + d[1]) * u + d[2]) * u + d[3]) * u + 1);
    }
    const u = q - 0.5;
    const t = u * u;
    return (((((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4]) * t + a[5]) * u /
      (((((b[0] * t + b[1]) * t + b[2]) * t + b[3]) * t + b[4]) * t + 1);
  }

  /* Post hoc power of the two sided two proportion z test at the observed rates.
     Pooled variance under the null for the critical value, unpooled under the
     alternative, which is exactly what statsmodels power_proportions_2indep
     computes. The second phi term is rejection on the wrong side, negligible
     except at tiny effects but kept so the formula is exact under the
     normal approximation. */
  function posthocPower(p1, p2, n1, n2, alpha) {
    if (!(n1 > 0) || !(n2 > 0)) { return NaN; }
    if (!(p1 > 0) || !(p1 < 1) || !(p2 > 0) || !(p2 < 1)) { return NaN; }
    if (!(alpha > 0) || !(alpha < 1)) { return NaN; }
    const za = invPhi(1 - alpha / 2);
    const diff = Math.abs(p1 - p2);
    const pbar = (n1 * p1 + n2 * p2) / (n1 + n2);
    const se0 = Math.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2));
    const se1 = Math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2);
    if (!(se1 > 0)) { return NaN; }
    return phi((diff - za * se0) / se1) + phi((-diff - za * se0) / se1);
  }

  function twoSidedP(p1, p2, n1, n2) {
    if (!(n1 > 0) || !(n2 > 0)) { return NaN; }
    if (!(p1 >= 0) || !(p2 >= 0)) { return NaN; }
    const pbar = (n1 * p1 + n2 * p2) / (n1 + n2);
    const se0 = Math.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2));
    if (!(se0 > 0)) { return NaN; }
    return 2 * (1 - phi(Math.abs(p2 - p1) / se0));
  }

  function readInt(id, min, max) {
    const n = el(id);
    if (!n) { return NaN; }
    const v = parseInt(n.value, 10);
    if (!isFinite(v) || v < min || v > max) { return NaN; }
    return v;
  }

  function readFloat(id, min, max) {
    const n = el(id);
    if (!n) { return NaN; }
    const v = parseFloat(n.value);
    if (!isFinite(v) || v < min || v > max) { return NaN; }
    return v;
  }

  /* DOM builders. Everything user visible goes through textContent, so no
     input, however malformed, can reach the page as markup. */
  function node(tag, cls, text) {
    if (typeof tag !== 'string' || tag.length === 0) { return null; }
    const e = document.createElement(tag);
    if (cls) { e.className = cls; }
    if (text != null) { e.textContent = String(text); }
    return e;
  }

  function statCard(value, cls, label) {
    const card = node('div', 'stat-card', null);
    card.appendChild(node('div', cls ? 'stat-value ' + cls : 'stat-value', value));
    card.appendChild(node('div', 'stat-label', label));
    return card;
  }

  function mutMessage(box, text) {
    if (!box || typeof text !== 'string') { return; }
    box.replaceChildren(node('p', 'mut', text));
  }

  function powerReadingText(power, p, alpha) {
    if (!isFinite(power)) { return ''; }
    let s = 'Observed power ' + fmt(power, 4) + ' is another way of writing p = '
      + fmt(p, 4) + '. The two numbers are locked together by the formula on this '
      + 'page, so the power figure adds nothing the p value did not already say.';
    if (p >= alpha) {
      s += ' Low observed power here does not show the test was underpowered, it '
        + 'restates that the result was not significant. To judge the design, read the '
        + 'confidence interval against the smallest lift that would matter to you.';
    } else {
      s += ' The result is significant at your alpha, and a significant result always '
        + 'reports observed power above one half, whatever the design deserved.';
    }
    return s;
  }

  function armsTable(o) {
    if (!o || !(o.n1 > 0) || !(o.n2 > 0)) { return null; }
    const wrap = node('div', 'tbl-wrap', null);
    const table = node('table', 'data', null);
    const thead = node('thead', null, null);
    const hrow = node('tr', null, null);
    ['Arm', 'Visitors', 'Conversions', 'Rate'].forEach(function (h) {
      hrow.appendChild(node('th', null, h));
    });
    thead.appendChild(hrow);
    const tbody = node('tbody', null, null);
    [['A', o.n1, o.c1, fmt(o.p1 * 100, 2) + '%'],
     ['B', o.n2, o.c2, fmt(o.p2 * 100, 2) + '%']].forEach(function (cells) {
      const row = node('tr', null, null);
      cells.forEach(function (c) { row.appendChild(node('td', null, c)); });
      tbody.appendChild(row);
    });
    table.appendChild(thead);
    table.appendChild(tbody);
    wrap.appendChild(table);
    return wrap;
  }

  function renderResult(box, o) {
    if (!box || !o) { return; }
    const sigCls = o.p < o.alpha ? 'sig-green' : 'sig-red';
    const grid = node('div', 'stat-grid', null);
    grid.appendChild(statCard(fmt(o.power, 4), sigCls, 'observed power at the observed effect'));
    grid.appendChild(statCard(fmt(o.p, 4), '', 'two sided p value'));
    grid.appendChild(statCard(fmt(o.z, 4), '', 'z statistic, pooled'));
    grid.appendChild(statCard(fmt(o.lo * 100, 2) + ' to ' + fmt(o.hi * 100, 2) + ' pp',
      '', '95% CI for the difference'));
    box.replaceChildren(grid, armsTable(o), node('p', null, powerReadingText(o.power, o.p, o.alpha)));
  }

  /* Returns {error} or the validated inputs. Kept separate from compute so each
     function stays simple enough to reason about at a glance. */
  function readTestInputs() {
    const n1 = readInt('na', 1, 1e9);
    const n2 = readInt('nb', 1, 1e9);
    const c1 = readInt('ca', 0, 1e9);
    const c2 = readInt('cb', 0, 1e9);
    const alpha = readFloat('alpha', 0.0001, 0.5);
    if (![n1, n2, c1, c2].every(isFinite)) {
      return { error: 'Each arm needs a visitor count and a conversion count, whole '
        + 'numbers, visitors at least 1.' };
    }
    if (!isFinite(alpha)) {
      return { error: 'Alpha must be a number between 0.0001 and 0.5.' };
    }
    if (c1 > n1 || c2 > n2) {
      return { error: 'Conversions cannot exceed visitors in an arm.' };
    }
    const p1 = c1 / n1;
    const p2 = c2 / n2;
    const inOpenUnit = function (r) { return r > 0 && r < 1; };
    if (!inOpenUnit(p1) || !inOpenUnit(p2)) {
      return { error: 'A rate of exactly 0% or 100% breaks the normal approximation '
        + 'this test relies on. Post hoc power is undefined there, and so is the z '
        + 'test itself. An exact method such as the Fisher exact test is the right '
        + 'tool for those data.' };
    }
    return { n1: n1, n2: n2, c1: c1, c2: c2, p1: p1, p2: p2, alpha: alpha, error: null };
  }

  function compute() {
    const box = el('out');
    if (!box) { return; }
    const inp = readTestInputs();
    if (inp.error) {
      mutMessage(box, inp.error);
      return;
    }
    const n1 = inp.n1;
    const n2 = inp.n2;
    const c1 = inp.c1;
    const c2 = inp.c2;
    const p1 = inp.p1;
    const p2 = inp.p2;
    const alpha = inp.alpha;
    const diff = p2 - p1;
    const pbar = (c1 + c2) / (n1 + n2);
    const se0 = Math.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2));
    const se1 = Math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2);
    renderResult(box, {
      n1: n1, n2: n2, c1: c1, c2: c2, p1: p1, p2: p2, alpha: alpha,
      z: diff / se0,
      p: twoSidedP(p1, p2, n1, n2),
      power: posthocPower(p1, p2, n1, n2, alpha),
      lo: diff - Z975 * se1,
      hi: diff + Z975 * se1
    });
  }

  /* The p value to observed power converter, the Hoenig and Heisey identity
     made into a control. power = phi(z_p - z_a2) + phi(-z_p - z_a2). */
  function convert() {
    const box = el('pout');
    if (!box) { return; }
    const p = readFloat('pin', 0.00001, 0.99999);
    const alpha = readFloat('palpha', 0.0001, 0.5);
    if (!isFinite(p) || !isFinite(alpha)) {
      mutMessage(box, 'The p value must sit between 0.00001 and 0.99999, alpha '
        + 'between 0.0001 and 0.5.');
      return;
    }
    const zp = invPhi(1 - p / 2);
    const za = invPhi(1 - alpha / 2);
    const power = phi(zp - za) + phi(-zp - za);
    box.replaceChildren(node('p', null, 'A two sided p value of ' + fmt(p, 4)
      + ' at alpha ' + fmt(alpha, 3) + ' always converts to an observed power of '
      + fmt(power, 4) + ', whatever the sample size, the metric or the platform. '
      + 'No data beyond the p value went into that number, which is the whole '
      + 'problem with it.'));
  }

  ['na', 'nb', 'ca', 'cb', 'alpha'].forEach(function (id) {
    const n = el(id);
    if (n) { n.addEventListener('input', compute); n.addEventListener('change', compute); }
  });
  ['pin', 'palpha'].forEach(function (id) {
    const n = el(id);
    if (n) { n.addEventListener('input', convert); n.addEventListener('change', convert); }
  });

  Array.prototype.forEach.call(document.querySelectorAll('.faq-q'), function (q) {
    q.addEventListener('click', function () {
      const item = q.parentElement;
      if (item) { item.classList.toggle('open'); }
    });
  });

  compute();
  convert();
}());

/* False discovery rate calculator.
   Bonferroni, Holm, Benjamini Hochberg and Benjamini Yekutieli adjusted p values,
   plus a browser replication of the published Monte Carlo cells.

   Written to NASA Power of 10 habits. Every loop has a fixed upper bound, there is
   no recursion, every allocation happens up front, and every return value that can
   fail is checked before use. */
(function () {
  'use strict';

  var MAXM = 200;              /* hard cap on metrics per run */
  var INV_NORM_STEPS = 120;    /* fixed bisection budget */
  var MAX_LINES = 5000;        /* hard cap on parsed input lines */
  var MIN_REPLICATE = 200;
  var MAX_REPLICATE = 5000;
  var NS = 'http://www.w3.org/2000/svg';
  var LABEL = { none: 'No correction', bonferroni: 'Bonferroni', holm: 'Holm', bh: 'Benjamini Hochberg' };

  /* Chart colours read from the stylesheet so the tool tracks the site theme. */
  var CSS = getComputedStyle(document.documentElement);
  function token(name, fallback) {
    var v = CSS.getPropertyValue(name);
    return (v && v.trim()) ? v.trim() : fallback;
  }
  var C_AXIS = token('--border', '#2a2a2a');
  var C_DIM = token('--text-dim', '#888');
  var C_ACCENT = token('--accent', '#2DD4BF');
  var C_GREEN = token('--green', '#22C55E');
  var C_RED = token('--red', '#EF4444');

  var SIM = null;
  var node = document.getElementById('sim-data');
  if (node) {
    try { SIM = JSON.parse(node.textContent); } catch (e) { SIM = null; }
  }
  if (SIM && (!SIM.cells || !SIM.cells.length)) { SIM = null; }

  function el(id) { return document.getElementById(id); }
  function esc(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }
  function pct(x, d) { return (100 * x).toFixed(d == null ? 1 : d) + '%'; }
  function fmtP(x) {
    if (!isFinite(x)) { return 'n/a'; }
    if (x > 0 && x < 0.0001) { return x.toExponential(2); }
    return x.toFixed(4);
  }
  function lab(i) { return i < 26 ? 'Metric ' + String.fromCharCode(65 + i) : 'Metric ' + (i + 1); }

  function normCdf(z) {
    var t = 1 / (1 + 0.2316419 * Math.abs(z));
    var d = 0.3989422804014327 * Math.exp(-z * z / 2);
    var p = d * t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    return z > 0 ? 1 - p : p;
  }
  function invNorm(q) {
    var lo = -12, hi = 12, mid = 0, i;
    for (i = 0; i < INV_NORM_STEPS; i += 1) {
      mid = (lo + hi) / 2;
      if (normCdf(mid) < q) { lo = mid; } else { hi = mid; }
    }
    return (lo + hi) / 2;
  }

  function orderAsc(ps) {
    var idx = ps.map(function (v, i) { return i; });
    idx.sort(function (a, b) { return ps[a] - ps[b]; });
    return idx;
  }
  function harmonic(m) {
    var s = 0, i;
    for (i = 1; i <= m; i += 1) { s += 1 / i; }
    return s;
  }
  function bonf(ps) {
    var m = ps.length;
    return ps.map(function (v) { return Math.min(1, v * m); });
  }
  function holm(ps) {
    var m = ps.length, idx = orderAsc(ps), adj = new Array(m), prev = 0, k, i, v;
    for (k = 0; k < m; k += 1) {
      i = idx[k];
      v = Math.min(1, Math.max(prev, (m - k) * ps[i]));
      adj[i] = v; prev = v;
    }
    return adj;
  }
  /* Step-up with an optional inflation factor. factor 1 gives Benjamini Hochberg,
     factor sum(1/i) gives Benjamini Yekutieli, matching R's p.adjust. */
  function stepUp(ps, factor) {
    var m = ps.length, idx = orderAsc(ps), adj = new Array(m), prev = 1, k, i, v;
    for (k = m - 1; k >= 0; k -= 1) {
      i = idx[k];
      v = Math.min(prev, ps[i] * m * factor / (k + 1));
      adj[i] = Math.min(1, v); prev = v;
    }
    return adj;
  }
  function bh(ps) { return stepUp(ps, 1); }
  function by(ps) { return stepUp(ps, harmonic(ps.length)); }

  function adjustAll(ps) {
    return { none: ps.slice(), bonferroni: bonf(ps), holm: holm(ps), bh: bh(ps), by: by(ps) };
  }

  function parseInput(raw) {
    var lines = String(raw || '').split('\n');
    var rows = [], dropped = 0, k = 0, i, j;
    var limit = Math.min(lines.length, MAX_LINES);
    for (i = 0; i < limit; i += 1) {
      var line = lines[i].trim();
      if (!line) { continue; }
      var nums = line.match(/-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?/g);
      if (!nums) { dropped += 1; continue; }
      var last = nums[nums.length - 1];
      var head = line.slice(0, line.lastIndexOf(last));
      var named = /[A-Za-z]/.test(head.replace(/\d[eE][+-]?\d/g, ' '));
      if (named) {
        k += 1;
        var pv = parseFloat(last);
        var nm = head.replace(/[,;=|]+ *$/, '').trim();
        if (isFinite(pv) && pv >= 0 && pv <= 1) { rows.push({ name: nm || lab(k - 1), p: pv }); } else { dropped += 1; }
      } else {
        for (j = 0; j < nums.length; j += 1) {
          var v2 = parseFloat(nums[j]);
          k += 1;
          if (isFinite(v2) && v2 >= 0 && v2 <= 1) { rows.push({ name: lab(k - 1), p: v2 }); } else { dropped += 1; }
        }
      }
    }
    var truncated = 0;
    if (rows.length > MAXM) { truncated = rows.length - MAXM; rows = rows.slice(0, MAXM); }
    return { rows: rows, dropped: dropped, truncated: truncated };
  }

  /* One adjusted value as a single cell. Colour carries the survives verdict so the
     table needs one column per correction instead of two. The title and the screen
     reader text keep the verdict available without colour. */
  function cell(v, alpha) {
    var lives = v <= alpha;
    var word = lives ? 'survives' : 'does not survive';
    return '<td class="' + (lives ? 'liv' : 'die') + '" title="' + fmtP(v) + ', ' + word + ' at alpha ' + alpha + '">'
      + fmtP(v) + '<span class="sr">, ' + word + '</span></td>';
  }

  function stepChart(sorted, alpha, m) {
    var W = 640, H = 220, L = 52, R = 16, T = 16, B = 34;
    var maxP = Math.max.apply(null, sorted.concat([alpha]));
    var yMax = Math.min(1, Math.max(maxP, alpha) * 1.15);
    if (!(yMax > 0)) { yMax = alpha > 0 ? alpha : 1; }
    var px = function (i) { return m === 1 ? L + (W - L - R) / 2 : L + i * (W - L - R) / (m - 1); };
    var py = function (v) { return T + (1 - Math.min(v, yMax) / yMax) * (H - T - B); };
    var svg = document.createElementNS(NS, 'svg');
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', 'Raw p values by rank against the Benjamini Hochberg line');
    var add = function (tag, attrs, text) {
      var n = document.createElementNS(NS, tag);
      Object.keys(attrs).forEach(function (a) { n.setAttribute(a, attrs[a]); });
      if (text != null) { n.textContent = text; }
      svg.appendChild(n);
    };
    add('line', { x1: L, y1: py(0), x2: W - R, y2: py(0), stroke: C_AXIS });
    add('line', { x1: L, y1: T, x2: L, y2: py(0), stroke: C_AXIS });
    add('text', { x: 6, y: T + 10, 'font-size': 11, fill: C_DIM }, yMax.toFixed(3));
    add('text', { x: 6, y: py(0), 'font-size': 11, fill: C_DIM }, '0');
    add('line', {
      x1: px(0), y1: py(alpha / m), x2: px(m - 1), y2: py(alpha),
      stroke: C_ACCENT, 'stroke-width': 1.5, 'stroke-dasharray': '5 4'
    });
    var i;
    for (i = 0; i < m; i += 1) {
      add('circle', { cx: px(i), cy: py(sorted[i]), r: 4, fill: sorted[i] <= alpha * (i + 1) / m ? C_GREEN : C_RED });
    }
    add('text', { x: L, y: H - 8, 'font-size': 11, fill: C_DIM }, 'rank 1');
    add('text', { x: W - R - 56, y: H - 8, 'font-size': 11, fill: C_DIM }, 'rank ' + m);
    return svg;
  }

  /* How many metrics each correction keeps at this alpha. */
  function survivorCounts(adj, m, alpha) {
    var counts = { bh: 0, by: 0, holm: 0, bonferroni: 0, none: 0 };
    Object.keys(counts).forEach(function (k) {
      var i;
      for (i = 0; i < m; i += 1) { if (adj[k][i] <= alpha) { counts[k] += 1; } }
    });
    return counts;
  }

  /* Largest rank still under the Benjamini Hochberg line, zero if none is. */
  function bhCutoff(ps, idx, m, alpha) {
    var r;
    for (r = m; r >= 1; r -= 1) {
      if (ps[idx[r - 1]] <= alpha * r / m) { return r; }
    }
    return 0;
  }

  function summaryCard(parsed, counts, m, alpha) {
    var h = '<div class="card"><h3 style="margin-top:0;">What survives</h3>';
    h += '<p>' + m + ' metric' + (m === 1 ? '' : 's') + ' in, alpha ' + alpha + ', p values used as supplied. '
      + 'Benjamini Hochberg keeps ' + counts.bh + ', Benjamini Yekutieli keeps ' + counts.by
      + ', Holm keeps ' + counts.holm + ', Bonferroni keeps ' + counts.bonferroni
      + '. With no correction at all you would have kept ' + counts.none + '.</p>';
    h += '<p>Benjamini Hochberg holds the false discovery rate at ' + alpha + ', the alpha you set, under independence or '
      + 'positive regression dependence. Read it as a level, not a count. On average at most that proportion of what it keeps '
      + 'is a false discovery, which scaled across the ' + counts.bh + ' kept here is ' + (alpha * counts.bh).toFixed(3) + '.</p>';
    if (parsed.dropped) {
      h += '<p class="mut">' + parsed.dropped + ' value' + (parsed.dropped === 1 ? '' : 's')
        + ' ignored, a p value has to sit between zero and one.</p>';
    }
    if (parsed.truncated) {
      h += '<p class="mut">' + parsed.truncated + ' value' + (parsed.truncated === 1 ? '' : 's')
        + ' past the limit of ' + MAXM + ' were not used.</p>';
    }
    return h;
  }

  function adjustedTable(rows, ps, adj, idx, m, alpha) {
    var h = '<div class="tbl-wrap"><table class="data adj"><caption>Adjusted p values, smallest raw p first. '
      + 'A value in green survives at alpha ' + alpha + '.</caption>'
      + '<thead><tr><th>Metric</th><th>Raw p</th><th>BH</th><th>BY</th><th>Holm</th><th>Bonferroni</th></tr></thead><tbody>';
    var a, i;
    for (a = 0; a < m; a += 1) {
      i = idx[a];
      h += '<tr><td>' + esc(rows[i].name) + '</td><td class="rawp">' + fmtP(ps[i]) + '</td>'
        + cell(adj.bh[i], alpha) + cell(adj.by[i], alpha)
        + cell(adj.holm[i], alpha) + cell(adj.bonferroni[i], alpha) + '</tr>';
    }
    return h + '</tbody></table></div></div>';
  }

  function stepUpCard(rows, ps, adj, idx, m, alpha, cut) {
    var h = '<div class="card"><h3 style="margin-top:0;">The step up, shown</h3>';
    h += '<p>' + (cut
      ? 'The largest rank still under alpha times rank over m is rank ' + cut + ', so Benjamini Hochberg keeps every metric '
        + 'at or below it, including any whose own raw p sits above its line.'
      : 'No rank came in under alpha times rank over m, so Benjamini Hochberg keeps nothing at this alpha.') + '</p>';
    h += '<div class="tbl-wrap"><table class="data"><thead><tr><th>Rank</th><th>Metric</th><th>Raw p</th>'
      + '<th>BH line, alpha x rank / m</th><th>Under its line</th><th>BH adjusted</th></tr></thead><tbody>';
    var b, i, rank, line;
    for (b = 0; b < m; b += 1) {
      i = idx[b]; rank = b + 1; line = alpha * rank / m;
      h += '<tr><td>' + rank + '</td><td>' + esc(rows[i].name) + '</td><td>' + fmtP(ps[i]) + '</td>'
        + '<td>' + fmtP(line) + '</td>'
        + '<td>' + (ps[i] <= line ? '<span class="go">yes</span>' : '<span class="stop">no</span>') + '</td>'
        + '<td>' + fmtP(adj.bh[i]) + '</td></tr>';
    }
    return h + '</tbody></table></div><div class="chart-container" id="chartslot"></div>'
      + '<p class="chart-note">Dots are your raw p values by rank. The dashed line is alpha times rank over m, the '
      + 'Benjamini Hochberg threshold. A green dot sits on or under its own line.</p></div>';
  }

  function render(parsed, alpha) {
    var out = el('out');
    if (!out) { return; }
    var rows = parsed.rows, m = rows.length;
    if (!m) {
      out.innerHTML = '<p class="mut">Nothing to adjust yet. Paste p values above, or draw an example.</p>';
      return;
    }
    var ps = rows.map(function (r) { return r.p; });
    var adj = adjustAll(ps), idx = orderAsc(ps);
    var counts = survivorCounts(adj, m, alpha);
    var cut = bhCutoff(ps, idx, m, alpha);
    out.innerHTML = summaryCard(parsed, counts, m, alpha)
      + adjustedTable(rows, ps, adj, idx, m, alpha)
      + stepUpCard(rows, ps, adj, idx, m, alpha, cut);
    var slot = el('chartslot');
    if (slot) { slot.appendChild(stepChart(idx.map(function (i) { return ps[i]; }), alpha, m)); }
  }

  function currentAlpha() {
    var box = el('alpha');
    var a = box ? parseFloat(box.value) : NaN;
    if (!isFinite(a) || a <= 0 || a >= 1) {
      a = (SIM && isFinite(SIM.alpha)) ? SIM.alpha : 0.05;
      if (box) { box.value = a; }
    }
    return a;
  }
  function compute() {
    var box = el('pvals');
    render(parseInput(box ? box.value : ''), currentAlpha());
  }

  function mulberry32(a) {
    return function () {
      a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function makeNormal(rnd) {
    var spare = null;
    return function () {
      if (spare !== null) { var s = spare; spare = null; return s; }
      var u = Math.max(rnd(), 1e-12), v = rnd(), r = Math.sqrt(-2 * Math.log(u));
      spare = r * Math.sin(2 * Math.PI * v);
      return r * Math.cos(2 * Math.PI * v);
    };
  }
  function drawExperiment(nrm, m, rho, m1, mu) {
    var c = nrm(), ps = [], truth = [], i, z;
    for (i = 0; i < m; i += 1) {
      z = Math.sqrt(rho) * c + Math.sqrt(1 - rho) * nrm();
      if (i < m1) { z += mu; }
      ps.push(2 * (1 - normCdf(Math.abs(z))));
      truth.push(i < m1);
    }
    return { p: ps, truth: truth };
  }

  function uniqSorted(vals) {
    var out = [];
    vals.forEach(function (v) { if (out.indexOf(v) < 0) { out.push(v); } });
    return out.sort(function (a, b) { return a - b; });
  }
  function fillSelect(sel, vals, fmt) {
    sel.innerHTML = '';
    vals.forEach(function (v) {
      var o = document.createElement('option');
      o.value = v; o.textContent = fmt(v);
      sel.appendChild(o);
    });
  }
  function currentCell() {
    if (!SIM) { return null; }
    var m = parseFloat(el('selM').value), rho = parseFloat(el('selRho').value), pi0 = parseFloat(el('selPi').value);
    var hit = SIM.cells.filter(function (c) { return c.m === m && c.rho === rho && c.pi0 === pi0; });
    return hit.length ? hit[0] : null;
  }

  function renderCell() {
    var c = currentCell(), box = el('cellout');
    if (!box) { return; }
    if (!c) { box.innerHTML = '<p class="mut">No cell was simulated at that combination.</p>'; return; }
    var h = '<p>' + c.m + ' metrics at correlation ' + c.rho + ', ' + c.m0 + ' of them true nulls and ' + c.m1
      + ' carrying a real effect. ';
    if (c.m1 > 0) {
      h += 'Each real effect is sized so a single uncorrected metric would find it '
        + pct(c.uncorrected_power, 0) + ' of the time, so the power column has something to measure. ';
    } else {
      h += 'With nothing real to find there is no power to report, so every rejection in this cell is a false one. ';
    }
    h += 'Every row below is ' + SIM.simulations_per_cell + ' simulated experiments at alpha ' + SIM.alpha + '.</p>';
    h += '<div class="tbl-wrap"><table class="data"><thead><tr><th>Correction</th><th>At least one false winner</th>'
      + '<th>Mean false discovery proportion</th><th>Power, true effects found</th></tr></thead><tbody>';
    SIM.corrections.forEach(function (k) {
      var r = c.results[k];
      h += '<tr><td>' + LABEL[k] + '</td><td>' + pct(r.any_false_win) + '</td><td>' + pct(r.fdr_mean) + '</td><td>'
        + (r.power_mean == null ? 'no true effects in this cell' : pct(r.power_mean)) + '</td></tr>';
    });
    h += '</tbody></table></div>';
    box.innerHTML = h;
    var rep = el('repout');
    if (rep) { rep.innerHTML = ''; }
  }

  function replicate() {
    var c = currentCell();
    if (!c) { return; }
    var n = Math.round(SIM.simulations_per_cell / 100);
    if (n < MIN_REPLICATE) { n = MIN_REPLICATE; }
    if (n > MAX_REPLICATE) { n = MAX_REPLICATE; }
    var rnd = mulberry32(SIM.seed >>> 0), nrm = makeNormal(rnd);
    var alpha = SIM.alpha, mu = invNorm(1 - alpha / 2) + invNorm(c.uncorrected_power);
    var acc = {};
    SIM.corrections.forEach(function (k) { acc[k] = { fw: 0, fdr: 0, pow: 0 }; });
    var t;
    for (t = 0; t < n; t += 1) {
      var e = drawExperiment(nrm, c.m, c.rho, c.m1, mu);
      var adj = adjustAll(e.p);
      SIM.corrections.forEach(function (k) {
        var V = 0, R = 0, S = 0, i;
        for (i = 0; i < c.m; i += 1) {
          if (adj[k][i] <= alpha) { R += 1; if (e.truth[i]) { S += 1; } else { V += 1; } }
        }
        if (V > 0) { acc[k].fw += 1; }
        acc[k].fdr += R > 0 ? V / R : 0;
        if (c.m1 > 0) { acc[k].pow += S / c.m1; }
      });
    }
    var h = '<p class="mut">Your browser just ran ' + n + ' experiments on the same recipe, seed ' + SIM.seed
      + '. A run this size wobbles, so read it as landing near the published cell rather than on it.</p>';
    h += '<div class="tbl-wrap"><table class="data"><thead><tr><th>Correction</th><th>False winner, your run</th>'
      + '<th>False winner, published</th><th>Power, your run</th><th>Power, published</th></tr></thead><tbody>';
    SIM.corrections.forEach(function (k) {
      var r = c.results[k];
      h += '<tr><td>' + LABEL[k] + '</td><td>' + pct(acc[k].fw / n) + '</td><td>' + pct(r.any_false_win) + '</td><td>'
        + (c.m1 > 0 ? pct(acc[k].pow / n) : 'n/a') + '</td><td>'
        + (r.power_mean == null ? 'n/a' : pct(r.power_mean)) + '</td></tr>';
    });
    h += '</tbody></table></div>';
    var rep = el('repout');
    if (rep) { rep.innerHTML = h; }
  }

  function drawIntoBox() {
    if (!SIM) { return; }
    var ms = uniqSorted(SIM.cells.map(function (c) { return c.m; }));
    var rhos = uniqSorted(SIM.cells.map(function (c) { return c.rho; }));
    var pis = uniqSorted(SIM.cells.map(function (c) { return c.pi0; }));
    var m = ms[Math.min(1, ms.length - 1)];
    var rho = rhos[Math.min(1, rhos.length - 1)];
    var pi0 = pis[Math.min(1, pis.length - 1)];
    var hit = SIM.cells.filter(function (c) { return c.m === m && c.rho === rho && c.pi0 === pi0; });
    var c = hit.length ? hit[0] : SIM.cells[0];
    var rnd = mulberry32(SIM.seed >>> 0), nrm = makeNormal(rnd);
    var mu = invNorm(1 - SIM.alpha / 2) + invNorm(c.uncorrected_power);
    var e = drawExperiment(nrm, c.m, c.rho, c.m1, mu);
    var box = el('pvals');
    if (box) {
      box.value = e.p.map(function (p, i) { return lab(i) + ', ' + p.toPrecision(4); }).join('\n');
    }
    var who = c.m1 === 1 ? lab(0) + ' carries the real effect' : lab(0) + ' through ' + lab(c.m1 - 1) + ' carry a real effect';
    var note = el('drawnote');
    if (note) {
      note.textContent = 'Drawn, not typed. ' + c.m + ' metrics at correlation ' + c.rho + ', seed ' + SIM.seed
        + ', alpha ' + SIM.alpha + ', two sided p values from the same equicorrelated Gaussian recipe the study uses. '
        + who + ' and the other ' + c.m0 + ' are true nulls, so you can watch which corrections keep the right ones.';
    }
    compute();
  }

  /* Wiring. Every handler checks its node exists first. */
  var form = el('calc');
  if (form) { form.addEventListener('submit', function (ev) { ev.preventDefault(); compute(); }); }
  var drawBtn = el('draw');
  if (drawBtn) { drawBtn.addEventListener('click', drawIntoBox); }
  var clearBtn = el('clear');
  if (clearBtn) {
    clearBtn.addEventListener('click', function () {
      var box = el('pvals'), note = el('drawnote'), out = el('out');
      if (box) { box.value = ''; }
      if (note) { note.textContent = ''; }
      if (out) { out.innerHTML = ''; }
    });
  }
  var pv = el('pvals');
  if (pv) { pv.addEventListener('input', compute); }
  var av = el('alpha');
  if (av) { av.addEventListener('input', compute); }

  if (SIM) {
    var ms2 = uniqSorted(SIM.cells.map(function (c) { return c.m; }));
    fillSelect(el('selM'), ms2, function (v) { return v + ' metrics'; });
    fillSelect(el('selRho'), uniqSorted(SIM.cells.map(function (c) { return c.rho; })),
      function (v) { return v === 0 ? 'independent, rho 0' : 'rho ' + v; });
    fillSelect(el('selPi'), uniqSorted(SIM.cells.map(function (c) { return c.pi0; })),
      function (v) { return pct(v, 0) + ' true nulls'; });
    el('selM').value = ms2[ms2.length - 1];
    el('selRho').value = 0;
    el('selPi').value = 1;
    el('selM').addEventListener('change', renderCell);
    el('selRho').addEventListener('change', renderCell);
    el('selPi').addEventListener('change', renderCell);
    var rp = el('replicate');
    if (rp) { rp.addEventListener('click', replicate); }
    renderCell();
    drawIntoBox();
  } else {
    compute();
  }

  /* FAQ accordion, matching the site pattern. */
  Array.prototype.forEach.call(document.querySelectorAll('.faq-q'), function (q) {
    q.addEventListener('click', function () {
      var item = q.parentElement;
      if (item) { item.classList.toggle('open'); }
    });
  });

  /* Mobile nav toggle, matching the site pattern. */
  var mt = document.querySelector('.mobile-toggle');
  var nv = document.querySelector('header nav');
  if (mt && nv && mt.dataset.navWired !== '1') {
    mt.dataset.navWired = '1';
    mt.setAttribute('aria-expanded', 'false');
    mt.addEventListener('click', function () {
      var open = nv.classList.toggle('open');
      mt.setAttribute('aria-expanded', open ? 'true' : 'false');
    });
  }
}());

/* Pooled standard deviation calculator.
   Everything runs in the browser. No input leaves the page and nothing is sent to a server.

   The page's whole argument is that a pooled SD is only meaningful under an equal variance
   assumption, so this file always computes BOTH the pooled Student t and the Welch t and
   shows them side by side. It never silently picks one.

   Written to the same bounds as the rest of this site: every loop is bounded, every function
   validates its arguments, and no function exceeds the size limit. */
(function () {
  'use strict';

  var MAX_VALUES = 5000;      /* hard cap on pasted raw values, keeps every loop bounded */
  var LANCZOS = [676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012,
    9.9843695780195716e-6, 1.5056327351493116e-7];

  function el(id) { return document.getElementById(id); }
  function num(id) { var n = el(id); return n ? parseFloat(n.value) : NaN; }
  function fmt(x, d) {
    if (!isFinite(x)) { return 'n/a'; }
    return x.toFixed(d == null ? 4 : d);
  }
  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c];
    });
  }

  /* ---- distribution maths -------------------------------------------------- */

  function logGamma(z) {
    if (!isFinite(z) || z <= 0) { return NaN; }
    if (z < 0.5) { return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z); }
    var x = 0.99999999999980993, i;
    var zz = z - 1;
    for (i = 0; i < LANCZOS.length; i += 1) { x += LANCZOS[i] / (zz + i + 1); }
    var t = zz + LANCZOS.length - 0.5;
    return 0.5 * Math.log(2 * Math.PI) + (zz + 0.5) * Math.log(t) - t + Math.log(x);
  }

  /* Continued fraction for the regularized incomplete beta, bounded at 300 iterations. */
  function betaCf(x, a, b) {
    var tiny = 1e-30, qab = a + b, qap = a + 1, qam = a - 1;
    var c = 1, d = 1 - qab * x / qap, m, m2, aa, del;
    if (Math.abs(d) < tiny) { d = tiny; }
    d = 1 / d;
    var h = d;
    for (m = 1; m <= 300; m += 1) {
      m2 = 2 * m;
      aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < tiny) { d = tiny; }
      c = 1 + aa / c; if (Math.abs(c) < tiny) { c = tiny; }
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < tiny) { d = tiny; }
      c = 1 + aa / c; if (Math.abs(c) < tiny) { c = tiny; }
      d = 1 / d; del = d * c; h *= del;
      if (Math.abs(del - 1) < 3e-12) { break; }
    }
    return h;
  }

  function betaInc(x, a, b) {
    if (!isFinite(x) || !isFinite(a) || !isFinite(b)) { return NaN; }
    if (x <= 0) { return 0; }
    if (x >= 1) { return 1; }
    var lb = logGamma(a + b) - logGamma(a) - logGamma(b)
      + a * Math.log(x) + b * Math.log(1 - x);
    var front = Math.exp(lb);
    if (x < (a + 1) / (a + b + 2)) { return front * betaCf(x, a, b) / a; }
    return 1 - Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b)
      + b * Math.log(1 - x) + a * Math.log(x)) * betaCf(1 - x, b, a) / b;
  }

  /* Two sided p value for Student t with df degrees of freedom. */
  function tTailTwoSided(t, df) {
    if (!isFinite(t) || !isFinite(df) || df <= 0) { return NaN; }
    var x = df / (df + t * t);
    return betaInc(x, df / 2, 0.5);
  }

  /* Critical t by bisection on the tail, bounded at 200 steps. */
  function tCritical(df, alpha) {
    if (!isFinite(df) || df <= 0 || !isFinite(alpha) || alpha <= 0 || alpha >= 1) { return NaN; }
    var lo = 0, hi = 500, mid, p, i;
    for (i = 0; i < 200; i += 1) {
      mid = (lo + hi) / 2;
      p = tTailTwoSided(mid, df);
      if (p > alpha) { lo = mid; } else { hi = mid; }
    }
    return (lo + hi) / 2;
  }

  /* ---- input handling ------------------------------------------------------ */

  /* Parse pasted raw values. Anything that is not a finite number is dropped and counted,
     so a stray word or a header row can never reach the arithmetic. */
  function parseValues(raw) {
    var out = [], dropped = 0, i;
    var toks = String(raw || '').split(/[\s,;|]+/);
    var limit = Math.min(toks.length, MAX_VALUES);
    for (i = 0; i < limit; i += 1) {
      var t = toks[i].trim();
      if (t) {
        var v = parseFloat(t);
        if (isFinite(v)) { out.push(v); } else { dropped += 1; }
      }
    }
    return { values: out, dropped: dropped };
  }

  function describe(values) {
    var n = values.length, i, sum = 0, ss = 0;
    if (n < 2) { return null; }
    for (i = 0; i < n; i += 1) { sum += values[i]; }
    var mean = sum / n;
    for (i = 0; i < n; i += 1) { ss += (values[i] - mean) * (values[i] - mean); }
    return { n: n, mean: mean, sd: Math.sqrt(ss / (n - 1)) };
  }

  /* A group is either summary statistics or pasted raw values. Raw values win when present,
     because a reader who pasted data expects the data to be used. */
  function readGroup(nId, mId, sId, rId) {
    var parsed = parseValues(el(rId) ? el(rId).value : '');
    if (parsed.values.length >= 2) {
      var d = describe(parsed.values);
      d.source = 'raw'; d.dropped = parsed.dropped;
      return d;
    }
    return { n: num(nId), mean: num(mId), sd: num(sId), source: 'summary', dropped: 0 };
  }

  function groupValid(g) {
    if (!g) { return false; }
    return isFinite(g.n) && g.n >= 2 && isFinite(g.mean) && isFinite(g.sd) && g.sd >= 0;
  }

  /* ---- the statistics ------------------------------------------------------ */

  function pooledStats(a, b, alpha) {
    var dfP = a.n + b.n - 2;
    var sp2 = ((a.n - 1) * a.sd * a.sd + (b.n - 1) * b.sd * b.sd) / dfP;
    var sp = Math.sqrt(sp2);
    var se = sp * Math.sqrt(1 / a.n + 1 / b.n);
    var diff = a.mean - b.mean;
    var t = se > 0 ? diff / se : NaN;
    var tc = tCritical(dfP, alpha);
    return {
      sp: sp, df: dfP, se: se, t: t, p: tTailTwoSided(t, dfP),
      lo: diff - tc * se, hi: diff + tc * se,
      d: sp > 0 ? diff / sp : NaN
    };
  }

  function welchStats(a, b, alpha) {
    var va = a.sd * a.sd / a.n, vb = b.sd * b.sd / b.n;
    var se = Math.sqrt(va + vb);
    var diff = a.mean - b.mean;
    var num1 = (va + vb) * (va + vb);
    var den = (va * va) / (a.n - 1) + (vb * vb) / (b.n - 1);
    var df = den > 0 ? num1 / den : NaN;
    var t = se > 0 ? diff / se : NaN;
    var tc = tCritical(df, alpha);
    return { df: df, se: se, t: t, p: tTailTwoSided(t, df), lo: diff - tc * se, hi: diff + tc * se };
  }

  /* ---- rendering ----------------------------------------------------------- */

  function summaryRow(label, g) {
    var note = g.source === 'raw'
      ? 'computed from ' + g.n + ' pasted values'
      : 'entered as summary statistics';
    if (g.dropped > 0) { note += ', ' + g.dropped + ' non numeric value(s) ignored'; }
    return '<tr><td>' + esc(label) + '</td><td>' + g.n + '</td><td>' + fmt(g.mean, 4)
      + '</td><td>' + fmt(g.sd, 4) + '</td><td class="mut">' + esc(note) + '</td></tr>';
  }

  function verdictText(pooled, welch, ratio) {
    var agree = (pooled.p < 0.05) === (welch.p < 0.05);
    if (ratio > 4) {
      return 'The larger variance is ' + fmt(ratio, 2) + ' times the smaller, which is well past '
        + 'the point where pooling is safe. Report the Welch result.';
    }
    if (!agree) {
      return 'The two tests disagree at the 0.05 level, so the equal variance assumption is doing '
        + 'the work rather than your data. Report the Welch result.';
    }
    return 'The variance ratio is ' + fmt(ratio, 2) + ' and both tests agree, so pooling did no '
      + 'harm here. Welch is still the safer one to report.';
  }

  function render(a, b, alpha) {
    var pooled = pooledStats(a, b, alpha);
    var welch = welchStats(a, b, alpha);
    var hi = Math.max(a.sd, b.sd), lo = Math.min(a.sd, b.sd);
    var ratio = lo > 0 ? (hi * hi) / (lo * lo) : Infinity;

    var h = '<div class="stat-grid">'
      + '<div class="stat-card"><div class="stat-value">' + fmt(pooled.sp, 4) + '</div>'
      + '<div class="stat-label">pooled standard deviation</div></div>'
      + '<div class="stat-card"><div class="stat-value">' + fmt(pooled.d, 3) + '</div>'
      + '<div class="stat-label">Cohen&#39;s d</div></div>'
      + '<div class="stat-card"><div class="stat-value">' + fmt(ratio, 2) + '</div>'
      + '<div class="stat-label">variance ratio, larger over smaller</div></div>'
      + '<div class="stat-card"><div class="stat-value">' + fmt(a.mean - b.mean, 4) + '</div>'
      + '<div class="stat-label">difference in means</div></div>'
      + '</div>';

    h += '<div class="tbl-wrap"><table class="data"><thead><tr><th>Group</th><th>n</th>'
      + '<th>Mean</th><th>SD</th><th>Source</th></tr></thead><tbody>'
      + summaryRow('A', a) + summaryRow('B', b) + '</tbody></table></div>';

    h += '<div class="tbl-wrap"><table class="data"><thead><tr><th>Test</th><th>t</th>'
      + '<th>df</th><th>p, two sided</th><th>' + fmt((1 - alpha) * 100, 0)
      + '% interval for the difference</th></tr></thead><tbody>'
      + '<tr><td>Student, pooled</td><td>' + fmt(pooled.t, 4) + '</td><td>' + fmt(pooled.df, 0)
      + '</td><td>' + fmt(pooled.p, 5) + '</td><td>' + fmt(pooled.lo, 4) + ' to ' + fmt(pooled.hi, 4) + '</td></tr>'
      + '<tr><td>Welch, not pooled</td><td>' + fmt(welch.t, 4) + '</td><td>' + fmt(welch.df, 2)
      + '</td><td>' + fmt(welch.p, 5) + '</td><td>' + fmt(welch.lo, 4) + ' to ' + fmt(welch.hi, 4) + '</td></tr>'
      + '</tbody></table></div>';

    h += '<p>' + esc(verdictText(pooled, welch, ratio)) + '</p>';
    h += '<p class="mut">Pooled from <code class="formula">s_p = sqrt( ((' + a.n + ' - 1) &middot; '
      + fmt(a.sd, 4) + '^2 + (' + b.n + ' - 1) &middot; ' + fmt(b.sd, 4) + '^2) / '
      + fmt(pooled.df, 0) + ' )</code></p>';
    return h;
  }

  function compute() {
    var box = el('out');
    if (!box) { return; }
    var a = readGroup('nA', 'mA', 'sA', 'rA');
    var b = readGroup('nB', 'mB', 'sB', 'rB');
    var alpha = num('alpha');
    if (!isFinite(alpha) || alpha <= 0 || alpha >= 1) { alpha = 0.05; }
    if (!groupValid(a) || !groupValid(b)) {
      box.innerHTML = '<p class="mut">Enter a size of at least 2, a mean and a standard '
        + 'deviation for both groups, or paste at least two raw values into each.</p>';
      return;
    }
    box.innerHTML = render(a, b, alpha);
  }

  /* ---- wiring, every handler checks its node exists first -------------------- */

  var ids = ['nA', 'mA', 'sA', 'rA', 'nB', 'mB', 'sB', 'rB', 'alpha'];
  ids.forEach(function (id) {
    var node = el(id);
    if (node) { node.addEventListener('input', compute); }
  });

  var fill = el('fill');
  if (fill) {
    fill.addEventListener('click', function () {
      var ra = el('rA'), rb = el('rB');
      if (ra) { ra.value = '12.1 13.4 11.8 12.9 14.0 11.2 13.1 12.5 12.8 13.6'; }
      if (rb) { rb.value = '14.9 13.2 15.5 16.1 12.8 15.9 14.4 13.9 16.8 15.1'; }
      compute();
    });
  }

  var clear = el('clear');
  if (clear) {
    clear.addEventListener('click', function () {
      ['rA', 'rB'].forEach(function (id) { var n = el(id); if (n) { n.value = ''; } });
      compute();
    });
  }

  Array.prototype.forEach.call(document.querySelectorAll('.faq-q'), function (q) {
    q.addEventListener('click', function () {
      var item = q.parentElement;
      if (item) { item.classList.toggle('open'); }
    });
  });

  compute();
}());

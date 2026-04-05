/* ABWex — A/B Test Calculator */
/* Frequentist + Bayesian + Multi-variant + Revenue + Sample Size */
/* All functions under 60 lines. Loops bounded to max 1000. */

(function() {
  'use strict';

  var state = { mode: 'frequentist', activeTab: 'calculator' };

  /* ===== MATH UTILITIES ===== */

  /* Normal CDF using Abramowitz & Stegun rational approximation (error < 7.5e-8) */
  function normalCDF(x) {
    if (x < -8) return 0;
    if (x > 8) return 1;
    var neg = x < 0;
    if (neg) x = -x;
    var t = 1 / (1 + 0.2316419 * x);
    var d = 0.3989422804014327;
    var p = d * Math.exp(-0.5 * x * x) *
      (t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))));
    return neg ? p : 1 - p;
  }

  /* Inverse normal CDF (rational approximation) */
  function normalInv(p) {
    if (p <= 0) return -8;
    if (p >= 1) return 8;
    if (p < 0.5) return -normalInv(1 - p);
    var t = Math.sqrt(-2 * Math.log(1 - p));
    var c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
    var d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
    return t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t);
  }

  /* Log-gamma using Stirling's approximation (Lanczos) */
  function logGamma(z) {
    if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z);
    z -= 1;
    var g = 7;
    var coef = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
      771.32342877765313, -176.61502916214059, 12.507343278686905,
      -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    var x = coef[0];
    for (var i = 1; i < g + 2 && i < 9; i++) x += coef[i] / (z + i);
    var t = z + g + 0.5;
    return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
  }

  /* Beta function log */
  function logBeta(a, b) {
    return logGamma(a) + logGamma(b) - logGamma(a + b);
  }

  /* Beta PDF */
  function betaPDF(x, a, b) {
    if (x <= 0 || x >= 1) return 0;
    return Math.exp((a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) - logBeta(a, b));
  }

  /* P(B > A) via Simpson's rule numerical integration */
  function probBbeatsA(alphaA, betaA, alphaB, betaB) {
    var n = 500;
    var h = 1.0 / n;
    var sum = 0;
    for (var i = 0; i <= n; i++) {
      var x = i * h;
      if (x <= 0 || x >= 1) continue;
      var pdfB = betaPDF(x, alphaB, betaB);
      var cdfA = betaRegularized(x, alphaA, betaA);
      var w = (i === 0 || i === n) ? 1 : (i % 2 === 0) ? 2 : 4;
      sum += w * pdfB * cdfA;
    }
    return (h / 3) * sum;
  }

  /* Regularized incomplete beta via continued fraction (bounded iterations) */
  function betaRegularized(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    if (x > (a + 1) / (a + b + 2)) return 1 - betaRegularized(1 - x, b, a);
    var lnpfx = a * Math.log(x) + b * Math.log(1 - x) - logBeta(a, b) - Math.log(a);
    var pfx = Math.exp(lnpfx);
    var m = 1, prev = 0, curr = 1;
    for (var i = 1; i <= 200; i++) {
      var idx = i;
      var d;
      if (idx % 2 === 0) {
        var k = idx / 2;
        d = (k * (b - k) * x) / ((a + 2 * k - 1) * (a + 2 * k));
      } else {
        var k2 = (idx - 1) / 2;
        d = -((a + k2) * (a + b + k2) * x) / ((a + 2 * k2) * (a + 2 * k2 + 1));
      }
      var next = curr + d * prev;
      prev = curr;
      curr = next;
      if (Math.abs(curr - prev) < 1e-10) break;
    }
    return pfx / curr;
  }

  /* ===== FREQUENTIST CALCULATION ===== */

  function calcFrequentist(n1, c1, n2, c2) {
    var p1 = c1 / n1, p2 = c2 / n2;
    var pPool = (c1 + c2) / (n1 + n2);
    var se = Math.sqrt(pPool * (1 - pPool) * (1 / n1 + 1 / n2));
    if (se === 0) return { error: 'Standard error is zero. Need more data.' };
    var z = (p2 - p1) / se;
    var pValue = 2 * (1 - normalCDF(Math.abs(z)));
    var seDiff = Math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2);
    var zCrit = 1.96;
    var ciLow = (p2 - p1) - zCrit * seDiff;
    var ciHigh = (p2 - p1) + zCrit * seDiff;
    var relImprove = p1 > 0 ? ((p2 - p1) / p1) * 100 : 0;
    return {
      p1: p1, p2: p2, z: z, pValue: pValue,
      ci: [ciLow, ciHigh], relImprove: relImprove,
      formula: 'z = (' + pct(p2) + ' - ' + pct(p1) + ') / sqrt(' + pct(pPool) + ' * ' + pct(1 - pPool) + ' * (1/' + n1 + ' + 1/' + n2 + ')) = ' + z.toFixed(4)
    };
  }

  /* ===== BAYESIAN CALCULATION ===== */

  function calcBayesian(n1, c1, n2, c2) {
    var alphaA = c1 + 1, betaA = n1 - c1 + 1;
    var alphaB = c2 + 1, betaB = n2 - c2 + 1;
    var pBbeatsA = probBbeatsA(alphaA, betaA, alphaB, betaB);
    var meanA = alphaA / (alphaA + betaA);
    var meanB = alphaB / (alphaB + betaB);
    var expectedLoss = Math.abs(meanB - meanA) * (1 - Math.max(pBbeatsA, 1 - pBbeatsA));
    return {
      alphaA: alphaA, betaA: betaA, alphaB: alphaB, betaB: betaB,
      pBbeatsA: pBbeatsA, meanA: meanA, meanB: meanB,
      expectedLoss: expectedLoss
    };
  }

  /* ===== SAMPLE SIZE CALCULATOR ===== */

  function calcSampleSize(baseline, mde, power, alpha) {
    var p1 = baseline;
    var p2 = p1 * (1 + mde / 100);
    var zAlpha = normalInv(1 - alpha / 2);
    var zBeta = normalInv(power);
    var delta = p2 - p1;
    if (delta === 0) return { error: 'MDE cannot be zero.' };
    var n = Math.pow(zAlpha + zBeta, 2) * (p1 * (1 - p1) + p2 * (1 - p2)) / (delta * delta);
    return { n: Math.ceil(n), p1: p1, p2: p2, delta: delta };
  }

  /* ===== RENDERING HELPERS ===== */

  function pct(v) { return (v * 100).toFixed(2) + '%'; }
  function pctShort(v) { return (v * 100).toFixed(1) + '%'; }

  function sigClass(pVal) {
    if (pVal < 0.05) return 'sig-green';
    if (pVal < 0.1) return 'sig-yellow';
    return 'sig-red';
  }

  function sigLabel(pVal) {
    if (pVal < 0.01) return 'Highly Significant';
    if (pVal < 0.05) return 'Statistically Significant';
    if (pVal < 0.1) return 'Marginally Significant';
    return 'Not Significant';
  }

  /* ===== SVG CHART BUILDERS ===== */

  function buildBarChart(variants) {
    var w = 500, h = 200, pad = 50;
    var maxRate = 0;
    for (var i = 0; i < variants.length && i < 10; i++) {
      var r = variants[i].conv / variants[i].visitors;
      if (r > maxRate) maxRate = r;
    }
    maxRate = maxRate * 1.3;
    if (maxRate === 0) maxRate = 0.1;
    var barW = Math.min(80, (w - pad * 2) / variants.length - 20);
    var colors = ['#2DD4BF', '#3B82F6', '#A855F7', '#F97316'];
    var svg = '<svg viewBox="0 0 ' + w + ' ' + h + '" xmlns="http://www.w3.org/2000/svg">';
    svg += '<line x1="' + pad + '" y1="' + (h - 30) + '" x2="' + (w - 10) + '" y2="' + (h - 30) + '" stroke="#2a2a2a" stroke-width="1"/>';
    for (var j = 0; j < variants.length && j < 4; j++) {
      var rate = variants[j].conv / variants[j].visitors;
      var barH = (rate / maxRate) * (h - 60);
      var x = pad + j * ((w - pad * 2) / variants.length) + 10;
      var y = h - 30 - barH;
      svg += '<rect x="' + x + '" y="' + y + '" width="' + barW + '" height="' + barH + '" fill="' + colors[j] + '" rx="3" opacity="0.85"/>';
      svg += '<text x="' + (x + barW / 2) + '" y="' + (y - 8) + '" text-anchor="middle" fill="#e8e8e8" font-size="12" font-family="Space Mono">' + pctShort(rate) + '</text>';
      svg += '<text x="' + (x + barW / 2) + '" y="' + (h - 12) + '" text-anchor="middle" fill="#888" font-size="11" font-family="IBM Plex Sans">' + variants[j].name + '</text>';
    }
    svg += '</svg>';
    return svg;
  }

  function buildBetaCurves(aA, bA, aB, bB) {
    var w = 500, h = 200, pad = 40;
    var pts = 200;
    var maxY = 0;
    var dataA = [], dataB = [];
    for (var i = 0; i <= pts; i++) {
      var x = 0.001 + (i / pts) * 0.998;
      var yA = betaPDF(x, aA, bA);
      var yB = betaPDF(x, aB, bB);
      dataA.push({ x: x, y: yA });
      dataB.push({ x: x, y: yB });
      if (yA > maxY) maxY = yA;
      if (yB > maxY) maxY = yB;
    }
    maxY *= 1.1;
    if (maxY === 0) maxY = 1;

    function toPath(data, color) {
      var path = 'M';
      for (var k = 0; k < data.length && k < 300; k++) {
        var px = pad + (data[k].x) * (w - pad * 2);
        var py = h - 30 - (data[k].y / maxY) * (h - 50);
        path += (k === 0 ? '' : ' L') + px.toFixed(1) + ' ' + py.toFixed(1);
      }
      return '<path d="' + path + '" fill="none" stroke="' + color + '" stroke-width="2" opacity="0.85"/>';
    }

    var svg = '<svg viewBox="0 0 ' + w + ' ' + h + '" xmlns="http://www.w3.org/2000/svg">';
    svg += '<line x1="' + pad + '" y1="' + (h - 30) + '" x2="' + (w - 10) + '" y2="' + (h - 30) + '" stroke="#2a2a2a" stroke-width="1"/>';
    svg += toPath(dataA, '#2DD4BF');
    svg += toPath(dataB, '#3B82F6');
    svg += '<text x="' + (w - 80) + '" y="20" fill="#2DD4BF" font-size="11" font-family="Space Mono">Control A</text>';
    svg += '<text x="' + (w - 80) + '" y="36" fill="#3B82F6" font-size="11" font-family="Space Mono">Variant B</text>';
    svg += '</svg>';
    return svg;
  }

  function buildPowerCurve(baseline, alpha, samples) {
    var w = 500, h = 200, pad = 50;
    var pts = [];
    for (var i = 1; i <= 50; i++) {
      var mde = i * 0.5;
      var p2 = baseline * (1 + mde / 100);
      var delta = p2 - baseline;
      if (delta <= 0) continue;
      var zAlpha = normalInv(1 - alpha / 2);
      var se = Math.sqrt(baseline * (1 - baseline) / samples + p2 * (1 - p2) / samples);
      if (se === 0) continue;
      var zBeta = Math.abs(delta) / se - zAlpha;
      var power = normalCDF(zBeta);
      pts.push({ mde: mde, power: power });
    }
    if (pts.length === 0) return '';
    var svg = '<svg viewBox="0 0 ' + w + ' ' + h + '" xmlns="http://www.w3.org/2000/svg">';
    svg += '<line x1="' + pad + '" y1="' + (h - 30) + '" x2="' + (w - 10) + '" y2="' + (h - 30) + '" stroke="#2a2a2a"/>';
    svg += '<line x1="' + pad + '" y1="10" x2="' + pad + '" y2="' + (h - 30) + '" stroke="#2a2a2a"/>';
    svg += '<text x="' + (w / 2) + '" y="' + (h - 5) + '" text-anchor="middle" fill="#888" font-size="10">Minimum Detectable Effect (%)</text>';
    svg += '<text x="15" y="' + (h / 2) + '" fill="#888" font-size="10" transform="rotate(-90,15,' + (h / 2) + ')">Power</text>';
    var path = 'M';
    for (var j = 0; j < pts.length && j < 100; j++) {
      var px = pad + (pts[j].mde / 25) * (w - pad * 2);
      var py = h - 30 - pts[j].power * (h - 50);
      path += (j === 0 ? '' : ' L') + px.toFixed(1) + ' ' + py.toFixed(1);
    }
    svg += '<path d="' + path + '" fill="none" stroke="#2DD4BF" stroke-width="2"/>';
    svg += '<line x1="' + pad + '" y1="' + (h - 30 - 0.8 * (h - 50)) + '" x2="' + (w - 10) + '" y2="' + (h - 30 - 0.8 * (h - 50)) + '" stroke="#EAB308" stroke-dasharray="4" opacity="0.5"/>';
    svg += '<text x="' + (w - 50) + '" y="' + (h - 30 - 0.8 * (h - 50) - 5) + '" fill="#EAB308" font-size="9">80% power</text>';
    svg += '</svg>';
    return svg;
  }

  /* ===== CALCULATOR TAB ===== */

  function buildCalcUI() {
    var el = document.getElementById('calc-content');
    if (!el) return;
    var html = '<div class="mode-toggle">';
    html += '<span>Frequentist</span>';
    html += '<div class="toggle-switch" id="mode-toggle" onclick="window.AB.toggleMode()"></div>';
    html += '<span>Bayesian</span></div>';
    html += '<div class="card"><h3>Control (A)</h3>';
    html += '<div class="input-group">';
    html += '<div class="input-field"><label>Visitors</label><input type="number" id="n1" value="5000" min="1"></div>';
    html += '<div class="input-field"><label>Conversions</label><input type="number" id="c1" value="200" min="0"></div>';
    html += '</div></div>';
    html += '<div class="card"><h3>Variant (B)</h3>';
    html += '<div class="input-group">';
    html += '<div class="input-field"><label>Visitors</label><input type="number" id="n2" value="5000" min="1"></div>';
    html += '<div class="input-field"><label>Conversions</label><input type="number" id="c2" value="240" min="0"></div>';
    html += '</div></div>';
    html += '<div class="card"><h3>Revenue Impact (optional)</h3>';
    html += '<div class="input-group">';
    html += '<div class="input-field"><label>Revenue per conversion ($)</label><input type="number" id="rev" value="50" min="0" step="0.01"></div>';
    html += '<div class="input-field"><label>Monthly visitors</label><input type="number" id="monthly" value="100000" min="0"></div>';
    html += '</div></div>';
    html += '<button class="btn btn-primary" onclick="window.AB.calculate()">Calculate Results</button>';
    html += '<div id="calc-results"></div>';
    html += '<div id="calc-chart" class="chart-container"></div>';
    html += '<div id="calc-formula" class="result-box" style="display:none;"></div>';
    el.innerHTML = html;
  }

  function calculate() {
    var n1 = parseInt(document.getElementById('n1').value, 10) || 1;
    var c1 = parseInt(document.getElementById('c1').value, 10) || 0;
    var n2 = parseInt(document.getElementById('n2').value, 10) || 1;
    var c2 = parseInt(document.getElementById('c2').value, 10) || 0;
    var rev = parseFloat(document.getElementById('rev').value) || 0;
    var monthly = parseInt(document.getElementById('monthly').value, 10) || 0;

    if (c1 > n1 || c2 > n2) {
      document.getElementById('calc-results').innerHTML = '<div class="result-box" style="color:var(--red);">Conversions cannot exceed visitors.</div>';
      return;
    }

    var html = '';
    if (state.mode === 'frequentist') {
      var r = calcFrequentist(n1, c1, n2, c2);
      if (r.error) {
        document.getElementById('calc-results').innerHTML = '<div class="result-box" style="color:var(--red);">' + r.error + '</div>';
        return;
      }
      html += '<div class="stat-grid">';
      html += statCard(pctShort(r.p1), 'Control Rate');
      html += statCard(pctShort(r.p2), 'Variant Rate');
      html += statCard('<span class="' + sigClass(r.pValue) + '">' + r.pValue.toFixed(4) + '</span>', 'P-Value');
      html += statCard((r.relImprove >= 0 ? '+' : '') + r.relImprove.toFixed(1) + '%', 'Relative Change');
      html += '</div>';
      html += '<div class="stat-grid">';
      html += statCard(r.z.toFixed(3), 'Z-Score');
      html += statCard(pctShort(r.ci[0]) + ' to ' + pctShort(r.ci[1]), '95% CI (diff)');
      html += statCard('<span class="' + sigClass(r.pValue) + '">' + sigLabel(r.pValue) + '</span>', 'Verdict');
      html += '</div>';
      var formulaEl = document.getElementById('calc-formula');
      formulaEl.style.display = 'block';
      formulaEl.innerHTML = '<span style="color:var(--text-dim);font-size:0.8rem;">' + r.formula + '</span>';
    } else {
      var b = calcBayesian(n1, c1, n2, c2);
      html += '<div class="stat-grid">';
      html += statCard(pctShort(b.meanA), 'Control Mean');
      html += statCard(pctShort(b.meanB), 'Variant Mean');
      html += statCard((b.pBbeatsA * 100).toFixed(1) + '%', 'P(B beats A)');
      html += statCard(pctShort(b.expectedLoss), 'Expected Loss');
      html += '</div>';
      document.getElementById('calc-chart').innerHTML = '<div class="card"><h3>Beta Distribution Curves</h3>' + buildBetaCurves(b.alphaA, b.betaA, b.alphaB, b.betaB) + '</div>';
      var formulaEl2 = document.getElementById('calc-formula');
      formulaEl2.style.display = 'block';
      formulaEl2.innerHTML = '<span style="color:var(--text-dim);font-size:0.8rem;">Beta(A) ~ Beta(' + b.alphaA + ', ' + b.betaA + '), Beta(B) ~ Beta(' + b.alphaB + ', ' + b.betaB + ')</span>';
    }

    if (rev > 0 && monthly > 0) {
      var p1 = c1 / n1, p2v = c2 / n2;
      var monthlyGain = (p2v - p1) * monthly * rev;
      html += '<div class="card revenue-card"><h3>Revenue Impact</h3>';
      html += '<div class="stat-grid">';
      html += statCard('$' + monthlyGain.toFixed(0), 'Monthly Revenue Change');
      html += statCard('$' + (monthlyGain * 12).toFixed(0), 'Annual Projection');
      html += '</div></div>';
    }

    if (state.mode === 'frequentist') {
      document.getElementById('calc-chart').innerHTML = '<div class="card"><h3>Conversion Rate Comparison</h3>' + buildBarChart([
        { name: 'Control', visitors: n1, conv: c1 },
        { name: 'Variant', visitors: n2, conv: c2 }
      ]) + '</div>';
    }

    document.getElementById('calc-results').innerHTML = html;
  }

  function statCard(value, label) {
    return '<div class="stat-card"><div class="stat-value">' + value + '</div><div class="stat-label">' + label + '</div></div>';
  }

  /* ===== MULTI-VARIANT TAB ===== */

  function buildMultiUI() {
    var el = document.getElementById('multi-content');
    if (!el) return;
    var html = '<div class="card"><h3>Multi-Variant Comparison (up to 4)</h3>';
    var names = ['A (Control)', 'B', 'C', 'D'];
    var visitors = [5000, 5000, 5000, 5000];
    var convs = [200, 240, 220, 210];
    for (var i = 0; i < 4; i++) {
      html += '<div class="variant-row">';
      html += '<div class="variant-label variant-' + 'abcd'[i] + '">' + names[i] + '</div>';
      html += '<div class="input-field"><label>Visitors</label><input type="number" id="mv-n' + i + '" value="' + visitors[i] + '" min="1"></div>';
      html += '<div class="input-field"><label>Conversions</label><input type="number" id="mv-c' + i + '" value="' + convs[i] + '" min="0"></div>';
      html += '<div></div>';
      html += '</div>';
    }
    html += '<p style="color:var(--text-dim);font-size:0.8rem;margin-bottom:12px;">Leave visitors at 0 to exclude a variant.</p>';
    html += '<button class="btn btn-primary" onclick="window.AB.calcMulti()">Compare All</button></div>';
    html += '<div id="multi-results"></div>';
    html += '<div id="multi-chart" class="chart-container"></div>';
    el.innerHTML = html;
  }

  
function betaRandom(alpha, beta) {
  var u, v, s;
  for (var i = 0; i < 200; i++) {
    u = Math.pow(Math.random(), 1.0 / alpha);
    v = Math.pow(Math.random(), 1.0 / beta);
    s = u + v;
    if (s <= 1.0) return u / s;
  }
  return alpha / (alpha + beta);
}

function calcMulti() {
    var variants = [];
    var names = ['A', 'B', 'C', 'D'];
    for (var i = 0; i < 4; i++) {
      var n = parseInt(document.getElementById('mv-n' + i).value, 10) || 0;
      var c = parseInt(document.getElementById('mv-c' + i).value, 10) || 0;
      if (n > 0) variants.push({ name: names[i], visitors: n, conv: c, rate: c / n });
    }
    if (variants.length < 2) {
      document.getElementById('multi-results').innerHTML = '<div class="result-box" style="color:var(--red);">Need at least 2 variants.</div>';
      return;
    }

    var html = '<div class="card"><h3>Pairwise Comparisons</h3>';
    var bestIdx = 0;
    for (var b = 1; b < variants.length; b++) {
      if (variants[b].rate > variants[bestIdx].rate) bestIdx = b;
    }

    for (var x = 0; x < variants.length && x < 10; x++) {
      for (var y = x + 1; y < variants.length && y < 10; y++) {
        var r = calcFrequentist(variants[x].visitors, variants[x].conv, variants[y].visitors, variants[y].conv);
        if (r.error) continue;
        html += '<div style="padding:8px 0;border-bottom:1px solid var(--border);font-size:0.9rem;">';
        html += '<strong>' + variants[x].name + ' vs ' + variants[y].name + ':</strong> ';
        html += 'p=' + r.pValue.toFixed(4) + ' <span class="' + sigClass(r.pValue) + '">(' + sigLabel(r.pValue) + ')</span>';
        html += ' | ' + variants[y].name + ' is ' + (r.relImprove >= 0 ? '+' : '') + r.relImprove.toFixed(1) + '% vs ' + variants[x].name;
        html += '</div>';
      }
    }

    html += '<div style="padding:12px 0;"><span class="winner-badge">Winner: Variant ' + variants[bestIdx].name + ' (' + pctShort(variants[bestIdx].rate) + ')</span></div>';
    html += '</div>';

    if (state.mode === 'bayesian') {
      html += '<div class="card"><h3>Bayesian: Probability Each Variant Is Best</h3>';
      html += '<div class="stat-grid">';
      for (var v = 0; v < variants.length && v < 4; v++) {
        var wins = 0, trials = 500;
        for (var t = 0; t < trials; t++) {
          var isBest = true;
      var sampleV = betaRandom(variants[v].conv + 1, variants[v].visitors - variants[v].conv + 1);
          for (var o = 0; o < variants.length && o < 4; o++) {
            if (o === v) continue;
            var sampleO = betaRandom(variants[o].conv + 1, variants[o].visitors - variants[o].conv + 1);
        if (sampleO >= sampleV) { isBest = false; break; }
          }
          if (isBest) wins++;
        }
        html += statCard((wins / trials * 100).toFixed(0) + '%', 'P(' + variants[v].name + ' best)');
      }
      html += '</div></div>';
    }

    document.getElementById('multi-results').innerHTML = html;
    document.getElementById('multi-chart').innerHTML = '<div class="card"><h3>Conversion Rates</h3>' + buildBarChart(variants) + '</div>';
  }

  /* ===== SAMPLE SIZE TAB ===== */

  function buildSampleUI() {
    var el = document.getElementById('sample-content');
    if (!el) return;
    var html = '<div class="card"><h3>Sample Size Pre-Calculator</h3>';
    html += '<div class="input-group">';
    html += '<div class="input-field"><label>Baseline conversion rate (%)</label><input type="number" id="ss-baseline" value="4" min="0.1" max="100" step="0.1"></div>';
    html += '<div class="input-field"><label>Min. detectable effect (%)</label><input type="number" id="ss-mde" value="10" min="0.1" step="0.1"></div>';
    html += '<div class="input-field"><label>Power (%)</label><select id="ss-power"><option value="0.8">80%</option><option value="0.9">90%</option><option value="0.95">95%</option></select></div>';
    html += '<div class="input-field"><label>Significance level</label><select id="ss-alpha"><option value="0.05">0.05 (95%)</option><option value="0.01">0.01 (99%)</option><option value="0.1">0.10 (90%)</option></select></div>';
    html += '</div>';
    html += '<div class="input-group">';
    html += '<div class="input-field"><label>Daily visitors (for duration estimate)</label><input type="number" id="ss-daily" value="1000" min="1"></div>';
    html += '</div>';
    html += '<button class="btn btn-primary" onclick="window.AB.calcSample()">Calculate Sample Size</button></div>';
    html += '<div id="sample-results"></div>';
    html += '<div id="sample-chart" class="chart-container"></div>';
    el.innerHTML = html;
  }

  function calcSample() {
    var baseline = parseFloat(document.getElementById('ss-baseline').value) / 100;
    var mde = parseFloat(document.getElementById('ss-mde').value);
    var power = parseFloat(document.getElementById('ss-power').value);
    var alpha = parseFloat(document.getElementById('ss-alpha').value);
    var daily = parseInt(document.getElementById('ss-daily').value, 10) || 1;

    var r = calcSampleSize(baseline, mde, power, alpha);
    if (r.error) {
      document.getElementById('sample-results').innerHTML = '<div class="result-box" style="color:var(--red);">' + r.error + '</div>';
      return;
    }

    var daysPerVariant = Math.ceil(r.n / (daily / 2));
    var html = '<div class="stat-grid">';
    html += statCard(r.n.toLocaleString(), 'Sample per Variant');
    html += statCard((r.n * 2).toLocaleString(), 'Total Sample Needed');
    html += statCard(daysPerVariant + ' days', 'Est. Duration');
    html += statCard(pctShort(r.p1) + ' -> ' + pctShort(r.p2), 'Baseline -> Target');
    html += '</div>';

    html += '<div class="card"><h3>When to Stop Advisor</h3>';
    html += '<p style="color:var(--text-dim);font-size:0.85rem;">Do not peek at results before collecting ' + r.n.toLocaleString() + ' visitors per variant. Early stopping inflates false positive rates due to the multiple testing problem. If you must peek, use sequential testing methods (alpha spending functions) that adjust significance thresholds at each interim analysis.</p>';
    var currentN = parseInt(document.getElementById('n1').value, 10) || 0;
    if (currentN > 0 && currentN < r.n) {
      var remaining = r.n - currentN;
      html += '<p style="color:var(--yellow);font-size:0.85rem;margin-top:8px;">Based on your calculator inputs, you need approximately ' + remaining.toLocaleString() + ' more visitors per variant to reach significance at this effect size.</p>';
    }
    html += '</div>';

    document.getElementById('sample-results').innerHTML = html;
    document.getElementById('sample-chart').innerHTML = '<div class="card"><h3>Power Curve (at current sample size)</h3>' +
      '<p style="color:var(--text-dim);font-size:0.8rem;margin-bottom:8px;">Shows probability of detecting various effect sizes with ' + r.n.toLocaleString() + ' visitors per variant.</p>' +
      buildPowerCurve(baseline, alpha, r.n) + '</div>';
  }

  /* ===== TAB / MODE SWITCHING ===== */

  function switchTab(tab) {
    state.activeTab = tab;
    var btns = document.querySelectorAll('.tab-btn');
    var panes = document.querySelectorAll('.tab-pane');
    for (var i = 0; i < btns.length && i < 20; i++) {
      btns[i].classList.toggle('active', btns[i].dataset.tab === tab);
    }
    for (var j = 0; j < panes.length && j < 20; j++) {
      panes[j].classList.toggle('active', panes[j].id === tab + '-content');
    }
  }

  function toggleMode() {
    state.mode = state.mode === 'frequentist' ? 'bayesian' : 'frequentist';
    var toggle = document.getElementById('mode-toggle');
    if (toggle) toggle.classList.toggle('active', state.mode === 'bayesian');
  }

  /* ===== FAQ ===== */

  function initFAQ() {
    var items = document.querySelectorAll('.faq-item');
    for (var i = 0; i < items.length && i < 50; i++) {
      items[i].querySelector('.faq-q').addEventListener('click', (function(item) {
        return function() { item.classList.toggle('open'); };
      })(items[i]));
    }
  }

  function initNav() {
    var toggle = document.querySelector('.mobile-toggle');
    var nav = document.querySelector('nav');
    if (toggle && nav) {
      toggle.addEventListener('click', function() { nav.classList.toggle('open'); });
    }
  }

  /* ===== INIT ===== */

  function init() {
    initNav();
    buildCalcUI();
    buildMultiUI();
    buildSampleUI();
    initFAQ();
    var tabBtns = document.querySelectorAll('.tab-btn');
    for (var i = 0; i < tabBtns.length && i < 20; i++) {
      tabBtns[i].addEventListener('click', (function(btn) {
        return function() { switchTab(btn.dataset.tab); };
      })(tabBtns[i]));
    }
  }

  window.AB = {
    calculate: calculate,
    calcMulti: calcMulti,
    calcSample: calcSample,
    toggleMode: toggleMode,
    switchTab: switchTab
  };

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();

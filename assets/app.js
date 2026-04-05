/* A/B Test Calculator - abwex.com */
(function() {
  'use strict';

  var bayesianMode = false;

  // Normal CDF approximation (Abramowitz and Stegun)
  function normalCDF(x) {
    var t = 1 / (1 + 0.2316419 * Math.abs(x));
    var d = 0.3989422804014327;
    var p = d * Math.exp(-x * x / 2) *
      (t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.8212560 + t * 1.3302744)))));
    return x > 0 ? 1 - p : p;
  }

  // Inverse normal (approximation)
  function normalInv(p) {
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    if (p < 0.5) return -normalInv(1 - p);
    var t = Math.sqrt(-2 * Math.log(1 - p));
    var c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
    var d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
    return t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t);
  }

  // Log-beta function for Bayesian
  function logBeta(a, b) {
    return logGamma(a) + logGamma(b) - logGamma(a + b);
  }

  // Stirling log-gamma approximation
  function logGamma(z) {
    if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z);
    z -= 1;
    var x = 0.99999999999980993;
    var p = [676.5203681218851, -1259.1392167224028, 771.32342877765313,
      -176.61502916214059, 12.507343278686905, -0.13857109526572012,
      9.9843695780195716e-6, 1.5056327351493116e-7];
    for (var i = 0; i < p.length; i++) x += p[i] / (z + i + 1);
    var t = z + p.length - 0.5;
    return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
  }

  // P(B > A) using Beta distribution Monte Carlo-free closed form approx
  function bayesianProbBBeatsA(alphaA, betaA, alphaB, betaB) {
    var total = 0;
    var steps = 200;
    for (var i = 0; i < steps; i++) {
      var x = (i + 0.5) / steps;
      var logPdfB = (alphaB - 1) * Math.log(x) + (betaB - 1) * Math.log(1 - x) - logBeta(alphaB, betaB);
      // CDF of A at x using regularized incomplete beta (numerical integration)
      var cdfA = 0;
      var innerSteps = 200;
      for (var j = 0; j < innerSteps; j++) {
        var t = (j + 0.5) / innerSteps * x;
        var logPdfA = (alphaA - 1) * Math.log(Math.max(t, 1e-15)) +
          (betaA - 1) * Math.log(Math.max(1 - t, 1e-15)) - logBeta(alphaA, betaA);
        cdfA += Math.exp(logPdfA) * (x / innerSteps);
      }
      total += Math.exp(logPdfB) * cdfA * (1 / steps);
    }
    return Math.min(1, Math.max(0, total));
  }

  function calculate() {
    var cVisitors = parseInt(document.getElementById('control-visitors').value, 10);
    var cConv = parseInt(document.getElementById('control-conversions').value, 10);
    var vVisitors = parseInt(document.getElementById('variant-visitors').value, 10);
    var vConv = parseInt(document.getElementById('variant-conversions').value, 10);

    if (isNaN(cVisitors) || isNaN(cConv) || isNaN(vVisitors) || isNaN(vConv)) return;
    if (cVisitors <= 0 || vVisitors <= 0) return;
    if (cConv > cVisitors || vConv > vVisitors) return;

    var p1 = cConv / cVisitors;
    var p2 = vConv / vVisitors;
    var pPooled = (cConv + vConv) / (cVisitors + vVisitors);
    var se = Math.sqrt(pPooled * (1 - pPooled) * (1 / cVisitors + 1 / vVisitors));
    var z = se > 0 ? (p2 - p1) / se : 0;
    var pValue = 2 * (1 - normalCDF(Math.abs(z)));
    var improvement = p1 > 0 ? ((p2 - p1) / p1) * 100 : 0;

    var sigLevel, sigText;
    if (pValue < 0.01) { sigLevel = '99%'; sigText = 'Highly significant'; }
    else if (pValue < 0.05) { sigLevel = '95%'; sigText = 'Significant'; }
    else if (pValue < 0.1) { sigLevel = '90%'; sigText = 'Marginally significant'; }
    else { sigLevel = '<90%'; sigText = 'Not significant'; }

    var isWinner = pValue < 0.05;
    var winner = isWinner ? (p2 > p1 ? 'Variant' : 'Control') : null;

    // Standard error for error bars
    var se1 = Math.sqrt(p1 * (1 - p1) / cVisitors) * 1.96;
    var se2 = Math.sqrt(p2 * (1 - p2) / vVisitors) * 1.96;

    renderResults(p1, p2, improvement, pValue, sigLevel, sigText, winner, se1, se2);

    // Bayesian
    if (bayesianMode) {
      var alphaA = 1 + cConv, betaA = 1 + cVisitors - cConv;
      var alphaB = 1 + vConv, betaB = 1 + vVisitors - vConv;
      var probBBeatsA = bayesianProbBBeatsA(alphaA, betaA, alphaB, betaB);
      renderBayesian(probBBeatsA);
    }
  }

  function renderResults(p1, p2, improvement, pValue, sigLevel, sigText, winner, se1, se2) {
    var el = document.getElementById('freq-results');
    el.classList.add('visible');

    var improvColor = improvement > 0 ? '#34D399' : improvement < 0 ? '#F87171' : 'var(--text)';
    var improvSign = improvement > 0 ? '+' : '';

    el.querySelector('.results-grid').innerHTML =
      '<div class="result-item"><div class="value">' + (p1 * 100).toFixed(2) + '%</div><div class="label">Control Rate</div></div>' +
      '<div class="result-item"><div class="value">' + (p2 * 100).toFixed(2) + '%</div><div class="label">Variant Rate</div></div>' +
      '<div class="result-item"><div class="value" style="color:' + improvColor + '">' + improvSign + improvement.toFixed(1) + '%</div><div class="label">Relative Change</div></div>' +
      '<div class="result-item"><div class="value">' + pValue.toFixed(4) + '</div><div class="label">p-value</div></div>' +
      '<div class="result-item"><div class="value">' + sigLevel + '</div><div class="label">Confidence</div><div class="sub">' + sigText + '</div></div>';

    // Verdict
    var verdictEl = el.querySelector('.verdict');
    if (winner) {
      verdictEl.className = 'verdict verdict-winner';
      verdictEl.textContent = winner + ' wins with ' + sigLevel + ' confidence (' + improvSign + improvement.toFixed(1) + '% change)';
    } else {
      verdictEl.className = 'verdict verdict-inconclusive';
      verdictEl.textContent = 'Not enough data to declare a winner (p = ' + pValue.toFixed(3) + ')';
    }

    // Bar chart
    var maxRate = Math.max(p1 + se1, p2 + se2, 0.01);
    var scale = 130 / maxRate;
    var h1 = Math.max(4, p1 * scale);
    var h2 = Math.max(4, p2 * scale);
    var errH1 = se1 * scale;
    var errH2 = se2 * scale;

    el.querySelector('.bar-chart').innerHTML =
      '<div class="bar-group">' +
        '<div class="bar-value">' + (p1 * 100).toFixed(2) + '%</div>' +
        '<div class="bar-wrapper"><div class="bar" style="height:' + h1 + 'px;background:var(--text-muted);position:relative;">' +
          '<div class="error-line" style="height:' + errH1 + 'px;top:-' + errH1 + 'px;"><div class="error-cap error-cap-top"></div><div class="error-cap" style="bottom:0;position:absolute;left:50%;transform:translateX(-50%);width:14px;height:2px;background:var(--text-muted);"></div></div>' +
        '</div></div>' +
        '<div class="bar-label">Control</div></div>' +
      '<div class="bar-group">' +
        '<div class="bar-value" style="color:var(--accent)">' + (p2 * 100).toFixed(2) + '%</div>' +
        '<div class="bar-wrapper"><div class="bar" style="height:' + h2 + 'px;background:var(--accent);position:relative;">' +
          '<div class="error-line" style="height:' + errH2 + 'px;top:-' + errH2 + 'px;"><div class="error-cap error-cap-top"></div><div class="error-cap" style="bottom:0;position:absolute;left:50%;transform:translateX(-50%);width:14px;height:2px;background:var(--text-muted);"></div></div>' +
        '</div></div>' +
        '<div class="bar-label">Variant</div></div>';
  }

  function renderBayesian(prob) {
    var el = document.getElementById('bayesian-result');
    el.style.display = 'block';
    el.innerHTML = '<div class="value">' + (prob * 100).toFixed(1) + '%</div>' +
      '<div class="label">Probability that Variant beats Control (Bayesian)</div>';
  }

  function calculateSampleSize() {
    var baseline = parseFloat(document.getElementById('ss-baseline').value) / 100;
    var mde = parseFloat(document.getElementById('ss-mde').value) / 100;
    var power = parseFloat(document.getElementById('ss-power').value) / 100;
    var alpha = parseFloat(document.getElementById('ss-alpha').value);

    if (isNaN(baseline) || isNaN(mde) || isNaN(power) || isNaN(alpha)) return;
    if (baseline <= 0 || baseline >= 1 || mde <= 0 || power <= 0 || power >= 1 || alpha <= 0 || alpha >= 1) return;

    var p1 = baseline;
    var p2 = baseline + baseline * mde;
    var zAlpha = normalInv(1 - alpha / 2);
    var zBeta = normalInv(power);

    var n = Math.ceil(
      Math.pow(zAlpha * Math.sqrt(2 * p1 * (1 - p1)) + zBeta * Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2)), 2) /
      Math.pow(p2 - p1, 2)
    );

    var el = document.getElementById('ss-results');
    el.classList.add('visible');
    el.innerHTML = '<div class="results-grid">' +
      '<div class="result-item"><div class="value">' + n.toLocaleString() + '</div><div class="label">Per Variation</div></div>' +
      '<div class="result-item"><div class="value">' + (n * 2).toLocaleString() + '</div><div class="label">Total Required</div></div>' +
      '<div class="result-item"><div class="value">' + (p2 * 100).toFixed(2) + '%</div><div class="label">Target Rate</div><div class="sub">(' + (baseline * 100).toFixed(1) + '% + ' + (mde * 100).toFixed(0) + '% MDE)</div></div>' +
      '</div>';
  }

  window.abwex = {
    calculate: calculate,
    calculateSampleSize: calculateSampleSize,
    switchTab: function(tab, btn) {
      var tabs = document.querySelectorAll('.tab-btn');
      var contents = document.querySelectorAll('.tab-content');
      for (var i = 0; i < tabs.length; i++) tabs[i].classList.remove('active');
      for (var j = 0; j < contents.length; j++) contents[j].classList.remove('active');
      btn.classList.add('active');
      document.getElementById('tab-' + tab).classList.add('active');
    },
    toggleBayesian: function() {
      var toggle = document.getElementById('bayesian-toggle');
      bayesianMode = !bayesianMode;
      if (bayesianMode) {
        toggle.classList.add('on');
      } else {
        toggle.classList.remove('on');
        document.getElementById('bayesian-result').style.display = 'none';
      }
    }
  };
})();

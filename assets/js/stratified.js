/* Stratified sampling calculator.
   Runs entirely in the browser. No input leaves the page.

   Two allocations, proportional and Neyman, plus the design effect against simple random
   sampling, which is the number that says whether stratifying bought anything at all.

   Bounds and validation follow the same rules as the rest of this site: every loop is
   bounded, every function validates its arguments, no function exceeds the size limit. */
(function () {
  'use strict';

  var MAX_STRATA = 200;   /* hard cap, keeps every loop bounded */

  function el(id) { return document.getElementById(id); }
  function fmt(x, d) { return isFinite(x) ? x.toFixed(d == null ? 2 : d) : 'n/a'; }
  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c];
    });
  }

  /* One line is name, population size, standard deviation, and optionally a stratum mean.
     The mean is what makes the design effect honest. Without it the only comparison
     available is against the within stratum variance, and proportional allocation scores
     exactly 1.000 against that by construction, which would read as "stratifying is
     useless" when it is really "you have not told me the thing that decides".
     A line that does not yield a usable size and spread is dropped and counted. */
  function parseLine(line, index) {
    var parts = String(line).split(/[,;\t|]+/);
    if (parts.length < 3) { return null; }
    var name = parts[0].trim() || ('Stratum ' + (index + 1));
    var N = parseFloat(parts[1]);
    var S = parseFloat(parts[2]);
    if (!isFinite(N) || N <= 0) { return null; }
    if (!isFinite(S) || S < 0) { return null; }
    var mean = parts.length > 3 ? parseFloat(parts[3]) : NaN;
    return { name: name, N: N, S: S, mean: isFinite(mean) ? mean : null };
  }

  function parseStrata(raw) {
    var lines = String(raw || '').split('\n');
    var out = [], dropped = 0, i;
    var limit = Math.min(lines.length, MAX_STRATA);
    for (i = 0; i < limit; i += 1) {
      var t = lines[i].trim();
      if (t) {
        var row = parseLine(t, out.length);
        if (row) { out.push(row); } else { dropped += 1; }
      }
    }
    return { strata: out, dropped: dropped, truncated: Math.max(0, lines.length - limit) };
  }

  /* Largest remainder rounding, so the per stratum counts always sum to exactly n. */
  function roundToTotal(shares, total) {
    var out = [], i, used = 0;
    for (i = 0; i < shares.length; i += 1) {
      var v = Math.floor(shares[i]);
      out.push(v); used += v;
    }
    var rema = [];
    for (i = 0; i < shares.length; i += 1) { rema.push({ i: i, r: shares[i] - Math.floor(shares[i]) }); }
    rema.sort(function (a, b) { return b.r - a.r; });
    var k = 0;
    while (used < total && k < rema.length * 2) {
      out[rema[k % rema.length].i] += 1;
      used += 1; k += 1;
    }
    return out;
  }

  /* A stratum can never give up more units than it contains. Where the cap binds we say so. */
  function capToPopulation(counts, strata) {
    var capped = 0, i, out = [];
    for (i = 0; i < counts.length; i += 1) {
      var c = counts[i];
      if (c > strata[i].N) { c = Math.floor(strata[i].N); capped += 1; }
      if (c < 0) { c = 0; }
      out.push(c);
    }
    return { counts: out, capped: capped };
  }

  function proportional(strata, total) {
    var i, N = 0;
    for (i = 0; i < strata.length; i += 1) { N += strata[i].N; }
    var shares = [];
    for (i = 0; i < strata.length; i += 1) { shares.push(total * strata[i].N / N); }
    return capToPopulation(roundToTotal(shares, total), strata);
  }

  function neyman(strata, total) {
    var i, denom = 0;
    for (i = 0; i < strata.length; i += 1) { denom += strata[i].N * strata[i].S; }
    if (!(denom > 0)) { return proportional(strata, total); }
    var shares = [];
    for (i = 0; i < strata.length; i += 1) {
      shares.push(total * (strata[i].N * strata[i].S) / denom);
    }
    return capToPopulation(roundToTotal(shares, total), strata);
  }

  /* Variance of the stratified mean, sum of W_h^2 S_h^2 / n_h. */
  function stratifiedVariance(strata, counts) {
    var i, N = 0, v = 0;
    for (i = 0; i < strata.length; i += 1) { N += strata[i].N; }
    for (i = 0; i < strata.length; i += 1) {
      if (!(counts[i] > 0)) { return Infinity; }
      var W = strata[i].N / N;
      v += W * W * strata[i].S * strata[i].S / counts[i];
    }
    return v;
  }

  function pooledWithin(strata) {
    var i, N = 0, v = 0;
    for (i = 0; i < strata.length; i += 1) { N += strata[i].N; }
    for (i = 0; i < strata.length; i += 1) {
      v += (strata[i].N / N) * strata[i].S * strata[i].S;
    }
    return v;
  }

  function haveAllMeans(strata) {
    var i;
    for (i = 0; i < strata.length; i += 1) {
      if (strata[i].mean === null) { return false; }
    }
    return strata.length > 0;
  }

  /* Between stratum variance, sum of W_h (mu_h - mu)^2. This is the entire source of the
     gain from stratifying. Total population variance is within plus between, and a simple
     random sample works against the total while a stratified sample only fights the within
     part. If every stratum mean is identical this term is zero and stratification really
     does buy nothing, which is a genuine finding rather than an artefact. */
  function betweenVariance(strata) {
    var i, N = 0, mu = 0, v = 0;
    for (i = 0; i < strata.length; i += 1) { N += strata[i].N; }
    for (i = 0; i < strata.length; i += 1) { mu += (strata[i].N / N) * strata[i].mean; }
    for (i = 0; i < strata.length; i += 1) {
      var d = strata[i].mean - mu;
      v += (strata[i].N / N) * d * d;
    }
    return v;
  }

  /* The variance a simple random sample of the same total size would deliver. With means
     supplied this uses the true population variance. Without them it falls back to the
     within stratum part only, which makes the reported gain conservative. */
  function simpleRandomVariance(strata, total) {
    var within = pooledWithin(strata);
    var between = haveAllMeans(strata) ? betweenVariance(strata) : 0;
    return { v: (within + between) / total, within: within, between: between,
      exact: haveAllMeans(strata) };
  }

  function allocationRows(strata, prop, ney, useNeyman) {
    var i, N = 0, h = '';
    for (i = 0; i < strata.length; i += 1) { N += strata[i].N; }
    for (i = 0; i < strata.length; i += 1) {
      var chosen = useNeyman ? ney.counts[i] : prop.counts[i];
      h += '<tr><td>' + esc(strata[i].name) + '</td>'
        + '<td>' + fmt(strata[i].N, 0) + '</td>'
        + '<td>' + fmt(100 * strata[i].N / N, 1) + '%</td>'
        + '<td>' + fmt(strata[i].S, 3) + '</td>'
        + '<td>' + prop.counts[i] + '</td>'
        + '<td>' + ney.counts[i] + '</td>'
        + '<td>' + chosen + '</td></tr>';
    }
    return h;
  }

  function render(strata, total, useNeyman, notes) {
    var prop = proportional(strata, total);
    var ney = neyman(strata, total);
    var vProp = stratifiedVariance(strata, prop.counts);
    var vNey = stratifiedVariance(strata, ney.counts);
    var srs = simpleRandomVariance(strata, total);
    var vSrs = srs.v;
    var chosenV = useNeyman ? vNey : vProp;

    var deffProp = vSrs > 0 ? vProp / vSrs : NaN;
    var deffNey = vSrs > 0 ? vNey / vSrs : NaN;
    var deff = useNeyman ? deffNey : deffProp;
    var gain = isFinite(deff) && deff > 0 ? (1 / deff) : NaN;

    var h = '<div class="stat-grid">'
      + '<div class="stat-card"><div class="stat-value">' + fmt(Math.sqrt(chosenV), 4) + '</div>'
      + '<div class="stat-label">standard error of the mean</div></div>'
      + '<div class="stat-card"><div class="stat-value">' + fmt(deff, 3) + '</div>'
      + '<div class="stat-label">design effect, below 1 is a gain</div></div>'
      + '<div class="stat-card"><div class="stat-value">' + fmt(gain, 2) + 'x</div>'
      + '<div class="stat-label">effective sample size multiplier</div></div>'
      + '<div class="stat-card"><div class="stat-value">' + strata.length + '</div>'
      + '<div class="stat-label">strata</div></div>'
      + '</div>';

    h += '<div class="tbl-wrap"><table class="data"><thead><tr><th>Stratum</th><th>Population</th>'
      + '<th>Share</th><th>SD</th><th>Proportional</th><th>Neyman</th><th>Chosen</th>'
      + '</tr></thead><tbody>' + allocationRows(strata, prop, ney, useNeyman)
      + '</tbody></table></div>';

    h += allocationNote(strata, prop, ney, deffProp, deffNey);
    h += deffNote(srs, deff, gain);
    h += cautionNote(useNeyman ? ney.capped : prop.capped, notes.dropped, srs.exact);
    return h;
  }

  /* Whether the two allocations actually differ, and by how much. */
  function allocationNote(strata, prop, ney, deffProp, deffNey) {
    var diff = 0, i;
    for (i = 0; i < strata.length; i += 1) { diff += Math.abs(prop.counts[i] - ney.counts[i]); }
    if (diff === 0) {
      return '<p>Proportional and Neyman give the same split here, so the spreads are close '
        + 'enough that weighting by them changes nothing. Use proportional, it needs less '
        + 'information.</p>';
    }
    return '<p>The two allocations differ by ' + fmt(diff / 2, 0) + ' units moved between strata. '
      + 'Proportional gives a design effect of ' + fmt(deffProp, 3) + ' and Neyman '
      + fmt(deffNey, 3) + '.</p>';
  }

  /* The design effect needs its own explanation, because a 1.000 means two completely
     different things depending on whether stratum means were supplied. */
  function deffNote(srs, deff, gain) {
    if (!srs.exact) {
      return '<p>No stratum means were supplied, so the comparison above uses only the within '
        + 'stratum variance. Proportional allocation scores exactly 1.000 against that by '
        + 'construction, so read a 1.000 here as not yet measured rather than as no gain. '
        + 'Add a fourth number on each line, the stratum mean, and the real design effect '
        + 'is computed from the between stratum spread, which is where the whole gain lives.</p>';
    }
    if (isFinite(deff) && deff >= 0.98) {
      return '<p>Your stratum means are close together, so the between stratum variance is small '
        + 'and stratification genuinely buys almost nothing here. That is a real result, not a '
        + 'missing input. Stratify on something that separates the metric more.</p>';
    }
    if (!isFinite(deff)) { return ''; }
    return '<p>Between stratum variance is ' + fmt(srs.between, 4) + ' against a within stratum '
      + 'variance of ' + fmt(srs.within, 4) + '. Stratifying removes the between part from '
      + 'the error, which is where the ' + fmt(gain, 2) + 'x comes from.</p>';
  }

  function cautionNote(capped, dropped, exact) {
    var h = '';
    if (capped > 0) {
      h += '<p class="mut">' + capped + ' stratum was capped at its own population size, so the '
        + 'allocation shown is the largest one that is actually possible.</p>';
    }
    if (dropped > 0) {
      h += '<p class="mut">' + dropped + ' line(s) ignored, each needs a name, a population '
        + 'size and a standard deviation separated by commas.</p>';
    }
    h += '<p class="mut">Comparison is against a simple random sample of the same total size. '
      + (exact
        ? 'Stratum means were supplied, so this uses the full population variance, within plus between.'
        : 'Without stratum means this uses the within stratum variance only, which understates the gain.')
      + '</p>';
    return h;
  }

  function compute() {
    var box = el('out');
    if (!box) { return; }
    var parsed = parseStrata(el('strata') ? el('strata').value : '');
    var total = el('total') ? parseInt(el('total').value, 10) : NaN;
    var useNeyman = el('cost') ? el('cost').value === '1' : true;
    if (parsed.strata.length < 2) {
      box.innerHTML = '<p class="mut">Enter at least two strata, one per line, as name, '
        + 'population size, standard deviation.</p>';
      return;
    }
    if (!isFinite(total) || total < parsed.strata.length) {
      box.innerHTML = '<p class="mut">Total sample size must be a number at least as large as '
        + 'the number of strata, so every stratum can receive at least one unit.</p>';
      return;
    }
    box.innerHTML = render(parsed.strata, total, useNeyman, parsed);
  }

  ['strata', 'total', 'cost'].forEach(function (id) {
    var n = el(id);
    if (n) { n.addEventListener('input', compute); n.addEventListener('change', compute); }
  });

  var clear = el('clear');
  if (clear) {
    clear.addEventListener('click', function () {
      var n = el('strata');
      if (n) { n.value = ''; }
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

'use strict';
{
  const byId = id => document.getElementById(id);
  const scriptUrl = document.currentScript ? document.currentScript.src : new URL('./ffh-ui.js', document.baseURI).href;
  const examples = {
    sparse: [[1,4,0],[3,0,2]], moderate: [[8,2,1],[1,5,4]], balanced: [[2,2,2],[2,2,2]],
    diagonal: [[3,0,0],[0,3,0],[0,0,3]], mixed: [[2,1,0],[0,2,1],[1,0,2]],
    asymmetric: [[4,1,0],[1,2,1],[0,1,2]], compatibility: [[1,9],[11,3]]
  };
  const counts = byId('counts');
  let worker = null, runId = 0, certificate = null, timeout = null;
  function el(tag, text, className) {
    const node = document.createElement(tag);
    if (text !== undefined) node.textContent = text;
    if (className) node.className = className;
    return node;
  }
  function fmt(value) {
    if (!Number.isFinite(value)) return 'Not available';
    if (value === 0) return '0';
    return Math.abs(value) < 0.0001 ? value.toExponential(6) : Number(value.toPrecision(8)).toString();
  }
  function probability(value) {
    return value === 0 ? 'Below numerical display range' : fmt(value);
  }
  function parseTable() {
    if (counts.value.length > 4096) throw new Error('Limit the counts table to 4,096 characters.');
    const raw = counts.value.trim();
    if (!raw) throw new Error('Enter at least two rows of counts.');
    const lines = raw.split(/\r?\n/);
    if (lines.length < 2 || lines.length > 4) throw new Error('Enter 2–4 rows, with one row per group.');
    const table = lines.map((line, i) => {
      const tokens = line.trim().includes(',') ? line.trim().split(',').map(s => s.trim()) : line.trim().split(/\s+/);
      if (tokens.length < 2 || tokens.length > 4) throw new Error('Row ' + (i+1) + ' needs 2–4 counts.');
      return tokens.map((token, j) => {
        if (!/^\d+$/.test(token) || !Number.isSafeInteger(Number(token))) throw new Error('Row ' + (i+1) + ', column ' + (j+1) + ': enter a nonnegative whole count.');
        return Number(token);
      });
    });
    if (table.some(row => row.length !== table[0].length)) throw new Error('Every row must have the same number of columns.');
    const total = table.flat().reduce((a,b) => a+b,0);
    if (!Number.isSafeInteger(total)) throw new Error('The total exceeds the supported safe integer range.');
    if (!total) throw new Error('The table must contain observations.');
    table.forEach((row,i) => { if (!row.some(Boolean)) throw new Error('Row ' + (i+1) + ' has a zero total. Remove this empty category.'); });
    table[0].forEach((_,j) => { if (!table.some(row => row[j] > 0)) throw new Error('Column ' + (j+1) + ' has a zero total. Remove this empty category.'); });
    return table;
  }
  function margins(table) {
    return { rows: table.map(row => row.reduce((a,b) => a+b,0)), cols: table[0].map((_,j) => table.reduce((a,row) => a+row[j],0)), total: table.flat().reduce((a,b) => a+b,0) };
  }
  function matrix(parent, table, title, showMargins, expected) {
    parent.replaceChildren();
    const wrap = el('div', undefined, 'table-scroll'); wrap.tabIndex = 0; wrap.setAttribute('role','region'); wrap.setAttribute('aria-label', title + '; scroll horizontally if needed');
    const grid = el('table'); grid.append(el('caption',title));
    const head = el('thead'), header = el('tr');
    const corner = el('th','Group / Outcome'); corner.scope = 'col'; header.append(corner);
    table[0].forEach((_,j) => { const cell = el('th','Outcome ' + (j+1)); cell.scope = 'col'; header.append(cell); });
    if (showMargins) { const cell = el('th','Total'); cell.scope = 'col'; header.append(cell); }
    head.append(header); grid.append(head);
    const body = el('tbody'), totals = margins(table);
    table.forEach((row,i) => {
      const tr = el('tr'), th = el('th','Group ' + (i+1)); th.scope = 'row'; tr.append(th);
      row.forEach(value => {
        const text = expected ? value.toFixed(3) + (value < 1 ? ' ‡' : value < 5 ? ' †' : '') : String(value);
        tr.append(el('td',text));
      });
      if (showMargins) tr.append(el('td',totals.rows[i]));
      body.append(tr);
    }); grid.append(body);
    if (showMargins) {
      const foot = el('tfoot'), tr = el('tr'), th = el('th','Total'); th.scope = 'row'; tr.append(th);
      totals.cols.forEach(value => tr.append(el('td',value))); tr.append(el('td',totals.total)); foot.append(tr); grid.append(foot);
    }
    wrap.append(grid); parent.append(wrap);
  }
  function rCode(table) {
    return 'observed <- matrix(c(' + table.map(row => row.join(', ')).join(',\n                     ') + '),\n                   nrow = ' + table.length + ', byrow = TRUE)\nfisher.test(observed, hybrid = FALSE, simulate.p.value = FALSE)';
  }
  function preview(announceError) {
    try {
      const table = parseTable();
      matrix(byId('preview'),table,'Observed counts and fixed margins',true,false);
      byId('r-code').textContent = rCode(table);
      counts.removeAttribute('aria-invalid'); byId('input-error').textContent = '';
      return table;
    } catch (error) {
      byId('preview').replaceChildren(); byId('r-code').textContent = 'Enter a valid counts table to generate the R command.';
      if (announceError) { counts.setAttribute('aria-invalid','true'); byId('input-error').textContent = error.message; }
      return null;
    }
  }
  function stop() {
    runId += 1;
    if (worker) worker.terminate(); worker = null;
    if (timeout) clearTimeout(timeout); timeout = null;
    byId('run').disabled = false; byId('cancel').hidden = true;
  }
  function invalidate(message) {
    stop(); certificate = null; byId('certificate').hidden = true; byId('results').hidden = true;
    byId('status').textContent = message;
  }
  function stat(label,value) {
    const card = el('div',undefined,'stat-card'); card.append(el('div',value,'stat-value'),el('div',label,'stat-label')); byId('metrics').append(card);
  }
  function renderCompletionStatus(result, complete, numericalError) {
    byId('completion').textContent = complete ? 'Enumeration completed. All feasible tables were included in the reference distribution.' : numericalError ? 'Enumeration finished, but numerical verification failed. No exact p-value is available. ' + (result.error || '') : 'Exact calculation not completed. ' + (result.error || result.status || 'The calculation stopped before completion.') + ' No exact p-value is available.';
    byId('status').textContent = complete ? 'Exact calculation completed. Results and audit are ready.' : numericalError ? 'Enumeration finished, but numerical verification failed. No exact result is available.' : 'Calculation stopped without an exact result. See the explanation and R command below.';
  }
  function renderSummary(result, alphaText, complete, numericalError) {
    renderCompletionStatus(result, complete, numericalError);
    if (complete) {
      stat('Two-sided exact p-value',probability(result.pValue));
      stat('Observed-table probability',probability(result.observedProbability));
      stat('Feasible tables enumerated',fmt(result.tableCount));
      stat('Tables included in the tail',fmt(result.tailCount));
      const association=result.exactDecision==='association';
      byId('interpretation').textContent = association ? 'At α = ' + alphaText + ', the exact test supports an association between group and outcome because p < α. It does not identify a winning variant or a significant individual cell.' : 'At α = ' + alphaText + ', the exact test does not provide sufficient evidence of an association because p is not below α. This does not establish independence or equivalence.';
    } else {
      stat('Exact p-value',numericalError ? 'Withheld' : 'Not computed');
      byId('interpretation').textContent = numericalError ? 'Enumeration finished, but its numerical checks failed. The exact p-value is withheld. Expected-count diagnostics remain available; use the R command for an independent calculation.' : 'Partial probability sums are not reported as an exact p-value. The expected-count diagnostics remain available. Use the R command for an independent calculation.';
    }
  }
  function comparisonRows(table, result, complete, numericalError, approx) {
    return [[table.length === 2 && table[0].length === 2 ? 'Fisher exact, two-sided' : 'Fisher–Freeman–Halton exact',complete?probability(result.pValue):(numericalError?'Withheld':'Not computed'),complete?'Completed enumeration':(numericalError?'Enumeration finished; numerical verification failed':'Incomplete')],['Pearson chi-square',approx===null?'Not available':probability(approx),'Approximation; no continuity correction']];
  }
  function renderComparison(table, result, complete, numericalError, approx) {
    const method = el('div',undefined,'table-scroll'); method.tabIndex=0; method.setAttribute('role','region'); method.setAttribute('aria-label','Method comparison');
    const comparison = el('table'); comparison.append(el('caption','Exact and approximate methods'));
    const head=el('thead'), hr=el('tr'); ['Method','p-value','Status'].forEach(label => { const th=el('th',label);th.scope='col';hr.append(th); });head.append(hr);comparison.append(head);
    const body=el('tbody');
    comparisonRows(table, result, complete, numericalError, approx).forEach(row => { const tr=el('tr');row.forEach((v,i)=>{const cell=el(i===0?'th':'td',v);if(i===0)cell.scope='row';tr.append(cell);});body.append(tr); });
    comparison.append(body);method.append(comparison);byId('comparison').append(method);
    if (Number.isFinite(result.chiSquare)) byId('comparison').append(el('p','Pearson χ² = ' + fmt(result.chiSquare) + '; degrees of freedom = ' + fmt(result.df) + '.'));
    if (complete && approx !== null) {
      byId('comparison').append(el('p','Absolute p-value difference: ' + fmt(Math.abs(result.pValue-approx)) + '. Pearson is an approximation; this calculator does not use it for the exact-test decision.'));
    }
  }
  function renderDiagnostics(expected) {
    matrix(byId('expected'),expected,'Expected counts under independence',false,true);
    const values=expected.flat(), small=values.filter(v=>v<5).length, tiny=values.filter(v=>v<1).length;
    byId('diagnostics').textContent='Minimum expected count: ' + fmt(Math.min(...values)) + '. ' + small + '/' + values.length + ' cells (' + fmt(100*small/values.length) + '%) below 5; ' + tiny + ' below 1.';
  }
  function renderAudit(table, result, alpha, alphaText, complete, numericalError, total) {
    if (complete) {
      byId('audit').append(el('p','Probability ordering includes tables no more probable than the observed table, including ties. Total probability mass across the complete reference set: ' + fmt(result.massSum) + '.'));
      if (result.method) byId('audit').append(el('p', result.method + '.'));
      if (result.tieMethod || result.tieTolerance !== undefined) byId('audit').append(el('p','Tie handling: ' + (result.tieMethod || ('numerical tolerance ' + result.tieTolerance)) + '.'));
      const tailExamples=Array.isArray(result.tailExamples)?result.tailExamples:[];
      tailExamples.forEach((item,i)=>{const div=el('div',undefined,'tail-item');div.append(el('p','Tail example '+(i+1)+' · probability '+probability(item.probability)+(item.tie?' · ties the observed probability':'')),el('pre',item.table.map(row=>row.join('\t')).join('\n')));byId('tail-examples').append(div);});
      if (!tailExamples.length) byId('tail-examples').append(el('p','No tail examples were returned by the engine.'));
      certificate={formatVersion:1,calculator:'ABWex Fisher–Freeman–Halton',calculatorRelease:'2026-09-12.1',route:'https://abwex.com/tools/fisher-freeman-halton-calculator/',generatedAt:new Date().toISOString(),observed:table,margins:total,alpha,alphaText,decisionRule:'association only when exact rational p < entered decimal alpha',method:'Fixed-margin conditional exact test; probability ordering including ties',result};
      byId('certificate').hidden=false;
    } else {
      byId('audit').append(el('p',numericalError ? 'Enumeration finished, but numerical verification failed. No verified exact certificate is available. Enumerated tables: '+fmt(result.tableCount)+'.' : 'Enumeration is incomplete; there is no completed exact certificate. Tables visited before stopping: '+fmt(result.tableCount)+'.'));
      certificate=null;byId('certificate').hidden=true;
    }
  }
  function renderResult(table, result, alpha, alphaText) {
    byId('results').hidden = false;
    ['metrics','comparison','audit','tail-examples'].forEach(id => byId(id).replaceChildren());
    const numericalError = result.status === 'numerical-error';
    const complete = !numericalError && result.complete === true && Number.isFinite(result.pValue) && (result.exactDecision==='association'||result.exactDecision==='insufficient');
    const approx = Number.isFinite(result.chiSquarePValue) ? result.chiSquarePValue : Number.isFinite(result.chiSquareP) ? result.chiSquareP : null;
    const total = margins(table);
    const expected = Array.isArray(result.expected) ? result.expected : table.map((row,i) => row.map((_,j) => total.rows[i]*total.cols[j]/total.total));
    renderSummary(result, alphaText, complete, numericalError);
    renderComparison(table, result, complete, numericalError, approx);
    renderDiagnostics(expected);
    renderAudit(table, result, alpha, alphaText, complete, numericalError, total);
  }
  counts.addEventListener('input',()=>{invalidate('Counts changed. Run the test to calculate new results.');preview(false);});
  byId('alpha').addEventListener('input',()=>{invalidate('Threshold changed. Run the test to update the interpretation.');byId('alpha').removeAttribute('aria-invalid');});
  byId('load-example').addEventListener('click',()=>{invalidate('Synthetic example loaded. Run the test to calculate results.');counts.value=examples[byId('example').value].map(row=>row.join(', ')).join('\n');preview(false);});
  byId('reset').addEventListener('click',()=>{invalidate('Cleared. Enter a counts table to begin.');counts.value='';byId('alpha').value='0.05';byId('input-error').textContent='';counts.removeAttribute('aria-invalid');byId('alpha').removeAttribute('aria-invalid');preview(false);counts.focus();});
  byId('cancel').addEventListener('click',()=>{invalidate('Exact calculation cancelled. No exact p-value was produced.');});
  byId('ffh-form').addEventListener('submit',event=>{
    event.preventDefault();invalidate('Validating counts.');
    const table=preview(true);if(!table){byId('status').textContent='Fix the counts table before running the test.';counts.focus();return;}
    const alphaText=byId('alpha').value.trim(),alpha=Number(alphaText);
    if(!alphaText||alphaText.length>64||!Number.isFinite(alpha)||alpha<=0||alpha>=1){byId('alpha').setAttribute('aria-invalid','true');byId('input-error').textContent='Enter a decision threshold strictly between 0 and 1.';byId('status').textContent='Fix the decision threshold before running the test.';byId('alpha').focus();return;}
    const current=runId;byId('run').disabled=true;byId('cancel').hidden=false;byId('status').textContent='Enumerating tables with these fixed margins. You can cancel while the calculation runs.';
    function finish(result){if(current!==runId)return;stop();renderResult(table,result,alpha,alphaText);}
    try {
      worker=new Worker(new URL('./ffh-worker.js',scriptUrl));
      worker.onmessage=event=>{if(event.data.runId===current)finish(event.data.result);};
      worker.onerror=()=>finish({complete:false,status:'error',error:'The calculation worker could not run. Reload the page or use the R command.'});
      timeout=setTimeout(()=>finish({complete:false,status:'limit',error:'The browser time limit of 20 seconds was reached.'}),20000);
      worker.postMessage({runId:current,table,alphaText});
    } catch(error){finish({complete:false,status:'error',error:error.message});}
  });
  byId('certificate').addEventListener('click',()=>{
    if(!certificate)return;
    const blob=new Blob([JSON.stringify(certificate,null,2)+'\n'],{type:'application/json'}), url=URL.createObjectURL(blob), link=el('a');
    link.href=url;link.download='abwex-ffh-certificate.json';document.body.append(link);link.click();link.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
  preview(false);
}

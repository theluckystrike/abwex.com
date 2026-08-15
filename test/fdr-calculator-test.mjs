import fs from 'fs';
const src = fs.readFileSync('/Users/mike/Desktop/sites/abwex-repo/assets/js/fdr-calculator.js','utf8');
// Pull the pure math out of the IIFE and evaluate it standalone.
const names = ['orderAsc','harmonic','bonf','holm','stepUp','bh','by','adjustAll','normCdf','invNorm'];
let code = '';
for (const n of names) {
  const m = src.match(new RegExp('\\n  function '+n+'\\s*\\([^)]*\\)\\s*\\{'));
  if (!m) { console.log('MISSING', n); continue; }
  let i = m.index + m[0].length - 1, depth = 0, start = m.index;
  while (i < src.length) { if (src[i]==='{') depth++; else if (src[i]==='}') { depth--; if (!depth) break; } i++; }
  code += src.slice(start, i+1) + '\n';
}
const consts = 'var MAXM=200, INV_NORM_STEPS=120, MAX_LINES=5000;\n';
const F = new Function(consts + code + '\nreturn {'+names.join(',')+'};')();
let fails = 0;
function eq(name, a, b, tol=1e-12){ const ok=Math.abs(a-b)<=tol; if(!ok){fails++;console.log('  FAIL',name,a,'!=',b);} return ok; }

// known vector, cross-checked against statsmodels earlier in this session
const p = [0.001,0.008,0.039,0.041,0.042,0.06,0.074,0.205,0.212,0.216];
const adj = F.adjustAll(p);
const EXP_BH_first = 0.0100, EXP_BY_first = 0.0293;
eq('BH[0]', +adj.bh[0].toFixed(4), EXP_BH_first, 1e-9);
eq('BY[0]', +adj.by[0].toFixed(4), EXP_BY_first, 1e-9);
const keeps = k => adj[k].filter(v=>v<=0.05).length;
eq('BH keeps', keeps('bh'), 2); eq('BY keeps', keeps('by'), 1);
eq('Holm keeps', keeps('holm'), 1); eq('Bonf keeps', keeps('bonferroni'), 1);
// monotonicity of every adjusted vector when sorted by raw p
for (const k of ['bh','by','holm','bonferroni']) {
  const idx = F.orderAsc(p); let prev = -1, mono = true;
  for (const i of idx) { if (adj[k][i] < prev - 1e-12) mono = false; prev = adj[k][i]; }
  if (!mono) { fails++; console.log('  FAIL monotonic', k); }
}
// bounds: every adjusted value in [0,1] and >= raw p
for (const k of ['bh','by','holm','bonferroni']) {
  adj[k].forEach((v,i)=>{ if(!(v>=0&&v<=1)){fails++;console.log('  FAIL range',k,v);} if(v < p[i]-1e-12){fails++;console.log('  FAIL below-raw',k,v,p[i]);} });
}
// BY is always >= BH (more conservative)
adj.by.forEach((v,i)=>{ if(v < adj.bh[i]-1e-12){fails++;console.log('  FAIL BY<BH',v,adj.bh[i]);} });
// harmonic
eq('harmonic(1)',F.harmonic(1),1); eq('harmonic(4)',F.harmonic(4),1+1/2+1/3+1/4);
// m=1 edge
const one = F.adjustAll([0.02]);
eq('m=1 bh',one.bh[0],0.02,1e-12); eq('m=1 by',one.by[0],0.02,1e-12);
// invNorm/normCdf round trip
eq('normCdf(0)',F.normCdf(0),0.5,1e-7);
eq('invNorm(0.975)',F.invNorm(0.975),1.959964,1e-4);
console.log(fails===0 ? `REGRESSION: all checks passed (${names.length} functions extracted)` : `REGRESSION: ${fails} FAILURES`);
process.exit(fails?1:0);

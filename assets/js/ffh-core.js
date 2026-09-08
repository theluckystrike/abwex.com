/* Local, bounded Fisher–Freeman–Halton enumeration, probability-ordered ties. */
'use strict';
{
  const LIMITS=Object.freeze({maxTotal:200,maxNodes:250000,maxTables:100000});
  function invalid(error){return {status:'invalid',complete:false,error};}
  function validateShape(input){
    if(!Array.isArray(input)||input.length<2||input.length>4) return invalid('Use 2 to 4 rows.');
    const C=Array.isArray(input[0])?input[0].length:0;
    if(C<2||C>4||Array.from(input).some(r=>!Array.isArray(r)||r.length!==C))
      return invalid('Use a rectangular table with 2 to 4 columns.');
    if(input.some(r=>Array.from(r).some(x=>!Number.isSafeInteger(x)||x<0||x>200)))
      return invalid('Counts must be nonnegative integers at most 200.');
    return {R:input.length,C};
  }
  function prepareCaps(options){
    const caps={};
    for(const key of ['maxNodes','maxTables']){
      const v=options[key]===undefined?LIMITS[key]:options[key];
      if(!Number.isSafeInteger(v)||v<1||v>LIMITS[key])
        return invalid(key+' must be a positive integer no greater than '+LIMITS[key]+'.');
      caps[key]=v;
    }
    return caps;
  }
  function prepare(input,options){
    const shape=validateShape(input);
    if(shape.status==='invalid') return shape;
    const {R,C}=shape;
    const rows=input.map(r=>r.reduce((a,b)=>a+b,0));
    const cols=Array.from({length:C},(_,j)=>input.reduce((s,r)=>s+r[j],0));
    const N=rows.reduce((a,b)=>a+b,0);
    if(!N||N>200) return invalid('Grand total must be from 1 to 200.');
    if(rows.some(x=>!x)||cols.some(x=>!x))
      return invalid('Every row and column needs a positive margin; no categories were removed.');
    const caps=prepareCaps(options);
    return caps.status==='invalid'?caps:{R,C,N,rows,cols,caps};
  }
  function initialize(input,p){
    const {R,C,N,rows,cols,caps}=p,fact=[1n],logFact=[0];
    for(let i=1;i<=N;i++){fact[i]=fact[i-1]*BigInt(i);logFact[i]=logFact[i-1]+Math.log(i);}
    const constant=rows.reduce((s,x)=>s+logFact[x],0)+cols.reduce((s,x)=>s+logFact[x],0)-logFact[N];
    let observedDenominator=1n,observedLog=0,chiSquare=0;
    for(const row of input) for(const x of row){observedDenominator*=fact[x];observedLog+=logFact[x];}
    const expected=rows.map(r=>cols.map(c=>r*c/N));
    for(let i=0;i<R;i++) for(let j=0;j<C;j++) chiSquare+=(input[i][j]-expected[i][j])**2/expected[i][j];
    const base={method:'Fixed-margins exact enumeration; probability ordering with exact factorial-product ties',
      total:N,rowMargins:rows,columnMargins:cols,observedProbability:Math.exp(constant-observedLog),
      expected,chiSquare,df:(R-1)*(C-1),chiSquarePValue:chiSquareSurvival(chiSquare,(R-1)*(C-1)),limits:caps};
    return {...p,base,fact,logFact,constant,observedDenominator,rr:rows.slice(),cc:cols.slice(),
      cells:Array(R*C).fill(0),stack:[],pos:0,nodes:0,tableCount:0,tailCount:0,
      mass:0,massComp:0,tailMass:0,tailComp:0,tailExamples:[]};
  }
  function validProbabilities(s){
    return Number.isFinite(s.mass)&&Math.abs(s.mass-1)<=1e-10&&
      Number.isFinite(s.tailMass)&&s.tailMass>=0&&s.tailMass<=1+1e-10&&
      Number.isFinite(s.base.observedProbability)&&s.base.observedProbability>0&&s.base.observedProbability<=1+1e-10;
  }
  function finish(s,complete,reason){
    const out={...s.base,status:complete?'complete':'incomplete',complete,nodeCount:s.nodes,
      tableCount:s.tableCount,tailCount:s.tailCount,massSum:s.mass,tailExamples:s.tailExamples,
      examplesDescription:'First at most 12 included tables in deterministic lexicographic enumeration; not the full tail.'};
    if(!complete){out.error=reason+' No exact p-value is reported because enumeration is incomplete.';return out;}
    const valid=validProbabilities(s);
    if(!valid){
      out.status='numerical-error';out.complete=false;out.enumerationComplete=true;
      out.error='Enumeration finished, but the probability-mass or numeric-range check failed. No p-value is reported.';
    }else out.pValue=Math.min(1,s.tailMass);
    return out;
  }
  function addProbability(s,value,tail){
    const key=tail?'tailMass':'mass',comp=tail?'tailComp':'massComp';
    const y=value-s[comp],t=s[key]+y;s[comp]=(t-s[key])-y;s[key]=t;
  }
  function visitTable(s){
    s.tableCount++;
    let denominator=1n,logDenominator=0;
    for(const x of s.cells){denominator*=s.fact[x];logDenominator+=s.logFact[x];}
    const probability=Math.exp(s.constant-logDenominator);addProbability(s,probability,false);
    // All probabilities have the same numerator; a greater denominator means a lower probability.
    if(denominator>=s.observedDenominator){
      s.tailCount++;addProbability(s,probability,true);
      if(s.tailExamples.length<12) s.tailExamples.push({
        table:Array.from({length:s.R},(_,i)=>s.cells.slice(i*s.C,(i+1)*s.C)),
        probability,tie:denominator===s.observedDenominator});
    }
  }
  function nextFrame(s){
    const i=Math.floor(s.pos/s.C),j=s.pos%s.C;
    let low=0,high=Math.min(s.rr[i],s.cc[j]);
    if(i===s.R-1) low=high=s.cc[j];
    if(j===s.C-1){if(i===s.R-1&&s.rr[i]!==s.cc[j]) high=-1;else low=high=s.rr[i];}
    if(low>s.rr[i]||low>s.cc[j]) high=-1;
    return {i,j,next:low,high,assigned:null};
  }
  function enumerate(s){
    // Each depth holds one cell's next candidate; no recursion or unbounded allocations.
    while(s.pos>=0){
      if(s.pos===s.R*s.C){
        if(s.tableCount>=s.caps.maxTables) return finish(s,false,'Table budget exceeded.');
        visitTable(s);s.pos--;
      }else{
        const frame=s.stack[s.pos]||(s.stack[s.pos]=nextFrame(s));
        if(frame.assigned!==null){
          s.rr[frame.i]+=frame.assigned;s.cc[frame.j]+=frame.assigned;frame.assigned=null;
        }
        if(frame.next>frame.high){
          s.stack.pop();s.pos--;
        }else{
          if(s.nodes>=s.caps.maxNodes) return finish(s,false,'Node budget exceeded.');
          s.nodes++;
          const value=frame.next++;frame.assigned=value;
          s.rr[frame.i]-=value;s.cc[frame.j]-=value;s.cells[s.pos]=value;s.pos++;
        }
      }
    }
    return finish(s,true);
  }
  function gammaSeries(a,x,scale){
    let term=1/a,sum=term;
    for(let n=1;n<=10000;n++){
      term*=x/(a+n);sum+=term;
      if(Math.abs(term)<Math.abs(sum)*1e-15) return Math.max(0,Math.min(1,1-scale*sum));
    }
    return null;
  }
  function gammaFraction(a,x,scale){
    const tiny=1e-300;
    let b=x+1-a,c=1/tiny,d=1/b,h=d;
    for(let n=1;n<=10000;n++){
      const an=-n*(n-a);b+=2;d=an*d+b;
      if(Math.abs(d)<tiny) d=tiny;
      c=b+an/c;if(Math.abs(c)<tiny) c=tiny;
      d=1/d;const delta=d*c;h*=delta;
      if(Math.abs(delta-1)<1e-15) return Math.max(0,Math.min(1,scale*h));
    }
    return null;
  }
  function halfIntegerLogGamma(a){
    let result=a%1===0?0:Math.log(Math.PI)/2;
    for(let v=a%1===0?1:0.5;v<a;v++) result+=Math.log(v);
    return result;
  }
  function chiSquareSurvival(value,df){
    if(!Number.isFinite(value)||value<0||!Number.isInteger(df)||df<1||df>9) return null;
    if(value===0) return 1;
    const a=df/2,x=value/2;
    const scale=Math.exp(a*Math.log(x)-x-halfIntegerLogGamma(a));
    return x<a+1?gammaSeries(a,x,scale):gammaFraction(a,x,scale);
  }
  function calculate(input,options){
    const prepared=prepare(input,options||{});
    return prepared.status==='invalid'?prepared:enumerate(initialize(input,prepared));
  }
  const api={calculate,chiSquareSurvival,limits:LIMITS};
  if(typeof module==='object'&&module.exports) module.exports=api;
  else globalThis.FFH=api;
}

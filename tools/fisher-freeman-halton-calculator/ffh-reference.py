"""Exact rational fixed-margin reference; independently checked with SciPy PMFs.
Run: python3 ffh-reference.py > ffh-fixtures.json
Not a production implementation. SciPy >=1.10 required only for cross-checks.
"""
import json
from fractions import Fraction
from math import factorial,prod,fsum,isclose
import numpy as np
import scipy
from scipy.stats import random_table, fisher_exact

def enumerate_tables(rows, cols):
    if len(rows)==1:
        if sum(cols)==rows[0]: yield [cols[:]]
        return
    def row_options(j,left,row):
        if j==len(cols)-1:
            if 0<=left<=cols[j]: yield row+[left]
            return
        lo=max(0,left-sum(cols[j+1:])); hi=min(left,cols[j])
        for v in range(lo,hi+1):
            yield from row_options(j+1,left-v,row+[v])
    for row in row_options(0,rows[0],[]):
        for rest in enumerate_tables(rows[1:],[c-v for c,v in zip(cols,row)]):
            yield [row]+rest

def fixture(name, observed):
    rows=list(map(sum,observed)); cols=list(map(sum,zip(*observed))); total=sum(rows)
    numerator=prod(factorial(v) for v in rows+cols)
    def probability(table):
        return Fraction(numerator,factorial(total)*prod(factorial(v) for row in table for v in row))
    obs=probability(observed); tables=list(enumerate_tables(rows,cols)); probs=[probability(a) for a in tables]
    selected=[p for p in probs if p<=obs]; pvalue=sum(selected,Fraction())
    assert sum(probs,Fraction())==1
    scipy_probs=random_table.pmf(np.array(tables), rows,cols)
    assert all(isclose(float(a),float(b),rel_tol=2e-12,abs_tol=1e-15) for a,b in zip(probs,scipy_probs))
    scipy_p=fsum(float(b) for a,b in zip(probs,scipy_probs) if a<=obs)
    assert isclose(float(pvalue),scipy_p,rel_tol=2e-12,abs_tol=1e-15)
    if len(rows)==len(cols)==2:
        assert isclose(fisher_exact(observed).pvalue,float(pvalue),rel_tol=2e-12)
    return dict(name=name,table=observed,row_margins=rows,column_margins=cols,total=total,
        feasible_tables=len(tables),extreme_tables=len(selected),observed_probability=float(obs),
        observed_probability_fraction=str(obs),p_value=float(pvalue),p_value_fraction=str(pvalue),
        scipy_pmf_sum_p_value=scipy_p,probability_mass_sum_exact=str(sum(probs,Fraction())))

cases=[('2x3_sparse',[[1,4,0],[3,0,2]]),('2x3_moderate',[[8,2,1],[1,5,4]]),
       ('2x3_balanced',[[2,2,2],[2,2,2]]),('3x3_diagonal',[[3,0,0],[0,3,0],[0,0,3]]),
       ('3x3_mixed',[[2,1,0],[0,2,1],[1,0,2]]),('3x3_asymmetric',[[4,1,0],[1,2,1],[0,1,2]]),
       ('2x2_compatibility',[[1,9],[11,3]])]
print(json.dumps(dict(method='exact rational complete enumeration, probability ordering, fixed margins',scipy_version=scipy.__version__,fixtures=[fixture(*a) for a in cases]),indent=2))

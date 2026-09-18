"""Cross check of the Cochran-Armitage trend statistic and exact p values on synthetic fixtures.

For each synthetic table this script computes, independently:
  - the Armitage Z with variance using N (the form behind R prop.trend.test) and its p value,
  - the conditional Z with variance using N-1 and its p value,
  - statsmodels Table.test_ordinal_association zscore and pvalue for the same table,
  - the exact conditional p value and mid p value by brute force enumeration of every table
    with the same group sizes and the same total conversions, weighted by prod C(n_i, x_i).
It fails loudly if the N-1 form disagrees with statsmodels beyond 1e-10, or if the
enumerated conditional probabilities do not sum to 1. The fixtures are synthetic
illustrations, never observed customer data.

Run: python3 verify_ca_against_statsmodels.py   (writes ca_trend_fixtures.json beside this file)
"""
import json
import math
import os
import platform
import sys
from itertools import product

import numpy as np
import scipy
import statsmodels
from scipy.stats import norm
from statsmodels.stats.contingency_tables import Table

FIXTURES = [
    {"id": "three_arms_equal", "n": [200, 200, 200], "x": [6, 10, 17]},
    {"id": "four_arms_sparse", "n": [50, 50, 50, 50], "x": [1, 2, 4, 6]},
    {"id": "three_arms_very_sparse", "n": [30, 30, 30], "x": [0, 1, 4]},
    {"id": "five_arms_flat_then_rise", "n": [100, 100, 100, 100, 100], "x": [3, 3, 4, 5, 8]},
    {"id": "three_arms_unequal", "n": [120, 80, 60], "x": [3, 4, 6]},
]


def armitage(n, x, scores):
    n, x, t = np.array(n, float), np.array(x, float), np.array(scores, float)
    big_n, r = n.sum(), x.sum()
    pbar = r / big_n
    stat_t = (t * x).sum() - pbar * (t * n).sum()
    var_n = pbar * (1 - pbar) * ((n * t * t).sum() - (n * t).sum() ** 2 / big_n)
    var_n1 = var_n * big_n / (big_n - 1)
    z_n, z_n1 = stat_t / math.sqrt(var_n), stat_t / math.sqrt(var_n1)
    return {"T": stat_t, "z_n": z_n, "chi2_n": z_n ** 2, "p_asym_n": 2 * norm.sf(abs(z_n)),
            "z_n1": z_n1, "p_asym_n1": 2 * norm.sf(abs(z_n1))}


def exact_conditional(n, x, scores):
    r = sum(x)
    ranges = [range(0, min(ni, r) + 1) for ni in n]
    big_n = sum(n)
    mean_s = r * sum(t * ni for t, ni in zip(scores, n)) / big_n
    obs = abs(sum(t * xi for t, xi in zip(scores, x)) - mean_s)
    tol = 1e-9 * max(1.0, obs)
    logw, dists = [], []
    for xs in product(*ranges):
        if sum(xs) != r:
            continue
        logw.append(sum(math.lgamma(ni + 1) - math.lgamma(xi + 1) - math.lgamma(ni - xi + 1) for ni, xi in zip(n, xs)))
        dists.append(abs(sum(t * xi for t, xi in zip(scores, xs)) - mean_s))
    logw = np.array(logw)
    w = np.exp(logw - logw.max())
    w /= w.sum()
    dists = np.array(dists)
    p_ge = float(w[dists >= obs - tol].sum())
    p_eq = float(w[np.abs(dists - obs) <= tol].sum())
    return {"tables_enumerated": int(len(w)), "p_exact": p_ge, "p_midp": p_ge - 0.5 * p_eq,
            "mass": float(w.sum())}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = []
    worst = 0.0
    for f in FIXTURES:
        k = len(f["n"])
        scores = list(range(k))
        a = armitage(f["n"], f["x"], scores)
        tab = np.array([[ni - xi, xi] for ni, xi in zip(f["n"], f["x"])], float)
        sm = Table(tab, shift_zeros=False).test_ordinal_association(row_scores=np.array(scores, float), col_scores=np.array([0.0, 1.0]))
        diff = max(abs(abs(a["z_n1"]) - abs(sm.zscore)), abs(a["p_asym_n1"] - sm.pvalue))
        worst = max(worst, diff)
        e = exact_conditional(f["n"], f["x"], scores)
        if abs(e["mass"] - 1.0) > 1e-12:
            sys.exit(f"conditional mass {e['mass']} is not 1 for {f['id']}")
        out.append({
            **f, "scores": scores,
            "rates": [round(xi / ni, 6) for ni, xi in zip(f["n"], f["x"])],
            "z_armitage_n": round(a["z_n"], 6), "chi2_armitage_n": round(a["chi2_n"], 6),
            "p_asymptotic_n": round(a["p_asym_n"], 6),
            "z_conditional_n1": round(a["z_n1"], 6), "p_asymptotic_n1": round(a["p_asym_n1"], 6),
            "statsmodels_zscore": round(float(sm.zscore), 6), "statsmodels_pvalue": round(float(sm.pvalue), 6),
            "p_exact_conditional": round(e["p_exact"], 6), "p_midp": round(e["p_midp"], 6),
            "tables_enumerated": e["tables_enumerated"],
        })
    if worst > 1e-10:
        sys.exit(f"N-1 form disagrees with statsmodels by {worst}")
    doc = {
        "id": "ca_trend_fixtures",
        "title": "Cochran-Armitage trend test worked examples on synthetic tables",
        "kind": "synthetic illustration tables, computed, not observed data",
        "check": {"reference": "statsmodels Table.test_ordinal_association", "max_abs_diff": worst, "tolerance": 1e-10,
                  "note": "R was not available on the build machine, the N form is the textbook Armitage variance. statsmodels Table is built with shift_zeros=False, its default adds 0.5 to every cell of a table containing a zero"},
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                        "statsmodels": statsmodels.__version__},
        "generator": "datasets/verify_ca_against_statsmodels.py",
        "fixtures": out,
    }
    with open(os.path.join(here, "ca_trend_fixtures.json"), "w") as fh:
        json.dump(doc, fh, indent=1)
    print(f"wrote {len(out)} fixtures, max |diff| vs statsmodels {worst:.3e}")
    for o in out:
        print(o["id"], o["x"], "z_n", o["z_armitage_n"], "p_n", o["p_asymptotic_n"], "p_n1", o["p_asymptotic_n1"],
              "sm_p", o["statsmodels_pvalue"], "exact", o["p_exact_conditional"], "midp", o["p_midp"], "tables", o["tables_enumerated"])


if __name__ == "__main__":
    main()

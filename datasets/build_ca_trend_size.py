"""Exact size study of the Cochran-Armitage test for trend in proportions.

Question answered: when every group converts at the same true rate p (the null is true),
how often does each version of the Cochran-Armitage trend test reject at nominal two sided
alpha 0.05? A correct test rejects at most 5 percent of the time.

Method, exact enumeration and no simulation:
  - k ordered groups, equally spaced scores 0..k-1, n visitors per group.
  - Each group's conversions are independent Binomial(n, p).
  - The joint distribution of R = total conversions and S = sum of score times conversions
    is built exactly by dynamic programming over the groups. The rejection decision of every
    test version depends only on (R, S), so the size is an exact finite sum.
  - A brute force enumeration over every outcome vector checks the dynamic programme on
    small cells before anything is written.

Test versions:
  asym_n    Armitage statistic, variance pbar(1-pbar)[sum n t^2 - (sum n t)^2 / N].
            Z squared equals the chi-square reported by R prop.trend.test.
  asym_n1   the same statistic with the conditional variance, larger by N/(N-1). This is the
            linear by linear association form reported by statsmodels
            Table.test_ordinal_association.
  exact     exact conditional (permutation) test given R, two sided by distance from the
            conditional mean, ties included.
  midp      mid p version of the exact test, half the probability of the observed |T|.

Outcomes with R = 0 or R = N have zero variance, no test is defined, and they count as
non rejections. Their probability is reported per cell as p_undefined.

Run: python3 build_ca_trend_size.py   (writes ca_trend_size.json beside this file)
"""
import json
import math
import os
import platform
import sys
from itertools import product

import numpy as np
import scipy
from scipy.stats import binom, norm

ALPHA = 0.05
K_VALUES = (3, 4, 5)
N_VALUES = (10, 20, 50, 100)
P_VALUES = (0.01, 0.02, 0.05, 0.10, 0.20, 0.50)
REL_TOL = 1e-9


def joint_rs(k, n, p):
    """Exact joint pmf of (R, S) for k independent Binomial(n, p) groups, scores 0..k-1."""
    r_max = k * n
    s_max = n * (k - 1) * k // 2
    pmf = binom.pmf(np.arange(n + 1), n, p)
    joint = np.zeros((r_max + 1, s_max + 1))
    joint[0, 0] = 1.0
    r_hi, s_hi = 0, 0
    for t in range(k):
        new = np.zeros_like(joint)
        for x in range(n + 1):
            if pmf[x] == 0.0:
                continue
            new[x:x + r_hi + 1, t * x:t * x + s_hi + 1] += pmf[x] * joint[:r_hi + 1, :s_hi + 1]
        joint = new
        r_hi += n
        s_hi += t * n
    return joint


def cond_weights(k, n):
    """Exact conditional pmf of S given R under the null, from prod C(n, x_i) weights."""
    r_max = k * n
    s_max = n * (k - 1) * k // 2
    logc = np.array([math.lgamma(n + 1) - math.lgamma(x + 1) - math.lgamma(n - x + 1) for x in range(n + 1)])
    w = np.full((r_max + 1, s_max + 1), -np.inf)
    w[0, 0] = 0.0
    r_hi, s_hi = 0, 0
    for t in range(k):
        new = np.full_like(w, -np.inf)
        for x in range(n + 1):
            blk = w[:r_hi + 1, :s_hi + 1] + logc[x]
            tgt = new[x:x + r_hi + 1, t * x:t * x + s_hi + 1]
            new[x:x + r_hi + 1, t * x:t * x + s_hi + 1] = np.logaddexp(tgt, blk)
        w = new
        r_hi += n
        s_hi += t * n
    out = np.zeros_like(w)
    for r in range(r_max + 1):
        row = w[r]
        m = row.max()
        if m == -np.inf:
            continue
        e = np.exp(row - m)
        out[r] = e / e.sum()
    return out


def decisions(k, n):
    """Rejection indicators for every (R, S) cell and every test version."""
    scores = np.arange(k, dtype=float)
    big_n = k * n
    sum_t = n * scores.sum()
    sum_t2 = n * (scores ** 2).sum()
    r_max, s_max = big_n, int(n * (k - 1) * k // 2)
    cond = cond_weights(k, n)
    rej = {v: np.zeros((r_max + 1, s_max + 1), dtype=bool) for v in ("asym_n", "asym_n1", "exact", "midp")}
    undefined = np.zeros((r_max + 1, s_max + 1), dtype=bool)
    s_axis = np.arange(s_max + 1, dtype=float)
    z_crit = norm.isf(ALPHA / 2)
    for r in range(r_max + 1):
        if r == 0 or r == big_n:
            undefined[r, :] = True
            continue
        pbar = r / big_n
        t_stat = s_axis - r * sum_t / big_n
        var_n = pbar * (1 - pbar) * (sum_t2 - sum_t ** 2 / big_n)
        var_n1 = var_n * big_n / (big_n - 1)
        rej["asym_n"][r] = np.abs(t_stat) / math.sqrt(var_n) > z_crit
        rej["asym_n1"][r] = np.abs(t_stat) / math.sqrt(var_n1) > z_crit
        dist = np.abs(t_stat)
        prob = cond[r]
        support = prob > 0
        d_s, p_s = dist[support], prob[support]
        order = np.argsort(-d_s, kind="stable")
        d_sorted, p_sorted = d_s[order], p_s[order]
        cum = np.cumsum(p_sorted)
        # p value for distance d: P(D >= d) with a relative tie tolerance on the distance
        tol = 1e-9 * max(1.0, float(d_sorted.max()))
        # -d_sorted is ascending. D >= d - tol  <=>  -D <= -d + tol ; D > d + tol  <=>  -D < -d - tol
        idx_ge = np.searchsorted(-d_sorted, -dist + tol, side="right") - 1
        idx_gt = np.searchsorted(-d_sorted, -dist - tol, side="left") - 1
        p_ge = np.where(idx_ge >= 0, cum[np.clip(idx_ge, 0, None)], 0.0)
        p_gt = np.where(idx_gt >= 0, cum[np.clip(idx_gt, 0, None)], 0.0)
        rej["exact"][r] = p_ge <= ALPHA
        rej["midp"][r] = (p_gt + 0.5 * (p_ge - p_gt)) <= ALPHA
    return rej, undefined


def brute_force_size(k, n, p):
    """Direct enumeration of every outcome vector, for checking the dynamic programme."""
    scores = list(range(k))
    big_n = k * n
    sum_t = n * sum(scores)
    sum_t2 = n * sum(t * t for t in scores)
    z_crit = norm.isf(ALPHA / 2)
    pmf = binom.pmf(np.arange(n + 1), n, p)
    size_n = size_n1 = undef = total = 0.0
    for xs in product(range(n + 1), repeat=k):
        pr = float(np.prod([pmf[x] for x in xs]))
        total += pr
        r = sum(xs)
        if r == 0 or r == big_n:
            undef += pr
            continue
        pbar = r / big_n
        t_stat = sum(t * x for t, x in zip(scores, xs)) - r * sum_t / big_n
        var_n = pbar * (1 - pbar) * (sum_t2 - sum_t ** 2 / big_n)
        if abs(t_stat) / math.sqrt(var_n) > z_crit:
            size_n += pr
        if abs(t_stat) / math.sqrt(var_n * big_n / (big_n - 1)) > z_crit:
            size_n1 += pr
    return {"asym_n": size_n, "asym_n1": size_n1, "p_undefined": undef, "total": total}


def cell(k, n, p, rej, undefined):
    joint = joint_rs(k, n, p)
    mass = float(joint.sum())
    out = {"k": k, "n_per_group": n, "p": p, "total_visitors": k * n}
    for v in ("asym_n", "asym_n1", "exact", "midp"):
        out[f"size_{v}"] = round(float(joint[rej[v]].sum()), 6)
    out["p_undefined"] = round(float(joint[undefined].sum()), 6)
    out["probability_mass"] = mass
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    checks = []
    for (k, n, p) in ((3, 10, 0.10), (3, 10, 0.02), (4, 8, 0.05), (5, 6, 0.20)):
        rej, undefined = decisions(k, n)
        joint = joint_rs(k, n, p)
        bf = brute_force_size(k, n, p)
        dp = {"asym_n": float(joint[rej["asym_n"]].sum()), "asym_n1": float(joint[rej["asym_n1"]].sum()),
              "p_undefined": float(joint[undefined].sum()), "total": float(joint.sum())}
        worst = max(abs(dp[key] - bf[key]) for key in dp)
        checks.append({"k": k, "n_per_group": n, "p": p, "max_abs_diff_dp_vs_brute_force": worst})
        if worst > 1e-12:
            sys.exit(f"dynamic programme disagrees with brute force at k={k} n={n} p={p}: {worst}")
    cells = []
    for k in K_VALUES:
        for n in N_VALUES:
            rej, undefined = decisions(k, n)
            for p in P_VALUES:
                c = cell(k, n, p, rej, undefined)
                if abs(c["probability_mass"] - 1.0) > REL_TOL:
                    sys.exit(f"probability mass {c['probability_mass']} is not 1 at k={k} n={n} p={p}")
                c["probability_mass"] = round(c["probability_mass"], 12)
                cells.append(c)
    doc = {
        "id": "ca_trend_size",
        "title": "Exact type I error of the Cochran-Armitage trend test at nominal alpha 0.05",
        "kind": "computed by exact enumeration, not simulated, not observed data",
        "alpha": ALPHA,
        "two_sided": True,
        "scores": "equally spaced 0 to k-1",
        "grid": {"k": list(K_VALUES), "n_per_group": list(N_VALUES), "p": list(P_VALUES)},
        "versions": {
            "asym_n": "asymptotic, Armitage variance with N, matches R prop.trend.test chi-square",
            "asym_n1": "asymptotic, conditional variance with N-1, matches statsmodels test_ordinal_association",
            "exact": "exact conditional permutation test given total conversions, ties included",
            "midp": "mid p version of the exact conditional test",
        },
        "undefined_rule": "R = 0 or R = N has zero variance and counts as no rejection, reported as p_undefined",
        "checks": {"dp_vs_brute_force": checks, "probability_mass_tolerance": REL_TOL},
        "generator": "datasets/build_ca_trend_size.py",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "cells": cells,
    }
    with open(os.path.join(here, "ca_trend_size.json"), "w") as fh:
        json.dump(doc, fh, indent=1)
    print(f"wrote {len(cells)} cells, dp checks {checks}")


if __name__ == "__main__":
    main()

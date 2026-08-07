#!/usr/bin/env python3
"""
Known-answer test for the correction routines used in build_fdr_simulation.py.

Checks this site's Benjamini-Hochberg, Holm and Bonferroni reject decisions against
statsmodels.stats.multitest.multipletests, which is the reference implementation the page
cites. Random p-value vectors, many shapes, exact agreement required.

Run: python3 verify_against_statsmodels.py
"""
import numpy as np
from statsmodels.stats.multitest import multipletests
from build_fdr_simulation import bh_reject, holm_reject, bonferroni_reject

ALPHA = 0.05
rng = np.random.default_rng(7)

fails = 0
checked = 0
for m in (2, 3, 5, 10, 20, 50):
    # mix of uniform p-values and deliberately small ones so rejections actually happen
    for trial in range(200):
        p = rng.uniform(size=m)
        if trial % 3 == 0:
            p[: max(1, m // 4)] *= 0.01
        row = p[None, :]

        mine = {
            "bh": bh_reject(row, ALPHA)[0],
            "holm": holm_reject(row, ALPHA)[0],
            "bonferroni": bonferroni_reject(row, ALPHA)[0],
        }
        ref = {
            "bh": multipletests(p, alpha=ALPHA, method="fdr_bh")[0],
            "holm": multipletests(p, alpha=ALPHA, method="holm")[0],
            "bonferroni": multipletests(p, alpha=ALPHA, method="bonferroni")[0],
        }
        for name in mine:
            checked += 1
            if not np.array_equal(mine[name], ref[name]):
                fails += 1
                print(f"MISMATCH {name} m={m} trial={trial}")
                print("  p    =", np.round(p, 5))
                print("  mine =", mine[name].astype(int))
                print("  ref  =", ref[name].astype(int))

print(f"\nchecked {checked} decisions across m in 2,3,5,10,20,50")
print("RESULT:", "all agree with statsmodels exactly" if fails == 0 else f"{fails} MISMATCHES")


# --- adjusted p-VALUE comparison -------------------------------------------------
# The page states that adjusted values agree with multipletests to within a tolerance, so the
# tolerance has to be measured, not asserted. Recompute the adjusted values the same way the
# page's calculator does and take the largest absolute difference against statsmodels.

def bh_adjusted(p):
    m = p.size
    order = np.argsort(p)
    ps = p[order]
    adj = ps * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]      # enforce monotonicity
    adj = np.clip(adj, 0, 1)
    out = np.empty_like(adj)
    out[order] = adj
    return out


def holm_adjusted(p):
    m = p.size
    order = np.argsort(p)
    ps = p[order]
    adj = ps * (m - np.arange(m))
    adj = np.maximum.accumulate(adj)                  # enforce monotonicity
    adj = np.clip(adj, 0, 1)
    out = np.empty_like(adj)
    out[order] = adj
    return out


def bonferroni_adjusted(p):
    return np.clip(p * p.size, 0, 1)


worst = {"bh": 0.0, "holm": 0.0, "bonferroni": 0.0}
n_vec = 0
for m in (2, 3, 5, 10, 20, 50):
    for trial in range(200):
        p = rng.uniform(size=m)
        if trial % 3 == 0:
            p[: max(1, m // 4)] *= 0.01
        n_vec += 1
        for name, fn, meth in (
            ("bh", bh_adjusted, "fdr_bh"),
            ("holm", holm_adjusted, "holm"),
            ("bonferroni", bonferroni_adjusted, "bonferroni"),
        ):
            ref = multipletests(p, alpha=ALPHA, method=meth)[1]
            worst[name] = max(worst[name], float(np.max(np.abs(fn(p) - ref))))

print(f"\nadjusted p-value agreement over {n_vec} vectors, max absolute difference:")
for k, v in worst.items():
    print(f"  {k:11s} {v:.3e}")
print(f"OVERALL max abs difference: {max(worst.values()):.3e}")

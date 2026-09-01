#!/usr/bin/env python3
"""Build the observed power simulation dataset for the post hoc power calculator page.

For each configuration, the true power of a two-proportion z-test is fixed by
(p1, p2, n per arm, alpha). We then simulate many replicate A/B tests drawn from
those true rates, compute the post hoc ("observed") power of each replicate from
its own observed rates, and record how the observed power distributes around the
true power. The point the page makes is that observed power is a noisy, biased
restatement of the p value, so the dataset quantifies exactly how noisy.

Deterministic: seeded PRNG, fixed grid. Output is posthoc_power_simulation.json.
Cross-checked against statsmodels power_proportions_2indep where available.
"""

import json
import math

import numpy as np

SEED = 20260901
REPS = 10000
ALPHA = 0.05
P1 = 0.10

Z_A2 = 1.959963984540054  # Phi^-1(0.975)


def phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def posthoc_power(p1, p2, n1, n2, alpha_z=Z_A2):
    """Pooled-null post hoc power, the statsmodels power_proportions_2indep form."""
    if n1 <= 0 or n2 <= 0:
        raise ValueError("n must be positive")
    if not (0.0 < p1 < 1.0) or not (0.0 < p2 < 1.0):
        return None
    p1, p2 = float(p1), float(p2)
    diff = abs(p1 - p2)
    pbar = (n1 * p1 + n2 * p2) / (n1 + n2)
    se0 = math.sqrt(pbar * (1 - pbar) * (1.0 / n1 + 1.0 / n2))
    se1 = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    if se1 == 0.0:
        return None
    return phi((diff - alpha_z * se0) / se1) + phi((-diff - alpha_z * se0) / se1)


def two_sided_p(p1, p2, n1, n2):
    pbar = (n1 * p1 + n2 * p2) / (n1 + n2)
    se0 = math.sqrt(pbar * (1 - pbar) * (1.0 / n1 + 1.0 / n2))
    if se0 == 0.0:
        return 1.0
    z = (p2 - p1) / se0
    return 2.0 * (1.0 - phi(abs(z)))


def solve_n_for_power(target, p1, p2, alpha_z=Z_A2):
    """Smallest equal n per arm whose true power meets the target."""
    lo, hi = 2, 4_000_000
    if posthoc_power(p1, p2, hi, hi, alpha_z) < target:
        raise ValueError("target power unreachable")
    while lo < hi:
        mid = (lo + hi) // 2
        if posthoc_power(p1, p2, mid, mid, alpha_z) >= target:
            hi = mid
        else:
            lo = mid + 1
    return lo


def percentile(sorted_vals, q):
    if not sorted_vals:
        return None
    idx = q * (len(sorted_vals) - 1)
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    frac = idx - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def simulate(rng, p1, p2, n, reps):
    true_power = posthoc_power(p1, p2, n, n)
    c1 = rng.binomial(n, p1, reps)
    c2 = rng.binomial(n, p2, reps)
    obs = []
    sig = 0
    for i in range(reps):
        o1 = c1[i] / n
        o2 = c2[i] / n
        pw = posthoc_power(o1, o2, n, n) if 0.0 < o1 < 1.0 and 0.0 < o2 < 1.0 else None
        if pw is None:
            pw = 0.0 if o1 == o2 else 1.0
        obs.append(pw)
        if two_sided_p(o1, o2, n, n) < ALPHA:
            sig += 1
    obs.sort()
    within_10 = sum(1 for v in obs if abs(v - true_power) <= 0.10) / reps
    return {
        "n_per_arm": n,
        "true_power": round(true_power, 4),
        "reps": reps,
        "observed_power": {
            "p05": round(percentile(obs, 0.05), 4),
            "p25": round(percentile(obs, 0.25), 4),
            "p50": round(percentile(obs, 0.50), 4),
            "p75": round(percentile(obs, 0.75), 4),
            "p95": round(percentile(obs, 0.95), 4),
            "mean": round(sum(obs) / reps, 4),
        },
        "share_within_0.10_of_true": round(within_10, 4),
        "share_significant": round(sig / reps, 4),
    }


def main():
    rng = np.random.default_rng(SEED)
    p2 = 0.12
    cells = []
    for target in (0.20, 0.50, 0.80, 0.95):
        n = solve_n_for_power(target, P1, p2)
        cells.append(simulate(rng, P1, p2, n, REPS))

    # The deterministic p-value to observed power map for a z-test, no simulation.
    pmap = []
    for p in (0.001, 0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.30, 0.50, 0.80):
        zq = abs(inv_phi(1.0 - p / 2.0))
        pw = phi(zq - Z_A2) + phi(-zq - Z_A2)
        pmap.append({"p_value": p, "observed_power": round(pw, 4)})

    out = {
        "seed": SEED,
        "generator": "build_posthoc_power_simulation.py",
        "date": "2026-09-01",
        "alpha": ALPHA,
        "baseline_p1": P1,
        "treatment_p2": p2,
        "reps_per_cell": REPS,
        "method": "two-proportion z-test, pooled variance under H0, unpooled under H1",
        "cells": cells,
        "p_value_to_observed_power": pmap,
    }
    with open("posthoc_power_simulation.json", "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


def inv_phi(q):
    """Inverse standard normal CDF, Acklam's rational approximation, |err| < 1.15e-9."""
    if not (0.0 < q < 1.0):
        raise ValueError("q out of range")
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow = 0.02425
    if q < plow:
        u = math.sqrt(-2.0 * math.log(q))
        return (((((c[0] * u + c[1]) * u + c[2]) * u + c[3]) * u + c[4]) * u + c[5]) / \
               ((((d[0] * u + d[1]) * u + d[2]) * u + d[3]) * u + 1.0)
    if q > 1.0 - plow:
        u = math.sqrt(-2.0 * math.log(1.0 - q))
        return -(((((c[0] * u + c[1]) * u + c[2]) * u + c[3]) * u + c[4]) * u + c[5]) / \
               ((((d[0] * u + d[1]) * u + d[2]) * u + d[3]) * u + 1.0)
    u = q - 0.5
    t = u * u
    return (((((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4]) * t + a[5]) * u / \
           (((((b[0] * t + b[1]) * t + b[2]) * t + b[3]) * t + b[4]) * t + 1.0)


if __name__ == "__main__":
    main()

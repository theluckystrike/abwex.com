#!/usr/bin/env python3
"""
First-party Monte Carlo study of multiple-testing corrections in A/B tests.

WHY THIS EXISTS
The false-discovery-rate calculators that rank for this query (tools.carbocation.com,
sdmproject.com, multipletesting.com, errickson.net, dreams.ucsc.edu) all adjust a list of
p-values and stop there. None of them answer the two questions an experimenter actually has:

  1. If I track m metrics on one A/B test, how often do I ship a false winner, and how much
     does each correction actually reduce that?
  2. What does the correction cost me in power, so how many real wins do I throw away?

Neither question has a closed form once the metrics are correlated, and A/B test metrics are
strongly correlated (add-to-cart, checkout-start, purchase all move together). So this
simulates it directly.

MODEL
One experiment tracks m metrics. Each metric produces a two-sided z-test of variant vs control.
Under the equicorrelated Gaussian model the m z-statistics are drawn from
    Z ~ Normal(mu, Sigma),  Sigma = (1 - rho) I + rho J
built as  Z_i = sqrt(rho) * C + sqrt(1 - rho) * E_i,  C and E_i iid standard normal.
This is the standard one-factor equicorrelation construction and it is PRDS, the dependence
condition under which Benjamini-Hochberg is proved to control FDR (Benjamini and Yekutieli 2001).

m0 of the m metrics are true nulls (mu = 0). The remaining m1 = m - m0 have a real effect,
parameterised as the per-metric power an UNCORRECTED two-sided test at alpha would have:
    mu = z_{1 - alpha/2} + z_{power}
so "power 0.80" means a single uncorrected metric would detect it 80 percent of the time.

MEASURED, per correction (none, Bonferroni, Holm, Benjamini-Hochberg):
  fwer          share of experiments with at least one false positive among the true nulls
  fdr_mean      mean of V/max(R,1), the realised false discovery proportion
  power_mean    mean share of the m1 true effects that were declared significant
  any_false_win share of experiments where at least one NULL metric was declared a winner
                (this is the practical "did I ship a fake win" number)

Deterministic: fixed seed per cell, so the published table is reproducible by anyone.

Run:  python3 build_fdr_simulation.py
Out:  fdr_simulation.json  (consumed as {{dataset:fdr_simulation}})
"""
import json
import itertools
import datetime as dt
import numpy as np
from scipy.stats import norm

SEED = 20260728
N_SIM = 200_000
ALPHA = 0.05

M_GRID = [3, 5, 10, 20]
RHO_GRID = [0.0, 0.3, 0.6]
PI0_GRID = [1.0, 0.8, 0.5]      # share of metrics that are true nulls
POWER_GRID = [0.80]             # uncorrected per-metric power for the true effects


def holm_reject(p, alpha):
    """Holm step-down. p: (n_sim, m). Returns boolean reject matrix."""
    n, m = p.shape
    order = np.argsort(p, axis=1)
    ps = np.take_along_axis(p, order, axis=1)
    thresh = alpha / (m - np.arange(m))            # alpha/m, alpha/(m-1), ...
    below = ps <= thresh
    # reject the first k where all previous held; stop at first failure
    keep = np.cumprod(below, axis=1).astype(bool)
    rej_sorted = keep
    rej = np.zeros_like(rej_sorted)
    np.put_along_axis(rej, order, rej_sorted, axis=1)
    return rej


def bh_reject(p, alpha):
    """Benjamini-Hochberg step-up. p: (n_sim, m). Returns boolean reject matrix."""
    n, m = p.shape
    order = np.argsort(p, axis=1)
    ps = np.take_along_axis(p, order, axis=1)
    thresh = alpha * (np.arange(1, m + 1) / m)
    below = ps <= thresh
    # largest k with p_(k) <= alpha*k/m ; reject all ranks <= k
    idx = np.where(below.any(axis=1),
                   m - 1 - np.argmax(below[:, ::-1], axis=1),
                   -1)
    ranks = np.arange(m)[None, :]
    rej_sorted = ranks <= idx[:, None]
    rej = np.zeros_like(rej_sorted)
    np.put_along_axis(rej, order, rej_sorted, axis=1)
    return rej


def bonferroni_reject(p, alpha):
    return p <= (alpha / p.shape[1])


def none_reject(p, alpha):
    return p <= alpha


CORRECTIONS = {
    "none": none_reject,
    "bonferroni": bonferroni_reject,
    "holm": holm_reject,
    "bh": bh_reject,
}


def run_cell(m, rho, pi0, power, rng):
    m0 = int(round(m * pi0))
    m0 = min(max(m0, 0), m)
    m1 = m - m0
    is_null = np.zeros(m, dtype=bool)
    is_null[:m0] = True

    mu = np.zeros(m)
    if m1 > 0:
        mu[m0:] = norm.ppf(1 - ALPHA / 2) + norm.ppf(power)

    # equicorrelated one-factor construction
    common = rng.standard_normal((N_SIM, 1))
    idio = rng.standard_normal((N_SIM, m))
    z = np.sqrt(rho) * common + np.sqrt(1 - rho) * idio + mu[None, :]
    p = 2 * norm.sf(np.abs(z))

    out = {}
    for name, fn in CORRECTIONS.items():
        rej = fn(p, ALPHA)
        V = rej[:, is_null].sum(axis=1)                       # false positives
        S = rej[:, ~is_null].sum(axis=1) if m1 > 0 else np.zeros(N_SIM)
        R = V + S
        fdp = np.where(R > 0, V / np.maximum(R, 1), 0.0)
        out[name] = {
            "fwer": round(float((V > 0).mean()), 4),
            "fdr_mean": round(float(fdp.mean()), 4),
            "power_mean": round(float((S / m1).mean()), 4) if m1 > 0 else None,
            "any_false_win": round(float((V > 0).mean()), 4),
        }
    return {
        "m": m, "rho": rho, "pi0": pi0, "m0": m0, "m1": m1,
        "uncorrected_power": power,
        "results": out,
    }


def main():
    rng = np.random.default_rng(SEED)
    cells = []
    for m, rho, pi0, power in itertools.product(M_GRID, RHO_GRID, PI0_GRID, POWER_GRID):
        cells.append(run_cell(m, rho, pi0, power, rng))

    doc = {
        "id": "fdr_simulation",
        "title": "Monte Carlo study of multiple-testing corrections in A/B tests",
        "generated_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d"),
        "method": (
            "Equicorrelated Gaussian z-statistics, Z_i = sqrt(rho)*C + sqrt(1-rho)*E_i with C and E_i "
            "iid standard normal, two-sided p-values, alpha 0.05. m0 true nulls have mu 0; the m1 true "
            "effects have mu = z_(1-alpha/2) + z_(power), so 'uncorrected_power' is the detection rate a "
            "single uncorrected metric would have."
        ),
        "simulations_per_cell": N_SIM,
        "alpha": ALPHA,
        "seed": SEED,
        "reproduce": "python3 build_fdr_simulation.py in sites/abwex/datasets/",
        "corrections": ["none", "bonferroni", "holm", "bh"],
        "metrics_defined": {
            "fwer": "share of experiments with at least one false positive among the true-null metrics",
            "fdr_mean": "mean realised false discovery proportion V/max(R,1)",
            "power_mean": "mean share of the true effects declared significant",
            "any_false_win": "share of experiments that declared at least one null metric a winner",
        },
        "cells": cells,
    }
    with open("fdr_simulation.json", "w") as f:
        json.dump(doc, f, indent=1)
    print(f"wrote fdr_simulation.json with {len(cells)} cells, {N_SIM} sims each")

    # console summary for the operator
    print("\nm  rho  pi0 | none fwer  bonf fwer  holm fwer   bh fdr | none pow  bonf pow  holm pow   bh pow")
    for c in cells:
        r = c["results"]
        pw = lambda k: "  n/a " if r[k]["power_mean"] is None else f"{r[k]['power_mean']:6.3f}"
        print(f"{c['m']:2d} {c['rho']:4.1f} {c['pi0']:4.2f} | "
              f"{r['none']['fwer']:9.4f} {r['bonferroni']['fwer']:10.4f} {r['holm']['fwer']:10.4f} "
              f"{r['bh']['fdr_mean']:7.4f} | {pw('none')} {pw('bonferroni')} {pw('holm')} {pw('bh')}")


if __name__ == "__main__":
    main()

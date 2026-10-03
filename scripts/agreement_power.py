"""
agreement_power.py -- how often could the registered agreement test have succeeded?

The registered criterion (ivygap/anatomic/agreement.py) is met when the Spearman correlation
between the ACS ranking and the accuracy ranking is >= 0.60 AND its percentile bootstrap CI over
methods excludes zero. With 12 comparable methods the observed rho was 0.081, and the paper calls
the result inconclusive -- underpowered -- rather than a refutation. This script puts a number on
"underpowered": the probability of meeting the criterion if the true rank correlation were
0.60, 0.70, 0.80 ..., at n = 12 and at larger panels.

It is a SENSITIVITY power analysis for effect sizes fixed in advance (the registered bar and
above), not "observed power" computed from the observed effect, which would be circular.
Data-generating model: (ACS, accuracy) pairs bivariate normal with Pearson r chosen so the
population Spearman correlation equals the target, rho_S = (6/pi) asin(r/2). The bootstrap is
the registered one (percentile, resampling methods, degenerate resamples skipped), vectorised;
`check_against_registered` compares the two on real-sized data.

    python3 scripts/agreement_power.py        # writes results/agreement_power.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config  # noqa: E402
from ivygap.anatomic import agreement  # noqa: E402

BAR = agreement.PREREGISTERED_RHO_THRESHOLD
TRUE_RHO = (0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9)
PANEL_N = (12, 20, 30, 50)
N_SIM, N_BOOT = 1000, 1000   # Monte Carlo SE of a power near 0.5: 1.6 points


def pearson_for_spearman(rho_s: float) -> float:
    return 2.0 * np.sin(np.pi * rho_s / 6.0)


def boot_ci(a: np.ndarray, t: np.ndarray, n_boot: int, rng) -> tuple[float, float]:
    """The registered percentile bootstrap over methods, vectorised."""
    n = a.size
    idx = rng.integers(0, n, size=(n_boot, n))
    ra = stats.rankdata(a[idx], axis=1)
    rt = stats.rankdata(t[idx], axis=1)
    ra -= ra.mean(1, keepdims=True)
    rt -= rt.mean(1, keepdims=True)
    den = np.sqrt((ra ** 2).sum(1) * (rt ** 2).sum(1))
    ok = den > 0                                  # a resample with no spread is skipped
    r = (ra[ok] * rt[ok]).sum(1) / den[ok]
    if r.size < 100:
        return float("nan"), float("nan")
    return float(np.quantile(r, 0.025)), float(np.quantile(r, 0.975))


def simulate(rho_s: float, n: int, n_sim: int, n_boot: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    r = pearson_for_spearman(rho_s)
    cov = np.array([[1.0, r], [r, 1.0]])
    met = above = excl = 0
    est = np.empty(n_sim)
    for i in range(n_sim):
        x = rng.multivariate_normal([0, 0], cov, size=n)
        rho = stats.spearmanr(x[:, 0], x[:, 1])[0]
        lo, _ = boot_ci(x[:, 0], x[:, 1], n_boot, rng)
        a, e = rho >= BAR, bool(np.isfinite(lo) and lo > 0)
        above += a
        excl += e
        met += a and e
        est[i] = rho
    return {"true_spearman": rho_s, "n_methods": n,
            "p_criterion_met": met / n_sim, "p_rho_at_or_above_bar": above / n_sim,
            "p_ci_excludes_zero": excl / n_sim,
            "median_estimated_rho": float(np.median(est))}


def check_against_registered(seed: int = 1) -> dict:
    """The vectorised bootstrap must agree with the registered loop to Monte Carlo error."""
    rng = np.random.default_rng(seed)
    x = rng.multivariate_normal([0, 0], [[1, .5], [.5, 1]], size=12)
    reg = agreement._bootstrap_rho(x[:, 0], x[:, 1], n_boot=5000, seed=seed)
    vec = boot_ci(x[:, 0], x[:, 1], 5000, np.random.default_rng(seed))
    return {"registered_ci": [round(v, 4) for v in reg], "vectorised_ci": [round(v, 4) for v in vec],
            "max_abs_difference": round(max(abs(a - b) for a, b in zip(reg, vec)), 4)}


def summarize(rows: list[dict]) -> dict:
    """Smallest simulated panel reaching 80% power, per true correlation. AT the bar it is
    unattainable at any size -- the criterion thresholds the estimate at the bar itself, so a
    true value equal to it is met about half the time -- and is reported as such, never as
    a number of methods."""
    need = {}
    for rs in (0.6, 0.7, 0.8):
        if rs == BAR:
            need[str(rs)] = "unattainable at any panel size (true value equals the bar)"
            continue
        ok = sorted(r["n_methods"] for r in rows if r["true_spearman"] == rs and r["p_criterion_met"] >= 0.8)
        lo = sorted(r["n_methods"] for r in rows if r["true_spearman"] == rs and r["p_criterion_met"] < 0.8)
        need[str(rs)] = (f"between {max(n for n in lo if n < ok[0])} and {ok[0]}" if ok and any(n < ok[0] for n in lo)
                         else (f"at most {ok[0]}" if ok else f"more than {max(PANEL_N)}"))
    return need


def main() -> int:
    if "--resummarize" in sys.argv:          # rewrite the summary from the stored rows only
        path = config.RESULTS_DIR / "agreement_power.json"
        out = json.loads(path.read_text())
        out["smallest_panel_with_80pct_power"] = summarize(out["rows"])
        path.write_text(json.dumps(out, indent=2))
        print("80% power needs:", out["smallest_panel_with_80pct_power"])
        return 0
    chk = check_against_registered()
    print("bootstrap check:", chk, flush=True)
    if chk["max_abs_difference"] > 0.05:
        print("ABORT: vectorised bootstrap does not reproduce the registered one")
        return 2
    rows = []
    for n in PANEL_N:
        for k, rs in enumerate(TRUE_RHO):
            row = simulate(rs, n, N_SIM, N_BOOT, seed=config.RANDOM_SEED + 1000 * n + k)
            rows.append(row)
            print(f"n={n:3d} true rho {rs:.1f}: criterion met {row['p_criterion_met']:.3f} "
                  f"(rho>=bar {row['p_rho_at_or_above_bar']:.3f}, CI>0 {row['p_ci_excludes_zero']:.3f})",
                  flush=True)
    need = summarize(rows)
    out = {"what_this_is": __doc__.split("\n")[1],
           "criterion": f"spearman >= {BAR} and percentile bootstrap CI over methods excludes zero",
           "model": "bivariate normal; Pearson r = 2 sin(pi rho_S / 6)",
           "n_sim": N_SIM, "n_boot": N_BOOT, "bootstrap_check": chk, "rows": rows,
           "smallest_panel_with_80pct_power": need,
           "kind": "sensitivity power analysis for pre-specified effect sizes (not observed power)"}
    (config.RESULTS_DIR / "agreement_power.json").write_text(json.dumps(out, indent=2))
    print("80% power needs:", need)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

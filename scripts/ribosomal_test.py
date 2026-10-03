"""Is the lymphoid inversion carried by RIBOSOMAL-PROTEIN genes admitted as B-cell markers? (WHY_B_OVER_T §7i)

HOW THIS HYPOTHESIS AROSE -- STATED BECAUSE IT IS POST HOC. scripts/b_column_diagnostics.py
asked which genes pull the fit toward B. Nine of the top twelve were ribosomal-protein genes
(RPL27A, RPL32, RPL13, ...). In the frozen signature the B_cell column carries 39.7% of its mass
on ribosomal-protein genes (T 14.0%), and in the TCGA arm's solved gene space all 13 RP genes
were admitted with B as their highest column -- 62.3% of the B column's mass there. Ribosomal
transcripts come from every cell in bulk tissue, so a ribosome-heavy column can absorb whatever
ribosomal signal the other columns leave unexplained.

Because the hypothesis was formed after seeing the data, it faces controls fixed before this
script first ran (2026-10-01):

  * TEST          remove ribosomal-protein genes -- the standard single-cell exclusion rule
                  ^RP[SL][0-9] / ^RPLP[0-9], not a list tuned to this result -- and refit;
  * MATCHED NULL  remove the same number of random NON-ribosomal genes whose highest column is
                  B (B markers that are not ribosomal), 50 draws. If these also restore T > B,
                  the effect is "losing B markers", not "losing ribosomes";
  * ANY-GENE NULL remove the same number of random genes of any kind, 50 draws;
  * IDENTIFIABILITY  the reduced signature must still recover a planted T:B of 7:3 in mixtures
                  built from it, with 25% multiplicative noise.

Frozen reference; the TCGA arm's own gene space; registered NNLS (all draws) and SVR (test, and
5 draws of each null). Methylation is a cohort-level truth only. Exploratory; selects nothing;
the registered results are unchanged whatever this shows. Writes results/ribosomal_test.json.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

sys.path.insert(0, str(Path(__file__).resolve().parent))
from b_column_diagnostics import cohort_inputs                 # noqa: E402
from lymphoid_ordering import _truth, k4s                      # noqa: E402

from ivygap import config                                      # noqa: E402
from ivygap.data.reference import load_frozen_reference        # noqa: E402
from ivygap.deconv.base import DeconvolutionInput              # noqa: E402
from ivygap.deconv.classical import NNLSDeconvolution, SVRDeconvolution  # noqa: E402

RP = re.compile(r"^(RP[SL]\d|RPLP\d)")
OUT = config.RESULTS_DIR / "ribosomal_test.json"
N_NULL_NNLS, N_NULL_SVR = 50, 0
# SVR NULL DRAWS DROPPED 2026-10-01, after the GBM SVR test came back: removing RP genes did not
# reduce B > T (90.2% -> 100%). The nulls exist to stop a FALSE POSITIVE -- an effect that random
# gene removal would also produce. With no effect to explain they cannot change the conclusion,
# and on this 2-core machine they would have cost ~4 h of CPU needed for BayesPrism. The NNLS nulls
# (cheap) and the identifiability check are kept.


def lymph_metrics(est: pd.DataFrame, truth: pd.DataFrame) -> dict:
    e = est.copy(); e.index = [k4s(i) for i in e.index]; e = e[~e.index.duplicated()]
    e = e.loc[[k for k in e.index if k in truth.index]]
    tot = e[["T_cell", "NK_cell", "B_cell"]].sum(axis=1); ok = tot > 0
    sh = e.loc[ok, ["T_cell", "NK_cell", "B_cell"]].div(tot[ok], axis=0)
    return {"n_with_lymphoid_signal": int(ok.sum()), "n_matched": int(len(e)),
            "frac_B_over_T": round(float((sh.B_cell > sh.T_cell).mean()), 4) if ok.any() else None,
            "frac_T_over_B": round(float((sh.T_cell > sh.B_cell).mean()), 4) if ok.any() else None,
            "mean_shares": {k: round(float(v), 4) for k, v in sh.mean().items()} if ok.any() else None}


def fit(mcls, ref, sub, man):
    return mcls().fit_predict(DeconvolutionInput(bulk=sub, references=(ref,), manifest=man,
                                                 cell_types=tuple(ref.cell_types)))


def planted_recovery(ref) -> dict:
    """Mixtures from the reduced signature with a known T:B = 7:3 among lymphocytes."""
    rng = np.random.default_rng(config.RANDOM_SEED)
    S = ref.profile.to_numpy(float); cts = list(ref.cell_types)
    iT, iB = cts.index("T_cell"), cts.index("B_cell")
    rec = []
    for _ in range(200):
        p = rng.dirichlet(np.ones(len(cts)))
        lym = p[iT] + p[iB]; p[iT], p[iB] = 0.7 * lym, 0.3 * lym
        b = (S @ p) * rng.lognormal(0, 0.25, S.shape[0])
        w = nnls(S, b)[0]
        if w[iT] + w[iB] > 0:
            rec.append(w[iT] / (w[iT] + w[iB]))
    rec = np.array(rec)
    return {"planted_T_share_of_T_plus_B": 0.7, "recovered_median": round(float(np.median(rec)), 4),
            "frac_recovering_T_over_B": round(float((rec > 0.5).mean()), 4), "n": int(len(rec))}


def main() -> int:
    base = load_frozen_reference()
    out = json.loads(OUT.read_text()) if OUT.exists() else {"cohorts": {}}
    out.update({"what_this_is": __doc__.split("\n")[0], "rule": "^RP[SL][0-9] | ^RPLP[0-9]",
                "n_null_draws": {"nnls": N_NULL_NNLS, "svr": N_NULL_SVR}})
    for c in ("gbm", "lgg"):
        sub, genes, man = cohort_inputs(c, base)
        base_c = "" if c == "gbm" else "_lgg"
        truth = _truth(config.RESULTS_DIR / f"methylation_celltypes{base_c}.csv").dropna()
        P = base.profile.loc[genes]
        rp = [g for g in genes if RP.match(g)]
        b_owned_nonrp = [g for g in genes if not RP.match(g) and P.loc[g].idxmax() == "B_cell"]
        nonrp = [g for g in genes if not RP.match(g)]
        rec = out["cohorts"].setdefault(c, {})
        rec.update({"n_genes": len(genes), "n_rp_removed": len(rp), "rp_genes": rp,
                    "n_b_owned_nonrp_pool": len(b_owned_nonrp),
                    "methylation_frac_T_over_B": round(float((truth.T_cell > truth.B_cell).mean()), 4)})

        def run(gene_list, mcls):
            ref = base.subset_genes(gene_list)
            return lymph_metrics(fit(mcls, ref, sub.loc[gene_list], man), truth)

        keep_all = list(genes)
        no_rp = [g for g in genes if g not in set(rp)]
        rng = np.random.default_rng(config.RANDOM_SEED)
        for mcls, n_null in ((NNLSDeconvolution, N_NULL_NNLS), (SVRDeconvolution, N_NULL_SVR)):
            m = mcls.name
            if f"{m}_baseline" not in rec:
                rec[f"{m}_baseline"] = run(keep_all, mcls)
                rec[f"{m}_rp_removed"] = run(no_rp, mcls)
                OUT.write_text(json.dumps(out, indent=2)); print(f"  {c}/{m}: baseline + rp_removed saved", flush=True)
            for null, pool in (("matched_null", b_owned_nonrp), ("anygene_null", nonrp)):
                key = f"{m}_{null}"
                if key in rec:
                    continue
                draws = []
                for _ in range(n_null):
                    drop = set(rng.choice(pool, size=len(rp), replace=False))
                    draws.append(run([g for g in genes if g not in drop], mcls)["frac_B_over_T"])
                draws = [d for d in draws if d is not None]
                rec[key] = {"frac_B_over_T_draws": draws,
                            "median": round(float(np.median(draws)), 4) if draws else None,
                            "min": round(float(np.min(draws)), 4) if draws else None}
                OUT.write_text(json.dumps(out, indent=2)); print(f"  {c}/{key} saved", flush=True)
        if "identifiability_rp_removed" not in rec:
            rec["identifiability_baseline"] = planted_recovery(base.subset_genes(keep_all))
            rec["identifiability_rp_removed"] = planted_recovery(base.subset_genes(no_rp))
            OUT.write_text(json.dumps(out, indent=2)); print(f"  {c}/identifiability saved", flush=True)

    for c, r in out["cohorts"].items():
        print(f"\n=== {c.upper()}  ({r['n_genes']} genes; {r['n_rp_removed']} RP removed; methylation T>B "
              f"{r['methylation_frac_T_over_B']:.1%})")
        for m in ("nnls", "svr"):
            b, t = r.get(f"{m}_baseline"), r.get(f"{m}_rp_removed")
            if b and t:
                print(f"   {m}: B>T baseline {b['frac_B_over_T']} ({b['n_with_lymphoid_signal']} w/ signal) -> "
                      f"RP removed {t['frac_B_over_T']} ({t['n_with_lymphoid_signal']} w/ signal); shares {t['mean_shares']}")
            for nl in ("matched_null", "anygene_null"):
                v = r.get(f"{m}_{nl}")
                if v:
                    print(f"      {nl}: B>T median {v['median']} min {v['min']} over {len(v['frac_B_over_T_draws'])} draws")
        for k in ("identifiability_baseline", "identifiability_rp_removed"):
            if k in r:
                print(f"   {k}: {r[k]}")
    print(f"\nwrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

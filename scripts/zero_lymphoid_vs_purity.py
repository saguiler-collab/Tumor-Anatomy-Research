"""
zero_lymphoid_vs_purity.py -- do methods return no lymphocytes preferentially in high-purity samples?

EXPLORATORY, not pre-registered. Until 2026-10-06 this analysis existed only as a one-off command in the
session transcript (line 19476); the verification re-run (scripts/verify_rerun.py) found that no script
produced the artefact the manuscript reads. The code below is that command, unchanged in logic; only
its file paths now come from ivygap.config. tests/ and verify_rerun.py check it reproduces the artefact.
"""
import pandas as pd, numpy as np, json
from ivygap import config
from scipy import stats
C=["T_cell","B_cell","NK_cell"]
k4s=lambda s:"-".join(str(s).split("-")[:3]+[str(s).split("-")[3][:2]]) if len(str(s).split("-"))>3 else str(s)
print("HYPOTHESIS (exploratory, NOT pre-registered): a method returns ZERO lymphoid signal")
print("more often in HIGH-purity samples -- losing the immune compartment exactly where the")
print("tumour dominates, which is where real tissue sits.")
print("Falsified if zero-lymphoid samples are not higher-purity than the rest.\n")
out={}
for tag,lab in [("","GBM"),("_lgg","LGG")]:
    est=pd.read_csv(config.RESULTS_DIR / f"estimates_full{tag}.csv")
    pur=pd.read_csv(config.RESULTS_DIR / f"absolute_purity_per_sample{tag}.csv",index_col=0)
    tcol="absolute_purity" if "absolute_purity" in pur else "purity"
    P=pd.Series(pur[tcol].values,index=[k4s(i) for i in pur.index]); P=P[~P.index.duplicated()]
    sc=next(c for c in ("sample","sample_id") if c in est.columns); est["k"]=est[sc].map(k4s)
    print(f"=== {lab} (n={len(P)} samples with purity) ===")
    print(f"{'method':16s} {'n zero':>7s} {'pur|zero':>9s} {'pur|nonzero':>12s} {'diff':>8s} {'MWU p':>10s}")
    print("-"*70)
    res={}; ps=[]
    for m,g in est.groupby("method"):
        g=g.drop_duplicates("k").set_index("k")
        idx=[i for i in g.index if i in P.index]
        if len(idx)<50: continue
        z=(g.loc[idx,C].sum(axis=1)==0).to_numpy()
        if z.sum()<5 or (~z).sum()<5:
            print(f"{m:16s} {int(z.sum()):7d}   (one group <5, not testable)")
            res[m]={"n_zero":int(z.sum()),"testable":False}; continue
        a,b=P.loc[idx][z].to_numpy(),P.loc[idx][~z].to_numpy()
        p=float(stats.mannwhitneyu(a,b,alternative="greater").pvalue); ps.append(p)
        res[m]={"n_zero":int(z.sum()),"purity_zero":round(float(a.mean()),4),
                "purity_nonzero":round(float(b.mean()),4),
                "diff":round(float(a.mean()-b.mean()),4),"mwu_p_greater":round(p,6),"testable":True}
        print(f"{m:16s} {int(z.sum()):7d} {a.mean():9.4f} {b.mean():12.4f} "
              f"{a.mean()-b.mean():+8.4f} {p:10.2e}")
    if ps:
        k=sum(p<0.05 for p in ps)
        print(f"  -> {k} of {len(ps)} testable methods: zero-lymphoid samples are "
              f"significantly HIGHER purity (one-sided MWU)")
        out[lab]={"methods":res,"n_significant":k,"n_testable":len(ps)}
    print()
json.dump({"EXPLORATORY":"not pre-registered; generated after the absence failure mode was found",
           "hypothesis":"methods return zero lymphoid signal preferentially in high-purity samples",
           "cohorts":out}, open(config.RESULTS_DIR / "zero_lymphoid_vs_purity.json","w"), indent=2)
print("wrote results/zero_lymphoid_vs_purity.json")

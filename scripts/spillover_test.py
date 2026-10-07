"""
spillover_test.py -- does the B column absorb macrophage/microglia spillover? (TCGA-LGG)

EXPLORATORY, not pre-registered. Until 2026-10-06 this analysis existed only as a one-off command in the
session transcript (line 18322); the verification re-run (scripts/verify_rerun.py) found that no script
produced the artefact the manuscript reads. The code below is that command, unchanged in logic; only
its file paths now come from ivygap.config. tests/ and verify_rerun.py check it reproduces the artefact.
"""
import numpy as np, pandas as pd, json
from scipy import stats
from ivygap import config
from ivygap.deconv.comparability import non_comparable

full = pd.read_csv(config.RESULTS_DIR / "estimates_full_lgg.csv")
scol = next(c for c in ("sample","sample_id") if c in full.columns)
M, T, B, NK = "Macrophage_Microglia", "T_cell", "B_cell", "NK_cell"
nc = non_comparable("frozen")
print("SPILLOVER HYPOTHESIS (exploratory, NOT pre-registered)")
print("  If B absorbs macrophage/microglia spillover, then across samples the B estimate")
print("  should track the Macrophage estimate MORE than the T estimate does.")
print("  Falsified if corr(B,M) <= corr(T,M) in most methods.\n")
print(f"{'method':18s} {'r(B,M)':>8s} {'r(T,M)':>8s} {'r(NK,M)':>8s}  {'B>T?':>6s}")
print("-"*58)
res={}; wins=0; n=0
for m,g in full.groupby("method"):
    g=g.dropna(subset=[M,T,B])
    if len(g)<50: continue
    def r(a,b):
        x,y=g[a].to_numpy(float),g[b].to_numpy(float)
        ok=np.isfinite(x)&np.isfinite(y)
        if ok.sum()<50 or np.std(x[ok])==0 or np.std(y[ok])==0: return float("nan")
        return float(stats.spearmanr(x[ok],y[ok]).statistic)
    rb,rt,rn=r(B,M),r(T,M),r(NK,M)
    w = (not np.isnan(rb)) and (not np.isnan(rt)) and rb>rt
    if not (np.isnan(rb) or np.isnan(rt)): n+=1; wins+=w
    res[m]={"r_B_M":None if np.isnan(rb) else round(rb,4),
            "r_T_M":None if np.isnan(rt) else round(rt,4),
            "r_NK_M":None if np.isnan(rn) else round(rn,4),"B_tracks_M_more":bool(w),
            "comparable": m not in nc}
    print(f"{m:18s} {rb:8.3f} {rt:8.3f} {rn:8.3f}  {str(w):>6s}"
          f"{'  [not comparable]' if m in nc else ''}")
p = stats.binomtest(wins,n,0.5,alternative="greater").pvalue if n else float("nan")
print(f"\n{wins} of {n} methods show r(B,M) > r(T,M).  sign-test p = {p:.4f}")
print("VERDICT:", "SUPPORTED" if p<0.05 else "NOT SUPPORTED -- the spillover story is not carried by these data")
json.dump({"hypothesis":"B_cell absorbs Macrophage_Microglia spillover",
           "status":"EXPLORATORY - not pre-registered, generated after seeing the inversion",
           "methods":res,"n_supporting":wins,"n_tested":n,
           "sign_test_p":None if np.isnan(p) else round(float(p),6),
           "verdict":"SUPPORTED" if p<0.05 else "NOT SUPPORTED"},
          open(config.RESULTS_DIR/"spillover_lgg.json","w"),indent=2)
print("wrote results/spillover_lgg.json")

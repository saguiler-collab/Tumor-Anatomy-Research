"""Does the lymphoid inversion persist under an INDEPENDENT atlas? -- prespecified/abdelfattah_cluster_mapping.md.

Refits TCGA GBM and LGG against the reference built from Abdelfattah et al. 2022 (GSE182109; not a
GBmap source study) with the registered NNLS and SVR, on that reference's own marker space chosen by
the same rule the frozen TCGA arm uses (`select_signature_genes`, SIGNATURE_GENES_PER_TYPE), and
compares the lymphoid ordering with DNA methylation exactly as scripts/lymphoid_ordering.py does.

READING, fixed in the pre-specification before any expression was seen: T > B under this atlas where
GBmap gives B > T means the inversion belongs to GBmap's lymphoid profiles; B > T persisting means it
does not depend on GBmap's annotation or donors.

Also records, for each atlas's B column, its ribosomal-protein mass share and its correlation with
the mean bulk profile -- the two properties measured for GBmap in WHY_B_OVER_T §7h.
Methylation is a cohort-level truth (results/lymphoid_tracking.json). Writes
results/independent_atlas_test.json. Exploratory; selects nothing.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from b_column_diagnostics import cohort_inputs                 # noqa: E402  (frozen-arm inputs, for samples)
from lymphoid_ordering import COMPARE, _truth, k4s, ordering_of, renormalise  # noqa: E402

from ivygap import config                                      # noqa: E402
from ivygap.data.reference import select_signature_genes      # noqa: E402
from ivygap.deconv.base import DeconvolutionInput, ReferenceBundle  # noqa: E402
from ivygap.deconv.classical import NNLSDeconvolution, SVRDeconvolution  # noqa: E402

RP = re.compile(r"^(RP[SL]\d|RPLP\d)")
IG = re.compile(r"^(IG[HKL][VDJC]|IGH[GAMDE]|IGKC|IGLC|IGLL|JCHAIN)")   # declared in the addendum, before the result
D = config.REFERENCE_DIR / "abdelfattah_2022"


def ref_dir(seed: int) -> Path:
    """The registered build, or a seed-sensitivity build (build_abdelfattah_reference.py --seeds)."""
    return D if seed == config.RANDOM_SEED else D.parent / f"{D.name}_seed{seed}"


def load() -> ReferenceBundle:
    prof = pd.read_csv(D / "profile.csv", index_col=0)
    prov = json.loads((D / "provenance.json").read_text())
    types = [t for t in config.CELL_TYPES if t in prof.columns]
    return ReferenceBundle(name="abdelfattah_2022", profile=prof[types],
                           sigma=pd.read_csv(D / "sigma.csv", index_col=0)[types],
                           cell_size=pd.read_csv(D / "cell_size.csv", index_col=0)["cell_size"].reindex(types),
                           n_donors=int(prov["n_patients"]), donor_profiles=None,
                           has_cross_donor_variance=True)


def main() -> int:
    import argparse
    import gzip
    ap = argparse.ArgumentParser()
    ap.add_argument("--drop-ig", action="store_true",
                    help="pre-declared sensitivity: remove immunoglobulin genes from the marker space")
    ap.add_argument("--seed", type=int, default=config.RANDOM_SEED,
                    help="use the reference built with this cell-sampling seed (seed sensitivity)")
    a = ap.parse_args()
    global D
    D = ref_dir(a.seed)
    if not (D / "profile.csv").exists():
        print(f"BLOCKED: no reference at {D} -- build it with build_abdelfattah_reference.py --seeds {a.seed}")
        return 2
    ref_all = load()
    out = {"what_this_is": __doc__.split("\n")[0], "reference": "abdelfattah_2022",
           "types": list(ref_all.cell_types),
           "provenance": json.loads((D / "provenance.json").read_text()), "cohorts": {}}
    from ivygap.data.reference import load_frozen_reference   # noqa: PLC0415
    frozen = load_frozen_reference()
    for c in ("gbm", "lgg"):
        base = "" if c == "gbm" else "_lgg"
        sub_frozen, _, man_frozen = cohort_inputs(c, frozen)          # the frozen arm's samples
        with gzip.open(config.PROCESSED_DIR / f"tcga_{c}_bulk_cpm.csv.gz", "rt") as fh:
            bulk = pd.read_csv(fh, index_col=0)
        cols = list(sub_frozen.columns)
        shared = [g for g in ref_all.profile.index if g in bulk.index]
        genes = [g for g in select_signature_genes(ref_all.subset_genes(shared),
                                                   n_per_type=config.SIGNATURE_GENES_PER_TYPE) if g in bulk.index]
        n_ig = sum(bool(IG.match(g)) for g in genes)
        if a.drop_ig:
            genes = [g for g in genes if not IG.match(g)]
        sub = bulk.loc[genes, cols]; sub = sub / sub.sum(axis=0) * 1e6
        ref = ref_all.subset_genes(genes)
        truth = _truth(config.RESULTS_DIR / f"methylation_celltypes{base}.csv").dropna()
        rec = {"n_samples": len(cols), "n_genes": len(genes), "n_ig_in_marker_space": n_ig, "ig_removed": a.drop_ig,
               "ig_mass_share_B_before_removal": round(float(ref_all.profile.loc[[g for g in ref_all.profile.index if g in set(genes) or (IG.match(g) and g in shared)], "B_cell"].pipe(lambda x: x[[g for g in x.index if IG.match(g)]].sum() / max(x.sum(), 1e-12))), 4),
               "methylation_ordering": ordering_of(truth.mean()),
               "methylation_frac_T_over_B": round(float((truth.T_cell > truth.B_cell).mean()), 4), "methods": {}}
        P = ref.profile; mb = sub.mean(axis=1)
        rec["profile_properties"] = {t: {"rp_mass_share": round(float(P.loc[[g for g in genes if RP.match(g)], t].sum() / P[t].sum()), 4),
                                         "r_with_mean_bulk": round(float(np.corrcoef(np.log1p(P[t]), np.log1p(mb))[0, 1]), 4)}
                                     for t in COMPARE if t in P.columns}
        for mcls in (NNLSDeconvolution, SVRDeconvolution):
            est = mcls().fit_predict(DeconvolutionInput(bulk=sub, references=(ref,), manifest=man_frozen.loc[cols],
                                                        cell_types=tuple(ref.cell_types)))
            e = est.reindex(columns=list(config.CELL_TYPES)).fillna(0.0)
            e.index = [k4s(i) for i in e.index]; e = e[~e.index.duplicated()]
            idx = [k for k in e.index if k in truth.index]
            ren = renormalise(e.loc[idx]); mu = ren.mean(); okr = ren.dropna()
            rec["methods"][mcls.name] = {
                "n_matched": len(idx), "n_with_lymphoid_signal": int(len(okr)),
                "lymphoid_mean": {k: round(float(mu[k]), 4) for k in COMPARE},
                "ordering": ordering_of(mu), "T_exceeds_B": bool(mu.T_cell > mu.B_cell),
                "frac_samples_B_over_T": round(float((okr.B_cell > okr.T_cell).mean()), 4) if len(okr) else None}
            reg = json.loads((config.RESULTS_DIR / f"lymphoid_ordering{base}.json").read_text())["methods"].get(mcls.name, {})
            rec["methods"][mcls.name]["gbmap_frozen_for_comparison"] = {
                "ordering": reg.get("ordering"), "frac_samples_B_over_T": reg.get("frac_samples_est_B_over_T"),
                "n_with_lymphoid_signal": reg.get("n_scored")}
        out["cohorts"][c] = rec
    out["ig_removed"] = a.drop_ig
    name = "independent_atlas_test_ig_removed.json" if a.drop_ig else "independent_atlas_test.json"
    if a.seed != config.RANDOM_SEED:
        name = name.replace(".json", f"_seed{a.seed}.json")
    out["seed"] = a.seed
    (config.RESULTS_DIR / name).write_text(json.dumps(out, indent=2))
    for c, r in out["cohorts"].items():
        print(f"\n=== {c.upper()} ({r['n_samples']} samples, {r['n_genes']} genes); methylation {r['methylation_ordering']}, "
              f"T>B in {r['methylation_frac_T_over_B']:.1%}")
        print("   profile properties:", r["profile_properties"])
        for m, v in r["methods"].items():
            gf = v["gbmap_frozen_for_comparison"]
            print(f"   {m:<5} Abdelfattah: {v['ordering']} (B>T in {v['frac_samples_B_over_T']}, signal {v['n_with_lymphoid_signal']}/{v['n_matched']})"
                  f"   | GBmap frozen: {gf['ordering']} (B>T in {gf['frac_samples_B_over_T']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

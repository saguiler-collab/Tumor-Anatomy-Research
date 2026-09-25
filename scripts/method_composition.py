"""What each estimator actually IS: input modality, reference type, and provenance.

WHY THIS EXISTS. The manuscript's alternative title says "bulk RNA deconvolution cannot
report the lymphoid compartment", and a reader is entitled to ask what "the methods" were.
Two facts have to be separable and were not written down anywhere:

1. **Input modality.** Every estimator here decomposes BULK RNA. None reads methylation,
   protein or imaging. That matters because the *truth* instruments deliberately do: DNA
   methylation (EpiDISH), DNA copy number (ABSOLUTE purity) and H&E histology. The whole
   design rests on the truth sharing no modality with the thing being tested, and saying
   "all our methods are bulk RNA" is a statement about the METHODS, never about the truth.

2. **Provenance.** They are NOT all published deconvolution tools. Of 15 estimators, 8 are
   published packages (`registry.PUBLISHED_TOOLS`), 4 are published algorithms reimplemented
   here, and 5 are classical regression baselines this project wrote. A claim about
   "deconvolution methods" that is really about NNLS and elastic net would be overbroad; a
   claim that the *published tools* fail is narrower and better supported.

Usage:  python3 scripts/method_composition.py [--write]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                                  # noqa: E402
from ivygap.deconv.registry import PUBLISHED_TOOLS                         # noqa: E402
from ivygap.deconv.r_bridge import (OWN_SIGNATURE_METHODS,                 # noqa: E402
                                    R_RETURNS_CELL_FRACTIONS,
                                    SIGNATURE_ONLY_METHODS)

#: What each estimator is, in one line. Sourced from the class docstrings in
#: `ivygap/deconv/registry.py`; kept here so the roster is one artefact rather than
#: fifteen docstrings a reader has to assemble.
WHAT_IT_IS = {
    "music":      "MuSiC — cross-donor-variance-weighted NNLS over individual cells.",
    "scdc":       "SCDC — ENSEMBLE-capable weighted regression over individual cells.",
    "scdc_ensemble": "SCDC in ENSEMBLE mode; with one reference it reduces exactly to SCDC.",
    "bisque":     "BisqueRNA — assay-transform regression; needs subjects assayed as BOTH "
                  "bulk and single cells.",
    "epic":       "EPIC — constrained least squares against a reference PROFILE, with "
                  "per-gene variance weighting.",
    "quantiseq":  "quanTIseq — constrained regression against its OWN built-in TIL10 "
                  "signature; ignores the reference it is handed.",
    "bayesprism": "BayesPrism — Bayesian hierarchical model over individual cells.",
    "dwls":       "DWLS — dampened weighted least squares with a condition-number gene search.",
    "cibersortx": "CIBERSORTx B-mode — nu-SVR plus bulk-mode batch correction.",
    "cibersortx_smode": "CIBERSORTx S-mode — nu-SVR plus batch correction defined against "
                        "the source cells.",
    "svr":        "nu-SVR -- the CIBERSORT core algorithm (nu swept over 0.25/0.5/0.75), "
                  "without CIBERSORTx batch correction. Named `svr`, not `cibersort`, because it is a reimplementation.",
    "nnls":       "Plain non-negative least squares, no regularisation. Baseline.",
    "elastic_net": "Non-negative elastic net with per-sample outcome-blind tuning. Baseline.",
    "bayesian":   "Per-sample collapsed Gibbs sampler over a Dirichlet-multinomial. Baseline.",
    "bayesian_hierarchical": "Hierarchical version of the above, pooling across samples. "
                             "Baseline.",
    "control_random": "NEGATIVE CONTROL — random cell proportions.",
    "control_shuffled_signature": "NEGATIVE CONTROL — signature matrix with gene labels "
                                  "shuffled.",
}

#: Published algorithms this project reimplemented rather than ran as the vendor's code.
#:
#: `svr` belongs here and NOT among the baselines. It is nu-support-vector regression on
#: standardised inputs with nu swept over {0.25, 0.5, 0.75} -- the CORE ALGORITHM OF
#: CIBERSORT (Newman 2015), without CIBERSORTx's batch correction. Filing it beside NNLS and
#: elastic net would understate the panel: it is a published method's engine, and it was the
#: predecessor project's frozen method. It runs under the name `svr` rather than `cibersort`
#: because CLAUDE.md forbids reporting a reimplementation under the published package's name.
REIMPLEMENTED = {"cibersortx", "cibersortx_smode", "bayesprism", "dwls", "svr"}

#: Estimators that are NOT any published deconvolution tool -- generic regression applied to
#: the deconvolution problem. Their job in the panel is to show that a failure is not an
#: artefact of one vendor's implementation.
BASELINES = {"nnls", "elastic_net", "bayesian", "bayesian_hierarchical"}

#: The tier label for BASELINES. Held as a constant because it contains the word
#: "published" ("not a published tool"), so a substring filter for published tools
#: silently matches it -- which it did, reporting 12 published tools where there
#: are 9.
BASELINE_TIER_NAME = "classical regression baseline, not a published tool"


def build() -> dict:
    lb = pd.read_csv(config.RESULTS_DIR / "anatomic/acs_leaderboard.csv").set_index("method")
    rows = []
    for m, r in lb.iterrows():
        impl = str(r["implementation"])
        is_ctl = bool(r["is_control"])
        if is_ctl:
            tier = "negative control"
        elif m in PUBLISHED_TOOLS and impl.startswith("R:"):
            tier = "published package, run as the vendor's code"
        elif m in PUBLISHED_TOOLS:
            tier = "published package, Python reimplementation used"
        elif m in REIMPLEMENTED:
            tier = "published algorithm, reimplemented here"
        elif m in BASELINES:
            tier = BASELINE_TIER_NAME
        else:
            tier = "UNCLASSIFIED -- add it to a tier before quoting this table"
        rows.append({
            "method": m,
            "what_it_is": WHAT_IT_IS.get(m, "—"),
            "input_modality": "bulk RNA",
            "reference_it_consumes":
                "its own built-in signature" if m in OWN_SIGNATURE_METHODS else
                "a reference profile (signature matrix)" if m in SIGNATURE_ONLY_METHODS else
                "individual cells with donor and cell-type labels"
                if m in PUBLISHED_TOOLS and m not in SIGNATURE_ONLY_METHODS else
                "a reference profile (signature matrix)",
            "tier": tier,
            "implementation_used": impl,
            "returns": "cell fractions" if m in R_RETURNS_CELL_FRACTIONS else "mRNA share",
            "acs": round(float(r["acs"]), 4),
            "is_control": is_ctl,
            "degenerate": bool(r["degenerate"]),
            "comparable": bool(r["comparable"]),
        })
    df = pd.DataFrame(rows)
    real = df[~df.is_control]
    return {
        "what_this_is":
            "Composition of the estimator panel: what each one is, what it reads, and "
            "whether it is a published tool or a baseline written here.",
        "input_modality_finding":
            f"All {len(df)} estimators decompose BULK RNA. None reads another modality. "
            f"The TRUTH instruments deliberately do — DNA methylation (EpiDISH), DNA copy "
            f"number (ABSOLUTE purity) and H&E histology — which is the point of the design "
            f"and must not be confused with the methods.",
        "n_total": int(len(df)),
        "n_real_methods": int(len(real)),
        "n_controls": int(df.is_control.sum()),
        "tiers": {k: int(v) for k, v in real.tier.value_counts().items()},
        "published_tools": sorted(PUBLISHED_TOOLS),
        "reimplemented_published_algorithms": sorted(REIMPLEMENTED),
        "baselines": sorted(BASELINES & set(real.method)),
        "unclassified": sorted(real.loc[real.tier.str.startswith("UNCLASSIFIED"), "method"]),
        "scope_caveat":
            "A claim about 'deconvolution methods' that rests on NNLS and elastic net is "
            "overbroad. The narrower and better-supported claim is that the PUBLISHED tools "
            "fail, with the baselines showing the failure is not an artefact of any one "
            "implementation.",
        "methods": df.to_dict(orient="records"),
        # DOES THE HEADLINE SURVIVE WITHOUT THE BASELINES? The lymphoid claim is the study's
        # strongest, and a reader is entitled to ask whether it leans on generic regressors.
        # Computed here rather than asserted, and beware the substring trap: the baseline
        # tier's own name contains the word "published".
        "lymphoid_by_tier": _lymphoid_by_tier(df),
    }


def _lymphoid_by_tier(df: pd.DataFrame) -> dict:
    out = {}
    for tag, cohort in (("", "GBM"), ("_lgg", "LGG")):
        f = config.RESULTS_DIR / f"lymphoid_ordering{tag}.json"
        if not f.exists():
            continue
        lo = json.loads(f.read_text())["methods"]
        tier = dict(zip(df.method, df.tier))
        rows = [(m, tier.get(m, "?"), bool(d.get("T_exceeds_B"))) for m, d in lo.items()]
        base = [r for r in rows if r[1] == BASELINE_TIER_NAME]
        pub = [r for r in rows if r[1] != BASELINE_TIER_NAME]
        out[cohort] = {
            "n_compared": len(rows),
            "n_correct_all": sum(r[2] for r in rows),
            "published_tools_and_algorithms": {
                "n": len(pub), "n_correct": sum(r[2] for r in pub),
                "methods": sorted(r[0] for r in pub)},
            "classical_baselines": {
                "n": len(base), "n_correct": sum(r[2] for r in base),
                "methods": sorted(r[0] for r in base)},
            "absent_no_estimate_returned": sorted(
                set(df.loc[~df.is_control, "method"]) - set(lo)),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    res = build()
    df = pd.DataFrame(res["methods"])
    pd.set_option("display.width", 200, "display.max_colwidth", 46)
    print(df[["method", "tier", "reference_it_consumes", "implementation_used",
              "returns", "acs"]].to_string(index=False))
    print()
    print(res["input_modality_finding"])
    print()
    for k, v in res["tiers"].items():
        print(f"  {v:>2}  {k}")
    print()
    print(res["scope_caveat"])
    if args.write:
        out = config.RESULTS_DIR / "method_composition.json"
        out.write_text(json.dumps(res, indent=2) + "\n")
        print(f"\nwrote {out.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

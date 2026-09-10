#!/usr/bin/env python3
"""
benchmark_concordance.py — does the ACS ranking agree with published benchmarks?

THE QUESTION
------------
This study's claim is that anatomic concordance can rank deconvolution methods without
ground truth. The internal test correlates the ACS ranking against accuracy on
donor-held-out pseudobulk built from the same atlas the reference comes from.

This asks the same question with **someone else's ground truth, on someone else's
tissue, ranked by someone else**. If ACS puts the methods an independent benchmark
called best near the top, that is external corroboration costing no new data. If it does
not, that is equally publishable and this project reports it.

WHAT IS PRE-SPECIFIED, AND WHY IT HAS TO BE
-------------------------------------------
The obvious way to fake this is to pick the published claim, the method subset, or the
statistic after seeing which combination agrees. So all three are fixed below, in this
file, before the comparison runs:

  * the published claim is a VERBATIM quotation naming a top tier — not a ranking read
    off a figure, and not a metric this project chose;
  * the method set is every method both this study and that paper evaluate, with no
    discretion;
  * the statistic is a one-sided Mann-Whitney U of ACS between the paper's top tier and
    the rest, decided before running.

Anything excluded is excluded here with its reason, in advance.

WHAT THIS CANNOT SHOW
---------------------
Agreement does not prove ACS is correct — both could be wrong in the same direction, and
both rest partly on simulated mixtures, which Nguyen et al. (2024) argue favours methods
sharing the simulator's assumptions. It shows that ACS, computed with no ground truth on
one tissue, orders methods compatibly with a benchmark that used ground truth on others.
That is worth something and it is not everything.

    python scripts/benchmark_concordance.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                        # noqa: E402

#: Published claims, quoted verbatim, each naming a top tier. Fixed before running.
PUBLISHED = {
    "avila_cobos_2020": {
        "citation": ("Avila Cobos F, et al. Benchmarking of cell type deconvolution "
                     "pipelines for transcriptomics data. Nat Commun 11:5650 (2020). "
                     "doi:10.1038/s41467-020-19015-1"),
        "quote": ("the five best bulk deconvolution methods (OLS, nnls, RLR, FARDEEP, "
                  "and CIBERSORT) and the three best methods that use scRNA-seq data as "
                  "reference (DWLS, MuSiC, SCDC) achieved median RMSE values lower "
                  "than 0.05"),
        "ground_truth": "simulated pseudobulk with known composition",
        # Their names -> ours. Only methods THIS study also runs.
        "top_tier": {"nnls": "nnls", "CIBERSORT": "svr", "DWLS": "dwls",
                     "MuSiC": "music", "SCDC": "scdc"},
        # Evaluated by them, run by us, NOT in their top tier. This is the comparison
        # group and it must be stated, or "top tier ranks high" is untestable.
        "evaluated_not_top": {"EPIC": "epic", "Bisque": "bisque",
                              "elastic net": "elastic_net"},
    },
    "sturm_2019": {
        "citation": ("Sturm G, et al. Comprehensive evaluation of transcriptome-based "
                     "cell-type quantification methods for immuno-oncology. "
                     "Bioinformatics 35:i436 (2019). doi:10.1093/bioinformatics/btz363"),
        "quote": ("due to a robust overall performance, we recommend EPIC and quanTIseq "
                  "for general purpose deconvolution"),
        "ground_truth": "flow cytometry and simulated immune mixtures",
        "top_tier": {"EPIC": "epic", "quanTIseq": "quantiseq"},
        "evaluated_not_top": {"CIBERSORT": "svr"},
    },
}

#: Excluded in advance, with reasons, so the exclusions are not a degree of freedom.
EXCLUSIONS = {
    "quantiseq": ("Excluded from the STATISTIC, reported separately. quanTIseq models "
                  "four of this project's eight roster types, so its ACS is computed "
                  "over 15 constraint-tumour pairs rather than 57 and carries "
                  "comparable=False. Including a score on a different denominator in a "
                  "rank test would compare two different quantities."),
    "dwls": ("Retained in the statistic but flagged. In the confirmatory run DWLS ran as "
             "a Python reimplementation in both anatomic stages, having exceeded its "
             "2,400s budget for the genuine R package. Avila Cobos evaluated the "
             "package. This is disclosed rather than dropped, because dropping the one "
             "method expected to disagree is exactly the move this design exists to "
             "prevent."),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out", default="results/benchmark_concordance.json")
    args = ap.parse_args()

    lb_path = Path(args.results_dir) / "anatomic" / "acs_leaderboard.csv"
    if not lb_path.exists():
        print(f"missing {lb_path}")
        return 1
    lb = pd.read_csv(lb_path).set_index("method")

    report = {
        "what_this_is": ("Does the ACS ranking agree with published benchmarks that used "
                         "real ground truth on other tissue? EXPLORATORY; pre-specified "
                         "in scripts/benchmark_concordance.py before running."),
        "acs_source": str(lb_path),
        "exclusions_declared_in_advance": EXCLUSIONS,
        "comparisons": {},
    }

    for key, spec in PUBLISHED.items():
        top = {ours: theirs for theirs, ours in spec["top_tier"].items()
               if ours in lb.index}
        rest = {ours: theirs for theirs, ours in spec["evaluated_not_top"].items()
                if ours in lb.index}

        # quanTIseq is reported but never enters the statistic.
        stat_top = {k: v for k, v in top.items() if k != "quantiseq"}
        stat_rest = {k: v for k, v in rest.items() if k != "quantiseq"}

        a = [float(lb.loc[m, "acs"]) for m in stat_top]
        b = [float(lb.loc[m, "acs"]) for m in stat_rest]

        res = {
            "citation": spec["citation"],
            "published_claim": spec["quote"],
            "their_ground_truth": spec["ground_truth"],
            "their_top_tier_our_acs": {m: round(float(lb.loc[m, "acs"]), 4)
                                       for m in top},
            "evaluated_not_top_our_acs": {m: round(float(lb.loc[m, "acs"]), 4)
                                          for m in rest},
            "n_overlapping_methods": len(top) + len(rest),
        }
        if len(a) >= 2 and len(b) >= 2:
            u, p = mannwhitneyu(a, b, alternative="greater")
            res["mann_whitney_u"] = float(u)
            res["p_one_sided_top_tier_higher"] = round(float(p), 5)
        else:
            res["mann_whitney_u"] = None
            res["p_one_sided_top_tier_higher"] = None
            res["why_no_statistic"] = (
                f"needs >=2 methods per group; have {len(a)} top and {len(b)} other. "
                f"With this few overlapping methods the comparison is descriptive.")

        # Rank position of each of their top-tier methods in our leaderboard.
        comparable = lb[(~lb["is_control"].astype(bool)) & lb["comparable"].astype(bool)]
        order = comparable["acs"].rank(ascending=False, method="min")
        res["our_rank_of_their_top_tier"] = {
            m: (int(order[m]) if m in order.index else None) for m in top}
        res["n_comparable_methods"] = int(len(comparable))
        report["comparisons"][key] = res

        print(f"\n=== {key}")
        print(f"  claim: \"{spec['quote'][:88]}...\"")
        print(f"  their top tier, our ACS and rank (of {len(comparable)}):")
        for m in top:
            r = res["our_rank_of_their_top_tier"][m]
            flag = "  [reimplementation]" if m == "dwls" else ""
            print(f"     {m:12s} ACS {lb.loc[m,'acs']:.4f}  rank "
                  f"{r if r else '-'}{flag}")
        if rest:
            print("  evaluated but not top tier:")
            for m in rest:
                print(f"     {m:12s} ACS {lb.loc[m,'acs']:.4f}")
        p = res["p_one_sided_top_tier_higher"]
        print(f"  one-sided Mann-Whitney p = {p if p is not None else 'not computed'}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out}")
    print("\nEXPLORATORY. Agreement does not prove ACS correct — both could err in the "
          "same direction. It shows ACS orders methods compatibly with benchmarks that "
          "used ground truth on other tissue.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

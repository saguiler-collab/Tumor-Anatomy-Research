#!/usr/bin/env python3
"""
reconstruct_gene_space.py — recover the leaderboard's exact 657-gene space, and prove it.

WHY THIS IS NEEDED
------------------
`docs/OPEN_DEFECTS.md` D10: the run persisted only the gene COUNT (657), never the gene
list. So `scripts/remeasure_method.py` had nothing to read, silently fell through to the
pipeline's FALLBACK gene space (1,591 genes), and both genuine-package re-measurements —
DWLS and BayesPrism — were made on a gene space no leaderboard method ever saw. The
comparisons drawn from them are withdrawn.

`run_all.py` now persists the list, but closing D10 by that route costs a full pipeline run,
whose benchmark stage alone took 34,721 s. This recovers the same list from the artefacts and
the atlas in minutes, because every step that produced it is deterministic:

    build_from_h5ad(atlas)                      seeded cell sampling
      -> split_donors(meta, seed=RANDOM_SEED)   seeded permutation
      -> build_reference(training cells)        deterministic
      -> select_signature_genes(ref, 80)        deterministic ranking + config markers
      -> intersect with the Ivy GAP bulk        deterministic

WHY THE RESULT CAN BE TRUSTED — TWO INDEPENDENT CHECKS
------------------------------------------------------
A reconstruction that cannot be checked is a guess, and a guess here is exactly what caused
D10. Two recorded facts gate it, and the script REFUSES to write if either fails:

  1. **The donor split must reproduce.** `method_selection_decision.json` records all 88
     training and 22 test donors by name. The split is a seeded permutation of the donor
     list, so reproducing it byte-for-byte confirms the atlas loaded to the same cells and
     the same seed took effect. This check is what makes the rest credible.
  2. **The gene count must be exactly 657.** Recorded independently in
     `equal_footing_certificate.json` (`n_genes`) and `method_selection_decision.json`
     (`n_signature_genes`). Eight types x 80 genes plus config markers, de-duplicated and
     intersected with the bulk, has no reason to land on 657 unless the chain matched.

The artefact it writes is labelled `provenance: RECONSTRUCTED`, not captured, with both
checks recorded, so nothing downstream can mistake it for a value the run itself emitted.

    python scripts/reconstruct_gene_space.py            # verify and write
    python scripts/reconstruct_gene_space.py --dry-run   # verify only
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                          # noqa: E402
from ivygap.bench.pseudobulk import split_donors                   # noqa: E402
from ivygap.data.load_ivygap import load_cached                    # noqa: E402
from ivygap.data.reference import (build_from_h5ad, build_reference,  # noqa: E402
                                   select_signature_genes)

#: What `run_all.py` passes. Not config.SIGNATURE_GENES_PER_TYPE, which is the default for
#: the frozen-reference path and would give a different space.
N_PER_TYPE_IN_RUN_ALL = 80


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="verify without writing")
    ap.add_argument("--out", default="results/benchmark/signature_genes.json")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    sel = json.loads((root / "results/benchmark/method_selection_decision.json").read_text())
    cert = json.loads((root / "results/anatomic/equal_footing_certificate.json").read_text())
    want_genes = int(cert["input_hashes"]["n_genes"])
    want_count_2 = int(sel["n_signature_genes"])
    if want_genes != want_count_2:
        print(f"ABORT: the two recorded gene counts disagree ({want_genes} vs "
              f"{want_count_2}). Resolve that before reconstructing anything.")
        return 1
    print(f"target: {want_genes} genes (agreed by both artefacts)")

    print("loading the Ivy GAP bulk ...")
    expr, meta_bulk = load_cached()

    print("loading the atlas (this is the slow step) ...")
    _, sc_expression, sc_meta = build_from_h5ad(config.REFERENCE_DIR / "gbmap_core.h5ad")
    print(f"  {sc_expression.shape[1]:,} cells x {sc_expression.shape[0]:,} genes")

    # --- CHECK 1: the donor split must reproduce the recorded one ------------
    train_donors, test_donors = split_donors(sc_meta, seed=config.RANDOM_SEED)
    rec_train = sorted(map(str, sel["train_donors"]))
    rec_test = sorted(map(str, sel["test_donors"]))
    ok_train = sorted(map(str, train_donors)) == rec_train
    ok_test = sorted(map(str, test_donors)) == rec_test
    print(f"\ncheck 1 — donor split reproduces: "
          f"train {len(train_donors)}/{len(rec_train)} {'OK' if ok_train else 'MISMATCH'}, "
          f"test {len(test_donors)}/{len(rec_test)} {'OK' if ok_test else 'MISMATCH'}")
    if not (ok_train and ok_test):
        extra = set(map(str, train_donors)) ^ set(rec_train)
        print(f"ABORT: the donor split does not reproduce. {len(extra)} donor(s) differ, "
              f"e.g. {sorted(extra)[:5]}.\nThe atlas did not load to the same cells, so "
              f"nothing downstream can be trusted. Nothing was written.")
        return 1

    # --- rebuild the reference from training donors only --------------------
    train_cells = sc_meta.index[sc_meta["donor"].astype(str).isin(train_donors)]
    print(f"\nbuilding the reference from {len(train_cells):,} training cells ...")
    reference = build_reference(sc_expression[train_cells], sc_meta.loc[train_cells],
                                name=config.PRIMARY_REFERENCE)

    genes = [g for g in select_signature_genes(reference, n_per_type=N_PER_TYPE_IN_RUN_ALL)
             if g in sc_expression.index]
    print(f"  select_signature_genes(n_per_type={N_PER_TYPE_IN_RUN_ALL}) -> {len(genes)} "
          f"genes before the bulk intersection")
    final = [g for g in genes if g in set(expr.index)]

    # --- CHECK 2: the count must match exactly ------------------------------
    print(f"\ncheck 2 — gene count: {len(final)} vs recorded {want_genes} "
          f"{'OK' if len(final) == want_genes else 'MISMATCH'}")
    if len(final) != want_genes:
        print(f"ABORT: reconstructed {len(final)} genes, the run recorded {want_genes}. "
              f"Off by {len(final) - want_genes}. A gene space that does not reproduce the "
              f"recorded count is not the leaderboard's, and using it would repeat D10 with "
              f"more confidence. Nothing was written.")
        return 1

    payload = {
        "what_this_is": ("The exact gene space every method on the leaderboard received, "
                         "RECONSTRUCTED from the atlas and the run's artefacts rather than "
                         "captured by the run itself. See docs/OPEN_DEFECTS.md D10."),
        "provenance": "RECONSTRUCTED",
        "reconstructed_utc": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc).isoformat(timespec="seconds"),
        "verified_by": [
            f"the donor split reproduces the recorded {len(rec_train)} training and "
            f"{len(rec_test)} test donors exactly, by name",
            f"the gene count reproduces the {want_genes} recorded independently in "
            f"equal_footing_certificate.json and method_selection_decision.json",
        ],
        "chain": ("build_from_h5ad(gbmap_core.h5ad) -> split_donors(seed="
                  f"{config.RANDOM_SEED}) -> build_reference(training cells) -> "
                  f"select_signature_genes(n_per_type={N_PER_TYPE_IN_RUN_ALL}) -> "
                  "intersect with the Ivy GAP bulk"),
        "n_genes": len(final),
        "sha256": config.sha256_strings(final),
        "genes": final,
    }
    if args.dry_run:
        print(f"\nDRY RUN — both checks passed. hash {payload['sha256'][:16]}…")
        return 0
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    print(f"\nboth checks passed. wrote {out}")
    print(f"  hash {payload['sha256'][:16]}…")
    print("\nD10 is now unblocked: scripts/remeasure_method.py will read this and run.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

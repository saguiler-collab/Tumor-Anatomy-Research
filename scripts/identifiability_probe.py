"""Is the B-over-T inversion a property of the SIGNATURE, or of applying it to real tissue?

WHY THIS EXISTS. `docs/IMMUNE_ARM.md` §6 established that all 12 methods place B cells above T
cells on TCGA tissue while DNA methylation places T above B in 93.2% of samples. It then
characterised that as "an identifiability failure ... a property of the problem as posed".

That characterisation is a CLAIM ABOUT THE SIGNATURE MATRIX, and it is testable without any
tissue at all: build mixtures FROM the reference, with a known T:B ratio, and see whether a
solver recovers it. If the signature genuinely cannot separate T from B, it will fail here too.

EXPLORATORY. Nothing in this script was pre-registered. It was written to test an explanation
this project had already published, and it refuted it.

    python scripts/identifiability_probe.py
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from ivygap import config
from ivygap.data.reference import load_frozen_reference, select_signature_genes

# Chosen to match the TISSUE, not to produce a result: methylation puts T:B at roughly 2.3
# within lymphocytes (results/methylation_celltypes_lgg.json), and glioma is macrophage-rich.
TRUTH = {"Tumor": 0.55, "Macrophage_Microglia": 0.25, "T_cell": 0.070, "NK_cell": 0.047,
         "B_cell": 0.030, "Oligodendrocyte": 0.03, "Endothelial": 0.015, "Astrocyte": 0.008}
NOISE_LEVELS = (0.0, 0.01, 0.05, 0.10, 0.25, 0.50, 1.00)
N_DRAWS = 200


def main() -> int:
    ref = load_frozen_reference()
    genes = select_signature_genes(ref, n_per_type=config.SIGNATURE_GENES_PER_TYPE)
    S = ref.profile.loc[genes]
    S = S / S.sum(axis=0)
    types = list(S.columns)
    A = S.to_numpy(dtype="float64")
    truth = pd.Series(TRUTH)[types]
    truth = truth / truth.sum()
    b0 = A @ truth.to_numpy()
    cond = float(np.linalg.cond(A))
    tb_truth = float(truth["T_cell"] / truth["B_cell"])

    print("IDENTIFIABILITY PROBE (exploratory, not pre-registered)")
    print("Mixtures are built FROM the reference, so any failure here is the signature's,")
    print("not the tissue's, not the solver implementation's.\n")
    print(f"planted truth: T {truth['T_cell']:.4f}  B {truth['B_cell']:.4f}  T:B {tb_truth:.2f}")
    print(f"condition number of the 8-type signature: {cond:.1f}\n")

    rng = np.random.default_rng(config.RANDOM_SEED)
    noise_rows = {}
    print(f"{'noise (CV)':>11s}  {'est T':>8s} {'est B':>8s} {'T:B':>7s}  verdict")
    print("-" * 54)
    for cv in NOISE_LEVELS:
        draws = []
        for _ in range(N_DRAWS if cv else 1):
            b = np.clip(b0 * (1 + rng.normal(0, cv, b0.size)), 0, None) if cv else b0
            x, _ = nnls(A, b)
            draws.append(x / x.sum())
        m = pd.DataFrame(draws, columns=types).mean()
        r = float(m["T_cell"] / m["B_cell"]) if m["B_cell"] > 1e-9 else float("inf")
        noise_rows[str(cv)] = {"T": round(float(m["T_cell"]), 4),
                               "B": round(float(m["B_cell"]), 4),
                               "ratio": round(r, 3), "T_exceeds_B": bool(m["T_cell"] > m["B_cell"])}
        print(f"{cv:11.2f}  {m['T_cell']:8.4f} {m['B_cell']:8.4f} {r:7.2f}  "
              f"{'T>B  recovered' if m['T_cell'] > m['B_cell'] else 'B>T  INVERTED'}")

    # A cell type present in tissue but absent from the reference is the standard explanation
    # for misassigned signal, so it is tested rather than invoked.
    print(f"\nMISSING-CELL-TYPE TEST — mixture has all 8; solver is given 7.")
    print(f"{'type withheld':24s} {'est T':>8s} {'est B':>8s} {'T:B':>7s}  verdict")
    print("-" * 62)
    drop_rows = {}
    for drop in types:
        if drop in ("T_cell", "B_cell"):
            continue
        keep = [t for t in types if t != drop]
        x, _ = nnls(S[keep].to_numpy(dtype="float64"), b0)
        m = pd.Series(x / x.sum(), index=keep)
        r = float(m["T_cell"] / m["B_cell"]) if m["B_cell"] > 1e-9 else float("inf")
        drop_rows[drop] = {"T": round(float(m["T_cell"]), 4), "B": round(float(m["B_cell"]), 4),
                           "ratio": round(r, 3), "T_exceeds_B": bool(m["T_cell"] > m["B_cell"])}
        print(f"{drop:24s} {m['T_cell']:8.4f} {m['B_cell']:8.4f} {r:7.2f}  "
              f"{'T>B' if m['T_cell'] > m['B_cell'] else 'B>T  INVERTED'}")

    held = all(v["T_exceeds_B"] for v in noise_rows.values()) and \
        all(v["T_exceeds_B"] for v in drop_rows.values())
    print("\nVERDICT:", "the signature separates T from B under the model in EVERY condition "
          "tested" if held else "the signature fails to separate T from B in at least one case")
    if held:
        print("  => the tissue inversion is NOT explained by the signature's conditioning,")
        print("     by noise, or by a missing cell type. Its cause is not identified.")

    out = {"EXPLORATORY": "not pre-registered; written to test an explanation this project "
                          "had already published, and it refuted it.",
           "planted_truth": {k: round(float(v), 4) for k, v in truth.items()},
           "planted_T_over_B": round(tb_truth, 3),
           "condition_number": round(cond, 2),
           "noise_sweep": noise_rows, "withheld_type": drop_rows,
           "signature_separates_T_from_B_in_all_conditions": bool(held)}
    (config.RESULTS_DIR / "identifiability_probe.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/identifiability_probe.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

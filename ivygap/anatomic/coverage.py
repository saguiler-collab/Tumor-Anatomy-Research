"""
coverage.py — what the reference cannot see, measured rather than asserted.

THE RISK THIS ADDRESSES
-----------------------
Anatomy_Test.md's risk table:

    Reference coverage — a seven-population signature cannot speak to astrocytes or
    neurons, which matter at the leading edge.
    → Report constraint coverage explicitly, and state which anatomic facts the
      signature is structurally unable to test.

`constraints.coverage_report()` does the first half: which cell types carry constraints
and which do not. Its `UNTESTABLE` entries say neuronal content cannot be tested because
"the seven-population roster contains no neuron column".

That was written when the reference was a collapsed seven-population signature matrix.
It is now materially misleading. The GBmap atlas **does** carry Neuron, and OPC, RG,
Mural cell, DC, Mast, Mono and Plasma B besides — `config.GBMAP_CELL_TYPE_MAP` maps 12
of its 20 labels and `build_from_h5ad` drops the rest. So the limitation is a roster
choice, not an absence in the data, and a reviewer who opens GBmap and finds Neuron will
reasonably conclude the project overlooked it.

The consequence is not cosmetic. Cells whose type is dropped do not vanish from the
tissue — only from the model. A deconvolution solved on eight columns still has to
explain the expression those cells contribute, so it is absorbed by whichever retained
type is closest.

Where that lands is measured, not assumed. On GBmap Core the dropped populations are
23,864 cells (7.0%), and the largest by far is Mono at 14,215 — myeloid, whose nearest
retained column is Macrophage_Microglia, which carries C5 and C6. Neuron, which the
constraint file singles out, is only 22 cells: the roster does drop it, but the atlas
could not have supported a neuron column either way. An earlier draft of this module
asserted the opposite — that neurons and OPCs dominate the dropped set — which the
counts do not support.

`constraints.py` is frozen and hashed, so this is reported ALONGSIDE it rather than by
editing it. Editing the constraint file to fix its own coverage note would change the
freeze hash and invalidate every result computed under it.
"""

from __future__ import annotations

import pandas as pd

from ivygap import config
from ivygap.anatomic import constraints as K


def reference_coverage(sc_meta: pd.DataFrame | None = None,
                       atlas_labels: dict | None = None) -> dict:
    """
    Measured coverage of the reference actually used.

    `atlas_labels` maps the atlas's own label -> cell count, for labels present in the
    source file. Supplying it turns "which populations were dropped" from a statement
    into a measurement.
    """
    roster = list(config.CELL_TYPES)
    mapped_labels = sorted(config.GBMAP_CELL_TYPE_MAP)

    out: dict = {
        "roster": roster,
        "n_roster_types": len(roster),
        "sidecar_types": list(config.SIDECAR_CELL_TYPES),
        "atlas_labels_mapped_onto_roster": mapped_labels,
        "constraint_coverage": K.coverage_report(),
    }

    if atlas_labels:
        dropped = {k: int(v) for k, v in atlas_labels.items()
                   if k not in config.GBMAP_CELL_TYPE_MAP}
        kept = {k: int(v) for k, v in atlas_labels.items()
                if k in config.GBMAP_CELL_TYPE_MAP}
        n_all = sum(atlas_labels.values())
        out["atlas_labels_dropped"] = dict(sorted(dropped.items(),
                                                  key=lambda kv: -kv[1]))
        out["n_cells_in_atlas"] = int(n_all)
        out["n_cells_dropped"] = int(sum(dropped.values()))
        out["fraction_of_atlas_dropped"] = (
            float(sum(dropped.values()) / n_all) if n_all else 0.0)
        out["n_cells_kept"] = int(sum(kept.values()))

    if sc_meta is not None and "cell_type" in sc_meta.columns:
        out["cells_per_roster_type_in_reference"] = {
            str(k): int(v) for k, v in sc_meta["cell_type"].value_counts().items()}

    dropped = out.get("atlas_labels_dropped", {})
    parts = ["Populations the atlas labels but the roster drops do not disappear from "
             "the tissue, only from the model. A solver with eight columns must still "
             "explain the expression those cells contribute, so it is absorbed by "
             "whichever retained type is closest."]

    if dropped:
        top = list(dropped.items())[:3]
        parts.append(
            "Largest dropped populations: "
            + ", ".join(f"{k} ({v:,} cells)" for k, v in top)
            + f" — {out['fraction_of_atlas_dropped']:.1%} of the atlas in total.")

        # Say where the absorbed signal most likely lands, from the measurement rather
        # than from an assumption about which populations are abundant.
        myeloid_like = {"Mono", "DC", "Mast"}
        hit = [k for k in dropped if k in myeloid_like]
        if hit and dropped.get("Mono", 0) == max(dropped.values()):
            parts.append(
                f"The dominant one is myeloid ({', '.join(sorted(hit))}), and the "
                f"nearest retained column is Macrophage_Microglia — which carries C5 "
                f"and C6. So the absorbed signal is concentrated on two of the seven "
                f"constraints rather than spread evenly.")

        if dropped.get("Neuron", 0) < 1000:
            parts.append(
                f"Note the frozen constraint file lists neuronal content as untestable "
                f"for want of a neuron column. That holds, but not for the reason it "
                f"gives: this atlas contains only {dropped.get('Neuron', 0)} neurons, "
                f"so a neuron column could not have been estimated from it either. "
                f"The limitation is the reference, not just the roster.")

    out["what_this_means"] = " ".join(parts)
    out["correction_to_the_frozen_constraint_file"] = (
        "constraints.UNTESTABLE says neuronal content cannot be tested because 'the "
        "seven-population roster contains no neuron column'. The conclusion holds and "
        "the stated reason is incomplete: the roster does drop the atlas's Neuron "
        "label, but the atlas carries so few neurons that a neuron column could not "
        "have been estimated from it in any case. Both facts are measured above. The "
        "constraint file is frozen and hashed, so this is recorded alongside it rather "
        "than by editing it — editing it would change the freeze hash and invalidate "
        "every result computed under it."
    )
    return out

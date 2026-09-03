"""
constraints.py — THE CONSTRAINT FILE.

This is step 2 of the build sequence in Anatomy_Test.md, and the step the entire design
rests on. Every constraint below is an ordinal prediction about cell composition across
anatomic structures, written from neuropathology BEFORE any deconvolution output was
looked at, and traceable to histology or in-situ hybridization — never to a
deconvolution result, which is exactly what would make the test circular.

`registration_payload()` emits the canonical JSON and `freeze_hash()` its SHA-256. Commit
the hash, register it publicly (OSF), and keep the registration timestamp ahead of every
result file. That timestamp is what converts "we predicted this" from a claim into a
checkable fact.

WHAT A CONSTRAINT MAY AND MAY NOT SAY
-------------------------------------
May: an ordering between two structures within one tumor, or that one structure holds
the maximum, or that three structures are monotone.

May not: any absolute fraction. Microdissected regions are still mixtures — the leading
edge is not pure brain and cellular tumor is not pure tumor — so a target value would be
indefensible. Ordering is what histology can honestly assert, and it is all that is
scored.

THE EXCLUSION IS PART OF THE PRE-REGISTRATION
---------------------------------------------
T cells are declared out of scope here, in the file, before any result. The pipeline's
own synthetic benchmark puts T cells at the detection floor, so a T-cell constraint would
score noise and would inflate or deflate ACS at random. Excluding it now is honest;
excluding it after seeing which methods it hurt would not be.

REFERENCE COVERAGE IS REPORTED, NOT HIDDEN
------------------------------------------
A seven-population signature cannot speak to neurons, and treats astrocytes as a
diagnostic sidecar known to be unidentifiable against tumor. Several true anatomic facts
about the leading edge are therefore structurally untestable here.
`coverage_report()` states which, so the constraint set's silence is visible rather than
mistaken for the biology being absent.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, asdict, field

from ivygap import config

#: Structures, periphery to core. MVP/PAN are the spec's names for the Ivy GAP
#: acronyms CTmvp/CTpan; `config.canonical_structure` maps between them so a constraint
#: written in the spec's vocabulary still matches the archive's labels.
LE, IT, CT, MVP, PAN = "LE", "IT", "CT", "MVP", "PAN"


@dataclass(frozen=True)
class Constraint:
    """
    One pre-registered ordinal prediction.

    kind:
        "pairwise"  — `cell_type` is higher in `structures[0]` than in `structures[1]`
        "maximum"   — `cell_type` is highest in `structures[0]` among `among`
        "monotone"  — `cell_type` increases along `structures` in order
    weight:
        Contribution to ACS. Not a confidence rating and not tunable after the fact:
        the monotone chain carries more because it is a conjunction of two steps and is
        correspondingly harder to satisfy by chance.
    """
    id: str
    kind: str
    cell_type: str
    structures: tuple[str, ...]
    rationale: str
    evidence: str
    weight: float = 1.0
    among: tuple[str, ...] = ()

    def describe(self) -> str:
        if self.kind == "pairwise":
            return f"{self.cell_type}: {self.structures[0]} > {self.structures[1]}"
        if self.kind == "maximum":
            return f"{self.cell_type}: {self.structures[0]} is the maximum"
        return f"{self.cell_type}: " + " < ".join(self.structures)


# =============================================================================
# THE CONSTRAINTS — frozen. Editing any line here changes freeze_hash() and
# invalidates every result computed under the old hash.
# =============================================================================

CONSTRAINTS: tuple[Constraint, ...] = (
    Constraint(
        id="C1", kind="pairwise", cell_type="Tumor", structures=(CT, LE),
        rationale="Near-definitional: the leading edge is brain infiltrated at low "
                  "tumour-cell density; cellular tumour is dense tumour. A method that "
                  "fails C1 is not measuring tumour.",
        evidence="Ivy GAP anatomic sampling criteria; H&E and ISH characterisation of "
                 "LE vs CT (Puchalski et al., Science 2018).",
        weight=1.0,
    ),
    Constraint(
        id="C2", kind="pairwise", cell_type="Oligodendrocyte", structures=(LE, CT),
        rationale="The leading edge retains normal white matter and cortex; cellular "
                  "tumour has displaced it.",
        evidence="Histological definition of the LE sampling region as infiltrated but "
                 "preserved brain parenchyma.",
        weight=1.0,
    ),
    Constraint(
        id="C3", kind="pairwise", cell_type="Endothelial", structures=(MVP, CT),
        rationale="Microvascular proliferation is defined by florid endothelial "
                  "proliferation. The single strongest constraint in the set.",
        evidence="WHO diagnostic criterion for MVP; Ivy GAP MVP dissection criterion.",
        weight=1.0,
    ),
    Constraint(
        id="C4", kind="maximum", cell_type="Endothelial", structures=(MVP,),
        among=(LE, IT, CT, MVP, PAN),
        rationale="Stricter form of C3: MVP must exceed all four other structures, not "
                  "just CT.",
        evidence="As C3. Scored separately because it is a materially harder claim and "
                 "a method can pass C3 while placing the endothelial maximum elsewhere.",
        weight=1.0,
    ),
    Constraint(
        id="C5", kind="pairwise", cell_type="Macrophage_Microglia", structures=(PAN, LE),
        rationale="The hypoxic peri-necrotic niche recruits myeloid cells; the leading "
                  "edge does not.",
        evidence="Hypoxia-driven myeloid recruitment in the pseudopalisading niche; "
                 "Ivy GAP ISH of myeloid markers by structure.",
        weight=1.0,
    ),
    Constraint(
        id="C6", kind="pairwise", cell_type="Macrophage_Microglia", structures=(MVP, CT),
        rationale="The perivascular niche is a documented myeloid reservoir.",
        evidence="Perivascular immune phenotype in GBM; CIBERSORTx analyses of Ivy GAP "
                 "reporting myeloid accumulation in the perivascular zone.",
        weight=1.0,
    ),
    Constraint(
        id="C7", kind="monotone", cell_type="Tumor", structures=(LE, IT, CT),
        rationale="A three-step monotone chain — much harder to satisfy by chance than "
                  "any single pairwise step, and scored as one unit.",
        evidence="Ivy GAP sampling definitions place IT between LE and CT in "
                 "tumour-cell density by construction.",
        weight=2.0,
    ),
)

#: Declared out of scope BEFORE any result. See the module docstring.
EXCLUSIONS: tuple[dict, ...] = (
    {
        "cell_type": "T_cell",
        "reason": "This pipeline's own synthetic benchmark puts T cells at the "
                  "detection floor, so a T-cell constraint would score noise.",
        "declared": "before any deconvolution output was examined",
    },
)

#: Anatomic facts the signature is structurally unable to test. Reported, not hidden.
UNTESTABLE: tuple[dict, ...] = (
    {"fact": "Neuronal content is highest at the leading edge",
     "why_untestable": "The seven-population roster contains no neuron column."},
    {"fact": "Non-neoplastic astrocytes are enriched at the leading edge",
     "why_untestable": "Astrocyte is a diagnostic-only sidecar, known to be "
                       "unidentifiable against Tumor; a constraint on it would score "
                       "the collinearity, not the biology."},
)


# =============================================================================
# ACCESS AND FREEZE
# =============================================================================

def by_id(cid: str) -> Constraint:
    for c in CONSTRAINTS:
        if c.id == cid:
            return c
    raise KeyError(f"no constraint {cid!r}; have {[c.id for c in CONSTRAINTS]}")


def total_weight() -> float:
    return float(sum(c.weight for c in CONSTRAINTS))


def as_records() -> list[dict]:
    return [asdict(c) for c in CONSTRAINTS]


def registration_payload() -> dict:
    """The canonical object to publish and timestamp."""
    return {
        "title": "The Anatomy Test — anatomic concordance constraint file",
        "constraints": as_records(),
        "exclusions": list(EXCLUSIONS),
        "untestable_by_this_signature": list(UNTESTABLE),
        "structures": [LE, IT, CT, MVP, PAN],
        "cell_type_roster": list(config.CELL_TYPES),
        "scoring": {
            "unit": "constraint x tumour pair",
            "aggregation": "weighted proportion satisfied, paired within tumour",
            "null": "structure labels permuted within each tumour",
            "interval": "bootstrap over tumours",
        },
    }


def freeze_hash() -> str:
    """SHA-256 of the registration payload. Register this before any result exists."""
    return hashlib.sha256(
        json.dumps(registration_payload(), sort_keys=True).encode()
    ).hexdigest()


def coverage_report() -> dict:
    """Which cell types the constraint set does and does not speak to."""
    covered = sorted({c.cell_type for c in CONSTRAINTS})
    return {
        "cell_types_with_constraints": covered,
        "cell_types_without_constraints": [
            t for t in config.CELL_TYPES if t not in covered],
        "explicitly_excluded": [e["cell_type"] for e in EXCLUSIONS],
        "structurally_untestable": list(UNTESTABLE),
        "n_constraints": len(CONSTRAINTS),
        "total_weight": total_weight(),
    }


def validate() -> None:
    """Structural checks. Run at import so a malformed constraint cannot reach scoring."""
    known = {LE, IT, CT, MVP, PAN}
    seen_ids = set()
    excluded = {e["cell_type"] for e in EXCLUSIONS}

    for c in CONSTRAINTS:
        if c.id in seen_ids:
            raise ValueError(f"duplicate constraint id {c.id!r}")
        seen_ids.add(c.id)
        if c.cell_type not in config.CELL_TYPES:
            raise ValueError(f"{c.id}: unknown cell type {c.cell_type!r}")
        if c.cell_type in excluded:
            raise ValueError(
                f"{c.id} scores {c.cell_type!r}, which EXCLUSIONS declares out of scope"
            )
        if c.cell_type in config.SIDECAR_CELL_TYPES:
            raise ValueError(
                f"{c.id} scores the sidecar type {c.cell_type!r}, which is known to be "
                f"unidentifiable; it would score collinearity, not biology"
            )
        bad = set(c.structures) - known
        if bad:
            raise ValueError(f"{c.id}: unknown structure(s) {sorted(bad)}")
        if c.weight <= 0:
            raise ValueError(f"{c.id}: weight must be positive")

        if c.kind == "pairwise":
            if len(c.structures) != 2 or c.structures[0] == c.structures[1]:
                raise ValueError(f"{c.id}: pairwise needs two distinct structures")
        elif c.kind == "maximum":
            if len(c.structures) != 1:
                raise ValueError(f"{c.id}: maximum names exactly one structure")
            if len(c.among) < 3:
                raise ValueError(f"{c.id}: maximum needs a comparison set of >= 3")
            if c.structures[0] not in c.among:
                raise ValueError(f"{c.id}: the maximum must be inside `among`")
            if set(c.among) - known:
                raise ValueError(f"{c.id}: unknown structure in `among`")
        elif c.kind == "monotone":
            if len(c.structures) < 3:
                raise ValueError(f"{c.id}: a monotone chain needs >= 3 structures")
            if len(set(c.structures)) != len(c.structures):
                raise ValueError(f"{c.id}: monotone chain repeats a structure")
        else:
            raise ValueError(f"{c.id}: unknown kind {c.kind!r}")

    if not CONSTRAINTS:
        raise ValueError("the constraint set is empty; ACS would be undefined")


validate()

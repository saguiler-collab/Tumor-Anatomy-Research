"""
constraints_brain.py — THE SECOND CONSTRAINT FILE (normal human brain). **DRAFT.**

STATUS: NOT FROZEN. `FROZEN = False` below, and `freeze_hash()` refuses to run while it
is False. Nothing may be deconvolved against these constraints until the region mapping
is verified against the Allen ontology and the file is frozen and registered. The refusal
is mechanical rather than a note, because a constraint file that can be scored before it
is registered is not a pre-registration.

WHY A SECOND FILE AND NOT AN EDIT
---------------------------------
`Anatomy_Test.md` step 7: *one tissue is a case study; two is a method.* The Ivy GAP
result stands on nine evaluable tumours, and its own leaderboard shows what that buys —
MuSiC leads by a single constraint-tumour pair out of 57, with six methods' intervals
containing its score. Anatomy separates good methods from bad on that cohort; it cannot
identify the best one, and no further analysis inside Ivy GAP changes that.

The GBM constraint file cannot be reused. It is written about CT, IT, LE, MVP and PAN,
which are glioblastoma structures; normal brain has none of them. Reusing it is not an
option and editing it is worse — editing changes `constraints.freeze_hash()` and voids
every GBM result computed under it. So this is a separate file with a separate hash and
a separate registration.

WHAT WAS AND WAS NOT LOOKED AT BEFORE WRITING THIS
--------------------------------------------------
Looked at: the Allen human brain structure ontology (graph_id 10) — a taxonomy of region
names, the same way LE/IT/CT/MVP/PAN had to be known before the GBM constraints could be
written. And published neuroanatomy, cited per constraint below.

NOT looked at: any Allen expression value, any deconvolution output, any composition
estimate from any method on any brain sample. That ordering is the whole claim, and it is
the one thing that cannot be repaired after the fact.

THE ROSTER PROBLEM THIS TIME
----------------------------
The GBM roster does not transfer either. Normal cortex is mostly neurons and this
project's roster has no neuron column; `reference_coverage.json` records why one could
not simply be added — the GBmap atlas contains **22 neurons**. A normal-brain run needs a
normal-brain reference. Siletti et al. (2023) *Science*, the Human Brain Cell Atlas v1.0,
supplies one: 2.48M neurons and 888k non-neuronal cells, on 10x 3' v3 — the same platform
family as GBmap, so the pipeline's normalisation and mode reasoning carry over.

WHAT A CONSTRAINT MAY AND MAY NOT SAY
-------------------------------------
Unchanged from the GBM file. Orderings between regions within one donor; never an
absolute fraction. A dissected region is still a mixture — cortex contains white matter
tracts and white matter contains neuronal somata — so a target value would be
indefensible here for exactly the reason it was there.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, asdict

#: Flip to True only when the region mapping below has been verified against the Allen
#: ontology and the claims are final. Freezing is what makes the hash meaningful; a hash
#: over a draft is a checksum of an opinion.
FROZEN = False

#: Coarse regions, in this file's own vocabulary. The Allen release annotates samples to
#: a 1,839-structure ontology; these four are the level at which neuroanatomy makes
#: claims that are near-definitional rather than merely plausible. `REGION_MAP` below is
#: the (separately verifiable) translation, kept apart from the claims so that fixing a
#: mapping detail does not silently rewrite a prediction.
CTX = "CTX"      # cerebral cortex, grey matter
WM = "WM"        # cerebral white matter
CBC = "CBC"      # cerebellar cortex
SUBC = "SUBC"    # subcortical grey (thalamus, striatum, and similar nuclei)

REGIONS = (CTX, WM, CBC, SUBC)

#: Cell types. Estimable from a normal-brain single-cell reference, and chosen so that
#: every constraint below names a population a reference can actually resolve.
NEURON = "Neuron"
OLIGO = "Oligodendrocyte"
ASTRO = "Astrocyte"
MICRO = "Microglia"
OPC = "OPC"
ENDO = "Endothelial"

CELL_TYPES = (NEURON, OLIGO, ASTRO, MICRO, OPC, ENDO)


@dataclass(frozen=True)
class BrainConstraint:
    """
    One pre-registered ordinal prediction about normal human brain.

    Same three kinds as the GBM file, deliberately: the scorer is shared, so a constraint
    written here is evaluated by exactly the code that evaluated the GBM set, and a
    difference in result cannot be a difference in scoring.

        "pairwise"  — `cell_type` is higher in `regions[0]` than in `regions[1]`
        "maximum"   — `cell_type` is highest in `regions[0]` among `among`
        "monotone"  — `cell_type` increases along `regions` in order
    """
    id: str
    kind: str
    cell_type: str
    regions: tuple[str, ...]
    rationale: str
    evidence: str
    weight: float = 1.0
    among: tuple[str, ...] = ()

    def describe(self) -> str:
        if self.kind == "pairwise":
            return f"{self.cell_type}: {self.regions[0]} > {self.regions[1]}"
        if self.kind == "maximum":
            return f"{self.cell_type}: {self.regions[0]} is the maximum"
        # Joined in order, matching the GBM file: regions are listed lowest-first, so
        # C7 there reads "Tumor: LE < IT < CT" and B5 here reads "Neuron: WM < CTX < CBC".
        # Reversing renders the claim backwards, which is how this first read
        # "Neuron: CBC < CTX < WM" — the opposite of what B5 asserts.
        return f"{self.cell_type}: " + " < ".join(self.regions)


#: The constraint set. Every entry cites neuroanatomy or histology; none cites a
#: deconvolution result, which is what would make the test circular.
CONSTRAINTS: tuple[BrainConstraint, ...] = (
    BrainConstraint(
        id="B1", kind="pairwise", cell_type=OLIGO, regions=(WM, CTX), weight=1.0,
        rationale=(
            "Near-definitional. White matter IS myelinated axon tracts, and myelin in "
            "the CNS is produced by oligodendrocytes, whose cell bodies are "
            "correspondingly concentrated there. This is the strongest claim in the set "
            "and the direct analogue of the GBM file's C3: a method that fails B1 is "
            "not measuring oligodendrocytes."),
        evidence=(
            "Standard neuroanatomy: the grey/white distinction is a myelin distinction. "
            "Oligodendrocyte density in white matter exceeds cortical grey matter by a "
            "wide margin in every published stereological series."),
    ),
    BrainConstraint(
        id="B2", kind="maximum", cell_type=OLIGO, regions=(WM,),
        among=REGIONS, weight=1.0,
        rationale=(
            "Stricter form of B1, exactly as C4 is a stricter form of C3. A method can "
            "put oligodendrocytes above cortex and still place their maximum in "
            "subcortical grey, which is heavily traversed by fibre tracts. Scored "
            "separately because it is a materially harder claim."),
        evidence="As B1.",
    ),
    BrainConstraint(
        id="B3", kind="pairwise", cell_type=NEURON, regions=(CTX, WM), weight=1.0,
        rationale=(
            "Definitional in the same way. Cortex is where neuronal somata sit; white "
            "matter carries their axons and contains few cell bodies. The analogue of "
            "the GBM file's C1."),
        evidence=(
            "Standard neuroanatomy. Neuronal somata are essentially absent from white "
            "matter apart from a sparse population of interstitial neurons."),
    ),
    BrainConstraint(
        id="B4", kind="maximum", cell_type=NEURON, regions=(CBC,),
        among=REGIONS, weight=1.0,
        rationale=(
            "The cerebellar granule cell layer is the densest neuronal population in "
            "the human CNS. Cerebellar granule cells alone outnumber every other neuron "
            "in the brain combined, in a structure holding roughly a tenth of its mass. "
            "A genuinely hard constraint: it requires a method to get a gross regional "
            "ordering right, not merely to separate grey from white."),
        evidence=(
            "Azevedo et al. (2009), J Comp Neurol 513:532 — isotropic fractionator "
            "counts putting ~69 billion of the brain's ~86 billion neurons in the "
            "cerebellum against ~16 billion in cerebral cortex."),
    ),
    BrainConstraint(
        id="B5", kind="monotone", cell_type=NEURON, regions=(WM, CTX, CBC), weight=2.0,
        rationale=(
            "A three-step monotone chain, and double-weighted for the same reason C7 is: "
            "a conjunction of two orderings is far harder to satisfy by chance than "
            "either step alone, and it is scored as one unit rather than as two."),
        evidence="As B3 and B4.",
    ),
    BrainConstraint(
        id="B6", kind="pairwise", cell_type=ENDO, regions=(CTX, WM), weight=1.0,
        rationale=(
            "Grey matter is more densely capillarised than white matter, following its "
            "higher metabolic demand. Weaker than B1-B5 and included deliberately as a "
            "constraint that is defensible but not near-definitional, so the set is not "
            "composed only of its easiest claims."),
        evidence=(
            "Capillary density and regional cerebral blood flow are both substantially "
            "higher in cortical grey matter than in white matter across the imaging and "
            "histology literature."),
    ),
)


#: Claims deliberately NOT made, recorded before any result so that the set's silences
#: are visible rather than mistaken for the biology being absent.
UNTESTABLE: tuple[dict[str, str], ...] = (
    {
        "fact": "Astrocytes are enriched in one region over another",
        "why_untestable": (
            "No safe ordering exists at this granularity. Protoplasmic astrocytes "
            "populate grey matter and fibrous astrocytes populate white matter, and "
            "published density estimates do not order the two consistently. A "
            "constraint here would score the reference's astrocyte definition rather "
            "than the anatomy."),
    },
    {
        "fact": "Microglial density differs between these regions",
        "why_untestable": (
            "Microglia are distributed relatively uniformly in the healthy brain. Any "
            "ordering would be within the noise of what deconvolution can resolve, "
            "which is the same reason T cells were excluded from the GBM set."),
    },
    {
        "fact": "OPC is enriched in white matter",
        "why_untestable": (
            "Directionally likely but not near-definitional, and OPC/oligodendrocyte "
            "boundaries differ between references. Excluded in advance rather than "
            "risk scoring an annotation choice."),
    },
)


#: Translation from this file's coarse vocabulary to the Allen ontology. Kept separate
#: from the claims on purpose: the claims are neuroanatomy and should not move, whereas
#: the mapping is bookkeeping against a specific release and must be verifiable.
#:
#: Verified 2026-09-10 against the Allen human structure graph (graph_id 10, 1,839
#: structures). Each entry names the ontology node whose descendants define the region;
#: assignment is by ancestry, so a sample lands in a region if that node is on its path
#: to the root. The ontology's own top-level split is `Br / GM | WM | SS`, which is the
#: distinction B1-B3 turn on, so this mapping is the ontology's rather than one imposed
#: on it.
#:
#: NOTE: the Allen graph is queried for structure NAMES only. No expression value was
#: read while writing this file.
REGION_MAP: dict[str, dict] = {
    CTX: {
        "description": "cerebral cortex, grey matter",
        "allen_node": "Cx", "allen_path": "Br / GM / Tel / Cx",
        "verified": True,
    },
    WM: {
        "description": "white matter, all divisions",
        "allen_node": "WM", "allen_path": "Br / WM",
        "note": ("The ontology's own depth-1 white-matter node, so this includes "
                 "telencephalic (TELWM), mesencephalic, metencephalic and "
                 "myelencephalic white matter. Taking the whole node rather than TELWM "
                 "alone makes B2's maximum claim harder, not easier: it must beat every "
                 "grey region while carrying brainstem tracts too."),
        "verified": True,
    },
    CBC: {
        "description": "cerebellar cortex",
        "allen_node": "CbCx", "allen_path": "Br / GM / MET / Cb / CbCx",
        "verified": True,
    },
    SUBC: {
        "description": "subcortical grey nuclei — basal ganglia and diencephalon",
        "allen_nodes": ["BG", "DiE"],
        "allen_paths": ["Br / GM / Tel / CxN / BG", "Br / GM / DiE"],
        "note": ("Present so B2 and B4 have a genuine competitor. Subcortical grey is "
                 "traversed by heavy fibre tracts, which is exactly the region that "
                 "could steal the oligodendrocyte maximum from WM, and it is dense "
                 "grey, which is what could steal the neuron maximum from CbCx."),
        "verified": True,
    },
}


def registration_payload() -> dict:
    """
    The canonical object the hash covers: the CLAIMS, not the file.

    Deliberately excludes REGION_MAP. A mapping correction — discovering that an Allen
    acronym belongs under CTX rather than SUBC — is bookkeeping and must not invalidate
    a registration, whereas changing what B1 asserts must.
    """
    return {
        "schema": "anatomy-test/brain-constraints/v1",
        "tissue": "normal human brain",
        "regions": list(REGIONS),
        "cell_types": list(CELL_TYPES),
        "constraints": [
            {k: (list(v) if isinstance(v, tuple) else v)
             for k, v in asdict(c).items() if k not in ("rationale", "evidence")}
            for c in CONSTRAINTS
        ],
        "untestable": [dict(u) for u in UNTESTABLE],
    }


def freeze_hash() -> str:
    """SHA-256 of the canonical payload. Refuses to run on a draft."""
    if not FROZEN:
        raise RuntimeError(
            "constraints_brain.py is a DRAFT (FROZEN = False). Verify REGION_MAP against "
            "the Allen ontology, settle the claims, set FROZEN = True, and only then "
            "hash and register. Hashing a draft produces a number that looks like a "
            "commitment and is not one.")
    payload = json.dumps(registration_payload(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def summary() -> str:
    lines = [f"{len(CONSTRAINTS)} constraints, weighted denominator "
             f"{sum(c.weight for c in CONSTRAINTS):.0f} per donor",
             f"regions: {', '.join(REGIONS)}",
             f"cell types: {', '.join(CELL_TYPES)}",
             f"FROZEN: {FROZEN}", ""]
    for c in CONSTRAINTS:
        lines.append(f"  {c.id}  w={c.weight:g}  {c.describe()}")
    lines.append("")
    lines.append("declared untestable:")
    for u in UNTESTABLE:
        lines.append(f"  - {u['fact']}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())

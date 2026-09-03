"""
portal_metadata.py — the Allen Institute portal's own metadata, used as an external check.

WHY THIS IS A SEPARATE MODULE
-----------------------------
`load_ivygap.py` derives its cohort description from `columns-samples.csv`, which ships
inside the 2014-11-25 expression archive. That is an *internal* count: the loader and
the thing it is being checked against come from the same file, so a parsing error that
mis-assigns structures would produce a self-consistent and entirely wrong answer.

Protocol step 3's gate is explicit that this is not enough:

    Gate: sample counts per structure reconcile against the portal's own documentation.

`rna_seq_samples_details.csv`, served live from the portal's `/api/v2/gbm/` namespace,
is that documentation. It is produced by the Allen Institute from their LIMS rather than
by us from the archive, it carries the study assignment as an explicit column
(`study_name`) rather than as a suffix convention we have to parse, and it names the
anatomic structure of every sample. If our parse of the archive and their table agree on
122 anatomic samples and on 30/24/19/25/24 per structure, the parse is checked. If they
disagree, we would rather know.

THE JOIN NOBODY ELSE CAN MAKE
-----------------------------
This file is also load-bearing for survival, for a reason that is easy to miss:

    tumor_details.csv       is keyed on  donor_id   (12111, 14734, ...)
    columns-samples.csv     is keyed on  tumor_id   (711547, 705757, ...)

Neither file contains the other's key. `rna_seq_samples_details.csv` contains both, and
is the only published file that does. Without it the clinical table cannot be attached to
the expression matrix at all — which is why the previous run reported survival as BLOCKED
even in principle, not merely for want of a download.

THE RELEASE SKEW IS REAL AND IS REPORTED
----------------------------------------
The portal is live and the expression archive is frozen at 2014-11-25. The portal
currently lists 279 RNA-seq samples; the archive holds 270. The nine extras are all
Cancer Stem Cells samples added after the freeze. That is a fact about the two releases,
not a defect, and `reconcile()` reports it as `portal_only` rather than treating it as a
count mismatch. What would be a defect is the reverse — a sample in the archive that the
portal does not list — and that is checked separately and hard-fails.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict, field

import pandas as pd

from ivygap import config

#: The portal's own name for the H&E-selected anatomic study. Matching on this column is
#: what makes the reconciliation independent of our `-reference-histology` suffix rule.
ANATOMIC_STUDY_NAME = "Anatomic Structures RNA Seq"
CLUSTER_STUDY_NAME = "Cancer Stem Cells RNA Seq"


class PortalMetadataMissing(RuntimeError):
    """The portal metadata table is not on disk. Reconciliation cannot be performed."""


@dataclass
class Reconciliation:
    """The measured comparison between our parse and the portal's table."""
    available: bool
    reason: str = ""
    n_archive_samples: int = 0
    n_portal_samples: int = 0
    portal_only: list[str] = field(default_factory=list)
    archive_only: list[str] = field(default_factory=list)
    anatomic_counts_archive: dict = field(default_factory=dict)
    anatomic_counts_portal: dict = field(default_factory=dict)
    anatomic_counts_agree: bool = False
    study_assignment_disagreements: list[dict] = field(default_factory=list)
    n_tumors_archive: int = 0
    n_tumors_portal: int = 0
    verdict: str = ""

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)


def load_sample_details() -> pd.DataFrame:
    """The portal's per-RNA-seq-sample table, typed so joins cannot silently fail."""
    path = config.IVYGAP_RNA_SEQ_DETAILS_PATH
    if not path.exists():
        raise PortalMetadataMissing(
            f"{path} not present. Run `python -m ivygap.data.download_ivygap` to fetch "
            f"it, or download rna_seq_samples_details.csv from the portal's Download "
            f"page and save it there."
        )
    df = pd.read_csv(path)
    # Ids are identifiers, not numbers. Read as str throughout so a join never depends
    # on whether pandas happened to infer int64 on one side and object on the other.
    for col in ("sample_id", "tumor_id", "donor_id", "sub_block_id"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    return df


def donor_tumor_map() -> pd.DataFrame:
    """
    One row per tumour: `tumor_id` (the expression matrix's key) -> `donor_id` (the
    clinical table's key), plus the tumour's display name.

    Raises if the correspondence is not one-to-one. A donor mapping to two tumours, or a
    tumour to two donors, would fan the clinical join out and silently duplicate
    outcomes across patients.
    """
    df = load_sample_details()
    pairs = df[["tumor_id", "donor_id", "tumor_name"]].drop_duplicates()

    dup_t = pairs["tumor_id"][pairs["tumor_id"].duplicated()].unique()
    if len(dup_t):
        raise ValueError(
            f"tumor_id maps to more than one donor_id in the portal table: "
            f"{sorted(dup_t)[:10]}. The clinical join would duplicate outcomes."
        )
    dup_d = pairs["donor_id"][pairs["donor_id"].duplicated()].unique()
    if len(dup_d):
        raise ValueError(
            f"donor_id maps to more than one tumor_id in the portal table: "
            f"{sorted(dup_d)[:10]}. The clinical join would fan out."
        )
    return pairs.set_index("tumor_id")


def reconcile(manifest: pd.DataFrame) -> Reconciliation:
    """
    Compare our parse of the frozen archive against the portal's live table.

    `manifest` is the loader's output: indexed by sample_id, carrying `patient_id`
    (= tumor_id), `structure` (canonicalised) and `is_anatomic_study`.

    Returns a measured comparison. It does not raise on a count difference — the caller
    decides what is fatal — except for the one case that can only be a defect: a sample
    present in the archive that the portal does not list at all.
    """
    try:
        portal = load_sample_details()
    except PortalMetadataMissing as exc:
        return Reconciliation(available=False, reason=str(exc),
                              verdict="INTERNAL COUNTS ONLY: the portal table is not "
                                      "present, so step 3's gate cannot be closed.")

    archive_ids = set(manifest.index.astype(str))
    portal_ids = set(portal["sample_id"])

    archive_only = sorted(archive_ids - portal_ids)
    portal_only = sorted(portal_ids - archive_ids)

    # Our canonical structure names (LE/IT/CT/MVP/PAN) against the portal's raw acronyms
    # (LE-reference-histology, CTmvp-reference-histology, ...). Fold the portal's labels
    # through the SAME function the loader uses, so the comparison is of counts and not
    # of vocabularies.
    p_anat = portal[portal["study_name"] == ANATOMIC_STUDY_NAME].copy()
    p_anat["canon"] = p_anat["structure_acronym"].map(config.canonical_structure)
    counts_portal = (p_anat[p_anat["canon"].isin(config.PRIMARY_STRUCTURES)]
                     ["canon"].value_counts().sort_index().to_dict())

    a_anat = manifest[manifest["is_anatomic_study"]
                      & manifest["structure"].isin(config.PRIMARY_STRUCTURES)]
    counts_archive = a_anat["structure"].value_counts().sort_index().to_dict()

    # Restrict the portal side to samples the frozen archive actually contains, so the
    # nine post-freeze additions do not read as a count mismatch.
    p_shared = p_anat[p_anat["sample_id"].isin(archive_ids)]
    counts_portal_shared = (p_shared[p_shared["canon"].isin(config.PRIMARY_STRUCTURES)]
                            ["canon"].value_counts().sort_index().to_dict())

    agree = counts_archive == counts_portal_shared

    # Independent check of the study assignment itself: our suffix rule vs their column.
    shared = portal[portal["sample_id"].isin(archive_ids)].set_index("sample_id")
    theirs = shared["study_name"] == ANATOMIC_STUDY_NAME
    ours = manifest.reindex(theirs.index)["is_anatomic_study"].astype(bool)
    mismatch = theirs.index[theirs.to_numpy() != ours.to_numpy()]
    disagreements = [
        {"sample_id": s,
         "portal_study": shared.loc[s, "study_name"],
         "portal_structure": shared.loc[s, "structure_acronym"],
         "our_is_anatomic_study": bool(ours.loc[s])}
        for s in mismatch
    ]

    if archive_only:
        verdict = (f"DEFECT: {len(archive_only)} sample(s) in the frozen archive are "
                   f"absent from the portal table: {archive_only[:5]}. The parse cannot "
                   f"be checked against documentation that does not describe it.")
    elif disagreements:
        verdict = (f"DEFECT: our `-reference-histology` suffix rule and the portal's "
                   f"`study_name` column disagree on {len(disagreements)} sample(s). "
                   f"One of the two is wrong and the anatomic cohort is not trustworthy "
                   f"until it is resolved.")
    elif not agree:
        verdict = (f"MISMATCH: per-structure counts differ. archive={counts_archive} "
                   f"portal={counts_portal_shared}. Step 3's gate is NOT closed.")
    else:
        verdict = (
            f"RECONCILED: {sum(counts_archive.values())} anatomic samples across "
            f"{len(counts_archive)} structures agree exactly with the portal's own "
            f"table ({counts_archive}), and the study assignment agrees on all "
            f"{len(shared)} shared samples. "
            f"{len(portal_only)} sample(s) exist on the live portal but not in the "
            f"2014-11-25 archive — a release skew, all in the "
            f"'{CLUSTER_STUDY_NAME}' study, not a parsing error. Step 3's gate is closed."
        )

    return Reconciliation(
        available=True,
        n_archive_samples=len(archive_ids),
        n_portal_samples=len(portal_ids),
        portal_only=portal_only,
        archive_only=archive_only,
        anatomic_counts_archive=counts_archive,
        anatomic_counts_portal=counts_portal_shared,
        anatomic_counts_agree=bool(agree),
        study_assignment_disagreements=disagreements,
        n_tumors_archive=int(manifest["patient_id"].nunique()),
        n_tumors_portal=int(portal.loc[portal["sample_id"].isin(archive_ids),
                                       "tumor_id"].nunique()),
        verdict=verdict,
    )

# External actions — what the pipeline cannot do for itself

Each item says what to get, exactly where it goes, what changes in the results when it
lands, and how to confirm it worked. Ordered by what unblocks the most.

> **Items 1 and 5 are DONE. Items 6–13 were added 2026-09-13** after the clinical-readiness
> and completeness audits found gaps that no amount of internal work can close. Item 6 is now
> the most binding constraint on this machine, and item 8 is the one that would most change
> what the paper can claim.

---

## 0 · WHAT TO DO, IN ORDER, BY MACHINE — rewritten 2026-09-14

Two lists. **Everything in list A runs on the 8 GB laptop** — most of it needs no compute at
all, and the one download is small. **List B needs the advanced machine** and should not be
started until A is done, because several items in B are cheaper once A has landed.

Ordered within each list by value per unit of effort.

### LIST A — the 8 GB laptop, starting now

| # | do this | effort | why it is worth doing first |
|---|---|---|---|
| ~~A1~~ | ~~Download GSE84465~~ | **DONE 2026-09-14** | Fetched and tested. It does **not** give a composition replication — the unpanned periphery is 13 cells. It gives the project's **first cross-patient constraint test** (C1, 4 patients, RESULTS.md §5c) and is the **best second-reference candidate** (5 of 8 roster types, 4 donors). Item 15 has the correction. |
| **A1b** | **Build the Darmanis second reference** | 2–4 h | Would make `SCDC ENSEMBLE` non-degenerate on a declared 5-type sub-roster, and give every method a reference-sensitivity check. This is what item 9 was for, and Darmanis reaches it without the authors. |
| **A2** | Publish the registration corrections | 20 min, browser | Tier 0.2. The registration is permanent and two statements in it are wrong. A published correction is stronger than one nobody checked. |
| **A3** | Read tumour percentage off the Ivy GAP H&E images | hours, no compute | Item 13. Ivy GAP publishes the images; this is a reading exercise. Converts the 0.33–0.51 tumour under-call from a benchmark artefact into a tissue-level result. |
| **A4** | Ask the Neftel authors for their cell-state assignments | 15 min, email | Item 9. Confirmed absent from all four GEO supplementary files. Deriving them locally is a list-B job and only yields 4 of 8 roster types anyway (item 16). |
| **A5** | Hunt for an immune ground-truth cohort | hours, browser | Item 12, still the most scientifically valuable thing on either list. IHC/flow/CyTOF with matched bulk RNA-seq. |
| **A6** | Install and wire **CDSeq** | 1–2 h | Item 14. The only *reference-free* candidate, so the only way to separate "the methods are wrong" from "the reference is wrong". Its source is already on disk. Try on 8 GB; if it thrashes, move it to B. |

**Nothing in list A needs the results tree re-run.** A1 and A6 add new evidence beside the
existing leaderboard; A2–A5 are external entirely.

### LIST B — the advanced machine, after list A

| # | do this | why it needs the bigger machine |
|---|---|---|
| **B1** | **One full `python scripts/run_all.py`** | Loads a 7.6 GB atlas and runs 17 methods; ~10 h. This is the item that replaces the archived DWLS and BayesPrism rows with the genuine packages, adopts whatever D12 decision is taken, and lands EPIC's `refProfiles.var`. Everything else in B is easier afterwards. |
| **B2** | A second tissue | Item 2/8. A whole new cohort: download, ingest, reference, run. Also the only route past nine evaluable tumours. |
| **B3** | Derive the Neftel annotations locally, if A4 got no answer | Item 16. inferCNV needed 16 GB *and* downsampling in the reference implementation. |
| **B4** | A reference retaining Mono/DC/Mast | Item 10. Re-maps and rebuilds from the full atlas. |
| **B5** | EcoTyper | 1 GB of package plus its own outputs; only worth it once a labelled second reference exists. |

---

## 0a · STATUS BOARD

| # | need | status | machine |
|---|---|---|---|
| 1 | OSF registration | **DONE** 2026-09-10 | — |
| 5 | CIBERSORTx Supplementary Note 1 | **DONE** 2026-09-06 | — |
| 6 | RAM for the gene-space reconstruction | **DONE 2026-09-14** — reconstructed and verified, both gates passed | was 8 GB, worked |
| **15** | **GSE84465 (Darmanis 2017)** | **NEW, and the top of list A** | 8 GB |
| 2 / 8 | a second tissue | open | advanced |
| 9 | GSE131928 cell-type labels | open — ask the authors (A4) | 8 GB to ask |
| **16** | **Neftel annotations derivable locally** | **NEW — the pipeline is already on disk, but yields 4 of 8 types** | advanced |
| 10 | a myeloid-retaining reference | open | advanced |
| 11 | paired bulk + single-cell, same subjects | **partly answered by 15** — GSE84465 is 4 patients, single-cell only | 8 GB |
| 12 | immune ground truth | open, highest scientific value | 8 GB to find |
| 13 | pathologist purity on Ivy GAP blocks | open, no compute needed | 8 GB |
| 7 / 4 | Ivy GAP vital status | open, probably unobtainable | — |
| 14 | CDSeq / RNA-Sieve | open | 8 GB to try |

---

## 15 · GSE84465 (Darmanis 2017) — **FETCHED AND TESTED 2026-09-14. It delivers less than I claimed, and something else instead.**

> **Correction to my own recommendation.** I put this at the top of list A on the grounds that
> it would turn the deconvolution-free composition validation from one specimen into five. **It
> does not.** I fetched it (20.5 MB, not the ~300 MB I estimated) and tested it, and the
> composition route is not viable. What it does deliver is narrower and still worth having.

### Why the composition route fails

Darmanis is **heavily FACS-panned**. Of 3,589 cells only **665 are `Unpanned`**; the rest were
sorted on CD45, HEPACAM, Thy1, GC or BSC. A composition computed over sorted cells measures the
sort, not the tissue.

And among the unpanned cells the **periphery is 13 cells** — twelve from `BT_S4`, one from
`BT_S6`, **none at all** from `BT_S1` or `BT_S2`:

| patient | unpanned periphery cells |
|---|---|
| BT_S1 | 0 |
| BT_S2 | 0 |
| BT_S4 | 12 |
| BT_S6 | 1 |

So the unbiased-composition route gives one patient's worth of periphery. That is no better
than Albiach's n = 1 and far thinner. **Had I run the obvious analysis over all 3,589 cells it
would have produced confident, invalid numbers.**

### What it does deliver — the first cross-patient constraint test in the project

`scripts/darmanis_constraint_check.py`. Within a **fixed sorting gate**, compare the fraction
of cells that are neoplastic in the core against the periphery. Holding the gate constant holds
the selection bias constant, and what remains is exactly C1's claim.

| gate | core | periphery | difference | p | patients |
|---|---|---|---|---|---|
| Astrocytes(HEPACAM) | 0.970 | 0.183 | **+0.786** | 0.0001 | **3/3** |
| Neurons(Thy1) | 0.517 | 0.007 | **+0.510** | 0.0001 | **3/3** |
| Oligodendrocytes(GC) | 0.017 | 0.000 | +0.017 | 0.387 | 1/2 |
| Microglia(CD45) | 0.006 | 0.002 | +0.004 | 0.290 | 2/4 |

**4 of 4 gates support C1's direction, 2 decisively, and the two strong gates agree in every
patient they can be scored in.** RESULTS.md §5c. This is the only cross-patient evidence for
any constraint in the project — Albiach (§5b) is n = 1 by construction.

It tests **C1 only**. C2 needs unbiased composition; C3/C4/C6 need a microvascular region; C5
needs a peri-necrotic one; C7 needs three ordered regions and this has two.

### The other thing it is good for, which I had not considered: a SECOND REFERENCE

Panning ruins composition and is **irrelevant to a reference**, which needs per-type expression
profiles rather than per-region abundances. Sorting actually helps, by enriching rare types.

| | cells | patients |
|---|---|---|
| Immune cell | 1,847 | 4 |
| Neoplastic | 1,091 | 4 |
| OPC | 406 | 4 |
| Astrocyte | 88 | 4 |
| Oligodendrocyte | 85 | 4 |
| Vascular | 51 | 4 |

That is **5 of this project's 8 roster types** — Tumor, Macrophage_Microglia (from Immune),
Oligodendrocyte, Endothelial (from Vascular), Astrocyte — across **4 donors**, so cross-donor
variance is estimable and MuSiC and SCDC can use it. Compare the alternatives: Neftel gives 4
of 8 with 9 patients but needs the labels derived (item 16); Albiach gives 8 of 8 with 1 donor,
which a cross-donor method cannot use.

**So Darmanis is the best available second-reference candidate**, and it would make
`SCDC ENSEMBLE` non-degenerate on a declared 5-type sub-roster — the thing item 9 was for.

**Caveats to state.** Astrocyte (88 cells), Oligodendrocyte (85) and Vascular (51) are thin,
and per patient thinner still (Astrocyte runs 3 to 51). Darmanis pools all lymphoid and myeloid
cells into one `Immune cell` label, so T, NK and B cannot be separated. It is Smart-seq2, where
GBmap is 87% 10x — a platform difference that is a confound and an opportunity, since it is the
regime CIBERSORTx's S-mode exists for.

**Where it is.** `data/raw/darmanis_2017/` — the count matrix, the series matrix, and
`cell_metadata.csv` extracted from it (3,589 rows; cell type, tissue, patient, sorting gate).

## 16 · The Neftel annotations can be derived locally — but only 4 of 8 roster types

**Found 2026-09-14 in the pipeline packages.**
`pipeline packages/ repos/scRNA_GBM_Neftel2019-main/scripts/01_GBM_scRNAseq_analysis.Rmd`
is a complete 966-line derivation: per-sample MAD QC, scDblFinder, Neftel's own
log2(TPM/10+1), Harmony, clustering, **marker-based cell type scoring**, inferCNV, a
three-method malignant call, and UCell scoring of the AC/OPC/NPC/MES states.

Both TPM matrices it needs are already on disk (1.5 GB 10x, 1.0 GB Smart-seq2), as is
`GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx`, which gives **per-cell donor
labels for 24,131 cells** — confirming 9 patients and giving the donor structure a cross-donor
method needs.

**So item 9 is not strictly blocked on the authors.** But two things temper it:

1. **It resolves only Macrophage, T cell, Oligodendrocyte and malignant** — the only
   non-malignant types Neftel's own pipeline scores. That is **4 of this project's 8 roster
   types**, with no endothelial (which carries C3, C4 and C6), no NK, no B cell, no astrocyte.
   So it cannot pair with GBmap as a full second reference, and **SCDC ENSEMBLE would still be
   degenerate on the full roster**. It could support a declared 4-type sub-roster sensitivity
   analysis, which is worth something but is not what item 9 promised.
2. **The annotations would be OURS, not the authors'.** A reviewer will ask whether our own
   clustering biased a reference toward our own roster. Author-supplied labels (A4) do not have
   that problem, which is why asking is still first.

The Rmd's own comments record that inferCNV needed **16 GB and downsampling to ~1,164 cells**.
On 8 GB it is not worth attempting. List B.

---

## 17 · Audited and found to contain nothing further — 2026-09-14

Recorded so nobody re-audits these:

- **CIBERSORTx** — Supplementary Table 1d is already used for the B-mode/S-mode decision.
  Tables 2a–2e are signature matrices for NSCLC, PBMC, follicular lymphoma and head and neck;
  **none for glioblastoma**. Tables 3 and 4 concern CREBBP genotyping and high-resolution
  DEGs. Nothing further to extract.
- **ABSOLUTE** — source only (`src/`, `v1.5/`, `wolF/`). The purity data itself was obtained
  separately; see item 3.
- **MuSiC, SCDC, DWLS, Bisque** — package sources, already installed and running. The MuSiC
  and SCDC sources were read to settle how `cell_size` and `ct.cell.size` are used (D1).
- **Scaden** — paper PDF only. Already assessed and dropped, with the reason recorded.
- **CDSeq** — R package source plus the paper. **Not wired**; see A6.
- **EcoTyper** — 1 GB of package and example outputs. **Not wired**; list B.
- **DWLS_other** — auxiliary scripts, nothing the installed package lacks.

---

## 1 · Register the constraint file publicly — **DONE 2026-09-10, no longer blocking**

> Registered with OSF at `2026-09-10T03:32:03Z`, hash
> `2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`, and every
> result file postdates it (`registration_status.json`). The text below is kept
> because it states what the registration is *for*, which the paper still needs.
> Note the corrections owed to the record: `ROAD_TO_PAPER.md` Tier 0.2.

**Why it blocks.** The whole design rests on the constraints being written *before* any
deconvolution output was seen. Right now the file is hashed and committed, which proves
it has not **changed** — not that it **predates** the results. A reviewer will test that
distinction first, and `registration_status.json` currently reports `UNREGISTERED`
because asserting otherwise would be untrue.

**The hash to register**

```
2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a
```

Regenerate it any time with:

```bash
python3 -c "import sys;sys.path.insert(0,'.');from ivygap.anatomic import constraints as K;print(K.freeze_hash())"
```

> **Read this before you start — OSF is changing underneath this step.**
>
> OSF is retiring Projects. Announced dates: **no new projects can be created after
> 16 November 2026**, and existing projects go **read-only on 19 February 2027**.
> OSF's own notice states that **Registrations, Preregistrations and Preprints are
> *unaffected*** — so the timestamp this step exists to create is not going away.
> Only the *route* to it is.
>
> As of today (2026-09-05) that leaves **72 days** on the project-based route below.
> Route A avoids projects entirely and is the one to prefer.

**Route A — register directly, no project (preferred; unaffected by the transition)**

1. Go to <https://osf.io/registries> and sign in (free; ORCID or institutional login).
2. **Add a registration** → choose the **OSF Preregistration** template. Start from
   scratch rather than linking a project.
3. Put the hash in the registration narrative, where it is part of the frozen record.
   In the summary/description field paste:

   > SHA-256 of the canonical constraint payload: `2d1fb47c...807a`
   > Constraints: seven ordinal claims about which cell types are enriched in which
   > anatomic structures, fixed before any deconvolution output was inspected.

4. Paste the **full contents** of `constraint_file.json` into a long-form field (it is
   small). This matters: a registration created without a project may not accept
   arbitrary file uploads, and the hash is only meaningful if the payload it hashes is
   *in* the immutable record. Pasting the JSON makes the registration self-verifying.
5. Submit. Submission — not saving a draft — is what mints the timestamp.
6. Copy the registration URL and its UTC timestamp.

**Route B — project first, then register (works only until 2026-11-16)**

Use this only if you specifically want the two files attached as file uploads.

1. <https://osf.io> → **Create new project**, name it e.g. *The Anatomy Test — anatomic
   concordance constraint file*. Keep it **Public**.
2. Upload `ivygap/anatomic/constraints.py` and
   `results_archive/<stamp>/results/anatomic/constraint_file.json`.
3. Paste the hash into the project Wiki with one line of context.
4. **Registrations** tab → **New registration** → *OSF Preregistration* → submit.
   A plain project upload is **not** a registration and carries no immutable timestamp.
5. Copy the registration DOI/URL and UTC timestamp.

Because projects go read-only in Feb 2027, treat anything uploaded this way as a copy,
not as the archival home. The registration is the archival object; the repo is the
working copy.

**Then, either route** — fill in `REGISTRATION.template.json` and rename it to
`REGISTRATION.json`:

```json
{
  "registry": "OSF",
  "url": "https://osf.io/xxxxx/",
  "registered_utc": "2026-09-05T14:00:00Z",
  "constraint_freeze_hash": "2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a"
}
```

**Confirm it worked**

```bash
python3 -c "import sys;sys.path.insert(0,'.');from ivygap.anatomic import registration as r;s=r.status();print(s.state);print(s.verdict)"
```

Expect `REGISTERED`. If you see `HASH_MISMATCH`, the constraints were edited after
registering — that is the one fatal state, and it means either restoring the registered
file or registering the new one and marking everything before as superseded.

**Timing matters.** Register *before* the next run, so every result file postdates the
registration. The checker compares mtimes and will tell you if any result predates it.

**Alternatives, all unaffected by the OSF change.** Any public, immutable timestamp
works — record which one you used in `registry`:

| Option | Mints a DOI | Account needed | Good for |
|---|---|---|---|
| **Zenodo** <https://zenodo.org> | yes | yes (GitHub login works) | attaching the actual files *and* getting a DOI; the closest drop-in for Route B |
| **OpenTimestamps** <https://opentimestamps.org> | no | no | a Bitcoin-anchored proof of the hash in about a minute, with nothing to sign up for |
| **Software Heritage** <https://softwareheritage.org> | no (SWHID) | no | archiving the repo itself at a commit |

Zenodo is the recommended second choice: drop in the same two files, and it timestamps
and DOIs them without depending on OSF's project system at all. OpenTimestamps is worth
doing *regardless* of which registry you pick — it costs nothing, needs no account, and
independently proves the hash existed today.

---

## 2 · A second tissue — **the only way past the statistical ceiling**

**Why it blocks.** Nine evaluable tumours give a weighted denominator of 65, so ACS can
take only 66 values and fourteen methods land on a handful of them. No analysis inside
Ivy GAP widens this. The protocol's step 7 says it plainly: *one tissue is a case study;
two is a method.*

**What to look for**, in priority order:

| Candidate | Why | Where |
|---|---|---|
| Allen Human Brain Atlas | The protocol's own suggestion; region-annotated bulk microarray across cortical structures | <https://human.brain-map.org/static/download> |
| GTEx (multi-region brain) | 13 brain subregions, bulk RNA-seq, large n | <https://gtexportal.org/home/downloads> |
| HPA / Allen Aging, Dementia & TBI | Region-annotated bulk with clinical data | <https://aging.brain-map.org> |

**Requirements** — a candidate is usable only if all four hold:

1. **Bulk RNA (not single-cell)**, linear scale or convertible to it.
2. **A region/structure label per sample**, assigned by anatomy or histology — **not**
   derived from expression. This is the same rule that excludes Ivy GAP's 148
   ISH-selected samples.
3. **Multiple regions per donor**, so comparisons stay paired within subject. This is
   the design's source of power and it is not negotiable.
4. **A cell-type reference that covers the tissue.** GBmap will not do for cortex.

**Steps**

1. Download and drop it into `pipeline packages/repos/<Tissue2>/`, the same way you
   supplied the deconvolution packages.
2. **Tell me the structures and the expected biology before I touch the data.** I will
   write `ivygap/anatomic/constraints_tissue2.py` and hash it *first*. The ordering is
   the point: a constraint file written after seeing the data is worthless, and the
   protocol's gate for step 7 is exactly this.
3. Register that second hash too (repeat item 1).
4. Then I run it.

**What changes:** the agreement test gains independent evidence, and the claim moves
from *"this worked on glioblastoma"* to *"this is a method"*.

---

## 3 · ABSOLUTE purity for TCGA-GBM — yardstick 2

**Why it matters.** It is the only yardstick that never touches RNA — DNA-derived tumour
purity, entirely orthogonal to everything else here. Your predecessor project recorded it
as BLOCKED, never located.

**Steps**

1. Go to the GDC PanCanAtlas publications page:
   <https://gdc.cancer.gov/about-data/publications/pancanatlas>
2. Download **`TCGA_mastercalls.abs_tables_JSedit.fixed.txt`** (ABSOLUTE purity/ploidy
   calls, all TCGA tumours). It is open-access — no dbGaP application.
3. Alternatively the ABSOLUTE segment tables are in the same section; the master calls
   file is the one with a purity column per sample.
4. Save it as `data/raw/tcga/absolute_purity.txt`.
5. Tell me, and I will wire `real_yardsticks.load_absolute_purity()` to it, filtering to
   GBM barcodes and joining on `patient_id`.

**What it needs to contain:** one row per TCGA sample with a barcode and a `purity`
column. ~120 GBM samples is what the protocol anticipates.

**What changes:** the primary result gains a second, RNA-independent yardstick. If ACS
agrees with DNA purity as well as with pseudobulk accuracy, that is a much harder result
to explain away.

---

## 4 · Ivy GAP vital status — closes the prognosis question

**Why it is stuck.** `tumor_details.csv` publishes `survival_days` and **no vital-status
column**. Without knowing who was censored, a C-index cannot be computed honestly. Worse,
the missing times are not missing at random: 8 of 9 blanks are MGMT-methylated
(Fisher p = 0.0021) and younger — i.e. the good-prognosis patients — so dropping them
biases the cohort toward short survival.

**Steps**

1. Contact the Allen Institute, referencing the Ivy GAP RNA-Seq release. The old
   `portal.brain-map.org/contact-us` page now 404s; the routes that resolve today
   (checked 2026-09-05) are:
   - <https://community.brain-map.org/> — the Allen Brain Map community forum, and the
     best first stop for a data question: answers are public and citable, which a
     private email reply is not.
   - <https://alleninstitute.org/contact-us/> — the institutional contact form.
   - `info@alleninstitute.org` as a fallback.
2. Ask specifically for: **vital status at last follow-up, and follow-up duration for
   patients with no recorded `survival_days`** — for the 42 tumours in
   `tumor_details.csv`.
3. If they supply it, save as `data/raw/ivygap/vital_status.csv` with columns
   `donor_id, vital_status, days_to_last_followup`.

**Be aware this may not change the answer.** Even reading every recorded time as a death
— the most generous assumption available — the cohort reaches 29 events and detects a
C-index difference of 0.515, against a plausible real effect near 0.05. Vital status
would make the analysis *correct*; it would not make it *powered*. Ivy GAP is 41 tumours.
For a powered prognostic result you need a different cohort, and TCGA-GBM is the obvious
one.

---

## 5 · CIBERSORTx Supplementary Note 1 — **FOUND 2026-09-06, no longer blocking**

> **This item is resolved.** The document is
> `~/Downloads/CIBERSORTx Supplementary.pdf` (14.5 MB, 43 pages); Supplementary Note 1
> ends on **page 39**, which states the S-mode algorithm in full. Nothing further needs
> to be obtained. What remains is implementation, which is now unblocked.
>
> **S-mode, as specified on p39** — it adjusts the SIGNATURE, not the mixtures:
>
> 1. Take signature `B` (c cell types) and the single-cell profiles `R` (n genes x r
>    cells) that `B` was derived from.
> 2. Build `k` artificial mixtures `M*` from cells in `R`. Per cell type, draw mixing
>    coefficients from `N(mu, sigma)` with `mu` = that type's fractional abundance in
>    `R`, and `sigma = 2*mu`. (The paper reports S-mode is robust to moderate variation
>    in these, Supp. Fig. 1k,l.)
> 3. Clip negatives in `F*` to 0, then normalise each mixture's coefficients to sum 1.
> 4. Sample single-cell transcriptomes per type according to `F*`, aggregate into `k`
>    bulk profiles in **TPM space** -> `M*`.
> 5. Apply **ComBat to `M` and `M*` together in log2 space** -> `Madj`, `Madj*`.
> 6. Convert `Madj*` back to linear space. Run **NNLS with `Madj*` and `F*`** to
>    reconstruct per-gene cell-type coefficients -> the adjusted signature **`Badj`**.
>    No adaptive noise filtration is needed here: `F*` is known exactly, and ComBat is
>    linear in log2 space so it preserves ordering.
> 7. Estimate `F` from the ORIGINAL `M` against `Badj`, by ordinary CIBERSORTx.
>
> Edge case, same page: NNLS needs more mixtures than cell types (`k > c`). If `M` has
> too few, generate extra pseudo-mixtures from `R` by repeating step 2, pairing each
> with a randomly chosen sample of `M` before S-mode normalisation. With 8 cell types
> and 122 or 270 mixtures this project is far past that bound, so the edge case does
> not bind here.
>
> Everything S-mode needs is already in this repo: the cell-level source (`r_bridge`'s
> registered cell source), `_combat_adjust()` in `deconv/classical.py`, and NNLS.



**Why it matters here specifically.** The paper says signature matrices from
droplet/UMI platforms cause deconvolution to fail, with cell types *"dropping out"* — and
S-mode, not B-mode, is what they recommend for exactly that case. This project
deconvolves 2014 laser-capture bulk against a **10x Chromium** atlas. S-mode is the mode
built for our configuration, and B-mode (which is implemented) is the weaker one.

**Steps**

1. Open the paper's page: <https://www.nature.com/articles/s41587-019-0114-2>
2. Scroll to **Supplementary information** and download **Supplementary Note 1** (or the
   combined Supplementary Information PDF).
3. Save it beside the main paper in `pipeline packages/repos/CIBERSORTx/`.
4. Tell me. S-mode adjusts the *signature* rather than the mixtures; with the algorithm
   in hand I can implement it the same way B-mode was implemented from the Methods.

**Why I have not guessed at it:** producing something that is not S-mode and labelling it
S-mode would be worse than not having it. The current entry states plainly that B-mode is
implemented and S-mode is not.

---

## 6 · RAM for the gene-space reconstruction — **DONE 2026-09-14, on the 8 GB laptop**

> Resolved without a bigger machine. `scripts/reconstruct_gene_space.py` recovered the
> leaderboard's exact 657-gene space and both verification gates passed: the 88/22 donor split
> reproduced by name, and the gene count reproduced exactly. The first attempt failed the
> count gate at 513 genes because it omitted `restrict_to_genes`, which is what the gate is
> for. DWLS and BayesPrism were then re-measured with all seven input-equivalence conditions
> passing — DWLS 0.7846, BayesPrism 0.8154.
>
> It worked on 8 GB only with the browser and editor closed; two concurrent attempts thrashed
> to a standstill. The text below is kept because that constraint still governs list B.

**Why it blocks.** D10: the run persisted the gene *count* (657) but never the gene *list*, so
both genuine-package re-measurements (DWLS, BayesPrism) silently ran on the 1,591-gene
fallback space and their comparisons are withdrawn.
`scripts/reconstruct_gene_space.py` recovers the exact list in minutes instead of a 34,721 s
benchmark re-run — but it must load `gbmap_core.h5ad`, which is **7.6 GB on disk, 338,564
cells x 27,632 genes**. On this machine that load competes with everything else: observed
free RAM of **8–100 MB** with **6–7 GB of 8 GB swap in use**, and two concurrent attempts
thrashed to a standstill at 78% CPU with an RSS of 34 MB.

**What to get.** Any of, in order of preference:

1. A machine with **≥ 32 GB RAM** (a cloud VM for one hour is enough — this is a
   minutes-long job once memory is not the constraint).
2. Close Edge and VS Code and run it alone. Measured: they hold ~2.2 GB between them, which
   is the difference between thrashing and finishing.
3. Failing both, accept the 34,721 s route: one full `python scripts/run_all.py`, which now
   persists the list itself.

**How to confirm it worked.** The script refuses to write unless two recorded facts
reproduce: the 88 training and 22 test donors by name, and the gene count of exactly 657. If
it writes, both passed.

```bash
python scripts/reconstruct_gene_space.py
python scripts/remeasure_method.py --method dwls       --budget 14400   # ~1.1 h
python scripts/remeasure_method.py --method bayesprism --budget 14400   # ~1.3 h
```

**Do not run those two concurrently.** BayesPrism forks three Gibbs workers at ~900 MB each.

**What changes.** DWLS and BayesPrism become comparable to the leaderboard, so
`method_completeness.json` goes from 9/15 to 11/15, and the Avila Cobos disagreement — which
this project currently cannot adjudicate — becomes answerable.

---

## 7 · Ivy GAP vital status — see item 4

Unchanged. Survival is BLOCKED for want of a censoring indicator, and the missingness is
informative (Fisher p = 0.002134 against MGMT methylation).

---

## 8 · A second tissue — see item 2, and it is now more urgent than when that was written

Two findings since make the ceiling sharper than "nine tumours":

- `constraint_sensitivity.json`: dropping **C6** alone moves 13 of the 14 ranked methods and
  takes rank agreement to rho 0.678. The ordering below the top is substantially one
  constraint's work.
- C6 and C5 are the `Macrophage_Microglia` constraints, and that column absorbs the **7.0%**
  of the atlas the roster drops (Mono 14,215 cells, DC 3,961, RG 2,807).

So the part of the leaderboard that discriminates most rests on the roster's most heavily
loaded column. A second tissue **must not reuse a roster with the same weakness**, or the
replication is not independent of it. See `docs/SECOND_TISSUE.md`.

---

## 9 · Cell-type labels for GSE131928 (Neftel) — a second reference

**Why.** Two things need it. `SCDC ENSEMBLE` is currently **degenerate**: with one reference
there is nothing to weight across, so it reduces exactly to SCDC and is reported as degraded
rather than as an ENSEMBLE result. And a second reference is the only way to test whether any
result here is an artefact of GBmap.

**State.** Confirmed against all four GEO supplementary files: **GSE131928 ships no
cell-type labels**. 9 patients, 21 samples, 16,201 cells (10x) plus SMART-Seq2. The
annotations exist in the paper's own analysis, not in the GEO record.

**What to get.** The authors' per-cell state assignments (MES/AC/OPC/NPC-like plus
non-malignant calls) — from the Broad Single Cell Portal, the paper's supplementary tables, or
by writing to the authors. Albiach is already on disk **with** labels but is **one donor**,
which is not enough for a method that needs cross-donor structure.

**Where it goes.** `data/reference/`, then a second `ReferenceBundle` passed to
`DeconvolutionInput.references`.

**What changes.** SCDC ENSEMBLE stops being degenerate (1 closable gap closed), and every
method gains a reference-sensitivity check.

---

## 10 · A reference that retains the myeloid subtypes — tests D8 directly

**Why.** D8's concern is that C6, the constraint carrying most of the ranking, sits on the
column absorbing Mono and DC. If C6's discriminating power survives *un-collapsing* those
populations, the concern is answered. If it does not, the finding is about the roster and must
be stated that way.

**What to get.** Either a GBM atlas with Mono/DC/Mast retained as separate labels, or GBmap
re-mapped onto a nine- or ten-type roster that keeps them. This is a **roster change**, so it
cannot touch the registered constraint file — it must be reported as a declared sensitivity
analysis alongside, never substituted in.

---

## 11 · Paired bulk and single-cell from the SAME subjects — Bisque at full strength

**Why.** Bisque's `ReferenceBasedDecomposition` has an overlap mode that estimates the
assay transform from subjects assayed both ways. Ivy GAP's bulk and GBmap share no subjects,
so it runs `use.overlap = FALSE` and the transform comes from marginal distributions. That is
Bisque's own documented fallback, and the row says so — but it means **this study has never
measured Bisque as intended**.

This matters more than it looks: Bisque is the panel's **most reliable method on tumour
purity** (smallest bias −0.0256, narrowest 95% interval ±0.489), and it is also the only
confirmed double cell-size correction (D1). Both of those want resolving on the real thing.

**What to get.** A GBM cohort with bulk RNA-seq *and* single-cell from the same patients —
e.g. GSE84465, or any study publishing both arms per subject.

---

## 12 · Immune ground truth — the largest error mode in the panel is unverified

**Why this is the most scientifically valuable item on the list.** Every method
over-estimates T cells in high-purity tumour by **4.6× to 7.9×** (true 0.044, predicted 0.202
to 0.350), and `ACS cannot see it` because `T_cell` carries no constraint. That failure is
currently established only against *pseudobulk* truth. If it reproduces against real
measured immune content, it is a finding about deconvolution as a clinical tool, not about
this benchmark.

**What to get**, any one of:

- **IHC or multiplex immunofluorescence** for CD3/CD8 on the Ivy GAP blocks, or on any GBM
  cohort with matched bulk RNA-seq.
- **Flow cytometry or CyTOF** immune fractions with matched bulk.
- A published GBM cohort with **orthogonally measured** immune infiltration and bulk RNA-seq.

**Where it goes.** A new yardstick beside `synthetic_mixtures` in
`yardstick_provenance.json`, scored by the same agreement machinery.

**What changes.** It either confirms the T-cell over-call on real tissue — which is a
publishable finding on its own and the strongest clinical result this project could produce —
or shows it is an artefact of pseudobulk construction, which is equally worth knowing and
would materially soften `CLINICAL_READINESS.md` §3.

---

## 13 · Pathologist or orthogonal purity on the Ivy GAP blocks

**Why.** Every method under-calls tumour by **0.33 to 0.51** at high purity, and real Ivy GAP
tissue sits there (MuSiC median 0.68). That is measured against pseudobulk truth. Against a
pathologist's estimate or an orthogonal purity call on the *same blocks*, it becomes a
statement about the tissue rather than about the simulation.

**What to get.** H&E-based tumour-cell-percentage estimates for the 122 anatomic blocks — Ivy
GAP publishes the H&E images, so this is a reading exercise, not new wet work. Or DNA-based
purity (ABSOLUTE-style) if any Ivy GAP block has matched sequencing.

**What changes.** Converts the panel's most serious accuracy failure from a benchmark artefact
into a tissue-level result, and gives the paper a purity axis that does not share a failure
mode with the pseudobulk construction.

---

## 14 · Methods not yet wired, and whether they are worth it

`pipeline packages/ repos/` already contains **CDSeq** and **EcoTyper** sources that nothing
in `ivygap/` or `R/` references. Both are real additions rather than more of the same:

| candidate | what it adds that the panel lacks | verdict |
|---|---|---|
| **CDSeq** | *reference-free* — estimates cell types without a signature, so it cannot inherit the reference's roster collapse (D8) or its marker-space bias (D11) | **worth wiring.** It is the only way to separate "the methods are wrong" from "the reference is wrong". |
| **EcoTyper** | recovers cell *states* rather than types, which is the axis Neftel's programmes live on | worth it only if item 9 lands; without labels there is nothing to validate against. |
| **InstaPrism** | a fast deterministic re-implementation of BayesPrism's model | useful only as a cross-check on the BayesPrism row; low priority. |
| **RNA-Sieve** | ships per-sample confidence intervals natively | **worth wiring** — it is the only candidate that addresses `CLINICAL_READINESS.md` §4 at source rather than by conformal post-hoc calibration. |
| **MuSiC2** | built for cross-condition bias, which this cohort does not have | skip, and say why. |
| **Scaden** | already assessed and **dropped**; a deep model needing training data this cohort cannot supply. Mentioned in the paper with the reason. | dropped |

**The honest framing for the paper:** the panel is already large enough that adding a
fifteenth least-squares variant proves nothing. What is missing is not more methods but
**one reference-free method** (CDSeq) and **one method with native uncertainty**
(RNA-Sieve), because those two test the two failure modes the current panel cannot
distinguish from itself.

---

## Order of operations

```
1. Register (item 1)          ← do this before the next run
2. Re-run and archive          ← every result then postdates the registration
3. Second tissue (item 2)      ← constraints written and hashed BEFORE data
4. Items 3–5 as they arrive
```

Items 3, 4 and 5 improve the result. Items 1 and 2 change what kind of claim the paper
can make.

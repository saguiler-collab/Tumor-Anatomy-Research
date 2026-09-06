# External actions — the five things the pipeline cannot do for itself

Each item says what to get, exactly where it goes, what changes in the results when it
lands, and how to confirm it worked. Ordered by what unblocks the most.

Nothing here is optional polish except items 3–5. Items 1 and 2 are what separate
"a strong internal result" from "a paper".

---

## 1 · Register the constraint file publicly — **blocking the central claim**

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

## Order of operations

```
1. Register (item 1)          ← do this before the next run
2. Re-run and archive          ← every result then postdates the registration
3. Second tissue (item 2)      ← constraints written and hashed BEFORE data
4. Items 3–5 as they arrive
```

Items 3, 4 and 5 improve the result. Items 1 and 2 change what kind of claim the paper
can make.

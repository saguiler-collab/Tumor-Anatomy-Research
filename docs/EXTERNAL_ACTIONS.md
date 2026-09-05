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

**Steps**

1. Go to <https://osf.io> and sign in (free; ORCID or institutional login works).
2. **Create new project** → name it e.g. *The Anatomy Test — anatomic concordance
   constraint file*. Keep it **Public**.
3. Upload two files from this repo:
   - `ivygap/anatomic/constraints.py` — the constraints themselves
   - `results_archive/<stamp>/results/anatomic/constraint_file.json` — the canonical
     JSON payload the hash is computed from
4. In the project **Wiki** or description, paste the hash above and one line:
   *"SHA-256 of the canonical constraint payload; see constraint_file.json."*
5. **Registrations** tab → **New registration** → *OSF Preregistration* template →
   submit. This is the step that creates the immutable timestamp. A plain project upload
   is not a registration.
6. Copy the registration's DOI/URL and its UTC timestamp.
7. Fill in `REGISTRATION.template.json` and rename it to `REGISTRATION.json`:

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

**Alternative if OSF is slow:** Zenodo (also DOI-minting) or an OpenTimestamps proof of
the hash. Any public, immutable timestamp works — record which one in `registry`.

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

1. Email the Allen Institute: <https://portal.brain-map.org/contact-us> (or
   `info@alleninstitute.org`), referencing the Ivy GAP RNA-Seq release.
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

## 5 · CIBERSORTx Supplementary Note 1 — enables S-mode

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

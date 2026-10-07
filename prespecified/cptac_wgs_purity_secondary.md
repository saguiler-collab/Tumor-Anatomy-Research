# Secondary analysis S2: whole-genome DNA purity as the per-sample truth for CPTAC tumour content

**Written 2026-10-06 22:34:47 EDT** (the file's birth time, copied from stat).

**Secondary, registered after the primary CPTAC analysis had been read.** Nothing here can change a
primary reading of `prespecified/cptac_per_sample_truth.md`.

## What had been seen at writing

- **The primary result.**
  - Its DNA control failed: snRNA tumour share against GIMiCC purity, Spearman 0.014, n 15 (S1 truth:
    -0.029). Every primary reading is therefore INCONCLUSIVE.
  - Its descriptive values (C1-C4) had been printed.
- **A post-hoc look at that control's instrument.**
  - GIMiCC's complete-case rule kept 769 of the 4,022 library CpGs.
  - Five samples cover only 54-67% of the library: C3N-01814, C3N-00662, C3N-02784, C3N-01815 and
    C3N-02188. The other 13 cover at least 95.7%.
  - Bulk methods' tumour estimates against GIMiCC purity: median Spearman +0.05 (frozen signature) and
    -0.05 (donor-level reference).
  - The same estimates against the snRNA tumour share: +0.24 and +0.35.
- **GDC metadata only.**
  - Every case has exactly one open AscatNGS whole-genome copy-number record, which carries the fields
    `tumor_purity` and `tumor_ploidy`.
  - Every case's tumour slides carry `percent_tumor_nuclei`.
  - Which tumour sample each file comes from:
    - For all 18 cases, the WGS, bulk RNA and methylation aliquots come from the same tumour
      sample(s), and the nuclei come from one of them.
    - Five cases pooled two or three tumour pieces for bulk, methylation and WGS: C3L-02705, C3N-01798,
      C3N-01814, C3N-02190 and C3N-02784.
- **No purity, ploidy or pathology value had been seen.**

## Why

The primary control failed, and it cannot say which side failed. Either the nuclei do not reflect
tumour content, or the DNA instrument failed: GIMiCC ran on a fifth of its library, and five samples
lack a third or more of it.

Whole-genome copy-number purity settles this:
- It is the same kind of truth as ABSOLUTE in TCGA.
- It is measured on the same aliquots as the bulk RNA.
- It answers the question directly, and it gives the TCGA-GBM purity result an independent cohort.

## The truth

- **AscatNGS tumour purity**, from the GDC harmonised whole-genome pipeline.
  - ASCAT: Van Loo et al., PNAS 2010. ascatNGS: Raine et al., Curr Protoc Bioinformatics 2016.
  - It is read from the `tumor_purity` field of each case's open "Allele-specific Copy Number Segment"
    file record. The file_id, data release and workflow version are recorded.
- **Used as GDC reports it.** A case without a value is excluded and counted, never imputed.
- **Values >= 0.99** may be fits to a genome without aberrations. They are counted and kept, and every
  S2 reading is also given without them.

## Questions and readings, fixed now

Every bar is one already registered: 0.40 (control 1, C1, E2) and 0.60 (the agreement bar). No new bar
is introduced. Permutation p values follow the primary rule's control 3: 10,000 shuffles of the
pairing, seed `config.RANDOM_SEED`.

### S2.1 Is the nuclei truth valid per sample?
- Spearman between the registered snRNA tumour share and WGS purity, over the 15 scored cases. The S1
  truth is reported beside it.
- **Pass if >= 0.40.**
  - On a pass, C1-C3 of the primary rule (already computed) may be reported, but only as **secondary
    readings** and always with this history: the registered DNA control failed, and a whole-genome
    control registered after that failure passed.
  - On a fail, the nuclei do not reflect tumour content here, and C1-C3 stay INCONCLUSIVE.

### S2.2 Was the methylation instrument the problem? (diagnostic)
- **(a)** GIMiCC purity as already computed (769 CpGs), against WGS purity, over the 18 cases.
- **(b)** GIMiCC re-run on the 13 samples that cover >= 90% of the library, using the CpGs complete
  across those 13. Scored against WGS purity over those 13.
  - Any threshold from 68% to 95% selects the same 13.
  - Coverage is a property of the input, and it was read before any purity value.
- **Reading:** if (b) reaches 0.40 and (a) does not, the failed primary control is attributed to
  methylation coverage. This never rescues the primary.

### S2.3 Do bulk methods recover tumour content per sample, against DNA?
- Per method and per reference build: Spearman between the bulk tumour estimate and WGS purity, over
  the 18 cases.
- Per reference, the median across methods. **Tracks if >= 0.40.**
- Set beside it: the TCGA-GBM median of the same statistic against ABSOLUTE
  (`absolute_purity_yardstick[_h5ad].json`). It is recomputed over the same methods (matched
  denominator).
- Descriptive:
  - recovery, the OLS slope of estimate on purity (as `check_recovery_definition` in
    `scripts/independent_verification.py`);
  - mean bias.

### S2.4 Does the methods' DNA-accuracy ranking transfer between cohorts?
- Per reference: Spearman between each method's S2.3 rho and its TCGA-GBM rho against ABSOLUTE, over
  the methods present in both.
- **Transfers if >= 0.60.**
- This decides whether a lab can choose a method from a published benchmark in a similar tissue
  (`docs/ACCURACY_FACTORS.md`, factor 3).

### S2.5 The registered question, with a whole-genome per-sample truth in a new cohort
- Spearman between ACS (the registered leaderboard) and each registered method's S2.3 rho. The methods
  are the registered 12 that are present, selected as C2 selects them (`yardstick_agreement.json`).
- **Bar 0.60.**
- The frozen-signature value corresponds to arm 2 of the registered test (rho 0.081, TCGA-GBM). The
  donor-level value is reported beside it.

### S2.6 Pathology (descriptive)
- Per case, the mean `percent_tumor_nuclei` over the tumour slides of the analysed sample(s).
- Spearman against WGS purity, and against the snRNA tumour share.
- CPTAC selected its samples on tumour content, so the range may be narrow. These values are reported,
  not read.

## Limits, stated now

- **Registered after the primary was read, so S2 is secondary throughout.**
- Only 18 cases, and 15 for anything that involves nuclei.
- ASCAT purity carries its own error: it is a model fit, not a count.
- Five cases pooled tumour pieces, and their nuclei come from one piece. Every S2.1 value is also given,
  descriptively, on the 10 single-piece cases scored.
- The cohort is glioblastoma only.

## Implementation

`scripts/cptac_wgs_purity.py`
- `--stage fetch`: writes `data/external/cptac_gbm/wgs_ascat_purity.tsv` and
  `pathology_tumor_nuclei.tsv`.
- `--stage gimicc`: runs S2.2(b).
- `--stage analyse`: writes `results/cptac_wgs_purity.json`.

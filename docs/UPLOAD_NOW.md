# Upload now, on the slow laptop. Tomorrow's machine is for one thing only.

**Written 2026-09-15 for a two-week deadline.** Everything here is browser work, no compute. The
project already has its result (`docs/ENDPOINT.md`); none of this changes it. These are the
steps that make it *citable and defensible*, which is what a judge or reviewer checks.

---

## TONIGHT — 4 items, about 90 minutes total, all browser

### 1. Post the corrections to OSF · **20 min** · *highest value of anything left*
Go to <https://osf.io/vuh64> (your project wiki) and add a link to
`docs/CORRECTIONS_REGISTRATION.md`.

**Why this one matters more than the rest.** Your registration is at <https://osf.io/dm2t8> and
is **retrospective** — that is a known weakness a reviewer will find. Ten corrections published
beside it, including **two retractions of your own claims** (C9, and D14's magnitude), turn that
weakness into the strongest thing about the project. A registration nobody re-examined is worth
less than one whose author publicly corrected it.

Post the section headed **"THE REGISTRATION UPDATE, IN SHORT"** — it is written to be pasted.

### 2. Push this repository · **5 min**
Six commits are sitting locally and unpushed. Everything below links to the repo, so nothing
else works until this is done. **Tell me and I'll push** — I don't push without being asked.

### 3. Zenodo DOI for the repository · **20 min**
<https://zenodo.org> → log in with GitHub → enable the repo → cut a release on GitHub. You get a
permanent DOI and a frozen snapshot.

**Why:** ISEF and STS both reward citable artefacts, and a DOI makes the pre-specifications
(`prespecified/`) verifiable by date to anyone. That is the thing that makes "I predicted this
before I ran it" a fact rather than a claim.

### 4. Archive the two external anchor PDFs · **10 min**
Save to the project folder, do not commit the PDFs:
- Ajaib et al., *Neuro-Oncology* 2023;25(7):1236–1248 (you have it) — **and download its
  Supplementary Tables S8 and S9.** Those carry the per-cell-type IMC correlations. Right now
  `docs/EXTERNAL_VALIDATION.md` cites the averages quoted in the body text; the tables would let
  the comparison be made per cell type, which is strictly better evidence.
- Avila Cobos et al., *Nat Commun* 2020;11:5650.

---

## TOMORROW, better machine — ONE job left, then stop

### B1 · One full `run_all.py` rebuilt from `raw/X` · **~10 h, unattended**
```
python scripts/build_gbmap_linear_reference.py     # ~10 min, already written
python scripts/run_all.py                          # overnight
```
Closes **D16** (the reference was log-transformed and treated as linear), lands the genuine DWLS
and BayesPrism rows, and produces the Problem 1 / Problem 2 numbers.

**It is not required for any claim the paper makes.** D16's direction is already measured
(+0.9161 — expression space barely moves the ordering), so this converts a stated limitation
into a measured one. Start it before bed and let it run.

### ~~A6 · CDSeq~~ · **DONE 2026-09-15, on the slow laptop — nothing to do**
Reference-free ACS **0.6829** with no atlas at all (p = 0.0009) and **0.7561** with the atlas
used only to name its output, against `control_random` 0.5122 which fails its own null. It ties
DWLS exactly and sits at the bottom of the real methods. **It cuts both ways**: the anatomy is
recoverable without an atlas, *and* the atlas earns its place, since 13 of 14 reference-using
methods beat it. `docs/ENDPOINT.md` §2e.

**Then stop.** Two weeks is not enough for a second tissue, and the ceiling of nine evaluable
tumours is a limitation to state, not to fix.

---

## DO NOT DO THESE

- **A second tissue (B2).** Right answer scientifically, impossible in two weeks. State the
  nine-tumour ceiling as a limitation; that is what it is.
- **Downloading the Ivy GAP BAMs.** 414 MB each, ~112 GB. RSEM already counted. There is no
  question you have that needs read-level data.
- **Emailing the Allen Institute** (old A22). The counts were public all along; you already have
  all 270.
- **Retuning ACS, the constraint file, or any weight** so the result looks better. It is
  forbidden by the project's own invariants and it is *the reason the negative result is
  credible*. The value of this paper is that you pre-registered, it failed, and you published
  the failure.
- **Chasing the IMC raw data.** It is the authors' own cohort and not public through the repo.
  The published verdict is a legitimate anchor and you have cited it correctly.

---

**Running here overnight, and nothing depends on it.** The per-arm gene-space comparison —
D16's one explicitly-unrun question, whether ACS prefers the log reference or the linear one when
each nominates its own markers. Slow (nu-SVR at 1,593 genes) and possibly hours. If it finishes
I will fold the number in; if it does not, the paper is unaffected, because the honest sentence
is already written: *if ACS were to prefer the model-violating reference, that is a finding about
ACS, not a reason to choose a reference.*

## Still genuinely open, and honest to leave open

**A24 — a second external ground truth with per-sample counts.** Your strongest negative result
(C10) rests on one study's published summary and nine tumours. A second anchor would make it
decisive. Realistically: not in two weeks, and the paper should say so in its limitations rather
than pretend otherwise.

**A3 — tumour percentage from the Ivy GAP H&E.** You said you'd get to this. It is worth doing
*if you have a spare afternoon*, because it converts the 0.33–0.51 tumour under-call from a
benchmark artefact into a tissue-level result — one of the paper's strongest clinical claims. It
is not required.

---

## What to say if someone asks "so did it work?"

> The constraints are real — two independent single-cell datasets confirm them without any
> deconvolution. Methods satisfy them far above chance and random controls don't, under every
> reference tested. So anatomy detects whether a method is working.
>
> It does not rank them. It fails against published benchmarks (p = 0.50), fails a
> pre-registered test against imaging mass cytometry (p = 0.0090, wrong direction), only about
> half the ordering survives a change of reference atlas, and it moves by 0.31 on marker-set
> size alone — three-quarters of the range it's ranking fifteen methods across. It can't even
> resolve its own top: four methods tie at 0.9846.
>
> I registered the hypothesis before testing it, it failed, and I published the failure along
> with two retractions of my own earlier claims.

**That last sentence is the paper.** Do not soften it.

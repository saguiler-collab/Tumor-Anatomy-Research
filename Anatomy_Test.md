CellAtlas GBM
Research direction
Proposed original research · computational biology
The Anatomy Test
Every method for estimating cell composition from bulk RNA needs a way to be checked. This one proposes using the shape of the tumor itself — and then asks whether that check actually works.


IVY Gap MetaData: https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fwww.cancerimagingarchive.net%2Fwp-content%2Fuploads%2FdoiJNLP-IvyGAP-nbia-digest-1.xlsx&wdOrigin=BROWSELINK

270
region-labeled RNA-seq samples, Ivy GAP
41
glioblastomas, each sampled at multiple sites
5
pathologist-defined anatomic structures
0
survival outcomes used to choose anything
01 — The deliverable
The name, and the question
The question
Can a tumor's own anatomy stand in for ground truth when choosing a cell-type deconvolution method — and does a method that gets the anatomy right also get the biology right?

Full title
Anatomic Concordance: Neuropathology as a Ground-Truth-Free Benchmark for Cell-Type Deconvolution in Glioblastoma
Short name
The Anatomy Test. The metric it produces is the Anatomic Concordance Score (ACS).
Hypothesis
Deconvolution methods that better reproduce the cell-composition gradients a pathologist would predict across a tumor's anatomic structures are also more accurate against independent ground truth — making anatomic concordance a usable substitute for ground truth in tissues that have none.
Null it can prove
That "the answer matches known biology" is not evidence that the answer is right. Which is a finding, because that check is exactly what most papers currently rely on.
Category
ISEF / Regeneron STS — Computational Biology & Bioinformatics (CBIO).
02 — Why it is open
The logic everyone runs, pointed backwards
Deconvolution estimates cell fractions from bulk RNA. To trust an estimate you need something to check it against, and the field has four options: synthetic mixtures you built yourself (circular — the mixture is made from the same reference the method solves with), flow cytometry on the same specimen (rarely collected), paired single-cell data (expensive, and rarely from the same tissue block), or tumor purity from DNA (one number, one cell type). For most tissues none of these exist. When they don't, the current best practice is a consensus — run several methods, average them, hope they aren't wrong together.

HOW IT IS DONE NOW
Deconvolution method
solves
Cell fractions per region
read as
Anatomy of the
microenvironment
anatomy is
the conclusion
THE ANATOMY TEST
Deconvolution method
Ordering constraints,
written in advance
Anatomy of the
microenvironment
derived from
scores
anatomy is
the test
One edge changes direction. That edge is the entire project.
Nothing is testable while the anatomy is an output; everything is, once it is an input.
The same three components in both rows. In current practice the anatomic pattern falls out of the deconvolution and is reported as a finding — so the method can never fail. Writing the anatomic ordering down first turns the same pattern into a prediction the method can miss.
Already published
CIBERSORTx has been run on Ivy GAP's anatomic structures — to report that immunosuppressive macrophages accumulate in the perivascular zone.
Benchmarks (CATD, omnideconv, Spotless) rank methods where ground truth exists, and fall back to method consensus where it doesn't.
Ivy GAP's regional biology is well characterized on its own terms, by in-situ hybridization and histology.
Not published
Anatomy used as a scoring function that ranks methods, rather than as a result to describe.
A ground-truth-free selection criterion built on hard biological constraints instead of on agreement between methods.
Any test of whether "looks biologically correct" predicts being numerically correct. This is assumed constantly and has never been measured.
That last gap is the one worth taking. A consensus can be confidently wrong in unison; a constraint derived from histology cannot be talked out of. And nobody has checked whether the sanity check everyone uses is a sanity check at all.

03 — The instrument
Five structures, and what each one must contain
The Ivy Glioblastoma Atlas Project laser-microdissected 270 RNA-seq samples from 41 tumors, each labeled with the anatomic structure it came from. Those labels are an independent physical fact about each sample — decided by a neuropathologist looking at tissue, with no reference to any RNA. That is what makes them usable as a test rather than as more data.

LE
IT
CT
MVP
PAN
leading edge
infiltrating
cellular tumor
microvascular
peri-necrotic
Tumor
Oligodendrocyte
Endothelial
Macrophage
HIGHER ON THE PAGE = LARGER PREDICTED FRACTION · LARGE DOT = PREDICTED MAXIMUM
No values are predicted. Only the ordering — which is all histology can honestly assert.
Four predicted profiles, written from neuropathology before any deconvolution is run. A method passes a constraint by reproducing the direction of a step, not its size. Endothelial fraction peaking in microvascular proliferation is close to definitional; oligodendrocyte fraction peaking at the leading edge follows from the leading edge being infiltrated brain parenchyma rather than dense tumor.
Each of those steps becomes a written, numbered constraint. They are ordinal, they are paired within a tumor (so a patient can never be compared against a different patient), and every one of them is traceable to histology or in-situ hybridization — never to a deconvolution result, which is what would make the whole thing circular.

ID	Prediction	Cell type	Why it must hold
C1	CT > LE	Tumor	Near-definitional: leading edge is brain infiltrated at low tumor-cell density; cellular tumor is dense tumor. A method that fails C1 is not measuring tumor.
C2	LE > CT	Oligodendrocyte	The leading edge retains normal white matter and cortex; cellular tumor has displaced it.
C3	MVP > CT	Endothelial	Microvascular proliferation is defined by florid endothelial proliferation. The single strongest constraint in the set.
C4	MVP is the max	Endothelial	Stricter form of C3: MVP must exceed all four other structures, not just CT.
C5	PAN > LE	Macrophage	The hypoxic peri-necrotic niche recruits myeloid cells; the leading edge does not.
C6	MVP > CT	Macrophage	The perivascular niche is a documented myeloid reservoir.
C7	LE < IT < CT	Tumor	A three-step monotone chain — much harder to satisfy by chance than any single pairwise step, and scored as one unit.
—	excluded	T cell	Declared out of scope in advance: this pipeline's own synthetic benchmark puts T cells at the detection floor, so a T-cell constraint would score noise. Excluding it before seeing data is the honest move; excluding it after is not.
04 — The measurement
Turning constraints into a score
Anatomic Concordance Score
For every constraint and every tumor that has both structures sampled, take the paired comparison and record whether it points the right way. ACS is the weighted proportion of satisfied constraint-tumor pairs, bootstrapped over tumors for a confidence interval.

The null comes from permuting structure labels within each tumor, ten thousand times — which preserves each patient's own composition and destroys only the anatomy. That is the comparison that matters.

What gets scored
Eight or more deconvolution methods on identical inputs and an identical signature: the project's frozen NNLS and SVR solvers, plus CIBERSORTx, BayesPrism, MuSiC, DWLS, EPIC and quanTIseq.

And critically, two negative controls — a random-fraction generator and a solver run against a shuffled signature matrix. If those score well, the constraint set is too easy and the result is reported as such rather than quietly retuned.

The experiment that makes this a study rather than a tool
A leaderboard on its own proves nothing — a method could satisfy every anatomic constraint and still be badly calibrated. So the real test is whether the ACS ranking agrees with the ranking from real ground truth, measured on three independent yardsticks where truth is actually known:

500 synthetic mixtures with known composition — already built and frozen in this project.
ABSOLUTE DNA purity on TCGA-GBM (n = 120), an orthogonal measurement that never touches RNA.
Pseudobulk from paired single-cell GBM data, where the composition is counted rather than estimated.
The primary result is one number: the rank correlation between ACS and true accuracy across all ten methods. Pre-specify the bar — say Spearman ρ ≥ 0.6 with a bootstrap CI excluding zero — before running it.

If the result is…	Then the claim is	Verdict
ACS rank tracks true accuracy	A validated way to choose a deconvolution method in any tissue with an anatomical atlas, using no ground truth and no patient outcomes. Transferable well past glioma.	Headline
ACS rank does not track true accuracy	"It matches known biology" is not evidence of correctness — a direct, quantified warning about the informal check the whole field leans on. Publishable, and arguably more useful.	Also a result
Negative controls score highly	The constraint set is too permissive. Report it, publish the constraint file anyway, and state what a harder set would need. Do not retune and re-report.	Report as-is
Every branch produces a defensible paper. That is the property to design for — it is what separates a project from a gamble, and judges recognize it immediately.

05 — Order of operations
Build sequence
Dependency-ordered. Each step is only safe once the one above it holds, and the pre-registration is deliberately placed before the first look at any deconvolution output.

Restore the frozen inputs
Vendor signature_matrix.tsv back into the repo with a recorded hash. The 13 parity tests currently fail on SIGNATURE_MISSING, and nothing below can run through a pipeline that cannot solve its own internal cohort.

Gate: 13/13 parity tests green.
Write and timestamp the constraint file
Every constraint, its direction, its weight, and a citation to the histology or in-situ evidence it came from. Commit it, hash it, and register it publicly — OSF takes minutes and is free. This is the step that makes the whole design credible.

Gate: registration timestamp precedes every result file.
Ingest Ivy GAP
270 RNA-seq samples with structure labels, from the Allen Institute portal or GEO (GSE107559), into a versioned bundle with a per-file SHA-256 manifest — the same ingestion pattern the project already uses for its release data.

Gate: sample counts per structure reconcile against the portal's own documentation.
Build the scorer, prove it on synthetic anatomy
Write the ACS implementation, then test it on simulated region data where you control the answer: composition that satisfies every constraint must score near 1, permuted labels near chance. Unit-test the metric before trusting it on real tissue.

Gate: known-answer fixtures pass; permutation null centered at chance.
Score every method, including the controls
Identical signature, identical samples, identical preprocessing. Report ACS per method with a CI, and per constraint — a constraint that no method satisfies is telling you something about the constraint.

Gate: no method silently drops samples; exclusions reported per step.
Run the agreement test
Rank the same ten methods on the three ground-truth yardsticks, correlate against the ACS ranking, and report against the bar you set in step 2 — whichever way it falls.

Gate: the pre-registered threshold is the one reported.
Prove it transfers
Repeat the construction on a second tissue with region-annotated bulk RNA — the Allen Human Brain Atlas, or a kidney atlas. One tissue is a case study; two is a method.

Gate: constraint file for tissue two written before its data is touched.
06 — Where it can break
Risks, and what neutralizes each
Risk	Why it matters	Neutralized by
Circularity	If a constraint was learned from deconvolution output, the test is rigged and a reviewer will find it.	Every constraint cites histology or in-situ hybridization, and the file is timestamped before any output exists.
A region is not a cell type	Microdissected regions are still mixtures. Absolute fraction targets would be indefensible.	Only ordinal, within-tumor comparisons are ever scored. No constraint asserts a value.
Reference coverage	A seven-population signature cannot speak to astrocytes or neurons, which matter at the leading edge.	Report constraint coverage explicitly, and state which anatomic facts the signature is structurally unable to test.
Sample size	The cleanly annotated anatomic-structures subset is 122 samples from 10 tumors.	The paired within-tumor design gets substantial power from repeated measures; the wider 270-sample set extends it. Report per-tumor results, not just the pool.
Method setup fairness	Any of these methods can be made to look bad by running it wrong.	One shared signature, published configs for every method, and defaults from each tool's own documentation — no per-method tuning.
07 — Fit
Why this competes
Fair judging rewards a question that can be answered wrong. Most computational-biology entries apply an existing tool to a new dataset and report what it found — which is a demonstration, not an experiment, and judges can tell. This design has a stated hypothesis, a pre-registered prediction set, negative controls, and a result that is informative in every direction it can fall.

Novel object
A metric that did not exist before — ACS — plus the constraint file and an open-source scorer that anyone can run on their own tissue.
Real stakes
Deconvolution results feed immunotherapy target selection and prognostic modeling across oncology. Which method you pick changes the answer, and right now that choice is frequently made by whichever one produced a nicer figure.
Independent
All data is open and unauthenticated — Ivy GAP, GEO, the GDC open tier. Compute is laptop-scale. Nothing here requires a lab, a mentor's dataset, or an institutional agreement.
Already half-built
This project already has the frozen solvers, the 500 synthetic mixtures, the ABSOLUTE purity comparison, a manifest-and-hash ingestion pattern, and a site to publish the leaderboard on. The new work is the constraints, the scorer, and the agreement test.
Rules it inherits
No outcome-driven method selection, forward-only application of frozen models, per-cohort before pooled, and the null stays reportable. Those constraints already govern this project — and here they are load-bearing rather than defensive, because the whole point is a selection criterion that never touches survival.
One more thing it fixes, quietly. The current frozen prognosis result is a null with an interval wide enough to hide anything smaller than ±0.063 Uno-C. The obvious objection to any follow-up is "you went looking for a method that worked." Choosing the method by anatomy — decided before survival is ever loaded — is the cleanest possible answer to that objection.

Cohort figures from the Allen Institute's Ivy GAP documentation: 270 laser-microdissected RNA-seq samples across 41 tumors, of which 122 samples from 10 tumors form the anatomic-structures study. Internal counts (153 / 92 / 70 events, the Uno-C panel, the E−M interval) read from the locked run deployable_rc_v1. Verify every external figure against the current portal release before citing it.

Sources — Ivy GAP RNA-Seq · Ivy GAP data documentation · GSE107559 · CATD benchmark · omnideconv · Spotless · Perivascular immune phenotype in GBM
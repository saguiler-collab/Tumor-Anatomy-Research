# Impact kit: what you can do with this research today

Written 2026-10-07. Everything below rests on results that are already verified (`docs/CLAIMS_LEDGER.md`,
`docs/VERIFICATION_RERUN.md`, `docs/REGISTRATION_AUDIT.md`). The drafts are yours to edit and send under your
name. Nothing here has been sent or posted.

---

## 1. What you can say today

In plain language, each with its evidence:

1. **Checking that a deconvolution result "matches the anatomy" catches broken methods, but cannot pick the
   accurate one.** Broken methods score 0.37-0.40 against 0.71-0.97 for real ones. Yet the anatomy ranking
   does not follow accuracy against DNA: rho 0.08 in TCGA, 0.12 in an independent cohort (CPTAC), 0.29 with
   17 methods. The bar was 0.60.
2. **No method gets the lymphocytes right in glioma.** About 20 methods tested; at most a few put T cells
   above B cells, and none robustly. Flow cytometry and two single-cell atlases put T cells 9-43 times above
   B. Repairing the software did not change this (0 of 10 repaired runs).
3. **Tumour and total immune content can be trusted for ranking samples.** The accuracy ranking of methods
   is the same whichever DNA test measures it (0.80-0.99). On the frozen signature it carries over to a new
   cohort (0.68); on the donor-level reference it does not (0.13).
4. **Software quality is as big a factor as the choice of method.** A flawed reimplementation turned an
   accurate method (BayesPrism, 0.72) into a useless one (0.16). A numerical bug made the DWLS package skip
   73 of 122 samples; a one-line fix removes it.
5. **The expensive measurement is not automatically the truth.** Single-nucleus RNA from the same tissue did
   not track tumour DNA (0.13), while a cheaper methylation measurement with complete probes did (0.92).

---

## 2. A fix the DWLS community can use: a comment for GitHub issue #12

**Where:** https://github.com/dtsoucas/DWLS/issues/12, which is open with no cause or fix posted. Posting it
is your call, under your account.

> **Cause and a solution-preserving fix for "constraints are inconsistent, no solution!"**
>
> We hit this error on 73 of 122 bulk samples (Ivy GAP glioblastoma, counts-per-million inputs) and traced it.
>
> **Cause.** `solveOLSInternal()` passes `D <- t(S) %*% S` and `d <- t(S) %*% B` to `quadprog::solve.QP`
> unscaled. On CPM-scale inputs, entries of `D` reach ~1e8-1e12. quadprog's absolute tolerances then report
> the problem as infeasible, although the only constraint is x >= 0 and the signature was well conditioned
> (kappa = 11.7). `solveDampenedWLSj()` already rescales by `norm(D, "2")` before calling `solve.QP`;
> `solveOLSInternal()` does not.
>
> **Fix (no package change).** Divide the signature and the bulk by the same constant before calling DWLS:
>
> ```r
> k <- max(signature)
> props <- solveDampenedWLS(signature / k, bulk / k)
> ```
>
> **Why the answer does not change.** Every DWLS solution is invariant to a common rescaling:
> - the least-squares minimiser is unchanged;
> - the dampening weights are rescaled by their minimum;
> - the `lm` coefficients in `findDampeningConstant()` are invariant;
> - the proportions are normalised at the end.
>
> **Measured.** On the 122 real samples, the first step failed on 73 as is and on 0 after rescaling. Where
> both versions solved, the full solutions agree to 3e-16.
>
> An equivalent in-package fix is to give `solveOLSInternal()` the same scaling as `solveDampenedWLSj()`:
> `sc <- norm(D, "2"); solve.QP(D/sc, d/sc, A, bzero)`.

---

## 3. A data request that would settle an open question: email to the CPTAC glioblastoma authors

**To:** the lead contact of Wang et al., *Cancer Cell* 39:509 (2021), doi:10.1016/j.ccell.2021.01.006, listed
in the paper as Li Ding (lding@wustl.edu). Check the paper for the current contact first.

> Subject: Request for snRNA-seq cell-type annotations, CPTAC glioblastoma (Wang et al., Cancer Cell 2021)
>
> Dear Dr Ding,
>
> I am a high-school researcher studying how to validate bulk RNA cell-type deconvolution without
> single-cell measurements. I use your CPTAC glioblastoma data from the GDC.
>
> Single-nucleus RNA-seq was made for 17 of the tumours from the same cryopulverised material as the bulk
> RNA. For the 15 I could score, I compared each tumour's single-nucleus composition with its whole-genome
> tumour purity. The cell labels came
> from GDC's automated per-sample clusters, named by marker panels. With those labels, the nucleus tumour
> fraction did not track DNA purity (Spearman 0.13). Bulk deconvolution did (median 0.45, 18 tumours).
>
> Your paper describes manually curated cell-type annotations on the merged dataset, which GDC does not
> distribute. Would you be willing to share the per-nucleus cell-type labels? They would show whether curated
> annotation recovers per-sample agreement with DNA. I would acknowledge the source in any write-up.
>
> Thank you for considering this.
>
> [Your name, school, contact]

---

## 4. A one-page guide for labs without single-cell or flow cytometry

Share it with a teacher, a mentor, or a lab that runs deconvolution. It is drawn from `docs/ACCURACY_FACTORS.md`
and the acceptance tests in `docs/CLAIMS_LEDGER.md`.

**How to use bulk-RNA deconvolution when you cannot afford to check it directly (glioma evidence):**
1. **Report tumour content and total immune content.** Treat lymphocyte subtypes as unmeasured.
2. **Use estimates to rank samples, not as exact percentages.** They compress the true range.
3. **Choose the method by its accuracy against a DNA measurement in a similar tissue, on the same reference
   build.** In glioma, on the frozen signature, a Bayesian regression, nu-SVR, a CIBERSORTx-style solver and DWLS (each
   this project's implementation) led in two cohorts.
4. **Run the published package.** If you use a reimplementation, check it against the package on real data.
5. **Check every sample got an estimate, and rescale inputs if a solver fails** (see the DWLS fix above).
6. **Use an anatomy or plausibility check only to throw out broken output,** never as proof of accuracy.
7. **Re-fit under several loss scales.** Do not trust a compartment that moves.
8. **Validate once against a DNA purity** (copy number, or methylation with complete probe coverage), not
   against single-cell composition.

---

## 5. Start writing now (the prose is yours)

The order, and where everything is:
- **The spine:** `docs/MANUSCRIPT_BLUEPRINT.md` §2 (R1-R4), with the CPTAC replication inside R2 and the
  repairs inside R3.
- **Facts and figures per section:** `docs/MANUSCRIPT.md`, where every number is generated from artefacts.
  Fill each `> WRITE:` beat.
- **Figures:** `docs/FIGURES.md` (16 figures with captions). Figure 1 is `Figure_design`.
- **Limitations:** blueprint §4b, 11 items.
- **What is verified, and what is new against prior work:** `docs/CLAIMS_LEDGER.md`.
- **Why each step was done:** blueprint §3 (the dated decision chain, 22 rows).

**Still running, and not blocking any of this:** the last verification steps (the LGG donor-level arm and
the extension panel), then the final report. They can only confirm or flag; every claim above is already
verified or explicitly labelled.

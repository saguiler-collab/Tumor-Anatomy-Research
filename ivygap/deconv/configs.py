"""
configs.py — the published configuration of every method, and every deviation from it.

THE RISK THIS ADDRESSES
-----------------------
Anatomy_Test.md's risk table, on method setup fairness:

    Any of these methods can be made to look bad by running it wrong.
    → One shared signature, published configs for every method, and defaults from each
      tool's own documentation — no per-method tuning.

The shared signature was enforced by the equal-footing certificate. The *configs* were
not published anywhere: every method carried its parameters on the object and no
artefact recorded them, so "we used each tool's defaults" was an assertion a reader had
to take on faith while reading a ranking those defaults produced.

It also stopped being true. Making DWLS run at all required two departures from its
published defaults — its differential-expression cutoffs are opened, because this
pipeline has already pre-filtered the gene space, and its signature build is capped per
cell type, because the published path costs half an hour per gene set here. Both are
defensible and both were decided before any DWLS score existed. Neither is defensible
if it is invisible.

So this module writes what each method was actually run with, marks anything that
differs from that method's own default, and states why. A reviewer should be able to
find every deviation in one file rather than by reading the source.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path

#: method name -> (what was changed, why, when the decision was made).
#: Only genuine departures from a tool's documented defaults belong here.
DECLARED_DEVIATIONS: dict[str, list[dict]] = {
    "bisque": [
        {
            "parameter": "the scale of the single-cell matrix this HARNESS exports, and "
                         "hence whether the package estimates cell size itself",
            "published_default": "raw per-cell counts",
            "used": ("cells normalised to 1e6 each before export, the same treatment every "
                     "cell-consuming method receives. NO cell-size factors are applied on "
                     "top: docs/OPEN_DEFECTS.md D12 measures the central conversion as the "
                     "IDENTITY, because run_benchmark rebuilds the reference from the "
                     "already-normalised matrix without passing cell_totals, so cell_size "
                     "comes out 1e6 for every type (spread 1.0000). An earlier version of "
                     "this entry said the factors were 'applied once and centrally'; they "
                     "are applied nowhere"),
            "why": ("Equal footing: one set of factors applied identically to all 16 "
                    "methods. Declared for this method even though its own convention is "
                    "UNMEASURED — the cell-size probe cannot read it, because on a roster "
                    "with six types absent from the mixture this package assigns ~75% of "
                    "its mass to the absent types and the remaining ratio measures nothing "
                    "about a cell-size convention. So this entry records what the harness "
                    "does, and explicitly does NOT claim the package would have corrected "
                    "cell size itself. An earlier version of docs/OPEN_DEFECTS.md D1 "
                    "asserted this package was double-corrected; that claim is withdrawn "
                    "as unmeasured. See scripts/verify_cell_size_semantics.py."),
            "decided": ("from the harness's own code path, which is a fact about this "
                        "pipeline and not about any score"),
        },
    ],
    "bayesprism": [
        {
            "parameter": "the scale of the single-cell matrix this HARNESS exports, and "
                         "hence whether the package estimates cell size itself",
            "published_default": "raw per-cell counts",
            "used": ("cells normalised to 1e6 each before export, the same treatment every "
                     "cell-consuming method receives. NO cell-size factors are applied on "
                     "top: docs/OPEN_DEFECTS.md D12 measures the central conversion as the "
                     "IDENTITY, because run_benchmark rebuilds the reference from the "
                     "already-normalised matrix without passing cell_totals, so cell_size "
                     "comes out 1e6 for every type (spread 1.0000). An earlier version of "
                     "this entry said the factors were 'applied once and centrally'; they "
                     "are applied nowhere"),
            "why": ("Equal footing, as above. This package's own convention is UNMEASURED "
                    "under the production export — the probe shows ~31% leakage onto roster "
                    "types absent from the mixture, which makes the reading unusable. Its "
                    "PYTHON reimplementation measures 0.7500 on the same probe, i.e. mRNA "
                    "share, which if it carries over to the R package would make this "
                    "project's conversion the first and correct one. Not assumed. See "
                    "scripts/verify_cell_size_semantics.py and docs/OPEN_DEFECTS.md D1."),
            "decided": ("from the harness's own code path, which is a fact about this "
                        "pipeline and not about any score"),
        },
    ],
    "scdc_ensemble": [
        {
            "parameter": "the scale of the single-cell matrix this HARNESS exports, and "
                         "hence whether the package estimates cell size itself",
            "published_default": ("raw per-cell counts, from which the package derives its "
                                  "own per-cell-type mRNA content — `music_basis` computes "
                                  "`M.S`, the mean library size per cell type, and "
                                  "`SCDC_basis` derives one the same way"),
            "used": ("cells normalised to 1e6 each before export, so every cell type's "
                     "mean library size is identical ACROSS THE FULL TRANSCRIPTOME. The "
                     "export is then restricted to the bulk's marker gene space, over which "
                     "per-type library sizes are NOT identical — measured at 2.17x spread on "
                     "the real atlas (Tumor 29,999 to Oligodendrocyte 65,003) — so whether "
                     "this package's own cell-size estimate is neutralised is UNRESOLVED. "
                     "CORRECTED 2026-09-14, docs/OPEN_DEFECTS.md D12: this project's central "
                     "conversion is the IDENTITY on the reference every score was computed "
                     "against, so it substitutes NOTHING and converts nothing. "
                     "run_benchmark rebuilds the reference from the already-normalised matrix "
                     "without passing cell_totals, so build_reference falls back to "
                     "expression.sum(axis=0) and cell_size comes out 1e6 for all eight types "
                     "(spread 1.0000). An earlier version of this entry claimed the factors "
                     "were 'applied once, centrally, to every method'. They are applied "
                     "nowhere."),
            "why": ("Equal footing WAS the intent: one set of factors, computed one way, "
                    "applied identically to all 16 methods rather than each package deriving "
                    "its own. D12 shows the intent was never realised — no factors are "
                    "applied. What the normalisation DOES still do is keep a deeply sequenced "
                    "cell from dominating its type's mean profile, which is a real effect and "
                    "the reason the deviation is still declared. MEASURED, not assumed: on a probe "
                    "where two types differ 3x in mRNA and are mixed 50/50 by cell count "
                    "(true cell fraction 0.500, true mRNA fraction 0.750), this package "
                    "returns 0.7501 when handed the probe's FULL gene space and 0.5336 when "
                    "handed raw per-cell counts — so it does convert cell size when its "
                    "input lets it estimate one. The probe's full-space export is NOT "
                    "what production writes: production exports a marker subset, over "
                    "which per-type library sizes are measured to span 2.17x on the real "
                    "atlas. So this deviation is DECLARED but its consequence is "
                    "UNRESOLVED, and no claim is made here that the package's own "
                    "conversion is neutralised. See "
                    "scripts/verify_cell_size_semantics.py, "
                    "results/cell_size_semantics_{normalized,raw}.json, and "
                    "docs/OPEN_DEFECTS.md D1."),
            "decided": ("from the measurement above, which is independent of any ACS or "
                        "benchmark score — the probe is synthetic and its truth is known "
                        "by construction"),
        },
    ],
    "scdc": [
        {
            "parameter": "the scale of the single-cell matrix this HARNESS exports, and "
                         "hence whether the package estimates cell size itself",
            "published_default": ("raw per-cell counts, from which the package derives its "
                                  "own per-cell-type mRNA content — `music_basis` computes "
                                  "`M.S`, the mean library size per cell type, and "
                                  "`SCDC_basis` derives one the same way"),
            "used": ("cells normalised to 1e6 each before export, so every cell type's "
                     "mean library size is identical ACROSS THE FULL TRANSCRIPTOME. The "
                     "export is then restricted to the bulk's marker gene space, over which "
                     "per-type library sizes are NOT identical — measured at 2.17x spread on "
                     "the real atlas (Tumor 29,999 to Oligodendrocyte 65,003) — so whether "
                     "this package's own cell-size estimate is neutralised is UNRESOLVED. "
                     "CORRECTED 2026-09-14, docs/OPEN_DEFECTS.md D12: this project's central "
                     "conversion is the IDENTITY on the reference every score was computed "
                     "against, so it substitutes NOTHING and converts nothing. "
                     "run_benchmark rebuilds the reference from the already-normalised matrix "
                     "without passing cell_totals, so build_reference falls back to "
                     "expression.sum(axis=0) and cell_size comes out 1e6 for all eight types "
                     "(spread 1.0000). An earlier version of this entry claimed the factors "
                     "were 'applied once, centrally, to every method'. They are applied "
                     "nowhere."),
            "why": ("Equal footing WAS the intent: one set of factors, computed one way, "
                    "applied identically to all 16 methods rather than each package deriving "
                    "its own. D12 shows the intent was never realised — no factors are "
                    "applied. What the normalisation DOES still do is keep a deeply sequenced "
                    "cell from dominating its type's mean profile, which is a real effect and "
                    "the reason the deviation is still declared. MEASURED, not assumed: on a probe "
                    "where two types differ 3x in mRNA and are mixed 50/50 by cell count "
                    "(true cell fraction 0.500, true mRNA fraction 0.750), this package "
                    "returns 0.7501 with the production export and 0.5336 with a "
                    "raw export, so it does convert when it can and production "
                    "prevents it ONLY when the export spans the full gene space, which "
                    "production's does not. UNRESOLVED for production. See "
                    "scripts/verify_cell_size_semantics.py, "
                    "results/cell_size_semantics_{normalized,raw}.json, and "
                    "docs/OPEN_DEFECTS.md D1."),
            "decided": ("from the measurement above, which is independent of any ACS or "
                        "benchmark score — the probe is synthetic and its truth is known "
                        "by construction"),
        },
    ],
    "music": [
        {
            "parameter": "the scale of the single-cell matrix this HARNESS exports, and "
                         "hence whether the package estimates cell size itself",
            "published_default": ("raw per-cell counts, from which the package derives its "
                                  "own per-cell-type mRNA content — `music_basis` computes "
                                  "`M.S`, the mean library size per cell type, and "
                                  "`SCDC_basis` derives one the same way"),
            "used": ("cells normalised to 1e6 each before export, so every cell type's "
                     "mean library size is identical ACROSS THE FULL TRANSCRIPTOME. The "
                     "export is then restricted to the bulk's marker gene space, over which "
                     "per-type library sizes are NOT identical — measured at 2.17x spread on "
                     "the real atlas (Tumor 29,999 to Oligodendrocyte 65,003) — so whether "
                     "this package's own cell-size estimate is neutralised is UNRESOLVED. "
                     "CORRECTED 2026-09-14, docs/OPEN_DEFECTS.md D12: this project's central "
                     "conversion is the IDENTITY on the reference every score was computed "
                     "against, so it substitutes NOTHING and converts nothing. "
                     "run_benchmark rebuilds the reference from the already-normalised matrix "
                     "without passing cell_totals, so build_reference falls back to "
                     "expression.sum(axis=0) and cell_size comes out 1e6 for all eight types "
                     "(spread 1.0000). An earlier version of this entry claimed the factors "
                     "were 'applied once, centrally, to every method'. They are applied "
                     "nowhere."),
            "why": ("Equal footing WAS the intent: one set of factors, computed one way, "
                    "applied identically to all 16 methods rather than each package deriving "
                    "its own. D12 shows the intent was never realised — no factors are "
                    "applied. What the normalisation DOES still do is keep a deeply sequenced "
                    "cell from dominating its type's mean profile, which is a real effect and "
                    "the reason the deviation is still declared. MEASURED, not assumed: on a probe "
                    "where two types differ 3x in mRNA and are mixed 50/50 by cell count "
                    "(true cell fraction 0.500, true mRNA fraction 0.750), this package "
                    "returns 0.7501 when handed the probe's FULL gene space and 0.5001 when "
                    "handed raw per-cell counts — so it does convert cell size when its "
                    "input lets it estimate one. The probe's full-space export is NOT "
                    "what production writes: production exports a marker subset, over "
                    "which per-type library sizes are measured to span 2.17x on the real "
                    "atlas. So this deviation is DECLARED but its consequence is "
                    "UNRESOLVED, and no claim is made here that the package's own "
                    "conversion is neutralised. See "
                    "scripts/verify_cell_size_semantics.py, "
                    "results/cell_size_semantics_{normalized,raw}.json, and "
                    "docs/OPEN_DEFECTS.md D1."),
            "decided": ("from the measurement above, which is independent of any ACS or "
                        "benchmark score — the probe is synthetic and its truth is known "
                        "by construction"),
        },
    ],
    "cibersortx": [
        {
            "parameter": "batch-correction mode (B-mode vs S-mode)",
            "published_default": ("no single default; Newman et al. (2019) choose per "
                                  "configuration. Their Supplementary Table 1d records "
                                  "the choice for every deconvolution in the paper."),
            "used": "B-mode",
            "why": ("Supplementary Table 1d shows a consistent split: every signature "
                    "derived from 10x Chromium and applied to a bulk RNA-seq mixture is "
                    "deconvolved in S-MODE, while B-mode is used for SMART-Seq2- and "
                    "microarray-derived signatures (LM22) against bulk. The split is "
                    "mechanistic — droplet 10x data carries a strong 3' bias and UMI "
                    "counting, so it sits further from bulk RNA-seq than full-length "
                    "SMART-Seq2 does. This project's reference is GBmap, measured at "
                    "87.1% 10x (214,284 cells 3' v2, 49,262 3' v3, 31,316 5' v1 of "
                    "338,564) against 2.7% Smart-seq2, and the mixtures are Ivy GAP "
                    "bulk RNA-seq. By the paper's own practice that is an S-mode "
                    "configuration, and this implementation runs B-mode. S-mode's "
                    "per-cell-type expression adjustment is not implemented here, so "
                    "the row is B-mode and says so rather than claiming to be the "
                    "mode the authors would have used."),
            "decided": ("from Supplementary Table 1d and the atlas assay counts, both "
                        "external to any score this project computed"),
        },
    ],
    "quantiseq": [
        {
            "parameter": "the bulk gene space this HARNESS hands the method",
            "published_default": "the full expression matrix",
            "used": "the full shared space, NOT the shared marker subset",
            "why": ("Every other method solves against this project's signature matrix, "
                    "so the marker subset — top-N genes per cell type ranked against "
                    "that signature, plus its markers — is the shared gene space and "
                    "equal footing holds. quanTIseq does not solve against it: it ships "
                    "TIL10 and ignores the signature it is handed. Ranking genes "
                    "against a reference it never reads kept 34 of TIL10's 138 "
                    "signature genes (24.6%), and quanTIseq answered with macrophages "
                    "at 100% of cells in the median sample and T cells identically zero "
                    "in all 122. On the full space it finds 136 of 138 (98.6%) and "
                    "returns a GBM-plausible 15.7% median myeloid fraction. The subset "
                    "was not equal footing for this method; it was mutilation. EPIC is "
                    "also signature-only but DOES read the supplied signature, so it "
                    "keeps the subset and is untouched by this."),
            "decided": ("from the signature-gene recovery rate, before any quanTIseq "
                        "ACS existed — the 2026-09-05T2154 run scored it NaN"),
        },
        {
            "parameter": "cell-size (mRNA) correction",
            "published_default": "scale_mRNA = TRUE, quanTIseq's own mRNA scaling",
            "used": "quanTIseq's own scaling only; this project's is not applied on top",
            "why": ("The pipeline's cell-size correction and quanTIseq's scale_mRNA are "
                    "the same correction. Applying both would double-correct, and this "
                    "project's version renormalises to sum 1, which on a method "
                    "covering four immune types would assert the tumour is entirely "
                    "immune. Keeping the published default and skipping ours leaves "
                    "quanTIseq on its own documented scale."),
            "decided": "from the two definitions, before any quanTIseq ACS existed",
        },
    ],
    "dwls": [
        {
            "parameter": "buildSignatureMatrixMAST(diff.cutoff, pval.cutoff)",
            "published_default": "diff.cutoff = 0.5, pval.cutoff = 0.01",
            "used": "diff.cutoff = 0, pval.cutoff = 1 (filter fully open)",
            "why": ("Equal footing requires every method to receive the same "
                    "pre-filtered gene space, so the informative genes have already "
                    "been chosen identically for all methods. DWLS's own differential-"
                    "expression pass is then a second filter, not a second opinion: at "
                    "the published cutoffs it left 31 genes for 8 cell types and the "
                    "driver refused to deconvolve on them."),
            "decided": "from the gene count, before any DWLS score existed",
        },
        {
            "parameter": "cells used for the signature build",
            "published_default": "all cells in the reference",
            "used": "<= 250 cells per cell type",
            "why": ("buildSignatureMatrixMAST costs ~30 minutes per gene set on this "
                    "hardware regardless of the DE method. A signature is a per-cell-"
                    "type summary, so type proportions are irrelevant to it and capping "
                    "per type spends the budget where it buys accuracy."),
            "decided": "from measured runtime, before any DWLS score existed",
        },
        {
            "parameter": "per-sample failure handling",
            "published_default": "an error aborts the call",
            "used": "a sample that fails to solve becomes NA and is counted",
            "why": ("solveDampenedWLS fails on particular mixtures with 'NA/NaN "
                    "argument'. One such sample was sending the whole cohort to the "
                    "Python fallback. NA is already what this pipeline means by a "
                    "failed sample and it is reported in n_failed_samples."),
            "decided": "from the failure mode, before any DWLS score existed",
        },
    ],
}


def _defaults_for(method) -> dict:
    """The method's own constructor defaults, read from its signature."""
    try:
        sig = inspect.signature(type(method).__init__)
    except (TypeError, ValueError):                          # pragma: no cover
        return {}
    return {name: p.default for name, p in sig.parameters.items()
            if p.default is not inspect.Parameter.empty and name != "self"}


def describe(methods: list, implementations: dict | None = None) -> list[dict]:
    """
    One record per method: what it ran as, what it was configured with, and which of
    those settings differ from the method's own defaults.
    """
    implementations = implementations or {}
    out = []
    for m in methods:
        inner = getattr(m, "fallback", m)
        defaults = _defaults_for(inner)
        params = dict(getattr(inner, "params", {}) or {})
        # RMethod carries wrapper bookkeeping in params; it is not a model setting.
        for k in ("r_method", "allow_fallback"):
            params.pop(k, None)

        differs = {k: {"used": v, "default": defaults.get(k)}
                   for k, v in params.items()
                   if k in defaults and v != defaults[k]}

        out.append({
            "method": m.name,
            "class": type(inner).__name__,
            "family": getattr(m, "family", None),
            "implementation": implementations.get(m.name, "python"),
            "parameters": {k: _jsonable(v) for k, v in params.items()},
            "differs_from_class_defaults": {
                k: {"used": _jsonable(d["used"]), "default": _jsonable(d["default"])}
                for k, d in differs.items()},
            "declared_deviations_from_published_defaults":
                DECLARED_DEVIATIONS.get(m.name, []),
        })
    return out


def _jsonable(v):
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    return str(v)


def write(methods: list, path: Path, implementations: dict | None = None) -> list[dict]:
    records = describe(methods, implementations)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({
        "what_this_is": (
            "Every method's configuration as actually run, with anything differing from "
            "its own defaults marked. Declared deviations from a tool's PUBLISHED "
            "defaults are listed per method with the reason and when the decision was "
            "made — the protocol's fairness rule is defaults with no per-method tuning, "
            "so a departure that is invisible is not defensible."),
        "n_methods": len(records),
        "methods_with_declared_deviations": sorted(
            r["method"] for r in records
            if r["declared_deviations_from_published_defaults"]),
        "configs": records,
    }, indent=2))
    return records

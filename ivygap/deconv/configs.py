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

# External validation: does ACS reproduce a protein ground truth?

**Run 2026-09-15 against `prespecified/imc_anchored_prediction.md`, which was committed before
any number below existed.**

## Why this is the only test of its kind in the project

Every other result here is internal. The registered primary outcome compares the ACS ranking
with a pseudobulk-accuracy ranking, and **both arms use GBmap** — the ACS arm deconvolves Ivy
GAP against it, the accuracy arm builds its mixtures from its cells. That shared dependence is
correction **C4**. So "anatomy agrees with truth" has, until now, meant "two GBmap-based
rankings agree with each other."

Ajaib et al. (**Neuro-Oncology** 2023;25(7):1236–1248) provide an outside anchor: single-cell
resolution **imaging mass cytometry**, 33 antibodies, on ten IDHwt GBM samples with matched bulk
RNA-seq on the same tissue. The core of their result is a near-controlled contrast of exactly
the kind this project uses — **one algorithm, MCPcounter, differing only in its marker set**:

| approach | r vs IMC (immune) | r vs IMC (neoplastic) |
|---|---|---|
| MCPcounter + GBM-specific markers (MCP_GBM) | **0.37** | **0.43** |
| MCPcounter + default markers | 0.27 | — |
| MCPcounter + **GBmap**-derived markers | **0.06** | 0.22 |
| CIBERSORTx | 0.05 | 0.02 |

GBmap is the atlas this entire leaderboard rests on, and against protein ground truth its
markers scored **0.06** — near chance, and worst of the four.

## The pre-specified test, and its result

> **Prediction, fixed in advance:** ACS(MCP_GBM) > ACS(MCP_GBmap) on the constraints both arms
> can express.

Genuine **MCPcounter 1.2.0**, one bulk matrix, one log transform applied identically to every
arm, the frozen constraint file, the same 122 anatomic samples. Only the marker set varies.
Scored on Tumor + Macrophage_Microglia (C1, C5, C6, C7) — the two cell types the IMC panel also
measured.

| arm | ACS | 95% CI | null p |
|---|---|---|---|
| MCP_GBM (Ajaib markers) | **0.5128** | [0.3571, 0.6857] | 0.0694 |
| MCP_GBmap | **0.6667** | [0.5250, 0.8333] | 0.0025 |

Paired over the same nine tumours: **delta +0.1389, 95% CI [+0.0444, +0.2500], p = 0.0090.**

> ### The prediction FAILED, in the direction named in advance as most damaging.

ACS ranked the GBmap-derived markers **above** the GBM-tissue-specific markers, and
significantly so. The pre-specification said: *"MCP_GBmap scoring at or near the top would be
the most damaging outcome. It is the marker set IMC ranks near chance (0.06), and it is derived
from the atlas this whole project uses. ACS preferring it would show ACS rewarding agreement
with its own reference rather than with tissue."*

It is worth adding that **MCP_GBM does not beat its own permutation null** (p = 0.0694). The
marker set that IMC ranks best achieves no significant anatomic concordance at all.

The three secondary comparisons point the same way. Full ACS ordering:
**MCP_default (0.9091) > MCP_GBmap (0.8485/0.8000) > MCP_GBM_moreno (0.7385) > MCP_GBM.**
IMC's ordering is **MCP_GBM > MCP_default > MCP_GBmap**. ACS puts the IMC winner last.

## The named alternative explanation, tested — and it is real

MCPcounter's estimate is a **mean over marker genes**, so more markers means a smoother score,
which could satisfy ordinal constraints more often for reasons unrelated to biology. MCP_GBmap
carries **200 markers per type**; Ajaib's set carries about **79**. That confound was named
before the run and is tested by rebuilding GBmap's markers at matched density and changing
nothing else:

| GBmap markers per type | n genes | ACS | 95% CI | null p | paired delta vs MCP_GBM |
|---|---|---|---|---|---|
| 200 *(the pre-specified arm)* | 1,600 | 0.6667 | [0.5250, 0.8333] | 0.0025 | **+0.1389**, p = 0.0090 |
| **79** *(matched to Ajaib)* | 632 | **0.3846** | [0.2955, 0.4872] | **0.4357** | **−0.1222**, p = 0.1056 |
| 15 | 120 | 0.6923 | [0.5000, 0.8864] | 0.0016 | +0.1722, p = 0.0560 |

**At matched marker density the direction reverses.** GBmap falls to 0.3846 — below MCP_GBM's
0.5128, and no longer distinguishable from its own permutation null (p = 0.4357). The
pre-specified comparison was not density-matched, so what it measured was partly marker count.

**This does not rescue the prediction.** The pre-specified arm was the 200-per-type one, and it
failed. Presenting the matched-density arm as the result would be choosing the analysis after
seeing it. The honest verdict on the registered question is **INCONCLUSIVE**: the test as
specified failed, the named confound is real, and at matched density the direction agrees with
IMC without reaching significance on nine tumours.

## The finding that is robust, and it is the more serious one

Look down that table again. **Same atlas. Same algorithm. Same bulk. Same constraints. Same
tumours. Only the number of markers per cell type changed** — and ACS moved:

**0.3846 → 0.6923**, with non-overlapping confidence intervals ([0.2955, 0.4872] against
[0.5000, 0.8864]), and from *indistinguishable from the null* (p = 0.4357) to *significant*
(p = 0.0016).

That swing is **0.31**. The entire published leaderboard spans roughly 0.60 to 1.00. So a
nuisance parameter — how many marker genes the reference contributes, with no change whatever to
the biology, the data or the constraints — moves ACS by about three-quarters of the range across
which it is being used to rank fifteen methods.

The movement is also **non-monotone** (0.667 at 200, 0.385 at 79, 0.692 at 15), which rules out
a simple "more markers is smoother" story and points instead at instability: nine evaluable
tumours and 31 constraint–tumour pairs is not enough to hold an ordinal score steady against a
change in gene space.

**For the registered question this is the result that matters most.** It is not about which
marker set wins. It is that ACS differences of the size that separate methods on the leaderboard
can be produced by a preprocessing choice that has nothing to do with the tissue. Combined with
D8 (one constraint, C6, carries most of the ranking) and D14/D16 (the ordering is substantially
a property of the reference atlas), the picture is consistent: **ACS is measuring something
real — every arm but one beats its permutation null — but it does not have the resolution to
rank methods, and it can be moved by choices the method itself did not make.**

## What is not claimed

- **Not** that GBM-specific markers are worse. At matched density they score higher.
- **Not** that the IMC study is wrong. It has protein ground truth on matched tissue; this
  project has ordinal anatomy on different tissue.
- **Not** that ACS is worthless. It separates real methods from negative controls under every
  reference tested, and most arms here beat their permutation null.
- **Not** a reason to change anything about ACS. The constraint file is frozen and hashed, and
  it has not been touched. Retuning it because a result came out badly is the specific failure
  this project's invariants forbid.

## Limits, stated as they were in advance

Nine evaluable tumours, 31 pairs in the primary comparison. Two cohorts, two modalities. My
GBmap markers are derived by this project's own selection rule; Ajaib's came from the GBmap
preprint's supplementary data — the same atlas, a different derivation. MCPcounter returns
abundance scores rather than proportions, so this class of method **can** be scored by ACS and
**cannot** be scored by the pseudobulk accuracy arm; the registered primary outcome could never
have contained one. The IMC data itself is not public through the vendored repository, so the
external anchor is the authors' published verdict, not a re-analysis.

Artefacts: `results/imc_anchored_test.json`, `results/imc_anchored_robustness.json`.
Scripts: `scripts/imc_anchored_test.py`, `scripts/imc_anchored_robustness.py`,
`scripts/mcp_marker_source.R`.

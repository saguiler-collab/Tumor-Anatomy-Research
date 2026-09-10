# The second tissue

`Anatomy_Test.md` step 7: *one tissue is a case study; two is a method.*

This is the single largest thing standing between what exists now and the claim the
protocol sets out to make. It is not a formality, and the Ivy GAP result itself is the
argument for it: **MuSiC leads the leaderboard by one constraint–tumour pair out of 57**,
and six methods' bootstrap intervals contain its score. Nine evaluable tumours can
separate good methods from bad; they cannot identify the best one, and no further
analysis inside Ivy GAP changes that. Only more tissue does.

---

## Status

| Piece | State |
|---|---|
| Constraint file | **Drafted, verified, NOT frozen** — `ivygap/anatomic/constraints_brain.py` |
| Region mapping | Verified against the Allen ontology (graph_id 10) 2026-09-10 |
| Bulk data | Identified, not downloaded — `scripts/fetch_allen.py` |
| Single-cell reference | Identified, not downloaded — Siletti et al. 2023 via CELLxGENE |
| Platform decision | **Open — yours to make.** See below. |
| Registration | Not submitted |

---

## 1 · What was blocking it, and what unblocked it

The download was never the blocker. Three things were, and two are now resolved.

**The constraints did not transfer.** `constraints.py` is about CT, IT, LE, MVP and PAN —
glioblastoma structures. Normal brain has none of them, and editing that file to fit
would change `constraints.freeze_hash()` and void every GBM result computed under it.
*Resolved:* `constraints_brain.py` is a separate file with its own hash and its own
registration.

**The roster did not transfer.** Normal cortex is mostly neurons; this project's roster
has no neuron column, and one could not simply be added because **the GBmap atlas
contains 22 neurons**. *Resolved:* Siletti et al. (2023) *Science*, the Human Brain Cell
Atlas v1.0, has 2.48M neurons and 888k non-neuronal cells on 10x 3′ v3 — the same
platform family as GBmap, so the pipeline's normalisation and its B-mode/S-mode reasoning
carry over unchanged. Per-supercluster files mean the 33 GB neuron archive is not
required.

**Platform.** Still open. See §4.

---

## 2 · The constraint set

Six constraints over four regions, deliberately mirroring the GBM set's shape.

| id | claim | kind | wt | why it is knowable in advance |
|---|---|---|---|---|
| B1 | Oligodendrocyte: WM > CTX | pairwise | 1 | White matter *is* myelinated axon tracts, and CNS myelin is made by oligodendrocytes. The analogue of C3. |
| B2 | Oligodendrocyte: WM is the maximum | maximum | 1 | Stricter form of B1, as C4 is of C3. Subcortical grey is heavily traversed by fibre tracts and could steal it. |
| B3 | Neuron: CTX > WM | pairwise | 1 | Cortex holds the somata; white matter carries their axons. The analogue of C1. |
| B4 | Neuron: CBC is the maximum | maximum | 1 | Cerebellar granule cells are the densest neuronal population in the CNS — ~69 of the brain's ~86 billion neurons (Azevedo 2009). |
| B5 | Neuron: WM < CTX < CBC | monotone | **2** | A conjunction of two orderings, double-weighted for the same reason C7 is. |
| B6 | Endothelial: CTX > WM | pairwise | 1 | Grey matter is more densely capillarised. Deliberately the weakest member, so the set is not composed only of its easiest claims. |

Weighted denominator: **7 per donor.**

**Declared untestable, in advance and with reasons:** astrocytes (grey/white density does
not order consistently across the literature — a constraint would score the reference's
astrocyte definition, not the anatomy), microglia (near-uniform in healthy brain, so any
ordering sits inside what deconvolution can resolve — the same reason T cells were
excluded from the GBM set), and OPC (directionally likely but not near-definitional, and
OPC/oligodendrocyte boundaries differ between references).

### Regions, verified against the Allen ontology

| id | Allen node | path |
|---|---|---|
| CTX | `Cx` | `Br / GM / Tel / Cx` |
| WM | `WM` | `Br / WM` (all divisions) |
| CBC | `CbCx` | `Br / GM / MET / Cb / CbCx` |
| SUBC | `BG`, `DiE` | `Br / GM / Tel / CxN / BG`, `Br / GM / DiE` |

The ontology's own depth-1 split is `Br / GM | WM | SS` — the grey/white distinction B1–B3
turn on is the ontology's, not one imposed on it. `SUBC` exists so B2 and B4 have a real
competitor rather than a walkover.

---

## 3 · What must happen in what order

The ordering is the whole scientific claim and cannot be repaired afterwards.

```
1. Settle the claims and freeze          FROZEN = True, then freeze_hash()
2. Register the hash publicly            a SECOND registration, separate from the GBM one
3. THEN download and build               Allen bulk + Siletti reference
4. THEN deconvolve and score             first time any brain output is looked at
```

Steps 1 and 2 must precede step 4. Nothing prevents step 3 from happening earlier —
downloading data is not looking at deconvolution output — but doing 1 and 2 first is
cleaner and costs nothing.

`freeze_hash()` raises while `FROZEN = False`, so step 4 cannot silently happen out of
order.

**Do not fold this into the GBM registration.** It needs its own constraint file and its
own hash, and two independent timestamps are stronger than one amended document.

---

## 4 · The platform decision — this one is yours

The Allen Human Brain Atlas offers two products and they trade off directly.

| | microarray | RNA-seq |
|---|---|---|
| Donors | **6** | 2 |
| Size | 1.57 GB | 44 MB |
| Samples | ~400–1000 sites per brain | 240 total |
| Platform | Microarray | RNA-seq |
| Weighted denominator at 7/donor | **42** | 14 |
| Compatible with this pipeline | **No** — CPM throughout, RNA-seq-derived signature | **Yes** |

Ivy GAP's weighted denominator is 65. So:

- **RNA-seq (2 donors, denominator 14)** is platform-clean and *badly* underpowered —
  materially weaker than the cohort you already have. It would not answer the question
  the second tissue exists to answer.
- **Microarray (6 donors, denominator 42)** is comparable in power to Ivy GAP but crosses
  a platform boundary.

**My reading:** the microarray path is the better study, and the platform gap is a
feature rather than only a cost. Cross-platform deconvolution is a known hard problem
that the field cares about; CIBERSORTx exists largely to address it, quanTIseq ships an
`is_arraydata` flag for it, and LM22 — the most widely used signature in the field — is
itself microarray-derived. A second tissue that also crosses platforms tests whether
anatomic concordance survives a harder setting, which is a stronger claim than repeating
the easy one.

It is a study-design decision with a real cost, so it is stated here rather than made
quietly by whichever file gets downloaded first.

---

## 5 · What will need building

None of this is blocked, and none of it should start before §3 steps 1–2.

- **A brain scorer path.** `acs.py` scores (tumour, structure); brain is (donor, region).
  The shape is identical — nested units, ordinal claims — so this is a parameter, not a
  rewrite.
- **A region assigner.** Map each Allen sample's fine structure to CTX/WM/CBC/SUBC by
  ancestry in the ontology. Deterministic and testable against the graph.
- **A Siletti reference builder.** `build_from_h5ad` already does donor-balanced
  subsampling; the work is a cell-type mapping onto the six-type roster, the same shape
  as `GBMAP_CELL_TYPE_MAP`.
- **A platform decision recorded as a declared deviation**, whichever way it goes.

---

## 6 · Honest expectations

Two constraints (B1, B3) are near-tautological and every competent method should satisfy
them — that is by design, exactly as C1 and C3 were. The discrimination will come from
B2, B4 and B5, the maximum and monotone claims.

If **the controls clear this set as easily as the real methods**, the finding is that the
constraint set does not discriminate in normal brain, and the protocol's instruction is
explicit: report it, publish the constraint file unchanged, and state what a harder set
would need. That branch is as publishable as the other and is registered in advance for
the same reason it was the first time.

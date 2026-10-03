"""
anatomic_auc.py -- the per-method AUC, computed so it can never be mistaken for ACS.

ACS and an AUC answer different questions. This script reports both side by side,
each under its own name, so neither is ever quoted as the other.

  ACS (registered, primary).  For each pre-registered constraint and each tumour,
      collapse to ONE composition per (tumour, structure), score 1 if the predicted
      ordering holds between those means and 0 if not (exact ties score 0), and take
      the weighted fraction satisfied. A threshold statistic on means.

  AUC_anat (exploratory, post-registration).  The same seven constraints, the same
      122 anatomic samples, the same tumour unit and the same null -- but graded, and at
      sample level. For each ordered structure pair (hi, lo) a constraint implies and
      each tumour that sampled both, the per-tumour AUC is the Mann-Whitney probability
      P(f(a) > f(b)) + 0.5 P(f(a) = f(b)) over every pair of that tumour's own samples,
      a in hi and b in lo. Pairs are NEVER formed across tumours: pooling would let
      between-tumour differences pass for anatomy and let a tumour with 18 samples
      outvote one with 2 (CLAUDE.md: never aggregate over nested samples).
        pairwise C1 C2 C3 C5 C6 -> the one pair
        maximum  C4             -> mean of MVP vs each other sampled structure (>= 2,
                                   the same evaluability rule ACS uses)
        monotone C7             -> mean of the two steps IT>LE and CT>IT (all three
                                   structures required, as in ACS)
      Method AUC_anat = sum over (constraint, tumour) of weight x AUC / sum of weights:
      ACS with its 0/1 indicator replaced by a graded probability. Chance is 0.5.

      Exactly where the two part, each pinned by a test:
        * ties -- ACS scores an exact tie 0, an AUC 0.5; so AUC_strict (ties = 0) and the
          tied-pair share are reported too;
        * conjunctions -- ACS scores C4 and C7 as all-or-nothing (MVP above EVERY other
          structure; BOTH steps of the chain), AUC_anat averages their comparisons and so
          gives partial credit (MVP above three of four structures: 0.75, not 0);
        * replication -- ACS compares means, AUC_anat every within-tumour sample pair.
      With one sample per (tumour, structure) and no ties, the five pairwise constraints
      (C1 C2 C3 C5 C6) coincide exactly; C4 and C7 do not.
      The chance levels therefore differ: 0.5 for AUC_anat by construction, about 0.37
      for ACS (its permutation-null mean), so the same number means different things on
      the two scales.

  C_purity (truth side).  Neither statistic above sees a ground truth. The accuracy side
      of this study is ABSOLUTE DNA purity; to put it on the same 0.5-is-chance scale it
      is re-expressed here as Harrell's concordance for a continuous outcome: the
      probability that a method orders two TCGA-GBM samples' tumour content the same way
      ABSOLUTE does (ties in the estimate 0.5). A re-expression of the registered
      Spearman rho, not a new yardstick.

Controls, both required before anything is read off:
  * positive: ACS recomputed here from the same estimate files must reproduce the
    registered leaderboard, and the registered agreement (rho = 0.081) must reproduce
    through the registered `agreement.test_agreement`. Either failing aborts.
  * negative: the two broken-input controls (random proportions, shuffled signature)
    are scored with AUC_anat exactly as the methods are. A graded score that rated them
    like real methods would be measuring nothing.

Exploratory. It selects nothing and changes no registered statistic.

    python3 scripts/anatomic_auc.py            # writes results/anatomic_auc.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config  # noqa: E402
from ivygap.anatomic import agreement  # noqa: E402
from ivygap.anatomic import constraints as K  # noqa: E402
from ivygap.anatomic.acs import _acs_from_table, _score_table, tumor_structure_means  # noqa: E402

LEADERBOARD = config.ANATOMIC_DIR / "acs_leaderboard.csv"
PURITY_REPORT = config.RESULTS_DIR / "absolute_purity_yardstick.json"
PURITY_PER_SAMPLE = config.RESULTS_DIR / "absolute_purity_per_sample.csv"
REGISTERED_AGREEMENT = config.RESULTS_DIR / "yardstick_agreement.json"
OUT = config.RESULTS_DIR / "anatomic_auc.json"


# --------------------------------------------------------------------------- pairs

def constraint_pairs(c: K.Constraint) -> list[tuple[str, str]]:
    """The ordered (hi, lo) structure pairs whose AUCs make up one constraint."""
    if c.kind == "pairwise":
        return [(c.structures[0], c.structures[1])]
    if c.kind == "maximum":
        return [(c.structures[0], s) for s in c.among if s != c.structures[0]]
    if c.kind == "monotone":
        s = c.structures
        return [(s[i + 1], s[i]) for i in range(len(s) - 1)]
    raise ValueError(f"unknown constraint kind {c.kind!r}")


def pair_auc(hi: np.ndarray, lo: np.ndarray) -> tuple[float, float, float]:
    """(AUC with ties 0.5, AUC with ties 0, share of tied pairs) for two samples."""
    hi = hi[np.isfinite(hi)]
    lo = lo[np.isfinite(lo)]
    if hi.size == 0 or lo.size == 0:
        return float("nan"), float("nan"), float("nan")
    d = hi[:, None] - lo[None, :]
    gt = float((d > 0).mean())
    eq = float((d == 0).mean())
    return gt + 0.5 * eq, gt, eq


def tumour_constraint_auc(vals: np.ndarray, labels: np.ndarray,
                          c: K.Constraint) -> tuple[float, float, float] | None:
    """
    One constraint, one tumour. None when the tumour did not sample what the constraint
    needs -- the same evaluability rule as `acs.evaluate_constraint`, so the two
    statistics are computed over the same (constraint, tumour) pairs.
    """
    present = {s for s in set(labels) if np.isfinite(vals[labels == s]).any()}
    pairs = constraint_pairs(c)
    if c.kind == "pairwise":
        if not set(c.structures) <= present:
            return None
    elif c.kind == "maximum":
        if c.structures[0] not in present:
            return None
        pairs = [(h, l) for h, l in pairs if l in present]
        if len(pairs) < 2:
            return None
    elif c.kind == "monotone":
        if not set(c.structures) <= present:
            return None
    got = [pair_auc(vals[labels == h], vals[labels == l]) for h, l in pairs]
    return tuple(float(np.mean([g[k] for g in got])) for k in range(3))


def auc_table(estimate: pd.DataFrame, tumours: np.ndarray, labels: np.ndarray) -> pd.DataFrame:
    """Long table, one row per evaluable (constraint, tumour): the graded twin of ACS's."""
    rows = []
    for c in K.CONSTRAINTS:
        v_all = estimate[c.cell_type].to_numpy(dtype="float64")
        for t in sorted(set(tumours)):
            m = tumours == t
            got = tumour_constraint_auc(v_all[m], labels[m], c)
            if got is None:
                continue
            rows.append({"constraint": c.id, "tumor_id": t, "weight": c.weight,
                         "auc": got[0], "auc_strict": got[1], "tied": got[2]})
    return pd.DataFrame(rows, columns=["constraint", "tumor_id", "weight", "auc",
                                       "auc_strict", "tied"])


def weighted(table: pd.DataFrame, col: str) -> float:
    if table.empty or table["weight"].sum() <= 0:
        return float("nan")
    return float((table[col] * table["weight"]).sum() / table["weight"].sum())


# --------------------------------------------------------------------------- null

class FastAUC:
    """
    AUC_anat for any relabelling, in matrix form, for the permutation loop.

    `auc_table` above is the readable definition and is ~100x too slow to call 10,000
    times per method. Here each tumour's pairwise comparison matrix G[i, j] = 1 / 0.5 / 0
    (f_i >, =, < f_j) is built once per cell type; a relabelling only changes the one-hot
    structure matrix L, and every structure-pair sum at once is L'GL. A test asserts this
    returns exactly what `auc_table` returns, on real estimates and on relabellings.
    """

    def __init__(self, estimate: pd.DataFrame, tumours: np.ndarray):
        self.structures = list(config.PRIMARY_STRUCTURES)
        self.code = {s: i for i, s in enumerate(self.structures)}
        self.cells = sorted({c.cell_type for c in K.CONSTRAINTS})
        self.idx = [np.where(tumours == t)[0] for t in sorted(set(tumours))]
        self.G, self.F = [], []
        for idx in self.idx:
            g_t, f_t = {}, {}
            for c in self.cells:
                v = estimate[c].to_numpy(dtype="float64")[idx]
                f = np.isfinite(v)
                d = np.where(f[:, None] & f[None, :], v[:, None] - v[None, :], np.nan)
                g_t[c] = np.where(np.isnan(d), 0.0, (d > 0) + 0.5 * (d == 0))
                f_t[c] = f.astype("float64")
            self.G.append(g_t)
            self.F.append(f_t)

    def score(self, labels: np.ndarray) -> float:
        num = den = 0.0
        for k, idx in enumerate(self.idx):
            L = np.zeros((idx.size, len(self.structures)))
            L[np.arange(idx.size), [self.code[s] for s in labels[idx]]] = 1.0
            S, N = {}, {}
            for c in self.cells:
                Lf = L * self.F[k][c][:, None]
                n = Lf.sum(0)
                S[c] = Lf.T @ self.G[k][c] @ Lf
                N[c] = n
            for con in K.CONSTRAINTS:
                n, s_ = N[con.cell_type], S[con.cell_type]
                pairs = constraint_pairs(con)
                if con.kind == "maximum":
                    if n[self.code[con.structures[0]]] == 0:
                        continue
                    pairs = [(h, l) for h, l in pairs if n[self.code[l]] > 0]
                    if len(pairs) < 2:
                        continue
                elif any(n[self.code[s]] == 0 for s in con.structures):
                    continue
                a = np.mean([s_[self.code[h], self.code[l]] /
                             (n[self.code[h]] * n[self.code[l]]) for h, l in pairs])
                num += con.weight * a
                den += con.weight
        return num / den if den > 0 else float("nan")

    def score_batch(self, codes: list[np.ndarray]) -> np.ndarray:
        """
        `score` for B relabellings at once. codes[k] is (B, n_k): structure codes of
        tumour k's samples under each relabelling. Same arithmetic as `score`, with the
        permutation axis carried through einsum; a test pins the two together.
        """
        B = codes[0].shape[0]
        num, den = np.zeros(B), np.zeros(B)
        n_s = len(self.structures)
        for k in range(len(self.idx)):
            L = np.eye(n_s)[codes[k]]                          # (B, n_k, n_s)
            S, N = {}, {}
            for c in self.cells:
                Lf = L * self.F[k][c][None, :, None]
                N[c] = Lf.sum(1)                                # (B, n_s)
                S[c] = np.einsum("bis,ij,bjt->bst", Lf, self.G[k][c], Lf)
            for con in K.CONSTRAINTS:
                n, s_ = N[con.cell_type], S[con.cell_type]
                pairs = constraint_pairs(con)
                if con.kind == "maximum":
                    h = self.code[con.structures[0]]
                    terms, ok_l = [], []
                    for _, l in pairs:
                        l = self.code[l]
                        ok = n[:, l] > 0
                        ok_l.append(ok)
                        terms.append(np.where(ok, s_[:, h, l] /
                                              np.maximum(n[:, h] * n[:, l], 1e-300), 0.0))
                    n_ok = np.sum(ok_l, axis=0)
                    use = (n[:, h] > 0) & (n_ok >= 2)
                    a = np.sum(terms, axis=0) / np.maximum(n_ok, 1)
                else:
                    use = np.all([n[:, self.code[s]] > 0 for s in con.structures], axis=0)
                    a = np.mean([s_[:, self.code[h], self.code[l]] /
                                 np.maximum(n[:, self.code[h]] * n[:, self.code[l]], 1e-300)
                                 for h, l in pairs], axis=0)
                num += np.where(use, con.weight * a, 0.0)
                den += np.where(use, con.weight, 0.0)
        return np.divide(num, den, out=np.full(B, np.nan), where=den > 0)


def auc_null(estimate: pd.DataFrame, tumours: np.ndarray, labels: np.ndarray,
             n_permutations: int, seed: int, batch: int = 1000) -> np.ndarray:
    """
    AUC_anat under within-tumour permutation of structure labels -- the ACS null:
    each tumour's labels are reshuffled among its own samples, preserving its
    compositions and its block count per structure.
    """
    rng = np.random.default_rng(seed)
    fast = FastAUC(estimate, tumours)
    base = [np.array([fast.code[s] for s in labels[idx]]) for idx in fast.idx]
    out = []
    for start in range(0, n_permutations, batch):
        b = min(batch, n_permutations - start)
        codes = [rng.permuted(np.tile(c, (b, 1)), axis=1) for c in base]
        out.append(fast.score_batch(codes))
    return np.concatenate(out)


# --------------------------------------------------------------------------- truth

def harrell_c(estimate: np.ndarray, truth: np.ndarray) -> tuple[float, int]:
    """P(estimate orders a pair as truth does), over pairs with distinct truth."""
    ok = np.isfinite(estimate) & np.isfinite(truth)
    e, t = estimate[ok], truth[ok]
    de = e[:, None] - e[None, :]
    dt = t[:, None] - t[None, :]
    usable = np.triu(dt != 0, k=1)
    s = np.sign(de[usable]) * np.sign(dt[usable])
    n = int(usable.sum())
    return (float(((s > 0).sum() + 0.5 * (s == 0).sum()) / n) if n else float("nan")), n


# --------------------------------------------------------------------------- main

def load_scoring_set() -> tuple[pd.DataFrame, list[str]]:
    man = pd.read_csv(config.SAMPLE_MANIFEST_PATH, sep="\t", dtype={"sample_id": str,
                                                                     "patient_id": str})
    man = man.set_index("sample_id")
    keep = man["structure"].isin(config.PRIMARY_STRUCTURES) & \
        man["is_anatomic_study"].astype(bool)
    return man[keep], list(man.index[keep])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-permutations", type=int, default=10_000)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()

    manifest, scored = load_scoring_set()
    board = pd.read_csv(LEADERBOARD).set_index("method")

    results, failures = {}, []
    for m in board.index:
        f = config.ESTIMATES_DIR / f"ivygap_{m}.csv"
        if not f.exists():
            failures.append(f"{m}: estimate file missing ({f.name})")
            continue
        est = pd.read_csv(f, index_col=0)
        est.index = est.index.astype(str)
        est = est.reindex([s for s in scored if s in est.index])
        if len(est) != len(scored):
            failures.append(f"{m}: {len(scored) - len(est)} of {len(scored)} scored "
                            f"samples absent from the estimate file")
            continue

        # Positive control: ACS from these exact rows must equal the registered value.
        acs_here = _acs_from_table(_score_table(tumor_structure_means(est, manifest)))
        acs_reg = float(board.loc[m, "acs"])
        if not np.isclose(acs_here, acs_reg, rtol=0, atol=1e-9):
            failures.append(f"{m}: ACS recomputed {acs_here!r} != registered {acs_reg!r}")
            continue

        tumours = manifest.loc[est.index, "patient_id"].astype(str).to_numpy()
        labels = manifest.loc[est.index, "structure"].astype(str).to_numpy()
        table = auc_table(est, tumours, labels)
        auc = weighted(table, "auc")
        null = auc_null(est, tumours, labels, args.n_permutations, config.RANDOM_SEED)
        null = null[np.isfinite(null)]
        p = float(((null >= auc).sum() + 1) / (null.size + 1))

        per_c = {}
        for c in K.CONSTRAINTS:
            sub = table[table["constraint"] == c.id]
            per_c[c.id] = {"description": c.describe(), "n_tumours": int(len(sub)),
                           "auc": float(sub["auc"].mean()) if len(sub) else float("nan"),
                           "auc_strict": (float(sub["auc_strict"].mean()) if len(sub)
                                          else float("nan"))}
        results[m] = {
            "acs": acs_reg, "acs_reproduced": True,
            "auc_anat": auc,
            "auc_anat_strict": weighted(table, "auc_strict"),
            "tied_pair_share": weighted(table, "tied"),
            "n_constraint_tumour_pairs": int(len(table)),
            "n_constraint_tumour_pairs_acs": int(board.loc[m, "n_constraint_tumor_pairs"]),
            "null_mean": float(null.mean()), "null_p": p,
            "is_control": bool(board.loc[m, "is_control"]),
            "degenerate": bool(board.loc[m, "degenerate"]),
            "comparable": bool(board.loc[m, "comparable"]),
            "per_constraint": per_c,
        }
        print(f"{m:28s} ACS {acs_reg:.4f}  AUC_anat {auc:.4f}  strict "
              f"{results[m]['auc_anat_strict']:.4f}  ties {results[m]['tied_pair_share']:.3f}"
              f"  null {null.mean():.3f}  p {p:.4g}", flush=True)

    if failures:
        print("ABORT -- a control failed or an input is missing:\n  " + "\n  ".join(failures))
        return 2

    # ---- truth side: concordance with ABSOLUTE, on the registered per-sample table
    ps = pd.read_csv(PURITY_PER_SAMPLE, index_col=0)
    rep = json.load(open(PURITY_REPORT))["methods"]
    c_purity, rho_purity = {}, {}
    for m in ps.columns.drop("absolute_purity"):
        c, n_pairs = harrell_c(ps[m].to_numpy(float), ps["absolute_purity"].to_numpy(float))
        c_purity[m] = {"c": c, "n_pairs": n_pairs}
        rho_purity[m] = float(rep[m]["spearman_vs_purity"])

    # ---- the registered agreement must reproduce before the AUC variant is read
    reg = json.load(open(REGISTERED_AGREEMENT))["arm_2_absolute_purity"]["all_methods"]
    acs_scores = {m: results[m]["acs"] for m in rho_purity if m in results}
    control = agreement.test_agreement(acs_scores, rho_purity, "absolute_purity",
                                       higher_truth_is_better=True)
    if round(control.rho, 3) != round(reg["spearman"], 3) or control.n_methods != reg["n_methods"]:
        print(f"ABORT -- registered agreement did not reproduce: rho {control.rho:.4f} "
              f"n {control.n_methods} vs registered {reg['spearman']} n {reg['n_methods']}")
        return 2

    auc_scores = {m: results[m]["auc_anat"] for m in acs_scores}
    variant = agreement.test_agreement(auc_scores, rho_purity, "absolute_purity",
                                       higher_truth_is_better=True)
    c_scores = {m: c_purity[m]["c"] for m in acs_scores}
    variant_c = agreement.test_agreement(auc_scores, c_scores, "absolute_purity_concordance",
                                         higher_truth_is_better=True)
    ms = sorted(acs_scores)
    rho_c_vs_rho = stats.spearmanr([c_scores[m] for m in ms], [rho_purity[m] for m in ms])[0]
    real = [m for m in results if not results[m]["is_control"] and results[m]["comparable"]]
    rho_acs_auc = stats.spearmanr([results[m]["acs"] for m in real],
                                  [results[m]["auc_anat"] for m in real])[0]
    best_control = max(results[m]["auc_anat"] for m in results if results[m]["is_control"])
    worst_real = min(results[m]["auc_anat"] for m in real)

    report = {
        "what_this_is": ("ACS and a per-method AUC, computed on identical inputs and kept "
                         "under separate names. Exploratory and post-registration: selects "
                         "nothing, changes no registered statistic."),
        "definitions": {
            "acs": "registered: weighted fraction of (constraint, tumour) pairs whose "
                   "ordering holds between per-(tumour, structure) MEANS; ties score 0",
            "auc_anat": "exploratory: the same pairs, graded -- per-tumour Mann-Whitney "
                        "probability over that tumour's own samples (ties 0.5), averaged "
                        "with ACS's weights; C4 and C7 averaged over their comparisons "
                        "where ACS requires all of them; chance 0.5; within-tumour "
                        "permutation null",
            "auc_anat_strict": "auc_anat with ties scored 0, ACS's tie convention",
            "c_purity": "truth side: Harrell's concordance of each method's TCGA-GBM "
                        "Tumor fraction with ABSOLUTE purity; chance 0.5",
            "what_none_of_them_is": ("ACS and AUC_anat measure agreement with anatomic "
                                     "expectation. Neither sees ground truth, so neither is "
                                     "accuracy. C_purity is accuracy against DNA."),
        },
        "n_samples_scored": len(scored),
        "n_permutations": args.n_permutations,
        "seed": config.RANDOM_SEED,
        "controls": {
            "acs_reproduced_for_every_method": True,
            "registered_agreement_reproduced": {"rho": round(control.rho, 4),
                                                "n": control.n_methods},
            "negative_controls_auc_anat": {m: round(results[m]["auc_anat"], 4)
                                           for m in results if results[m]["is_control"]},
            "margin_worst_real_minus_best_control": round(worst_real - best_control, 4),
        },
        "methods": results,
        "c_purity": c_purity,
        "rank_agreement": {
            "acs_vs_auc_anat_spearman_comparable_methods": round(float(rho_acs_auc), 4),
            "n_comparable_methods": len(real),
            "c_purity_vs_rho_purity_spearman": round(float(rho_c_vs_rho), 4),
        },
        "agreement_with_auc_in_place_of_acs": {
            "vs_rho_purity": variant.to_dict(),
            "vs_c_purity": variant_c.to_dict(),
            "registered_for_comparison": {"rho": reg["spearman"], "p": reg["p"],
                                          "ci": reg["ci"], "n": reg["n_methods"]},
        },
    }
    Path(args.out).write_text(json.dumps(report, indent=2, default=float))
    print(f"\nregistered agreement reproduced: rho {control.rho:+.3f} (n {control.n_methods})")
    print(f"AUC_anat in ACS's place:         rho {variant.rho:+.3f} p {variant.p_value:.3f} "
          f"CI [{variant.ci_low:+.3f}, {variant.ci_high:+.3f}]")
    print(f"ACS vs AUC_anat ranking (comparable methods): rho {rho_acs_auc:+.3f}")
    print(f"C_purity vs rho_purity ranking: rho {rho_c_vs_rho:+.3f}")
    print(f"margin worst real - best control (AUC_anat): {worst_real - best_control:+.4f}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

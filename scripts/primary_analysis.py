"""
primary_analysis.py

Primary held-out analysis for the frozen five-model RAMP score, computed
directly from the raw per-call response files in heldout_results/. This
is the script behind the manuscript's primary endpoint and its
robustness analyses:

  - claim extraction (GO:\\d{7} regex, set semantics within a response)
  - frozen score S = n_models + max_prompt_consistency + 5*rago_support
    (FROZEN_PROTOCOL.md Sec. 3 -- implemented as specified, not tuned)
  - precision at 10%-coverage under the whole-score-group tie rule
    (FROZEN_PROTOCOL.md Sec. 4), and the 100%-coverage baseline
  - paired protein-level bootstrap for (precision@10% - baseline)
    (B=2,000, seed 42; FROZEN_PROTOCOL.md Sec. 9)
  - precision by individual score and cumulatively by threshold
  - organism and namespace breakdowns
  - component ablations (ranking by sub-combinations of the components)
  - leave-one-model-out, including each reduced pool's size and baseline
  - concentration of correct claims across proteins (count, max, Gini)
  - exact-10% tie-breaking sensitivity value (descriptive only; not the
    primary endpoint)
  - consistency check of the read-only claim->score lookup used by every
    external-validation analysis
  - aggregate inference accounting (calls, tokens) from the raw records

Inputs (all already in this repository; nothing is fetched or generated):
  heldout_results/*.jsonl                 raw five-model responses
  heldout_benchmark_243.json              the 243-protein cohort
  scripts/heldout_sample_250.json         ground truth (experimental GO terms)
  heldout_rag_candidates.json             top-15 BLAST candidate terms per
                                          protein/namespace (run_blast_retrieval.py)
  openai_external_validation/analysis/frozen_claim_score_lookup_READONLY.json
                                          compared against, never modified

No LLM call, BLAST search, or network access is made.

Usage:
  python scripts/primary_analysis.py
  python scripts/primary_analysis.py --out primary_analysis/primary_analysis_report.json
"""

import argparse
import json
import re
from pathlib import Path
from collections import defaultdict

import numpy as np

HERE = Path(__file__).resolve().parent
RAGO = HERE.parent
GO_PATTERN = re.compile(r"GO:\d{7}")
NS_SHORT = {"molecular_function": "MF", "biological_process": "BP", "cellular_component": "CC"}
FROZEN_PROMPTS = {"P1_zeroshot", "P2_constrained", "P3_fewshot", "P4_cot"}
FROZEN_SCORE_LOOKUP_PATH = RAGO / "openai_external_validation" / "analysis" / "frozen_claim_score_lookup_READONLY.json"

COVERAGE = 0.10
RAGO_WEIGHT = 5
B_BOOTSTRAP = 2000
SEED = 42


def load_inputs():
    cohort = {r["accession"] for r in json.load(open(RAGO / "heldout_benchmark_243.json"))}
    heldout = json.load(open(RAGO / "scripts" / "heldout_sample_250.json"))
    org_of = {r["accession"]: r["organism"] for r in heldout}
    gt = {(r["accession"], ns): set(t["go_id"] for t in r["go_terms"].get(s, []))
          for r in heldout for ns, s in NS_SHORT.items()}
    rag = json.load(open(RAGO / "heldout_rag_candidates.json"))
    rag_terms = {(acc, ns): set(t["go_id"] for t in rag[acc].get(s, []))
                 for acc in rag for ns, s in NS_SHORT.items()}

    # Final successful record per (accession, namespace, prompt, model). The
    # 18 calls that failed once (HTTP 502) and succeeded on retry keep both
    # attempts in the raw file; only the successful attempt is used.
    latest = {}
    for f in sorted((RAGO / "heldout_results").glob("*.jsonl")):
        for line in open(f):
            r = json.loads(line)
            if r["accession"] not in cohort or r["prompt_id"] not in FROZEN_PROMPTS:
                continue
            if r.get("success"):
                latest[(r["accession"], r["namespace"], r["prompt_id"], r["model"])] = r
    return cohort, org_of, gt, rag_terms, latest


def build_claims(latest):
    """claim (acc, ns, go_id) -> {model: set(prompts that produced it)}"""
    claim_models = defaultdict(lambda: defaultdict(set))
    for (acc, ns, prompt, model), r in latest.items():
        for go_id in set(GO_PATTERN.findall(r.get("response") or "")):
            claim_models[(acc, ns, go_id)][model].add(prompt)
    return claim_models


def components(claims, claim_models, rag_terms, models):
    """n_models, max_prompt_consistency, rago_support restricted to `models`.
    Claims produced by none of `models` get n_models == 0 and are dropped by
    the caller (they do not exist in that configuration's claim pool)."""
    M = np.array([sum(1 for m in claim_models[c] if m in models) for c in claims])
    P = np.array([max([len(claim_models[c][m]) for m in claim_models[c] if m in models] or [0])
                  for c in claims])
    R = np.array([int(c[2] in rag_terms.get((c[0], c[1]), set())) for c in claims])
    return M, P, R


def gini(x):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def precision_at_coverage(score, y, coverage=COVERAGE):
    """Whole-score-group rule: walk down from the highest score, never
    splitting a tied group, until at least coverage*N claims are included."""
    target = coverage * len(score)
    for s in np.sort(np.unique(score))[::-1]:
        sel = score >= s
        if sel.sum() >= target:
            return {"threshold": int(s), "n_selected": int(sel.sum()), "n_correct": int(y[sel].sum()),
                    "precision": float(y[sel].mean()), "coverage": float(sel.sum() / len(score))}


def paired_bootstrap(score, y, prot, B=B_BOOTSTRAP, seed=SEED):
    accs = np.unique(prot)
    idx_of = {a: np.where(prot == a)[0] for a in accs}
    rng = np.random.default_rng(seed)
    p10, base, diff = [], [], []
    for _ in range(B):
        ii = np.concatenate([idx_of[a] for a in rng.choice(accs, size=len(accs), replace=True)])
        a = precision_at_coverage(score[ii], y[ii])["precision"]
        b = float(y[ii].mean())
        p10.append(a); base.append(b); diff.append(a - b)
    ci = lambda v: [float(x) for x in np.percentile(v, [2.5, 97.5])]
    return {"B": B, "seed": seed, "n_proteins": int(len(accs)),
            "precision_at_10_ci": ci(p10), "baseline_ci": ci(base), "paired_difference_ci": ci(diff)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RAGO / "primary_analysis" / "primary_analysis_report.json"))
    args = ap.parse_args()

    cohort, org_of, gt, rag_terms, latest = load_inputs()
    models = sorted({k[3] for k in latest})
    claim_models = build_claims(latest)
    claims = sorted(claim_models)
    y = np.array([int(c[2] in gt.get((c[0], c[1]), set())) for c in claims])
    prot = np.array([c[0] for c in claims])
    M, P, R = components(claims, claim_models, rag_terms, set(models))
    S = M + P + RAGO_WEIGHT * R

    report = {"inputs": {"n_proteins": len(cohort), "models": models,
                         "n_successful_calls": len(latest), "n_claims": len(claims)}}

    # ---------------- Primary endpoint ----------------
    report["primary"] = precision_at_coverage(S, y)
    report["baseline"] = {"n_correct": int(y.sum()), "n_claims": len(y), "precision": float(y.mean())}
    report["primary"]["bootstrap"] = paired_bootstrap(S, y, prot)

    # ---------------- Score--correctness table ----------------
    report["by_score"] = []
    for s in range(int(S.min()), int(S.max()) + 1):
        at, ge = S == s, S >= s
        report["by_score"].append({"score": s, "claims": int(at.sum()), "correct": int(y[at].sum()),
                                   "precision": float(y[at].mean()) if at.any() else None,
                                   "cum_claims": int(ge.sum()), "cum_precision": float(y[ge].mean())})

    # ---------------- Organism / namespace breakdowns ----------------
    # Each subgroup applies the same whole-group 10%-coverage rule to its own pooled claims.
    org = np.array([org_of[a] for a in prot])
    ns = np.array([c[1] for c in claims])
    report["by_organism"] = {o: {**precision_at_coverage(S[org == o], y[org == o]),
                                 "baseline": float(y[org == o].mean())} for o in sorted(set(org))}
    report["by_namespace"] = {n: {**precision_at_coverage(S[ns == n], y[ns == n]),
                                  "baseline": float(y[ns == n].mean())} for n in NS_SHORT}

    # ---------------- Concentration across proteins ----------------
    # Per-protein counts over all 243 cohort proteins (zeros included).
    sel = S >= report["primary"]["threshold"]
    accs = sorted(cohort)
    n_correct = defaultdict(int)
    n_sel = defaultdict(int)
    for a, s_, t in zip(prot, sel, y):
        n_sel[a] += int(s_)
        n_correct[a] += int(s_ and t)
    report["concentration"] = {
        "proteins_with_correct": sum(1 for a in accs if n_correct[a] > 0),
        "max_correct_per_protein": max(n_correct.values()),
        "gini_correct": gini([n_correct[a] for a in accs]),
        "gini_selected": gini([n_sel[a] for a in accs])}

    # ---------------- Exact-10% tie-breaking sensitivity ----------------
    # Descriptive sensitivity only (NOT the primary endpoint, NOT prespecified):
    # take every claim above the cut-off group, then fill to ceil(10% * N)
    # from the cut-off group in lexicographic (accession, full namespace
    # name, GO ID) order. The value depends on this arbitrary ordering.
    thr = report["primary"]["threshold"]
    n_target = int(np.ceil(COVERAGE * len(S)))
    above = np.where(S > thr)[0]
    fill = np.where(S == thr)[0][:n_target - len(above)]  # claims are sorted
    picked = np.concatenate([above, fill])
    report["exact_coverage_tiebreak_sensitivity"] = {
        "n_selected": int(len(picked)), "n_correct": int(y[picked].sum()),
        "precision": float(y[picked].mean()), "tie_order": "lexicographic (accession, full namespace name, GO ID)"}

    # ---------------- Component ablations ----------------
    ablations = {"n_models": M, "max_prompt_consistency": P, "n_models+max_prompt_consistency": M + P,
                 "n_models+rago_support": M + RAGO_WEIGHT * R,
                 "max_prompt_consistency+rago_support": P + RAGO_WEIGHT * R, "full_combined_score": S}
    report["ablations"] = {k: precision_at_coverage(v, y) for k, v in ablations.items()}

    # ---------------- Leave-one-model-out ----------------
    # Removing a model removes the claims only that model produced, so each
    # reduced configuration has its own (smaller) claim pool and baseline.
    report["leave_one_model_out"] = {}
    for drop in models:
        Mm, Pm, Rm = components(claims, claim_models, rag_terms, set(models) - {drop})
        keep = Mm > 0
        Sm = (Mm + Pm + RAGO_WEIGHT * Rm)[keep]
        report["leave_one_model_out"][drop] = {
            **precision_at_coverage(Sm, y[keep]),
            "pool_claims": int(keep.sum()), "pool_correct": int(y[keep].sum()),
            "pool_baseline": float(y[keep].mean())}

    # ---------------- Read-only lookup consistency ----------------
    lookup = json.load(open(FROZEN_SCORE_LOOKUP_PATH))
    rebuilt = {f"{c[0]}|{c[1]}|{c[2]}": int(s) for c, s in zip(claims, S)}
    report["lookup_check"] = {"lookup_entries": len(lookup), "rebuilt_entries": len(rebuilt),
                              "identical": lookup == rebuilt,
                              "n_key_or_score_mismatches": sum(1 for k in set(lookup) | set(rebuilt)
                                                               if lookup.get(k) != rebuilt.get(k))}

    # ---------------- Inference accounting (aggregate, from raw usage fields) ----------------
    # These are token counts reported by the inference gateway for the
    # final successful call per key. They are not a cost or runtime figure.
    acct = {}
    for m in models:
        recs = [r for k, r in latest.items() if k[3] == m]
        acct[m] = {"calls": len(recs),
                   "prompt_tokens": sum((r.get("usage") or {}).get("prompt_tokens", 0) or 0 for r in recs),
                   "completion_tokens": sum((r.get("usage") or {}).get("completion_tokens", 0) or 0 for r in recs)}
    tot = {k: sum(v[k] for v in acct.values()) for k in ("calls", "prompt_tokens", "completion_tokens")}
    tot["calls_per_protein"] = tot["calls"] / len(cohort)
    report["inference_accounting"] = {"per_model": acct, "total": tot}

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(report, open(out, "w"), indent=2)

    p, b, bs = report["primary"], report["baseline"], report["primary"]["bootstrap"]
    print(f"claims {len(claims)}  calls {len(latest)}")
    print(f"PRIMARY  {p['n_correct']}/{p['n_selected']} = {p['precision']:.2%}  "
          f"(coverage {p['coverage']:.2%}, score >= {p['threshold']})")
    print(f"BASELINE {b['n_correct']}/{b['n_claims']} = {b['precision']:.2%}")
    print(f"95% CI precision {[round(x*100, 2) for x in bs['precision_at_10_ci']]}  "
          f"baseline {[round(x*100, 2) for x in bs['baseline_ci']]}  "
          f"paired diff (pp) {[round(x*100, 2) for x in bs['paired_difference_ci']]}")
    for k, v in report["ablations"].items():
        print(f"ablation {k}: {v['precision']:.2%}")
    for k, v in report["leave_one_model_out"].items():
        print(f"LOMO -{k}: pool {v['pool_claims']} correct {v['pool_correct']} base {v['pool_baseline']:.2%} "
              f"selected {v['n_correct']}/{v['n_selected']} = {v['precision']:.2%}")
    for k, v in {**report["by_organism"], **report["by_namespace"]}.items():
        print(f"subgroup {k}: {v['precision']:.2%} (baseline {v['baseline']:.2%})")
    print("concentration:", report["concentration"])
    print("exact-10% tie-break sensitivity:", report["exact_coverage_tiebreak_sensitivity"])
    print("lookup check:", report["lookup_check"])
    print("inference totals:", tot)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

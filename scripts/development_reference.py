"""
development_reference.py

Development-set reference value for the frozen five-model score
(manuscript "Development-set reference" and FROZEN_PROTOCOL.md Sec. 8):
precision at 10%-coverage on PROBE's 1,000-protein development benchmark,
its protein-level bootstrap CI, the 100%-coverage baseline, and the
development-set leave-one-model-out check used to justify the model set.

This is the procedure that produced the value recorded in the frozen
protocol before the held-out evaluation (0.0595, 95% CI [0.0563, 0.0631]).
It used B = 1,000 bootstrap replicates (the prespecified minimum), not the
B = 2,000 used for the held-out analyses; pass --B 2000 to recompute with
the held-out setting ([0.0563, 0.0628]).

Inputs:
  <probe_results>/*.jsonl   PROBE development responses for the five frozen
                            models, from https://github.com/utpoln/PROBE
                            (probe_results/); ground truth is read from each
                            record's ground_truth field
  rag_candidates.json       top-15 BLAST candidate terms per development
                            protein/namespace (this repository)

No LLM call, BLAST search, or network access is made.

Usage:
  python scripts/development_reference.py --probe-results /path/to/PROBE/probe_results
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
GEN_PROMPTS = {"P1_zeroshot", "P2_constrained", "P3_fewshot", "P4_cot"}
RULE_MODELS = ["mistral-large-123b", "qwen2.5-72b", "google_gemma-4-31b", "mixtral-8x7b", "mistral-7b"]
NS_SHORT = {"molecular_function": "MF", "biological_process": "BP", "cellular_component": "CC"}
MAX_SCORE = 14
RAGO_WEIGHT = 5
SEED = 42


def load_support(probe_results):
    raw_support = defaultdict(set)  # (acc, ns, go_id) -> {(model, prompt)}
    gt = {}
    for model in RULE_MODELS:
        for line in open(Path(probe_results) / f"{model}.jsonl"):
            rec = json.loads(line)
            if rec.get("prompt_id") not in GEN_PROMPTS or not rec.get("success"):
                continue
            acc, ns = rec["accession"], rec["namespace"]
            gt[(acc, ns)] = set(GO_PATTERN.findall(rec.get("ground_truth", "") or ""))
            for go_id in set(GO_PATTERN.findall(rec.get("response", "") or "")):
                raw_support[(acc, ns, go_id)].add((model, rec["prompt_id"]))
    return raw_support, gt


def protein_matrix(raw_support, gt, rag, models):
    """Per-protein claim and correct counts by score (0..MAX_SCORE)."""
    n = defaultdict(lambda: np.zeros(MAX_SCORE + 1, dtype=np.int64))
    c = defaultdict(lambda: np.zeros(MAX_SCORE + 1, dtype=np.int64))
    for (acc, ns, go_id), support in raw_support.items():
        by_model = defaultdict(set)
        for m, p in support:
            if m in models:
                by_model[m].add(p)
        if not by_model:
            continue
        rago = int(any(t["go_id"] == go_id for t in rag.get(acc, {}).get(NS_SHORT[ns], [])))
        s = len(by_model) + max(len(v) for v in by_model.values()) + RAGO_WEIGHT * rago
        n[acc][s] += 1
        c[acc][s] += int(go_id in gt.get((acc, ns), set()))
    accs = sorted(n)
    return np.array([n[a] for a in accs]), np.array([c[a] for a in accs])


def precision_at_coverage(n_agg, c_agg, coverage=0.10):
    """Whole-score-group rule (FROZEN_PROTOCOL.md Sec. 4)."""
    target = n_agg.sum() * coverage
    cum_n = cum_c = 0
    for s in range(MAX_SCORE, -1, -1):
        if cum_n >= target:
            break
        cum_n += n_agg[s]
        cum_c += c_agg[s]
    return cum_c / cum_n


def bootstrap_ci(n_mat, c_mat, B, seed=SEED):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(B):
        idx = rng.choice(n_mat.shape[0], size=n_mat.shape[0], replace=True)
        vals.append(precision_at_coverage(n_mat[idx].sum(axis=0), c_mat[idx].sum(axis=0)))
    return [float(x) for x in np.percentile(vals, [2.5, 97.5])]


def summarize(n_mat, c_mat, B):
    n_agg, c_agg = n_mat.sum(axis=0), c_mat.sum(axis=0)
    return {"n_proteins": int(n_mat.shape[0]), "n_claims": int(n_agg.sum()),
            "precision_at_10": float(precision_at_coverage(n_agg, c_agg)),
            "ci_95": bootstrap_ci(n_mat, c_mat, B),
            "baseline": float(c_agg.sum() / n_agg.sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe-results", required=True)
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--out", default=str(RAGO / "primary_analysis" / "development_reference_report.json"))
    args = ap.parse_args()

    raw_support, gt = load_support(args.probe_results)
    rag = json.load(open(RAGO / "rag_candidates.json"))

    report = {"B": args.B, "seed": SEED,
              "all_five": summarize(*protein_matrix(raw_support, gt, rag, set(RULE_MODELS)), args.B),
              "leave_one_model_out": {m: summarize(*protein_matrix(raw_support, gt, rag, set(RULE_MODELS) - {m}),
                                                   args.B) for m in RULE_MODELS}}

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(report, open(out, "w"), indent=2)
    for label, v in [("all five", report["all_five"])] + [(f"minus {k}", v) for k, v in
                                                          report["leave_one_model_out"].items()]:
        print(f"{label:28s} P@10% {v['precision_at_10']:.4f} [{v['ci_95'][0]:.4f}, {v['ci_95'][1]:.4f}]  "
              f"baseline {v['baseline']:.4f}  ({v['n_proteins']} proteins, {v['n_claims']} claims)")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

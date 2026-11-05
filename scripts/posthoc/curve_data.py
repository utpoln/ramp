"""Real-data series for the proposed 'supported claims vs review budget' figure (post hoc). Writes CSV argv[1]."""
import json, sys, csv
from pathlib import Path
from collections import Counter
import numpy as np
from common import *  # frozen loader: claims, y, M, P, R, S, gt, rag_list, helpers
pool = json.load(open(RAGO / "pool_go_lookup.json"))
freq = Counter((s_, t["go_id"]) for d in pool.values() for s_, lst in d.items() for t in lst)
F = np.array([freq[(NS_SHORT[c[1]], c[2])] for c in claims], dtype=float)
rows = [(a, ns, g, rk) for (a, ns), lst in rag_list.items() for rk, g in enumerate(lst)]
hy = np.array([int(g in gt.get((a, ns), set())) for a, ns, g, _ in rows]); hr = -np.array([r for *_, r in rows], dtype=float)
series = {"full": S, "n_models": M, "n_models+rago": M + 5 * R, "n_models+prompt": M + P, "term_frequency": F}
Ks = sorted(set(np.unique(np.round(np.logspace(np.log10(50), np.log10(N), 60)).astype(int))) | {606, 3030, 6059, 7431, 12117, 30293, N})
with open(sys.argv[1], "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["budget_claims", "budget_pct"] + [f"exp_supported_{k}" for k in series] + ["exp_supported_homology_only", "exp_supported_random", "recall_full_of_423", "recall_full_of_933"])
    for K in Ks:
        vals = [expected_at_budget(v, y, K) for v in series.values()]
        h = expected_at_budget(hr, hy, K) if K <= len(rows) else ""
        w.writerow([K, round(100 * K / N, 3)] + [round(v, 2) for v in vals] + [round(h, 2) if h != "" else "", round(K * y.mean(), 2),
                   round(vals[0] / 423, 4), round(vals[0] / 933, 4)])
print("rows", len(Ks), "homology list size", len(rows))

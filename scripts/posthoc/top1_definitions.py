"""Per-protein top-1 under two aggregations (post-hoc clarification of an unspecified tie rule)."""
import json, sys
from pathlib import Path
import numpy as np
from common import *  # frozen loader: claims, y, M, P, R, S, gt, rag_list, helpers
accs = np.unique(prot); idx = {a: np.where(prot == a)[0] for a in accs}
per = {}
for a in accs:
    ii = idx[a]; sel = ii[S[ii] == S[ii].max()]
    per[a] = (len(sel), int(y[sel].sum()))
rng = np.random.default_rng(42)
pooled, mean_ = [], []
for _ in range(2000):
    b = rng.choice(accs, len(accs), replace=True)
    pooled.append(sum(per[a][1] for a in b) / sum(per[a][0] for a in b))
    mean_.append(np.mean([per[a][1] / per[a][0] for a in b]))
out = dict(proteins=len(accs), tied_top_claims=sum(v[0] for v in per.values()), correct=sum(v[1] for v in per.values()),
           pooled=sum(v[1] for v in per.values()) / sum(v[0] for v in per.values()), pooled_ci=list(np.percentile(pooled, [2.5, 97.5])),
           mean_over_proteins=float(np.mean([v[1] / v[0] for v in per.values()])), mean_ci=list(np.percentile(mean_, [2.5, 97.5])),
           random_pick_baseline_mean_over_proteins=float(np.mean([y[idx[a]].mean() for a in accs])))
json.dump(out, open(sys.argv[1], "w"), indent=1); print(json.dumps(out, indent=1))

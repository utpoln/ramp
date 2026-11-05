"""Completes analysis A over the protocol's coverage grid (1, 5, 10, 20, 50%).
10% and 20% were designated as the comparison points before any interval was computed;
the other grid points are reported for completeness, not selected. Same settings as paired_frequency.py."""
import json, sys
from pathlib import Path
from collections import Counter
import numpy as np
from common import *  # frozen loader: claims, y, M, P, R, S, gt, rag_list, helpers
pool = json.load(open(RAGO / "pool_go_lookup.json"))
freq = Counter((s_, t["go_id"]) for d in pool.values() for s_, lst in d.items() for t in lst)
F = np.array([freq[(NS_SHORT[c[1]], c[2])] for c in claims], dtype=float)
accs = np.unique(prot); idx = {a: np.where(prot == a)[0] for a in accs}
rng = np.random.default_rng(42)
boot_ids = [np.concatenate([idx[a] for a in rng.choice(accs, len(accs), replace=True)]) for _ in range(2000)]
out = {}
for f in (0.01, 0.05, 0.10, 0.20, 0.50):
    K = int(np.ceil(f * N)); er, ef = expected_at_budget(S, y, K), expected_at_budget(F, y, K)
    d = [(expected_at_budget(S[ii], y[ii], k) - expected_at_budget(F[ii], y[ii], k)) / k for ii in boot_ids for k in [int(np.ceil(f * len(ii)))]]
    out[f"{f:g}"] = dict(K=K, ramp=er, freq=ef, ramp_prec=er / K, freq_prec=ef / K, diff_pp=100 * (er - ef) / K,
                         diff_ci_pp=list(100 * np.percentile(d, [2.5, 97.5])), extra_claims=er - ef)
json.dump(out, open(sys.argv[1], "w"), indent=1); print(json.dumps(out, indent=1))

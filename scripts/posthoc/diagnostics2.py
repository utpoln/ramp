"""Second batch of exploratory diagnostics (read-only)."""
import json, re, sys
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np

from common import *  # claims, y, M, P, R, S, gt, rag_list, helpers
require_obo()
rep = {}

# term-frequency prior ("Naive" baseline): frequency of the GO term among experimental annotations in the frozen reference pool
pool = json.load(open(RAGO / "pool_go_lookup.json"))
freq = Counter()
for acc, d in pool.items():
    for s_, lst in d.items():
        for t in lst:
            freq[(s_, t["go_id"])] += 1
F = np.array([freq[(NS_SHORT[c[1]], c[2])] for c in claims], dtype=float)
K10 = int(np.ceil(.1 * N))
rep["naive_freq_prior"] = {b: expected_at_budget(F, y, int(np.ceil(f * N))) for b, f in (("1%", .01), ("5%", .05), ("10%", .1), ("20%", .2))}
rep["full_same_budgets"] = {b: expected_at_budget(S, y, int(np.ceil(f * N))) for b, f in (("1%", .01), ("5%", .05), ("10%", .1), ("20%", .2))}
# full score with frequency as a within-tie secondary key is NOT computed (would be a new score)
# combined check: does S add information beyond term frequency? precision of S-top vs F-top within same budget, and S within frequency strata
q = np.quantile(F, [.5, .9])
strata = {"F<=median": F <= q[0], "median<F<=p90": (F > q[0]) & (F <= q[1]), "F>p90": F > q[1]}
rep["S_within_frequency_strata"] = {k: dict(n=int(m.sum()), base=float(y[m].mean()),
                                            wg10=whole_group(S[m], y[m], .1)) for k, m in strata.items()}

# homology-only list at matched budgets
rows = [(a, ns, g, rk) for (a, ns), lst in rag_list.items() for rk, g in enumerate(lst)]
hy = np.array([int(g in gt.get((a, ns), set())) for a, ns, g, _ in rows]); hr = -np.array([r for *_, r in rows], dtype=float)
rep["homology_vs_full_matched"] = {K: dict(homology=expected_at_budget(hr, hy, K), full=expected_at_budget(S, y, K))
                                   for K in (495, 1362, 2042, 3964)}

# proteins with no BLAST hit at all (R==0 for all their claims)
has = {a for (a, ns), l in rag_list.items() if l}
nohit = np.array([a not in has for a in prot])
rep["no_homolog_proteins"] = dict(proteins=len(set(prot[nohit])), claims=int(nohit.sum()), base=float(y[nohit].mean()),
                                  wg10=whole_group(S[nohit], y[nohit], .1))

# valid-ID-only claim set
valid = status == "valid"
rep["valid_only"] = dict(n=int(valid.sum()), base=float(y[valid].mean()), wg10=whole_group(S[valid], y[valid], .1),
                         expected_exact10=expected_at_budget(S[valid], y[valid], int(np.ceil(.1 * valid.sum()))) / int(np.ceil(.1 * valid.sum())))
# alt_id remapping and wrong-namespace cross matches
alt2prim = {}
cur = None
for line in open(OBO):
    if line.startswith("id: GO:"): cur = line[4:14]
    elif line.startswith("alt_id: "): alt2prim[line[8:18]] = cur
alt_match = sum(1 for c, st in zip(claims, status) if st == "alt_id" and alt2prim.get(c[2]) in gt.get((c[0], c[1]), set()))
wrong_ns_match = sum(1 for c, st in zip(claims, status) if st == "wrong_namespace" and c[2] in gt.get((c[0], nsof[c[2]]), set()))
rep["alt_id_claims_matching_gt_after_remap"] = alt_match
rep["wrong_namespace_claims_matching_gt_in_true_namespace"] = wrong_ns_match

# few-shot echo: claims that are few-shot example IDs and produced only under P3
FEW = {"GO:0003700", "GO:0046872", "GO:0042802", "GO:0004714", "GO:0005515", "GO:0016301", "GO:0005524", "GO:0003779", "GO:0005198",
       "GO:0006915", "GO:0006351", "GO:0007050", "GO:0006281", "GO:0045786", "GO:0007049", "GO:0006468", "GO:0045736",
       "GO:0005654", "GO:0005829", "GO:0000785", "GO:0005925", "GO:0015629", "GO:0032587", "GO:0005737", "GO:0005874", "GO:0000779"}
p3only = np.array([all(v == {"P3_fewshot"} for v in cm[c].values()) for c in claims])
fs = np.array([c[2] in FEW for c in claims])
rep["fewshot_echo"] = dict(fewshot_claims=int(fs.sum()), fewshot_P3_only=int((fs & p3only).sum()),
                           P3_share_of_fewshot_supports=float(np.mean([sum(1 for v in cm[c].values() if "P3_fewshot" in v) / len(cm[c]) for c, f_ in zip(claims, fs) if f_])),
                           primary_without_fewshot_ids=dict(n=int((S >= 4)[~fs].sum()), c=int(y[(S >= 4) & ~fs].sum()), prec=float(y[(S >= 4) & ~fs].mean()), base=float(y[~fs].mean())))
# per-prompt and per-model raw precision
pp = defaultdict(lambda: [0, 0]); mm = defaultdict(lambda: [0, 0])
for c, t in zip(claims, y):
    for m_, ps in cm[c].items():
        mm[m_][0] += 1; mm[m_][1] += t
        for p_ in ps: pp[p_][0] += 0
for (a, ns, p_, m_), r in latest.items():
    for g in set(GO.findall(r.get("response") or "")):
        pp[p_][0] += 1; pp[p_][1] += int(g in gt.get((a, ns), set()))
rep["raw_precision_by_prompt"] = {k: [v[0], v[1], v[1] / v[0]] for k, v in sorted(pp.items())}
rep["claims_by_model"] = {k: [v[0], v[1], v[1] / v[0]] for k, v in mm.items()}
# claims per protein-namespace produced by P2 for BP/CC (affected by P2 wording bug) vs MF
p2 = {ns: [0, 0] for ns in NS_SHORT}
for (a, ns, p_, m_), r in latest.items():
    if p_ == "P2_constrained":
        gs = set(GO.findall(r.get("response") or ""))
        p2[ns][0] += 1; p2[ns][1] += len(gs)
rep["P2_mean_ids_per_call_by_ns"] = {k: v[1] / v[0] for k, v in p2.items()}
json.dump(rep, open(sys.argv[1], "w"), indent=1, default=float)
print(json.dumps(rep, indent=1, default=float))

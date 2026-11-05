"""Post-hoc analysis A: RAMP (frozen score) vs annotation-frequency ranking at matched 10% and 20% budgets.

Same claim population (60,585 held-out claims), same budget K = ceil(f * N) for both rankings,
expected-value handling of boundary ties for both. Paired protein-level bootstrap
(B = 2,000, numpy default_rng(42), all 243 proteins), K recomputed per replicate.
Also computes two small items needed for figures/tables:
  S5  GO-identifier validity class by score band
  S4  GPT-5.5 cumulative-precision CIs at S>=12 and S>=14 (same method as analyze_external_model.py)
Writes JSON to argv[1].
"""
import json, sys
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np

from common import *  # claims, y, M, P, R, S, prot, nsa, gt, rag_list, whole_group, expected_at_budget
require_obo()
rep = {}

# annotation-frequency score: number of reference-pool proteins (ADAPT pool, experimental codes) annotated with the term
pool = json.load(open(RAGO / "pool_go_lookup.json"))
freq = Counter((s_, t["go_id"]) for d in pool.values() for s_, lst in d.items() for t in lst)
F = np.array([freq[(NS_SHORT[c[1]], c[2])] for c in claims], dtype=float)

accs = np.unique(prot); idx = {a: np.where(prot == a)[0] for a in accs}
rng = np.random.default_rng(42)
boot_ids = [np.concatenate([idx[a] for a in rng.choice(accs, len(accs), replace=True)]) for _ in range(2000)]

A = {}
for f in (0.10, 0.20):
    K = int(np.ceil(f * N))
    er, ef = expected_at_budget(S, y, K), expected_at_budget(F, y, K)
    d, pr, pf = [], [], []
    for ii in boot_ids:
        k = int(np.ceil(f * len(ii)))
        a_, b_ = expected_at_budget(S[ii], y[ii], k) / k, expected_at_budget(F[ii], y[ii], k) / k
        pr.append(a_); pf.append(b_); d.append(a_ - b_)
    A[f"{int(f*100)}pct"] = dict(
        budget_claims=K, ramp_expected_supported=er, freq_expected_supported=ef,
        ramp_precision=er / K, freq_precision=ef / K,
        ramp_precision_ci=list(np.percentile(pr, [2.5, 97.5])), freq_precision_ci=list(np.percentile(pf, [2.5, 97.5])),
        diff_pp=100 * (er - ef) / K, diff_ci_pp=list(100 * np.percentile(d, [2.5, 97.5])),
        additional_supported_claims_ramp=er - ef,
        share_reps_ramp_higher=float(np.mean(np.array(d) > 0)))
rep["A_ramp_vs_frequency"] = A
rep["A_settings"] = dict(N_claims=N, proteins=len(accs), B=2000, seed=42, tie_rule="expected value at boundary group, both rankings",
                         frequency_source="pool_go_lookup.json (ADAPT 25,460-protein pool, experimental evidence codes)")

# S5: validity class by score band (same OBO parsing as audit diagnostics.py)
terms, nsof, obs = set(), {}, set(); cur = None
for line in open(OBO):
    line = line.rstrip("\n")
    if line == "[Term]": cur = {}; continue
    if line.startswith("[") and line != "[Term]": cur = None; continue
    if cur is None: continue
    if line.startswith("id: GO:"): cur["id"] = line[4:]; terms.add(cur["id"])
    elif line.startswith("namespace: "): nsof[cur["id"]] = line[11:]
    elif line.startswith("is_obsolete: true"): obs.add(cur["id"])
    elif line.startswith("alt_id: "): terms.add(line[8:]); nsof.setdefault(line[8:], None)
status = np.array(["nonexistent" if g not in terms else "obsolete" if g in obs else "alt_id" if nsof.get(g) is None
                   else "wrong_namespace" if nsof[g] != ns else "valid" for (a, ns, g) in claims])
bands = ((2, 2), (3, 4), (5, 8), (9, 14))
rep["S5_validity_by_band"] = {f"{lo}-{hi}": {c: int(((S >= lo) & (S <= hi) & (status == c)).sum()) for c in
                              ("valid", "wrong_namespace", "obsolete", "alt_id", "nonexistent")} for lo, hi in bands}
rep["S5_validity_totals"] = dict(Counter(status.tolist()))

# S4: GPT-5.5 cumulative precision CIs at S>=9/12/14, method of external_validation/analyze_external_model.py
frozen_S = {c: int(s_) for c, s_ in zip(claims, S)}
cl = defaultdict(set)
for line in open(RAGO / "openai_external_validation/production/raw/production_master_output.jsonl"):
    r = json.loads(line); a, ns, p = r["custom_id"].split("|")
    if a not in cohort: continue
    ch = ((r.get("response") or {}).get("body") or {}).get("choices") or [{}]
    for g in set(GO.findall((ch[0].get("message") or {}).get("content") or "")): cl[(a, ns, g)].add(p)
g55 = {}
for thr in (9, 12, 14):
    per = defaultdict(lambda: [0, 0])
    for k in cl:
        if frozen_S.get(k, -1) >= thr:
            per[k[0]][0] += 1; per[k[0]][1] += int(k[2] in gt.get((k[0], k[1]), set()))
    ac = list(per); r2 = np.random.default_rng(42); vals = []
    for _ in range(2000):
        ii = r2.choice(len(ac), size=len(ac), replace=True)
        nn = sum(per[ac[i]][0] for i in ii); cc = sum(per[ac[i]][1] for i in ii)
        if nn: vals.append(cc / nn)
    n = sum(v[0] for v in per.values()); c = sum(v[1] for v in per.values())
    g55[f"S>={thr}"] = dict(n=n, c=c, point=c / n, ci=list(np.percentile(vals, [2.5, 97.5])), n_proteins=len(ac))
rep["S4_gpt55_threshold_cis"] = g55

json.dump(rep, open(sys.argv[1], "w"), indent=1, default=float)
print(json.dumps(rep, indent=1, default=float))

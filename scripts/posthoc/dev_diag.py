import json, re, sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from paths import RAGO, PR, require_probe; require_probe()
GO = re.compile(r"GO:\d{7}"); NS = {"molecular_function":"MF","biological_process":"BP","cellular_component":"CC"}
PROMPTS = {"P1_zeroshot","P2_constrained","P3_fewshot","P4_cot"}
MODELS = ["mistral-large-123b","qwen2.5-72b","google_gemma-4-31b","mixtral-8x7b","mistral-7b"]
sup = defaultdict(lambda: defaultdict(set)); gt = {}
for m in MODELS:
    for l in open(PR/f"{m}.jsonl"):
        r = json.loads(l)
        if r.get("prompt_id") not in PROMPTS or not r.get("success"): continue
        k = (r["accession"], r["namespace"]); gt[k] = set(GO.findall(r.get("ground_truth","") or ""))
        for g in set(GO.findall(r.get("response","") or "")): sup[(k[0],k[1],g)][m].add(r["prompt_id"])
rag = json.load(open(RAGO/"rag_candidates.json"))
cl = sorted(sup); y = np.array([int(c[2] in gt.get((c[0],c[1]),set())) for c in cl])
M = np.array([len(sup[c]) for c in cl]); P = np.array([max(len(v) for v in sup[c].values()) for c in cl])
R = np.array([int(any(t["go_id"]==c[2] for t in rag.get(c[0],{}).get(NS[c[1]],[]))) for c in cl])
N = len(cl)
def wg(s, cov=.1):
    t = cov*len(s)
    for v in np.sort(np.unique(s))[::-1]:
        m = s>=v
        if m.sum()>=t: return dict(thr=int(v), n=int(m.sum()), prec=round(float(y[m].mean()),4), cov=round(float(m.mean()),4))
def exp_at(s, K):
    tk, c = 0, 0.
    for v in np.sort(np.unique(s))[::-1]:
        g = s==v; n = int(g.sum())
        if tk+n <= K: tk += n; c += y[g].sum()
        else: c += (K-tk)*y[g].mean(); break
    return c
K = int(np.ceil(.1*N))
out = {"n_claims": N, "baseline": float(y.mean()), "n_proteins": len({c[0] for c in cl})}
out["weights"] = {w: dict(wg=wg(M+P+w*R), exact10=round(exp_at(M+P+w*R,K)/K,4)) for w in (0,1,2,3,5,8)}
out["ablation_exact10"] = {k: round(exp_at(v,K)/K,4) for k,v in {"full":M+P+5*R,"M":M,"P":P,"R":R,"M+P":M+P,"M+5R":M+5*R,"P+5R":P+5*R}.items()}
out["ablation_wholegroup"] = {k: wg(v) for k,v in {"full":M+P+5*R,"M":M,"M+P":M+P,"M+5R":M+5*R}.items()}
print(json.dumps(out, indent=1))

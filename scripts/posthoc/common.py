import json, re, sys
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
from paths import RAGO, OBO, require_obo  # repository-relative input locations
GO = re.compile(r"GO:\d{7}")
NS_SHORT = {"molecular_function": "MF", "biological_process": "BP", "cellular_component": "CC"}
PROMPTS = {"P1_zeroshot", "P2_constrained", "P3_fewshot", "P4_cot"}
# ---------- inputs (same conventions as scripts/primary_analysis.py) ----------
cohort = {r["accession"] for r in json.load(open(RAGO / "heldout_benchmark_243.json"))}
held = json.load(open(RAGO / "scripts/heldout_sample_250.json"))
held = [r for r in held if r["accession"] in cohort]
func = {r["accession"]: (r.get("function_text") or "") for r in held}
gt = {(r["accession"], ns): {t["go_id"] for t in r["go_terms"].get(s, [])} for r in held for ns, s in NS_SHORT.items()}
rag = json.load(open(RAGO / "heldout_rag_candidates.json"))
rag_list = {(a, ns): [t["go_id"] for t in rag[a].get(s, [])] for a in rag for ns, s in NS_SHORT.items()}

latest = {}
for f in sorted((RAGO / "heldout_results").glob("*.jsonl")):
    for line in open(f):
        r = json.loads(line)
        if r["accession"] in cohort and r["prompt_id"] in PROMPTS and r.get("success"):
            latest[(r["accession"], r["namespace"], r["prompt_id"], r["model"])] = r
models = sorted({k[3] for k in latest})
cm = defaultdict(lambda: defaultdict(set))
for (a, ns, p, m), r in latest.items():
    for g in set(GO.findall(r.get("response") or "")):
        cm[(a, ns, g)][m].add(p)
claims = sorted(cm)
prot = np.array([c[0] for c in claims]); nsa = np.array([c[1] for c in claims])
y = np.array([int(c[2] in gt.get((c[0], c[1]), set())) for c in claims])
M = np.array([len(cm[c]) for c in claims])
P = np.array([max(len(v) for v in cm[c].values()) for c in claims])
R = np.array([int(c[2] in set(rag_list.get((c[0], c[1]), []))) for c in claims])
S = M + P + 5 * R
assert len(claims) == 60585 and y.sum() == 423
N = len(claims)


def whole_group(score, yy, cov):
    t = cov * len(score)
    for s in np.sort(np.unique(score))[::-1]:
        sel = score >= s
        if sel.sum() >= t:
            return dict(thr=float(s), n=int(sel.sum()), c=int(yy[sel].sum()), prec=float(yy[sel].mean()), cov=float(sel.mean()))


def expected_at_budget(score, yy, K):
    """Expected correct among the top-K claims when ties at the boundary are
    broken uniformly at random (no arbitrary ordering)."""
    above = score > 0  # placeholder
    order = np.sort(np.unique(score))[::-1]
    taken, corr = 0, 0.0
    for s in order:
        g = score == s
        ng = int(g.sum())
        if taken + ng <= K:
            taken += ng; corr += yy[g].sum()
        else:
            corr += (K - taken) * yy[g].mean(); taken = K
            break
    return corr



terms, parents, obs, nsof, name = set(), defaultdict(set), set(), {}, {}
cur = None
for line in open(OBO):
    line = line.rstrip("\n")
    if line == "[Term]":
        cur = {}; continue
    if line.startswith("[") and line != "[Term]":
        cur = None; continue
    if cur is None:
        continue
    if line.startswith("id: GO:"):
        cur["id"] = line[4:]; terms.add(cur["id"])
    elif line.startswith("namespace: "):
        nsof[cur["id"]] = line[11:]
    elif line.startswith("name: "):
        name[cur["id"]] = line[6:]
    elif line.startswith("is_obsolete: true"):
        obs.add(cur["id"])
    elif line.startswith("is_a: "):
        parents[cur["id"]].add(line[6:16])
    elif line.startswith("relationship: part_of "):
        parents[cur["id"]].add(line[22:32])
    elif line.startswith("alt_id: "):
        terms.add(line[8:]); nsof.setdefault(line[8:], None)
anc_cache = {}
def ancestors(t):
    if t in anc_cache: return anc_cache[t]
    out, st = set(), [t]
    while st:
        u = st.pop()
        for p in parents.get(u, ()):
            if p not in out:
                out.add(p); st.append(p)
    anc_cache[t] = out; return out
status = []
for (a, ns, g) in claims:
    if g not in terms: status.append("nonexistent")
    elif g in obs: status.append("obsolete")
    elif nsof.get(g) is None: status.append("alt_id")
    elif nsof[g] != ns: status.append("wrong_namespace")
    else: status.append("valid")
status = np.array(status)

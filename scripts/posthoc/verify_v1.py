"""Phase-2 verification analyses for the RAMP revision package (read-only, post hoc).

Reads only released files. Writes one JSON report given as argv[1].
Every quantity here is exploratory unless it is a re-derivation of a frozen value.
"""
import json, re, sys
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np

from common import *  # claims, y, M, P, R, S, gt, rag_list, latest, cm, helpers
from paths import PR, require_probe
require_probe()
rep = {}
B = 2000

accs = np.unique(prot); idx = {a: np.where(prot == a)[0] for a in accs}
rng = np.random.default_rng(42)
boot = [rng.choice(accs, len(accs), replace=True) for _ in range(B)]
boot_ids = [np.concatenate([idx[a] for a in b]) for b in boot]
sel4 = S >= 4

# ---- V1. frozen primary re-derivation (must equal 7,431 / 355) ----
rep["V1_primary"] = dict(n=int(sel4.sum()), c=int(y[sel4].sum()), N=N, base_c=int(y.sum()))

# ---- V2. per-namespace primary precision with protein-bootstrap CIs (same replicates) ----
ns_out = {}
for ns in NS_SHORT:
    m = nsa == ns
    bp = [y[ii][(sel4[ii]) & (nsa[ii] == ns)].mean() for ii in boot_ids]
    bb = [y[ii][nsa[ii] == ns].mean() for ii in boot_ids]
    ns_out[ns] = dict(sel_n=int((sel4 & m).sum()), sel_c=int(y[sel4 & m].sum()), prec=float(y[sel4 & m].mean()),
                      prec_ci=list(np.percentile(bp, [2.5, 97.5])), all_n=int(m.sum()), all_c=int(y[m].sum()),
                      base=float(y[m].mean()), base_ci=list(np.percentile(bb, [2.5, 97.5])))
rep["V2_namespace_primary"] = ns_out

# ---- V3. truncation-proxy sensitivity: drop frozen-five responses with completion_tokens >= 512 ----
def ctok(r):
    u = r.get("usage")
    if isinstance(u, str):
        mm = re.search(r"'completion_tokens':\s*(\d+)", u); return int(mm.group(1)) if mm else None
    return (u or {}).get("completion_tokens")
cap = {k for k, r in latest.items() if (ctok(r) or 0) >= 512}
by_model = Counter(k[3] for k in cap)
cm2 = defaultdict(lambda: defaultdict(set))
for (a, ns, p, m), r in latest.items():
    if (a, ns, p, m) in cap: continue
    for g in set(GO.findall(r.get("response") or "")): cm2[(a, ns, g)][m].add(p)
c2 = sorted(cm2)
y2 = np.array([int(c[2] in gt.get((c[0], c[1]), set())) for c in c2])
S2 = np.array([len(cm2[c]) + max(len(v) for v in cm2[c].values()) + 5 * int(c[2] in set(rag_list.get((c[0], c[1]), []))) for c in c2])
rep["V3_truncation_proxy"] = dict(responses=len(latest), capped=len(cap), capped_by_model=dict(sorted(by_model.items())),
                                  claims_without_capped=len(c2), base=float(y2.mean()), wg10=whole_group(S2, y2, .1))

# ---- V4. precision by (M,P,R) cell for score values 8-11 (non-monotonicity) ----
cells = {}
for s_ in (8, 9, 10, 11, 12, 13, 14):
    m = S == s_
    sub = Counter(zip(M[m], P[m], R[m]))
    cells[s_] = {f"M{a}P{b}R{c}": [int(n), int(y[m & (M == a) & (P == b) & (R == c)].sum())] for (a, b, c), n in sorted(sub.items())}
rep["V4_score_cells"] = cells

# ---- V5. external models: overlap-conditioned paired contrasts, all 243 proteins resampled ----
def load_ext(path, gpt=False):
    out = {}
    for line in open(path):
        r = json.loads(line)
        if gpt:
            a, ns, p = r["custom_id"].split("|")
            body = (r.get("response") or {}).get("body") or {}
            ch = body.get("choices") or [{}]
            out[(a, ns, p)] = (ch[0].get("message") or {}).get("content") or ""
        elif r.get("success") in (True, "True") and r["prompt_id"] in PROMPTS:
            out[(r["accession"], r["namespace"], r["prompt_id"])] = r.get("response") or ""
    return out
ext_files = {"gpt-5.5": (RAGO / "openai_external_validation/production/raw/production_master_output.jsonl", True),
             "llama3.3-70b": (RAGO / "external_validation/llama3_3_70b/llama3.3-70b.jsonl", False),
             "qwen3-32b": (RAGO / "external_validation/qwen3_32b/Qwen_Qwen3-32B.jsonl", False),
             "gemma3-12b": (RAGO / "external_validation/gemma3_12b/gemma3-12b.jsonl", False),
             "gpt-oss-120b": (RAGO / "external_validation/gpt_oss_120b/openai_gpt-oss-120b.jsonl", False)}
frozen_S = {c: int(s_) for c, s_ in zip(claims, S)}
ecm = defaultdict(lambda: defaultdict(set))
ext = {}
for mname, (pth, gpt) in ext_files.items():
    resp = load_ext(pth, gpt)
    cl = defaultdict(set)
    for (a, ns, p), txt in resp.items():
        if a not in cohort: continue
        for g in set(GO.findall(txt)): cl[(a, ns, g)].add(p); ecm[(a, ns, g)][mname].add(p)
    keys = sorted(cl)
    ey = np.array([int(k[2] in gt.get((k[0], k[1]), set())) for k in keys])
    ov = np.array([k in frozen_S for k in keys])
    fs = np.array([frozen_S.get(k, -1) for k in keys])
    ep = np.array([k[0] for k in keys]); eidx = defaultdict(list)
    for i_, a in enumerate(ep): eidx[a].append(i_)
    eidx = {a: np.array(v) for a, v in eidx.items()}
    res = dict(responses=len(resp), claims=len(keys), correct=int(ey.sum()), overlap_n=int(ov.sum()), overlap_c=int(ey[ov].sum()),
               nonoverlap_n=int((~ov).sum()), nonoverlap_c=int(ey[~ov].sum()), proteins_with_claims=len(eidx))
    for thr in (9, 14):
        hi = fs >= thr
        d1, d2 = [], []
        for b in boot:  # same 243-protein replicates as the primary; proteins without claims contribute nothing
            ii = np.concatenate([eidx[a] for a in b if a in eidx])
            o, h = ov[ii], hi[ii]; yy = ey[ii]
            if h.sum() == 0: continue
            d1.append(yy[h].mean() - yy[o].mean())   # vs overlap-subset baseline (the fair comparator)
            d2.append(yy[h].mean() - yy.mean())      # vs the model's overall precision (manuscript comparator)
        res[f"S>={thr}"] = dict(n=int(hi.sum()), c=int(ey[hi].sum()), prec=float(ey[hi].mean()),
                                diff_vs_overlap_base_pp=100 * float(ey[hi].mean() - ey[ov].mean()),
                                diff_vs_overlap_ci_pp=list(100 * np.percentile(d1, [2.5, 97.5])),
                                diff_vs_overall_pp=100 * float(ey[hi].mean() - ey.mean()),
                                diff_vs_overall_ci_pp=list(100 * np.percentile(d2, [2.5, 97.5])), reps_used=len(d1))
    ext[mname] = res
rep["V5_external_overlap_conditioned"] = ext

# ---- V6. external ensemble (frozen formula re-applied): CI construction details ----
ek = sorted(ecm); eyy = np.array([int(k[2] in gt.get((k[0], k[1]), set())) for k in ek])
eM = np.array([len(ecm[k]) for k in ek]); eP = np.array([max(len(v) for v in ecm[k].values()) for k in ek])
eR = np.array([int(k[2] in set(rag_list.get((k[0], k[1]), []))) for k in ek]); eS = eM + eP + 5 * eR
eprot = np.array([k[0] for k in ek]); eidx = {a: np.where(eprot == a)[0] for a in np.unique(eprot)}
pt = whole_group(eS, eyy, .1)
dd, pp, thr_ct, cov_ct = [], [], Counter(), []
for b in boot:
    ii = np.concatenate([eidx[a] for a in b if a in eidx])
    w = whole_group(eS[ii], eyy[ii], .1); thr_ct[w["thr"]] += 1; cov_ct.append(w["cov"])
    dd.append(w["prec"] - eyy[ii].mean()); pp.append(w["prec"])
# fixed-threshold variant (S >= point-estimate threshold held constant across replicates)
dfix = []
for b in boot:
    ii = np.concatenate([eidx[a] for a in b if a in eidx]); h = eS[ii] >= pt["thr"]
    dfix.append(eyy[ii][h].mean() - eyy[ii].mean())
rep["V6_external_ensemble"] = dict(claims=len(ek), correct=int(eyy.sum()), base=float(eyy.mean()), proteins=len(eidx), point=pt,
                                   diff_point_pp=100 * (pt["prec"] - float(eyy.mean())),
                                   diff_ci_pp_rethresholded=list(100 * np.percentile(dd, [2.5, 97.5])),
                                   prec_ci_rethresholded=list(100 * np.percentile(pp, [2.5, 97.5])),
                                   thresholds_chosen_across_reps=dict(thr_ct), coverage_range=[min(cov_ct), max(cov_ct)],
                                   diff_ci_pp_fixed_threshold=list(100 * np.percentile(dfix, [2.5, 97.5])),
                                   overlap_with_frozen_claims=int(sum(k in frozen_S for k in ek)),
                                   note="different seed stream from audit_v3 (same seed 42, replicates drawn over all 243 accessions)")
# external ensemble restricted to the two models never in the PROBE development benchmark (GPT-5.5, GPT-oss-120B)
for label, keep in (("never_in_dev_benchmark", {"gpt-5.5", "gpt-oss-120b"}), ("in_dev_benchmark", {"llama3.3-70b", "qwen3-32b", "gemma3-12b"})):
    sub = {k: {m: v for m, v in d.items() if m in keep} for k, d in ecm.items()}
    sk = sorted(k for k, d in sub.items() if d)
    sy = np.array([int(k[2] in gt.get((k[0], k[1]), set())) for k in sk])
    sS = np.array([len(sub[k]) + max(len(v) for v in sub[k].values()) + 5 * int(k[2] in set(rag_list.get((k[0], k[1]), []))) for k in sk])
    sp = np.array([k[0] for k in sk]); sidx = {a: np.where(sp == a)[0] for a in np.unique(sp)}
    w0 = whole_group(sS, sy, .1); dl = []
    for b in boot:
        ii = np.concatenate([sidx[a] for a in b if a in sidx]); w = whole_group(sS[ii], sy[ii], .1); dl.append(w["prec"] - sy[ii].mean())
    rep[f"V6b_subensemble_{label}"] = dict(models=sorted(keep), claims=len(sk), base=float(sy.mean()), wg10=w0,
                                          diff_point_pp=100 * (w0["prec"] - float(sy.mean())), diff_ci_pp=list(100 * np.percentile(dl, [2.5, 97.5])))

# ---- V7. original development-stage computation (10 PROBE models, P1-P4), as recorded in the PROBE-project session of 2026-10-05 ----
sup = defaultdict(lambda: defaultdict(set)); dgt = {}
for f in sorted(PR.glob("*.jsonl")):
    for line in open(f):
        r = json.loads(line)
        if r.get("prompt_id") not in PROMPTS or not r.get("success"): continue
        k = (r["accession"], r["namespace"]); dgt[k] = set(GO.findall(r.get("ground_truth", "") or ""))
        for g in set(GO.findall(r.get("response", "") or "")): sup[(k[0], k[1], g)][f.stem].add(r["prompt_id"])
drag = json.load(open(RAGO / "rag_candidates.json"))
dk = sorted(sup); dy = np.array([int(k[2] in dgt.get((k[0], k[1]), set())) for k in dk])
dS = np.array([len(sup[k]) + max(len(v) for v in sup[k].values()) + 5 * int(any(t["go_id"] == k[2] for t in drag.get(k[0], {}).get(NS_SHORT[k[1]], []))) for k in dk])
rep["V7_dev_10model"] = dict(models=sorted({m for d in sup.values() for m in d}), claims=len(dk), base=float(dy.mean()),
                             exact10_expected=expected_at_budget(dS, dy, int(len(dk) * .1)) / int(len(dk) * .1),
                             exact1_expected=expected_at_budget(dS, dy, int(len(dk) * .01)) / int(len(dk) * .01),
                             wg10=whole_group(dS, dy, .1))

json.dump(rep, open(sys.argv[1], "w"), indent=1, default=float)
print(json.dumps(rep, indent=1, default=float))

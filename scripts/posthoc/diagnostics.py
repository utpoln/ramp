"""Exploratory audit diagnostics on the frozen RAMP held-out data.

Read-only: reads repository files, writes only to the output path given.
Nothing here is a prespecified analysis unless labelled as such.
"""
import json, re, sys
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np

from paths import RAGO, OBO, require_obo  # repository-relative input locations
require_obo()
OUT = Path(sys.argv[1])
GO = re.compile(r"GO:\d{7}")
NS_SHORT = {"molecular_function": "MF", "biological_process": "BP", "cellular_component": "CC"}
PROMPTS = {"P1_zeroshot", "P2_constrained", "P3_fewshot", "P4_cot"}
rep = {}

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


# ---------- 1. prespecified-but-unreported secondary endpoints ----------
rng = np.random.default_rng(42)
accs = np.unique(prot); idx = {a: np.where(prot == a)[0] for a in accs}
boot_ids = [np.concatenate([idx[a] for a in rng.choice(accs, len(accs), replace=True)]) for _ in range(2000)]
sec = {}
for cov in (0.01, 0.05, 0.10, 0.20, 0.50):
    pt = whole_group(S, y, cov)
    bs = [whole_group(S[ii], y[ii], cov)["prec"] - y[ii].mean() for ii in boot_ids]
    bp = [whole_group(S[ii], y[ii], cov)["prec"] for ii in boot_ids]
    pt["prec_ci"] = list(np.percentile(bp, [2.5, 97.5])); pt["diff_ci"] = list(np.percentile(bs, [2.5, 97.5]))
    sec[f"cov_{cov}"] = pt
rep["secondary_coverage_levels_whole_group"] = sec
# per-protein top-1: all claims tied at each protein's own maximum score (whole-group rule within protein)
top1n = top1c = 0; anyhit = 0
per = {}
for a in accs:
    ii = idx[a]; mx = S[ii].max(); sel = ii[S[ii] == mx]
    per[a] = (len(sel), int(y[sel].sum()))
    top1n += len(sel); top1c += int(y[sel].sum()); anyhit += int(y[sel].sum() > 0)
pb = []
for ii_b in [rng.choice(accs, len(accs), replace=True) for _ in range(2000)]:
    n_ = sum(per[a][0] for a in ii_b); c_ = sum(per[a][1] for a in ii_b); pb.append(c_ / n_)
rep["per_protein_top1"] = dict(claims=top1n, correct=top1c, pooled_prec=top1c / top1n,
                               proteins_with_correct_top=anyhit, prec_ci=list(np.percentile(pb, [2.5, 97.5])),
                               mean_tied_top_claims_per_protein=top1n / len(accs))

# ---------- 2. matched-budget ablations ----------
configs = {"full": S, "M": M, "P": P, "R": R, "M+P": M + P, "M+5R": M + 5 * R, "P+5R": P + 5 * R}
budgets = {"1%": int(np.ceil(.01 * N)), "5%": int(np.ceil(.05 * N)), "10%": int(np.ceil(.10 * N)),
           "12.27%": 7431, "20%": int(np.ceil(.2 * N)), "50%": int(np.ceil(.5 * N))}
mb = {}
for name, sc in configs.items():
    mb[name] = {b: round(expected_at_budget(sc, y, K), 1) for b, K in budgets.items()}
    mb[name]["whole_group_10"] = whole_group(sc, y, .10)
mb["random"] = {b: round(K * y.mean(), 1) for b, K in budgets.items()}
rep["matched_budget_expected_correct"] = mb
# paired bootstrap of matched-budget differences (exact 10% budget, expected-value tie handling)
pairs = {}
for other in ("M", "M+5R", "M+P"):
    d = []
    for ii in boot_ids:
        K = int(np.ceil(.1 * len(ii)))
        d.append((expected_at_budget(S[ii], y[ii], K) - expected_at_budget(configs[other][ii], y[ii], K)) / K)
    K = budgets["10%"]
    pairs[f"full_minus_{other}"] = dict(point_pp=100 * (expected_at_budget(S, y, K) - expected_at_budget(configs[other], y, K)) / K,
                                        ci_pp=list(100 * np.percentile(d, [2.5, 97.5])))
rep["matched_budget_paired_bootstrap_10pct"] = pairs
# area under precision-recall-style curve: average precision over full ranking (ties expected)
def avg_prec_curve(sc):
    pts = []
    for K in np.unique(np.linspace(100, N, 300).astype(int)):
        pts.append(expected_at_budget(sc, y, K) / K)
    return float(np.mean(pts))
rep["mean_precision_over_budgets_100_to_N"] = {k: avg_prec_curve(v) for k, v in configs.items()}

# ---------- 3. weight sensitivity (descriptive; held-out) ----------
ws = {}
for w in (0, 1, 2, 3, 4, 5, 6, 8, 10):
    sc = M + P + w * R
    wg = whole_group(sc, y, .10)
    ws[w] = dict(whole_group=wg, expected_exact10=expected_at_budget(sc, y, budgets["10%"]) / budgets["10%"])
rep["weight_sensitivity_heldout"] = ws
sel4 = S >= 4
rep["primary_set_composition"] = dict(R1_total=int(R.sum()), R1_in_primary=int((sel4 & (R == 1)).sum()),
                                      R1_correct=int(y[R == 1].sum()), R1_prec=float(y[R == 1].mean()),
                                      R0_in_primary=int((sel4 & (R == 0)).sum()), R0_in_primary_correct=int(y[sel4 & (R == 0)].sum()))

# ---------- 4. recall of experimental annotations ----------
gt_total = sum(len(v) for v in gt.values())
gt_ns = {ns: sum(len(gt[(a, ns)]) for a in cohort) for ns in NS_SHORT}
rep["recall"] = dict(gt_terms_total=gt_total, gt_by_ns=gt_ns, recovered_any_claim=int(y.sum()),
                     recovered_primary=int(y[sel4].sum()), recall_any=y.sum() / gt_total, recall_primary=y[sel4].sum() / gt_total,
                     supported_per_1000_reviewed_primary=1000 * y[sel4].mean(), supported_per_1000_random=1000 * y.mean())

# ---------- 5. homology-only baseline (no LLM) ----------
hb = {"all_candidates": None}
rows = [(a, ns, g, rank) for (a, ns), lst in rag_list.items() if a in cohort for rank, g in enumerate(lst)]
hy = np.array([int(g in gt.get((a, ns), set())) for a, ns, g, _ in rows])
hrank = np.array([r for *_, r in rows])
hb["all_candidates"] = dict(n=len(rows), correct=int(hy.sum()), prec=float(hy.mean()) if len(rows) else None)
for k in (1, 3, 5):
    s = hrank < k
    hb[f"top{k}_per_protein_ns"] = dict(n=int(s.sum()), correct=int(hy[s].sum()), prec=float(hy[s].mean()))
llm_keys = set(claims)
hb["candidates_also_LLM_claims"] = sum(1 for a, ns, g, _ in rows if (a, ns, g) in llm_keys)
hb["correct_candidates_not_claimed_by_any_LLM"] = int(sum(1 for (a, ns, g, _), t in zip(rows, hy) if t and (a, ns, g) not in llm_keys))
hb["union_correct_LLM_or_homology"] = len({(a, ns, g) for (a, ns, g, _), t in zip(rows, hy) if t} | {c for c, t in zip(claims, y) if t})
hb["proteins_with_any_candidate"] = len({a for (a, ns), l in rag_list.items() if a in cohort and l})
rep["homology_only_baseline"] = hb
# BLAST top-hit identity
best = {}
for line in open(RAGO / "heldout_blast_hits.tsv"):
    q, s_, pid, ln, ev, bit = line.rstrip("\n").split("\t")
    if q not in best or float(bit) > best[q][1]:
        best[q] = (float(pid), float(bit))
ids = np.array([best[a][0] for a in cohort if a in best])
rep["blast_top_hit_identity"] = dict(queries_with_hit=len(ids), queries_without_hit=len(cohort) - len(ids),
                                     ge90=int((ids >= 90).sum()), ge70=int((ids >= 70).sum()), ge50=int((ids >= 50).sum()),
                                     lt30=int((ids < 30).sum()), median=float(np.median(ids)))
# precision of top set stratified by whether protein has a >=50% identity homolog
hi = np.array([best.get(a, (0, 0))[0] >= 50 for a in prot])
rep["primary_by_homolog_identity"] = {k: dict(n=int((sel4 & m_).sum()), c=int(y[sel4 & m_].sum()), prec=float(y[sel4 & m_].mean()),
                                              base=float(y[m_].mean()), proteins=len(set(prot[m_])))
                                      for k, m_ in (("top_hit_ge50pct", hi), ("top_hit_lt50pct_or_none", ~hi))}

# ---------- 6. GO validity and hierarchy (OBO hash-matching file) ----------
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
rep["claim_validity"] = {s_: dict(n=int((status == s_).sum()), pct=float((status == s_).mean() * 100),
                                  correct=int(y[status == s_].sum()),
                                  n_in_primary=int((sel4 & (status == s_)).sum())) for s_ in np.unique(status)}
for k in ("valid",):
    v = status == k
    rep["validity_by_score_band"] = {f"S={lo}-{hi_}": float(((status != "valid") & (S >= lo) & (S <= hi_)).sum() / max(1, ((S >= lo) & (S <= hi_)).sum()))
                                     for lo, hi_ in ((2, 2), (3, 4), (5, 8), (9, 14))}
# hierarchy-aware: claim counted supported if it equals OR is an ancestor of any experimental term (true-path rule)
gt_closure = {}
for key, s_ in gt.items():
    clo = set(s_)
    for t in s_:
        clo |= ancestors(t)
    gt_closure[key] = clo
yh = np.array([int(c[2] in gt_closure.get((c[0], c[1]), set())) for c in claims])
roots = {"GO:0003674", "GO:0008150", "GO:0005575"}
yh_nr = np.array([int(t and c[2] not in roots) for c, t in zip(claims, yh)])
rep["hierarchy_aware"] = dict(
    baseline_exact=float(y.mean()), baseline_propagated=float(yh.mean()), baseline_propagated_noroot=float(yh_nr.mean()),
    primary_exact=float(y[sel4].mean()), primary_propagated=float(yh[sel4].mean()), primary_propagated_noroot=float(yh_nr[sel4].mean()),
    by_score_propagated={int(s_): float(yh[S == s_].mean()) for s_ in np.unique(S)},
    whole_group_10_propagated=whole_group(S, yh_nr, .10))
pdiff = [yh_nr[ii][S[ii] >= whole_group(S[ii], y[ii], .1)["thr"]].mean() - yh_nr[ii].mean() for ii in boot_ids[:1000]]
rep["hierarchy_aware"]["paired_diff_ci_noroot_B1000"] = list(np.percentile(pdiff, [2.5, 97.5]))

# ---------- 7. generic terms, few-shot echo, namespace ----------
FEWSHOT = {"molecular_function": {"GO:0003700", "GO:0046872", "GO:0042802", "GO:0004714", "GO:0005515", "GO:0016301", "GO:0005524", "GO:0003779", "GO:0005198"},
           "biological_process": {"GO:0006915", "GO:0006351", "GO:0007050", "GO:0006281", "GO:0045786", "GO:0007049", "GO:0006468", "GO:0045736"},
           "cellular_component": {"GO:0005654", "GO:0005829", "GO:0000785", "GO:0005925", "GO:0015629", "GO:0032587", "GO:0005737", "GO:0005874", "GO:0000779"}}
fs = np.array([c[2] in FEWSHOT[c[1]] for c in claims])
rep["fewshot_example_ids"] = dict(n=int(fs.sum()), correct=int(y[fs].sum()), prec=float(y[fs].mean()),
                                  in_primary=int((sel4 & fs).sum()), correct_in_primary=int(y[sel4 & fs].sum()),
                                  share_of_primary_correct=float(y[sel4 & fs].sum() / y[sel4].sum()))
cnt = Counter(c[2] for c, t, s_ in zip(claims, y, sel4) if t and s_)
rep["top_correct_terms_in_primary"] = [(g, name.get(g, "?"), n_) for g, n_ in cnt.most_common(15)]
top10 = {g for g, _ in cnt.most_common(10)}
rep["share_primary_correct_from_10_terms"] = sum(n_ for g, n_ in cnt.most_common(10)) / y[sel4].sum()
excl = np.array([c[2] not in top10 for c in claims])
rep["primary_excluding_10_most_frequent_correct_terms"] = dict(n=int((sel4 & excl).sum()), c=int(y[sel4 & excl].sum()),
                                                               prec=float(y[sel4 & excl].mean()), base=float(y[excl].mean()))
rep["ns_share_of_primary_correct"] = {ns: int(y[sel4 & (nsa == ns)].sum()) for ns in NS_SHORT}

# ---------- 8. FUNCTION-text lexical overlap ----------
def norm(s): return re.sub(r"[^a-z0-9 ]", " ", s.lower())
lex = np.array([bool(name.get(c[2])) and len(name[c[2]]) >= 6 and f" {norm(name[c[2]]).strip()} " in f" {norm(func[c[0]])} " for c in claims])
rep["function_text_lexical_overlap"] = dict(
    claims_with_name_in_function_text=int(lex.sum()), pct=float(lex.mean() * 100),
    prec_lex=float(y[lex].mean()), prec_nonlex=float(y[~lex].mean()),
    correct_lex=int(y[lex].sum()), correct_total=int(y.sum()),
    primary_lex=int((sel4 & lex).sum()), primary_correct_lex=int(y[sel4 & lex].sum()),
    primary_prec_excluding_lex=float(y[sel4 & ~lex].mean()), base_excluding_lex=float(y[~lex].mean()),
    function_texts_with_GO_ids=sum(1 for t in func.values() if GO.search(t)))
gtl = [(a, ns, g) for (a, ns), s_ in gt.items() for g in s_]
rep["gt_terms_name_in_function_text"] = dict(n=sum(1 for a, ns, g in gtl if name.get(g) and len(name[g]) >= 6 and f" {norm(name[g]).strip()} " in f" {norm(func[a])} "), total=len(gtl))

# ---------- 9. external models: overlap selection + re-application of frozen formula to a new ensemble ----------
def load_ext(path, gpt=False):
    out = {}
    for line in open(path):
        r = json.loads(line)
        if gpt:
            a, ns, p = r["custom_id"].split("|")
            body = (r.get("response") or {}).get("body") or {}
            ch = body.get("choices") or [{}]
            txt = (ch[0].get("message") or {}).get("content") or ""
            out[(a, ns, p)] = txt
        else:
            if r.get("success") and r["prompt_id"] in PROMPTS:
                out[(r["accession"], r["namespace"], r["prompt_id"])] = r.get("response") or ""
    return out
ext_files = {"gpt-5.5": (RAGO / "openai_external_validation/production/raw/production_master_output.jsonl", True),
             "llama3.3-70b": (RAGO / "external_validation/llama3_3_70b/llama3.3-70b.jsonl", False),
             "qwen3-32b": (RAGO / "external_validation/qwen3_32b/Qwen_Qwen3-32B.jsonl", False),
             "gemma3-12b": (RAGO / "external_validation/gemma3_12b/gemma3-12b.jsonl", False),
             "gpt-oss-120b": (RAGO / "external_validation/gpt_oss_120b/openai_gpt-oss-120b.jsonl", False)}
ecm = defaultdict(lambda: defaultdict(set))
ext = {}
frozen_S = {c: s_ for c, s_ in zip(claims, S)}
for m, (pth, gpt) in ext_files.items():
    resp = load_ext(pth, gpt)
    cl = defaultdict(set)
    for (a, ns, p), txt in resp.items():
        if a not in cohort: continue
        for g in set(GO.findall(txt)):
            cl[(a, ns, g)].add(p); ecm[(a, ns, g)][m].add(p)
    keys = list(cl); ey = np.array([int(k[2] in gt.get((k[0], k[1]), set())) for k in keys])
    ov = np.array([k in frozen_S for k in keys])
    own_P = np.array([len(cl[k]) for k in keys]); own_R = np.array([int(k[2] in set(rag_list.get((k[0], k[1]), []))) for k in keys])
    own = own_P + 5 * own_R  # single-model analogue (exploratory; not the frozen score)
    ext[m] = dict(n=len(keys), c=int(ey.sum()), prec=float(ey.mean()), overlap_n=int(ov.sum()),
                  overlap_prec=float(ey[ov].mean()), nonoverlap_prec=float(ey[~ov].mean()),
                  share_correct_in_overlap=float(ey[ov].sum() / max(1, ey.sum())),
                  own_P_plus_5R_whole_group_10=whole_group(own, ey, .10))
rep["external_overlap"] = ext
# frozen formula re-applied to the five external models as a NEW ensemble (exploratory, not the frozen score)
ek = sorted(ecm); eyy = np.array([int(k[2] in gt.get((k[0], k[1]), set())) for k in ek])
eM = np.array([len(ecm[k]) for k in ek]); eP = np.array([max(len(v) for v in ecm[k].values()) for k in ek])
eR = np.array([int(k[2] in set(rag_list.get((k[0], k[1]), []))) for k in ek]); eS = eM + eP + 5 * eR
eprot = np.array([k[0] for k in ek]); eidx = {a: np.where(eprot == a)[0] for a in np.unique(eprot)}
rng2 = np.random.default_rng(42); ea = np.unique(eprot); dd = []
for _ in range(2000):
    ii = np.concatenate([eidx[a] for a in rng2.choice(ea, len(ea), replace=True)])
    dd.append(whole_group(eS[ii], eyy[ii], .1)["prec"] - eyy[ii].mean())
rep["external_ensemble_reapplied_formula"] = dict(n_claims=len(ek), baseline=float(eyy.mean()), wg10=whole_group(eS, eyy, .1),
                                                  diff_ci=list(np.percentile(dd, [2.5, 97.5])),
                                                  by_score={int(s_): [int((eS == s_).sum()), float(eyy[eS == s_].mean())] for s_ in np.unique(eS)})

json.dump(rep, open(OUT, "w"), indent=1, default=float)
print(json.dumps(rep, indent=1, default=float))

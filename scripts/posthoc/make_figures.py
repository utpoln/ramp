"""New v4 figures, drawn only from released data and saved analysis outputs.

fig_baselines.pdf          supported claims vs review budget, precision vs coverage (post hoc comparators)
fig_namespace_validity.pdf namespace precision vs baseline; GO-identifier validity by score band (post hoc)
fig_external.pdf           frozen-score ranking within external models' overlapping claims; panel re-applied as an ensemble

Every plotted number is recomputed or read from a JSON/CSV output and checked with an assert
against the value reported in the manuscript. Run: python3 scripts/posthoc/make_figures.py (after run_all.sh has produced output/)
"""
import csv, json, re
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8.5, "axes.spines.top": False,
                            "axes.spines.right": False, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
                            "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.grid": True,
                            "grid.color": "#e6e5e1", "grid.linewidth": 0.6, "axes.axisbelow": True})
import os
from paths import RAGO
AN = AUD = Path(os.environ.get("RAMP_POSTHOC_OUT", Path(__file__).parent / "output"))
HERE = AN / "figures"; HERE.mkdir(parents=True, exist_ok=True)
C1, C2, C3, C4, C5, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#8a8985"
TXT2 = "#52514e"
N, SUPPORTED = 60585, 423

def close(a, b, tol=0.06):
    assert abs(a - b) <= tol, (a, b)

# ------------------------------------------------------------------ Figure A: baselines
rows = list(csv.DictReader(open(AN / "curve_data.csv")))
K = np.array([int(r["budget_claims"]) for r in rows])
def col(name):
    return np.array([float(r[name]) if r[name] != "" else np.nan for r in rows])
series = [  # label, column, color, style, marker
    ("Frozen RAMP score", "exp_supported_full", C1, "-", "o"),
    ("n_models + rago_support", "exp_supported_n_models+rago", C3, (0, (5, 2)), "s"),
    ("Annotation-frequency ranking", "exp_supported_term_frequency", C2, (0, (1, 1.5)), "^"),
    ("n_models alone", "exp_supported_n_models", C4, (0, (6, 2, 1, 2)), "D"),
    ("Homology candidates only (no LLM)", "exp_supported_homology_only", C5, "-", "v"),
    ("Random order", "exp_supported_random", GREY, (0, (3, 3)), None),
]
# checks against file 02 / manuscript values
i10, i20, i1 = list(K).index(6059), list(K).index(12117), list(K).index(606)
close(col("exp_supported_full")[i10], 345.7); close(col("exp_supported_term_frequency")[i10], 329.0)
close(col("exp_supported_full")[i20], 370.3); close(col("exp_supported_term_frequency")[i20], 375.5)
close(col("exp_supported_n_models")[i10], 314.4); close(col("exp_supported_n_models+rago")[i10], 345.2)
close(col("exp_supported_full")[i1], 137.8)
hom = col("exp_supported_homology_only"); assert np.nanmax(hom) <= 304.01
cum_table = [(14, 128, 49), (13, 239, 74), (12, 350, 101), (11, 435, 114), (10, 540, 126), (9, 898, 190), (8, 1355, 232),
             (7, 2097, 285), (6, 2831, 312), (5, 4325, 334), (4, 7431, 355), (3, 15073, 380), (2, 60585, 423)]
# verify the frozen threshold points against the audit loader's whole-group outputs
sec = json.load(open(AUD / "diagnostics_out.json"))["secondary_coverage_levels_whole_group"]
for key, thr in (("cov_0.01", 9), ("cov_0.05", 5), ("cov_0.1", 4), ("cov_0.2", 3)):
    n_, c_ = [(n, c) for s, n, c in cum_table if s == thr][0]
    assert sec[key]["n"] == n_ and sec[key]["c"] == c_, key

fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3))
for lab, cname, col_, ls, mk in series:
    v = col(cname); m = ~np.isnan(v)
    ax[0].plot(K[m], v[m], color=col_, ls=ls, lw=2.2 if cname == "exp_supported_full" else 1.5,
               marker=mk, markevery=8, ms=4.5, label=lab, zorder=3 if cname == "exp_supported_full" else 2)
    ax[1].plot(100 * K[m] / N, 100 * v[m] / K[m], color=col_, ls=ls, lw=2.2 if cname == "exp_supported_full" else 1.5,
               marker=mk, markevery=8, ms=4.5, zorder=3 if cname == "exp_supported_full" else 2)
ax[0].plot(3964, 304, marker="v", color=C5, ms=6, zorder=4)
# homology list end point (3,964 candidates, 173 proteins) is described in the caption
ax[0].axhline(SUPPORTED, color=TXT2, lw=0.8, ls=":")
ax[0].text(45, SUPPORTED + 8, "423 supported claims in pool", fontsize=6.6, color=TXT2)
for a in ax:
    a.set_xscale("log")
for x, lab in ((6059, "10% budget"), (12117, "20%")):
    ax[0].axvline(x, color=TXT2, lw=0.7, ls="--"); ax[0].text(x * 1.05, 15, lab, fontsize=6.8, color=TXT2, rotation=90)
ax[0].set_xlabel("Review budget (claims reviewed, log scale)")
ax[0].set_ylabel("Expected supported claims found")
ax[0].set_ylim(0, 450)
ax[0].set_title("A  Supported claims recovered", loc="left", fontsize=9)
cov = [100 * n / N for _, n, _ in cum_table]; prec = [100 * c / n for _, n, c in cum_table]
ax[1].scatter(cov, prec, s=22, facecolor="white", edgecolor=C1, lw=1.3, zorder=5, label="Frozen whole-group thresholds")
p4 = [i for i, t in enumerate(cum_table) if t[0] == 4][0]
ax[1].scatter([cov[p4]], [prec[p4]], s=46, color=C1, zorder=6)
ax[1].annotate("prespecified operating point\n(score ≥ 4: 12.27%, 4.78%)", (cov[p4], prec[p4]), xytext=(-118, -30),
               textcoords="offset points", fontsize=6.6, color=TXT2, arrowprops=dict(arrowstyle="-", color=TXT2, lw=0.6))
ax[1].set_yscale("log")
ax[1].set_xlabel("Coverage (% of 60,585 claims, log scale)")
ax[1].set_ylabel("Precision (%, log scale)")
ax[1].set_title("B  Precision versus coverage", loc="left", fontsize=9)
h, l = ax[0].get_legend_handles_labels(); h2, l2 = ax[1].get_legend_handles_labels()
fig.legend(h + h2, l + l2, loc="lower center", ncol=3, fontsize=6.8, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(rect=(0, 0.14, 1, 1))
fig.savefig(HERE / "fig_baselines.pdf"); fig.savefig(HERE / "fig_baselines.png", dpi=220); plt.close(fig)

# ------------------------------------------------------------------ Figure B: namespace and validity
ver = json.load(open(AN / "verify_v1_out.json"))
pf = json.load(open(AN / "paired_frequency_out.json"))
ns = ver["V2_namespace_primary"]
order = [("molecular_function", "MF"), ("biological_process", "BP"), ("cellular_component", "CC")]
close(100 * ns["cellular_component"]["prec"], 12.02, 0.006); close(100 * ns["molecular_function"]["prec"], 1.81, 0.006)
close(100 * ns["biological_process"]["prec"], 1.52, 0.006)
band = pf["S5_validity_by_band"]
assert sum(sum(v.values()) for v in band.values()) == N
assert sum(v["valid"] for v in band.values()) == 36567
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
x = np.arange(3); w = 0.36
for j, (key, lab, colr) in enumerate((("prec", "Primary set (score ≥ 4)", C1), ("base", "All claims (baseline)", GREY))):
    vals = [100 * ns[k][key] for k, _ in order]
    lo = [100 * ns[k][key + "_ci"][0] for k, _ in order]; hi = [100 * ns[k][key + "_ci"][1] for k, _ in order]
    ax[0].bar(x + (j - 0.5) * w, vals, w * 0.92, color=colr, label=lab)
    ax[0].errorbar(x + (j - 0.5) * w, vals, yerr=[np.array(vals) - lo, np.array(hi) - vals], fmt="none", ecolor="#0b0b0b", lw=0.8, capsize=2)
ax[0].set_xticks(x, [f"{l}\n{ns[k]['sel_c']}/{ns[k]['sel_n']:,} vs\n{ns[k]['all_c']}/{ns[k]['all_n']:,}" for k, l in order], fontsize=7)
ax[0].set_ylabel("Precision (%)")
ax[0].legend(fontsize=6.8, frameon=False, loc="upper left")
ax[0].set_title("A  Precision by GO namespace", loc="left", fontsize=9)
cls = [("valid", "Valid in requested namespace", C1), ("wrong_namespace", "Valid ID, other namespace", C2),
       ("obsolete", "Obsolete", C4), ("alt_id", "Alternative (merged) ID", C5), ("nonexistent", "Non-existent", GREY)]
blabels = ["2", "3–4", "5–8", "9–14"]
left = np.zeros(4)
for k, lab, colr in cls:
    tot = np.array([sum(band[b].values()) for b in ("2-2", "3-4", "5-8", "9-14")])
    v = 100 * np.array([band[b][k] for b in ("2-2", "3-4", "5-8", "9-14")]) / tot
    ax[1].barh(blabels, v, left=left, color=colr, label=lab, edgecolor="white", linewidth=1)
    left += v
for i, b in enumerate(("2-2", "3-4", "5-8", "9-14")):
    t = sum(band[b].values()); ax[1].text(101, i, f"n={t:,}", va="center", fontsize=6.6, color=TXT2)
ax[1].set_xlim(0, 100); ax[1].set_xlabel("Share of claims (%)"); ax[1].set_ylabel("Frozen score band")
ax[1].grid(axis="y", visible=False)
ax[1].legend(fontsize=6.4, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=2)
ax[1].set_title("B  GO-identifier validity by score", loc="left", fontsize=9)
fig.tight_layout()
fig.savefig(HERE / "fig_namespace_validity.pdf", bbox_inches="tight"); fig.savefig(HERE / "fig_namespace_validity.png", dpi=220, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ Figure C: external models
from common import claims, S, gt, cohort, GO, PROMPTS  # frozen loader
frozen_S = {c: int(s_) for c, s_ in zip(claims, S)}
def load(path, gpt):
    out = {}
    for line in open(path):
        r = json.loads(line)
        if gpt:
            a, n_, p = r["custom_id"].split("|")
            ch = ((r.get("response") or {}).get("body") or {}).get("choices") or [{}]
            out[(a, n_, p)] = (ch[0].get("message") or {}).get("content") or ""
        elif r.get("success") in (True, "True") and r["prompt_id"] in PROMPTS:
            out[(r["accession"], r["namespace"], r["prompt_id"])] = r.get("response") or ""
    return out
ext = [("GPT-5.5", RAGO / "openai_external_validation/production/raw/production_master_output.jsonl", True),
       ("Llama 3.3 70B", RAGO / "external_validation/llama3_3_70b/llama3.3-70b.jsonl", False),
       ("Qwen3 32B", RAGO / "external_validation/qwen3_32b/Qwen_Qwen3-32B.jsonl", False),
       ("Gemma3 12B", RAGO / "external_validation/gemma3_12b/gemma3-12b.jsonl", False),
       ("GPT-oss-120B", RAGO / "external_validation/gpt_oss_120b/openai_gpt-oss-120b.jsonl", False)]
reported_ge9 = {"GPT-5.5": 28.46, "Llama 3.3 70B": 25.09, "Qwen3 32B": 23.85, "Gemma3 12B": 24.50, "GPT-oss-120B": 25.36}
reported_ge14 = {"GPT-5.5": 42.48, "Llama 3.3 70B": 37.17, "Qwen3 32B": 38.21, "Gemma3 12B": 36.84, "GPT-oss-120B": 40.00}
figC = {}
for name, pth, gpt in ext:
    keys = set()
    for (a, n_, p), txt in load(pth, gpt).items():
        if a in cohort:
            keys |= {(a, n_, g) for g in GO.findall(txt)}
    yk = {k: int(k[2] in gt.get((k[0], k[1]), set())) for k in keys}
    ov = [k for k in keys if k in frozen_S]; non = [k for k in keys if k not in frozen_S]
    pts = {t: 100 * np.mean([yk[k] for k in ov if frozen_S[k] >= t]) for t in (2, 9, 12, 14)}
    pts["non"] = 100 * np.mean([yk[k] for k in non])
    close(pts[9], reported_ge9[name], 0.006); close(pts[14], reported_ge14[name], 0.006)
    figC[name] = pts
close(figC["GPT-5.5"][2], 14.85, 0.006); close(figC["Gemma3 12B"]["non"], 0.06, 0.006)
ens = json.load(open(AUD / "diagnostics_out.json"))["external_ensemble_reapplied_formula"]
assert ens["wg10"]["n"] == 4878 and ens["wg10"]["c"] == 531
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.7), gridspec_kw={"width_ratios": [1.25, 1]})
xs = np.arange(5)
marks = [(2, "All overlapping claims", GREY, "o"), (9, "Score ≥ 9", C2, "s"), (12, "Score ≥ 12", C3, "^"), (14, "Score ≥ 14", C1, "D")]
for t, lab, colr, mk in marks:
    ax[0].plot(xs, [figC[n][t] for n, *_ in ext], ls="none", marker=mk, ms=6, color=colr, label=lab, zorder=3)
ax[0].plot(xs, [figC[n]["non"] for n, *_ in ext], ls="none", marker="x", ms=6, mew=1.4, color="#0b0b0b",
           label="Non-overlapping claims (cannot be scored)", zorder=3)
for i, (n, *_) in enumerate(ext):
    ax[0].vlines(i, figC[n]["non"], figC[n][14], color="#d6d5d0", lw=1, zorder=1)
ax[0].set_xticks(xs, [n for n, *_ in ext], rotation=20, ha="right")
ax[0].set_ylabel("Cumulative precision (%)")
ax[0].legend(fontsize=6.6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=3)
ax[0].set_ylim(-1, 50)
ax[0].set_title("A  Frozen score within overlapping claims", loc="left", fontsize=9)
bs = ens["by_score"]; sc = sorted(int(k) for k in bs)
ax[1].bar(sc, [100 * bs[str(s)][1] for s in sc], 0.8, color=[C1 if s >= 5 else "#b9c9e0" for s in sc])
ax[1].axhline(100 * ens["baseline"], color=TXT2, ls="--", lw=0.9)
ax[1].text(1.6, 100 * ens["baseline"] + 1.0, f"baseline {100*ens['baseline']:.2f}%", fontsize=6.6, color=TXT2, ha="left")
ax[1].text(1.6, 27, "whole-group 10% rule\nselects score ≥ 5:\n531/4,878 = 10.89%", fontsize=6.6, color=TXT2)
ax[1].set_xlabel("Score from the five panel models")
ax[1].set_ylabel("Precision at score (%)")
ax[1].set_xticks(sc, [str(s) for s in sc], fontsize=6.8)
ax[1].set_title("B  Formula re-applied (exploratory)", loc="left", fontsize=9)
fig.tight_layout()
fig.savefig(HERE / "fig_external.pdf", bbox_inches="tight"); fig.savefig(HERE / "fig_external.png", dpi=220, bbox_inches="tight"); plt.close(fig)
json.dump({k: {str(a): round(b, 4) for a, b in v.items()} for k, v in figC.items()}, open(HERE / "fig_external_values.json", "w"), indent=1)
print("All figure checks passed; figures written to", HERE)

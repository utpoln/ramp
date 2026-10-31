"""
cohort_comparison.py

Held-out cohort characteristics versus the development benchmark
(manuscript Table "Held-out cohort characteristics"):

  - sequence length (mean, Cohen's d, Mann-Whitney U)
  - total experimental GO terms per protein (mean, Cohen's d), with and
    without E. coli
  - per-namespace presence (% of proteins with >= 1 experimental term) and
    Cohen's d on per-namespace term counts
  - share of proteins at exactly the three-term eligibility floor

Cohen's d uses the pooled standard deviation (held-out minus development).

Inputs (all already in this repository):
  heldout_benchmark_243.json      the 243-protein cohort
  scripts/heldout_sample_250.json held-out lengths and experimental GO terms
  benchmark_master.json           development benchmark (1,000 proteins)
  probe_1000_lengths.json         development sequence lengths

Usage:
  python scripts/cohort_comparison.py
  python scripts/cohort_comparison.py --out primary_analysis/cohort_comparison_report.json
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
RAGO = HERE.parent
NS_SHORT = {"molecular_function": "MF", "biological_process": "BP", "cellular_component": "CC"}


def cohens_d(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / pooled)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RAGO / "primary_analysis" / "cohort_comparison_report.json"))
    args = ap.parse_args()

    cohort = {r["accession"] for r in json.load(open(RAGO / "heldout_benchmark_243.json"))}
    held = [r for r in json.load(open(RAGO / "scripts" / "heldout_sample_250.json")) if r["accession"] in cohort]
    dev = json.load(open(RAGO / "benchmark_master.json"))
    dev_len = json.load(open(RAGO / "probe_1000_lengths.json"))

    h_counts = [{"MF": int(r["n_exp_go_mf"]), "BP": int(r["n_exp_go_bp"]), "CC": int(r["n_exp_go_cc"])}
                for r in held]
    d_counts = []
    for r in dev:
        c = {"MF": 0, "BP": 0, "CC": 0}
        for t in json.loads(r["go_terms_json"]):
            c[NS_SHORT[t["namespace"]]] += 1
        d_counts.append(c)
    h_tot = [sum(c.values()) for c in h_counts]
    d_tot = [sum(c.values()) for c in d_counts]

    h_l = [int(r["sequence_length"]) for r in held]
    d_l = [dev_len[r["accession"]] for r in dev]

    h_noeco = [t for t, r in zip(h_tot, held) if r["organism"] != "ecoli"]
    d_noeco = [t for t, r in zip(d_tot, dev) if r["organism"] != "ecoli"]

    report = {
        "n_proteins": {"heldout": len(held), "development": len(dev)},
        "sequence_length": {"heldout_mean": float(np.mean(h_l)), "development_mean": float(np.mean(d_l)),
                            "cohens_d": cohens_d(h_l, d_l),
                            "mann_whitney_p": float(mannwhitneyu(h_l, d_l).pvalue)},
        "total_exp_go_terms": {"heldout_mean": float(np.mean(h_tot)), "development_mean": float(np.mean(d_tot)),
                               "cohens_d": cohens_d(h_tot, d_tot),
                               "cohens_d_excluding_ecoli": cohens_d(h_noeco, d_noeco)},
        "namespaces": {ns: {"heldout_presence_pct": 100 * float(np.mean([c[ns] > 0 for c in h_counts])),
                            "development_presence_pct": 100 * float(np.mean([c[ns] > 0 for c in d_counts])),
                            "cohens_d_term_counts": cohens_d([c[ns] for c in h_counts], [c[ns] for c in d_counts])}
                       for ns in ("MF", "BP", "CC")},
        "pct_at_three_term_floor": {"heldout": 100 * float(np.mean([t == 3 for t in h_tot])),
                                    "development": 100 * float(np.mean([t == 3 for t in d_tot]))},
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(report, open(out, "w"), indent=2)
    print(json.dumps(report, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

# Post-hoc exploratory analyses

**Status: post hoc.** Nothing in this directory was specified in `FROZEN_PROTOCOL.md`. These analyses were run after the held-out results had been examined. They use only released files, make no LLM calls, and change neither the frozen score nor the primary endpoint. They are exploratory and not corrected for multiplicity.

The prespecified primary analysis is `scripts/primary_analysis.py`.

## Requirements

**Software:** Python 3.9 or later, `numpy`, and `matplotlib` (figures only). Install with `pip install -r scripts/posthoc/requirements.txt`. No BLAST or network access is needed.

**Input files.** Everything is in this repository except the PROBE development responses.

| Input | Location | Used by |
|---|---|---|
| Held-out cohort and ground truth | `heldout_benchmark_243.json`, `scripts/heldout_sample_250.json` | all |
| Frozen five-model responses | `heldout_results/*.jsonl` | all |
| Frozen homology candidates | `heldout_rag_candidates.json`, `heldout_blast_hits.tsv` | all |
| Homology-reference GO lookup (term-frequency ranking) | `pool_go_lookup.json` | `diagnostics2.py`, `paired_frequency*.py`, `curve_data.py` |
| Development homology candidates | `rag_candidates.json` | `dev_diag.py`, `verify_v1.py` |
| External-model responses | `openai_external_validation/production/raw/production_master_output.jsonl`, `external_validation/*/*.jsonl` | `diagnostics.py`, `verify_v1.py`, `paired_frequency.py`, `make_figures.py` |
| GO release 2026-05-19 (matches the frozen 16-character SHA-256 prefix) | `reference_data/go-basic_2026-05-19.obo` (see `reference_data/GO_FILE_VERIFICATION.md`) | `diagnostics.py`, `diagnostics2.py`, `paired_frequency.py` |
| **PROBE development responses** | `probe_results/` of https://github.com/utpoln/PROBE | `dev_diag.py`, `verify_v1.py` |

**Environment variables** (all optional):

| Variable | Meaning | Default |
|---|---|---|
| `RAMP_PROBE_RESULTS` | Path to `probe_results/` | `../PROBE/probe_results`, i.e. a sibling clone of the PROBE repository |
| `RAMP_GO_OBO` | Path to the GO file. Scripts stop unless its full SHA-256 equals that of the distributed release file (`923645000afcb1df…`, which begins with the frozen 16-character prefix). | `reference_data/go-basic_2026-05-19.obo` |
| `RAMP_POSTHOC_OUT` | Output directory | `scripts/posthoc/output/` |

## Running

```bash
git clone https://github.com/utpoln/ramp.git
git clone https://github.com/utpoln/PROBE.git      # sibling of ramp/, for probe_results/
cd ramp
pip install -r scripts/posthoc/requirements.txt
bash scripts/posthoc/run_all.sh
```

`run_all.sh` runs every analysis, draws the figures, and compares each output with `expected/`. It ends with `ALL OUTPUTS REPRODUCED` if everything matches exactly. Runtime is a few minutes; the bootstraps (B = 2,000) dominate.

## Scripts and what they report

All bootstraps resample proteins (all 243), with B = 2,000 and `numpy.random.default_rng(42)`. The one exception is the true-path interval, which uses B = 1,000.

| Script | Output | Content (manuscript / Additional file 1 location) |
|---|---|---|
| `common.py` | — | Frozen loader. Rebuilds the 60,585 claims and the frozen score S = M + P + 5R, and asserts 60,585 claims / 423 supported. Also provides the whole-group and expected-value tie helpers. |
| `diagnostics.py` | `diagnostics_out.json` | Prespecified secondary endpoints (Table 8); matched-budget component comparisons (Table 9); held-out weight sensitivity (Table S6); recall; homology-only baseline (Table S11); GO-identifier validity (Table S7); true-path credit; few-shot, FUNCTION-text and generic-term checks (Tables S8, S9); external overlap; formula re-applied to the external models |
| `diagnostics2.py` | `diagnostics2_out.json` | Term-frequency ranking and within-stratum results; homology versus score at matched budgets; no-homolog proteins; valid-identifier-only analysis |
| `dev_diag.py` | `dev_diag_out.json` | Development-data weight sensitivity and matched-budget ablation (needs PROBE) |
| `verify_v1.py` | `verify_v1_out.json` | Namespace CIs (Table S2); token-cap sensitivity; score-cell composition (Table S10); overlap-conditioned external contrasts (Table 10); external-ensemble interval and sub-ensembles; reproduction of the original 10-model development curve (needs PROBE) |
| `top1_definitions.py` | `top1_definitions_out.json` | Per-protein top-1, pooled and protein-averaged (Table 8) |
| `paired_frequency.py` | `paired_frequency_out.json` | Frozen score versus term-frequency ranking at the 10% and 20% budgets (paired); validity by score band; GPT-5.5 CIs at S ≥ 9/12/14 |
| `paired_frequency_grid.py` | `paired_frequency_grid_out.json` | The same comparison at every protocol budget (Table S12) |
| `curve_data.py` | `curve_data.csv` | Series for Figure 4 |
| `make_figures.py` | `output/figures/*.pdf`, `fig_external_values.json` | Figures 4–6. Every plotted value is asserted against the manuscript before drawing. |
| `check_reproduction.py` | — | Exact comparison of `output/` with `expected/` |

The analysis code was moved here from the authors' working directories. Only the input-path lines were changed, so `expected/` holds the outputs exactly as reported in the manuscript.

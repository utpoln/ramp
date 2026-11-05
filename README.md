# RAMP: A Frozen Multi-Signal Reliability Framework for Prioritizing LLM-Generated Gene Ontology Annotations

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Status: manuscript in preparation](https://img.shields.io/badge/status-manuscript%20in%20preparation-orange.svg)](#)

> **Kallol Naha, Hasan Jamil**
> Department of Computer Science, University of Idaho
> naha7197@vandals.uidaho.edu | jamil@uidaho.edu

---

## Overview

This repository holds the code, frozen protocol, and raw results for a study of whether a **reliability score for LLM-generated Gene Ontology (GO) claims**, built and frozen using only existing development data, concentrates experimentally supported claims on proteins excluded from its construction. A panel of five further models, added after the primary result, examines how the frozen ranking behaves on claims that other models also produce.

The reliability score combines three signals already observable at inference time: cross-model agreement (`n_models`), cross-prompt consistency (`max_prompt_consistency`), and homology-based retrieval support (`rago_support`, from BLAST against a frozen candidate pool). All three components, the five models and four prompts used to compute them, and the primary endpoint were fixed and timestamped **before** any prediction was generated on the held-out cohort (see `FROZEN_PROTOCOL.md` / `frozen_experiment_config.json`).

**Primary held-out validation** (the only confirmatory result; five frozen models, 243 held-out proteins excluded from score construction): at the prespecified 10%-coverage operating point, the score selected 7,431 of 60,585 claims (12.27%; tied score groups are not split), with precision **4.78%** against a **0.70% all-claim (100%-coverage) baseline**, a protein-level paired bootstrap difference of **+4.08 percentage points** (95% CI [+3.71, +4.47]).

**Post-hoc analyses** (exploratory, run after the held-out results were known; `scripts/posthoc/`): at matched review budgets, cross-model agreement and homology support carry the ranking and prompt consistency adds no detectable benefit; a ranking of the same claims by GO-term annotation frequency performs comparably (difference at the 10% budget +0.28 pp, 95% CI [−0.05, +0.63]); most supported claims are common cellular-component terms.

**External-model panel** (post hoc; GPT-5.5, Llama 3.3 70B, Qwen3 32B, Gemma3 12B, GPT-oss-120B; none is among the five models used to calculate the frozen score. GPT-5.5 and GPT-oss-120B played no part in development. Llama 3.3 70B, Qwen3 32B and Gemma3 12B are PROBE development-benchmark models whose development responses were involved in exploratory score construction, and Llama 3.3 70B and Gemma3 12B also took part in the pre-freeze retrieval pilot): among each model's claims that the frozen five-model ensemble also produced, higher frozen scores went with higher precision (score ≥ 9: 13.6–18.2 pp above the overlapping subset's own precision). Overlapping claims are already enriched and non-overlapping claims cannot be scored, so this does not show transfer to other models' full output.

> **Note:** The manuscript and LaTeX source will be released upon journal acceptance.

---

## Key Findings

| # | Result |
|---|---|
| 1 | Primary 10%-coverage precision: **4.78%** vs. **0.70%** baseline (+4.08 pp, 95% CI [+3.71, +4.47]) on 243 held-out proteins |
| 2 | Precision increases with the frozen score across its full range (0.09% at score 2 → 38.28% at score 14), not strictly monotonic at every individual value |
| 3 | Prespecified secondary endpoints (reported after the primary result was known): 21.16% at the 1% coverage target, 7.72% at 5%, 2.52% at 20%; per-protein top-1 31.5% (pooled) |
| 4 | Leave-one-model-out: removing any one frozen model never drops precision below 4.78% (each reduced ensemble has a smaller pool and higher baseline). Post hoc, at matched budgets: full score − `n_models` +0.52 pp [+0.30, +0.75]; full score − (`n_models` + `rago_support`) +0.01 pp [−0.12, +0.13]. The whole-group component ablations are at different coverages and are not directly comparable |
| 5 | Post hoc: term-frequency ranking comparable to the score at 10% and 20% budgets, better at 50%; homology candidates alone close at their own budgets; 39.6% of extracted GO identifiers are invalid for their namespace and concentrate at low scores; 270 of 355 supported claims are cellular-component terms |
| 6 | External-model panel (post hoc, see above): within overlapping claims, cumulative precision 23.85%–28.46% at score ≥ 9 across the five models |
| 7 | Held-out cohort is substantially sparser in experimental GO annotation than the development benchmark (mean 3.84 vs. 7.67 terms/protein, Cohen's *d* = −0.645) — reported as a property of the validation design |

---

## Repository Structure

```
ramp/
├── FROZEN_PROTOCOL.md                  # Full frozen protocol (timestamped; prespecified, not lodged in a public registry; later edited only to rename RAGO -> RAMP, commit c657d0c)
├── frozen_experiment_config.json       # Machine-readable frozen config + 16-character SHA-256 prefixes for six governed files (BLAST FASTA, GO file, four scripts)
├── GO_RELEASE_CORRECTION.md            # Correction: GO release label in the two frozen files is 2026-05-19, not 2026-03-25
│
├── benchmark_master.json               # PROBE's 1,000-protein development set (reused, not modified)
├── pool_proteins.fasta                 # ADAPT's 25,460-protein candidate pool (BLAST retrieval source)
├── pool_blastdb.*                      # Frozen BLAST database built from the pool above
├── pool_go_lookup.json                 # GO-term lookup for the candidate pool
│
├── heldout_benchmark_243.json          # The 243-protein held-out cohort (human 63, mouse 63, yeast 62, E. coli 55)
├── heldout_query_proteins.fasta        # Held-out protein sequences (for BLAST)
├── heldout_rag_candidates.json         # BLAST retrieval candidates for the held-out cohort
├── heldout_blast_hits.tsv              # Raw BLAST hits, held-out cohort
│
├── reference_data/
│   ├── go-basic_2026-05-19.obo         # GO release matching the frozen 16-character SHA-256 hash prefix 923645000afcb1df
│   └── GO_FILE_VERIFICATION.md         # Source URL, size, full SHA-256, and check against the frozen 16-character prefix
│
├── scripts/
│   ├── draw_heldout_sample.py          # Held-out cohort sampling (seed=42, full audit log)
│   ├── simulate_eligibility_rules.py   # Eligibility-rule simulation used before freezing the design
│   ├── fetch_heldout_sequences.py      # Pull sequences for sampled accessions
│   ├── extract_sequences.py            # Sequence extraction utilities
│   ├── build_pool_go_lookup.py         # Build the GO lookup for the candidate pool
│   ├── run_blast_retrieval.py          # BLASTp retrieval against the frozen pool DB
│   ├── check_retrieval_quality.py      # Retrieval sanity checks
│   ├── evaluate_rag.py                 # RAG/BLAST-support evaluation
│   ├── decompose_rag_unmatched.py      # Diagnostic decomposition of unmatched RAG candidates
│   ├── temporal_holdout_rag.py         # Temporal-holdout style checks
│   ├── llm_prompting_heldout.py        # Primary 5-model inference on the held-out cohort (P1–P4)
│   ├── primary_analysis.py             # Frozen score, primary endpoint, paired bootstrap, ablations, leave-one-model-out (from raw responses)
│   ├── cohort_comparison.py            # Held-out vs. development cohort characteristics
│   ├── development_reference.py        # Development-set reference value (needs PROBE probe_results/)
│   ├── llm_prompting_rag.py            # RAG-augmented inference variant (pre-freeze pilot)
│   └── posthoc/                        # POST-HOC exploratory analyses (not in the frozen protocol); run_all.sh, README.md, expected/
│
├── primary_analysis/
│   ├── primary_analysis_report.json    # Output of scripts/primary_analysis.py
│   ├── cohort_comparison_report.json   # Output of scripts/cohort_comparison.py
│   └── development_reference_report.json # Output of scripts/development_reference.py
│
├── heldout_results/                    # Raw primary-model responses (JSONL, one file per model)
│   ├── mistral-large-123b.jsonl
│   ├── qwen2.5-72b.jsonl
│   ├── google_gemma-4-31b.jsonl
│   ├── mixtral-8x7b.jsonl
│   ├── mistral-7b.jsonl
│   └── run_manifest.json               # Full reproducibility manifest
│
├── rag_results/                        # Pre-freeze retrieval-selection pilot on 100 development proteins (Mistral Large, Llama 3.3 70B, Gemma3 12B)
├── rag_evaluation/
│   └── rag_leaderboard.json
│
├── openai_external_validation/         # External-model panel, first stage: GPT-5.5 (post hoc)
│   ├── README.md                       # Subfolder-specific documentation
│   ├── OPENAI_EXTERNAL_PROTOCOL.md     # Protocol for this stress test
│   ├── OPENAI_EXTERNAL_VALIDATION_REPORT.md
│   ├── openai_prompting_heldout.py     # GPT-5.5 inference on the identical 243-protein cohort
│   ├── config.json
│   ├── raw/gpt-5.5-2026-04-23.jsonl    # Raw GPT-5.5 API records
│   ├── parsed/openai_claims.json       # Deduplicated GPT-5.5 claims
│   ├── analysis/                       # Audit + final analysis scripts
│   ├── production/                     # Production batch submission scripts
│   └── pilot_cost_optimization/        # Pilot batch-cost study (preceding the full run)
│
├── external_validation/                # External-model panel, second stage: four open-weight models (post hoc)
│   ├── run_external_model.py           # Generic runner (reuses the frozen P1-P4 prompt builder)
│   ├── analyze_external_model.py       # Generic integrity audit + analysis pipeline
│   ├── qwen3_truncation_sensitivity.py # Qwen3 truncation-sensitivity analysis
│   ├── combined_external_validation_summary.json
│   ├── llama3_3_70b/                   # Llama 3.3 70B raw responses + analysis
│   ├── qwen3_32b/                      # Qwen3 32B raw responses + analysis
│   ├── gemma3_12b/                     # Gemma3 12B raw responses + analysis
│   ├── gpt_oss_120b/                   # GPT-oss-120B raw responses + analysis
│   └── command_a_111b/                 # Discontinued (pilot only, see DISCONTINUED.md)
│
├── reproducibility/                    # How to reproduce the analysis vs. regenerate responses
│   ├── README.md
│   └── ollama.md                       # Optional local-inference guidance (not required)
│
└── manuscript/                         # Manuscript source -- not included in this repository;
                                         # released upon journal acceptance (see note above)
```

---

## Documentation

Every markdown document in this repository, linked from here:

- [`FROZEN_PROTOCOL.md`](FROZEN_PROTOCOL.md) — the full frozen protocol for the primary five-model experiment (score definition, held-out eligibility rule, tie-handling rule, bootstrap procedure, timestamps)
- [`GO_RELEASE_CORRECTION.md`](GO_RELEASE_CORRECTION.md) — correction to the GO release label recorded in `FROZEN_PROTOCOL.md` and `frozen_experiment_config.json` (label only; no file, hash, or result changed)
- **External-model panel (post hoc): GPT-5.5**
  - [`openai_external_validation/README.md`](openai_external_validation/README.md) — subfolder overview
  - [`openai_external_validation/OPENAI_EXTERNAL_PROTOCOL.md`](openai_external_validation/OPENAI_EXTERNAL_PROTOCOL.md) — protocol for the GPT-5.5 stress test
  - [`openai_external_validation/OPENAI_EXTERNAL_VALIDATION_REPORT.md`](openai_external_validation/OPENAI_EXTERNAL_VALIDATION_REPORT.md) — results report
  - [`openai_external_validation/pilot_cost_optimization/PILOT_PROTOCOL.md`](openai_external_validation/pilot_cost_optimization/PILOT_PROTOCOL.md) — pilot batch-cost protocol
  - [`openai_external_validation/pilot_cost_optimization/pilot_report.md`](openai_external_validation/pilot_cost_optimization/pilot_report.md) — pilot results
- **External-model panel (post hoc): open-weight models**
  - [`external_validation/command_a_111b/DISCONTINUED.md`](external_validation/command_a_111b/DISCONTINUED.md) — why the Cohere Command-A 111B run was discontinued after a technical pilot (never used as scientific data)
- **Post-hoc analyses**
  - [`scripts/posthoc/README.md`](scripts/posthoc/README.md) — inputs, dependencies, and how to reproduce every post-hoc value reported in the manuscript
  - [`reference_data/GO_FILE_VERIFICATION.md`](reference_data/GO_FILE_VERIFICATION.md) — provenance and hash check of the GO file
- **Reproducibility**
  - [`reproducibility/README.md`](reproducibility/README.md) — reproducing the analysis from released data vs. regenerating LLM responses, plus the full frozen-artifact index
  - [`reproducibility/ollama.md`](reproducibility/ollama.md) — optional local-inference guidance (Ollama is not required; MindRouter is the backend actually used for every published result)

---

## The Frozen Reliability Score

```
combined_score = n_models + max_prompt_consistency + 5 x rago_support     (range: 2-14)
```

| Component | Definition | Range |
|---|---|---|
| `n_models` | Number of the five frozen models producing this exact GO identifier for this protein/namespace, pooled across P1–P4 | 1–5 |
| `max_prompt_consistency` | Maximum, over the five models, of the number of P1–P4 prompts under which that model produced this exact GO identifier | 1–4 |
| `rago_support` | 1 if the GO identifier appears among the top-15 BLAST-retrieved candidate terms for this protein/namespace, else 0 | 0–1 |

The score, the five models, the four prompts, the BLAST database/parameters, the GO release, the held-out eligibility rule, the tie-handling rule, and the bootstrap procedure were all fixed in `frozen_experiment_config.json` **before** any held-out prediction was generated, and are not retuned based on held-out or external-model results.

---

## Models Evaluated

### Primary (frozen score, five models, P1–P4)

Selected from PROBE's published ten-model leaderboard by a fixed, mechanical rule (ranks 1, 3, 5, 7, 9). The rule was adopted after exploratory development analyses that used all ten PROBE models, and recorded in the frozen protocol before any held-out prediction:

| Model | MindRouter identity |
|---|---|
| Mistral Large 123B | `mistral-large:123b` |
| Qwen2.5 72B | `qwen2.5:72b` |
| Gemma-4 31B | `google/gemma-4-31b` |
| Mixtral 8x7B | `mixtral:8x7b` |
| Mistral 7B | `mistral:7b` |

P5 (closed-set candidate selection) is excluded from the frozen score by design — it is structurally different from the open-ended generation task used by P1–P4, and none of the three score components depend on P5 data.

### External-model panel (post hoc; none used to compute the score)

The panel was not part of the frozen protocol; it was assembled after the primary result, in two stages. Cohere Command-A 111B was selected but discontinued for runtime after a technical pilot and 8 production calls (`external_validation/command_a_111b/DISCONTINUED.md`) and replaced by GPT-oss-120B.

No panel model is one of the five models whose responses are used to calculate the frozen score. Three of them are PROBE development-benchmark models, and their development-set responses were involved in the exploratory construction of the score. Their held-out results are therefore not independent of score development. Separately from that, two of the same three (Llama 3.3 70B and Gemma3 12B) took part in the pre-freeze retrieval pilot on 100 development proteins (`rag_results/`), which evaluated LLM selection from homology-retrieved candidate lists.

| Model | Identifier | Used to calculate the frozen score | PROBE development responses involved in exploratory score construction | Participant in the pre-freeze retrieval pilot (`rag_results/`) |
|---|---|---|---|---|
| GPT-5.5 | `gpt-5.5-2026-04-23` | No | No (not a PROBE development model) | No |
| Llama 3.3 70B | `llama3.3:70b` | No | **Yes** (PROBE development model, rank 2) | **Yes** |
| Qwen3 32B | `Qwen/Qwen3-32B` | No | **Yes** (PROBE development model, rank 4) | No |
| Gemma3 12B | `gemma3:12b` | No | **Yes** (PROBE development model, rank 10) | **Yes** |
| GPT-oss-120B | `openai/gpt-oss-120b` | No | No (not a PROBE development model; same vendor as GPT-5.5) | No |

---

## Held-Out Cohort

243 proteins, disjoint from both PROBE's 1,000-protein development set and ADAPT's 25,460-protein candidate pool (verified zero overlap):

| Organism | N |
|---|---|
| *Homo sapiens* | 63 |
| *Mus musculus* | 63 |
| *Saccharomyces cerevisiae* | 62 |
| *Escherichia coli* | 55 |
| **Total** | **243** |

Eligibility rule (identical to the rule used to build the development pool, not changed after observing the held-out distribution shift): ≥3 total experimental-evidence GO terms (evidence codes IDA, IMP, IPI, IGI, IEP, EXP, HDA, HMP, HGI, HEP) and a non-empty UniProt FUNCTION comment. Sampling seed 42.

---

## Reproduce Results

### Step 1 — Draw the held-out cohort (pre-built sample included)
```bash
python scripts/draw_heldout_sample.py
python scripts/fetch_heldout_sequences.py
```

### Step 2 — BLAST retrieval (RAGO support signal)
```bash
python scripts/run_blast_retrieval.py
python scripts/check_retrieval_quality.py
```

### Step 3 — Primary five-model inference (P1–P4, held-out cohort)
```bash
python scripts/llm_prompting_heldout.py
```
Requires a MindRouter-compatible endpoint; see `scripts/llm_prompting_heldout.py` for configuration. 14,580 calls (243 proteins x 3 namespaces x 5 models x 4 prompts).

### Step 4 — Primary analysis (frozen score, primary endpoint, robustness analyses)
```bash
python scripts/primary_analysis.py
```
Reads only files already in this repository (`heldout_results/*.jsonl`, `heldout_benchmark_243.json`, `scripts/heldout_sample_250.json`, `heldout_rag_candidates.json`); makes no LLM call, BLAST search, or network request. Computes the frozen score, precision at 10%-coverage under the whole-score-group rule, the 100%-coverage baseline, the paired protein-level bootstrap (B = 2,000, seed 42), the per-score table, organism and namespace breakdowns, component ablations, and leave-one-model-out (with each reduced pool's size and baseline). It also checks that the read-only claim→score lookup used by every external-validation analysis is identical to the score rebuilt from raw responses, and reports the concentration-across-proteins statistics, the exact-10% tie-breaking sensitivity value (descriptive only), and aggregate call and token counts. Output: `primary_analysis/primary_analysis_report.json`.

```bash
python scripts/cohort_comparison.py
```
Held-out vs. development cohort characteristics (manuscript Table 5; requires `scipy`). Output: `primary_analysis/cohort_comparison_report.json`.

```bash
python scripts/development_reference.py --probe-results /path/to/PROBE/probe_results
```
Development-set reference value (0.0595, 95% CI [0.0563, 0.0631]; B = 1,000, as recorded in `FROZEN_PROTOCOL.md`) and the development-set leave-one-model-out check, from the PROBE development responses (https://github.com/utpoln/PROBE). Output: `primary_analysis/development_reference_report.json`.

### Step 5 — RAG/BLAST-support evaluation
```bash
python scripts/evaluate_rag.py
```

### Step 6 — External-model panel, GPT-5.5 (post hoc)
```bash
cd openai_external_validation
python openai_prompting_heldout.py
python analysis/audit_and_analyze.py
python analysis/final_analysis.py
```
Requires `OPENAI_API_KEY` as an environment variable — never pass API keys on the command line or hardcode them. See `openai_external_validation/README.md` and `OPENAI_EXTERNAL_PROTOCOL.md` for the full protocol, including a methodological finding on reasoning-token budget exhaustion during early pilot runs.

### Step 7 — Post-hoc exploratory analyses
```bash
git clone https://github.com/utpoln/PROBE.git ../PROBE   # development responses, needed by two scripts
pip install -r scripts/posthoc/requirements.txt
bash scripts/posthoc/run_all.sh
```
Reproduces every post-hoc value in the manuscript (matched-budget comparisons, baselines, GO-identifier validity, sensitivity analyses, external-model contrasts, Figures 4–6) and checks each output against `scripts/posthoc/expected/`. See [`scripts/posthoc/README.md`](scripts/posthoc/README.md).

All raw results (`heldout_results/*.jsonl`, `openai_external_validation/raw/*.jsonl`, `external_validation/<model>/*.jsonl`) and the full frozen-protocol audit trail are included so the primary and external-model panel (post hoc) results can be re-derived directly from source without re-running inference. See [`reproducibility/README.md`](reproducibility/README.md) for exactly which reported values each script regenerates and which are not covered by a released script.

---

## License

MIT License, consistent with the companion PROBE repository.
Development data derived from UniProt Swiss-Prot ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).

---

## Contact

| | |
|---|---|
| **Kallol Naha** | naha7197@vandals.uidaho.edu |
| **Hasan Jamil** | jamil@uidaho.edu |
| **Institution** | Department of Computer Science, University of Idaho |
| **Issues** | Please open a GitHub issue |

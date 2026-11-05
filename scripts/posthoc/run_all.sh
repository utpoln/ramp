#!/usr/bin/env bash
# Runs every post-hoc analysis, then the figure script, then compares outputs with expected/.
# Usage (from the repository root):  bash scripts/posthoc/run_all.sh
# Optional: RAMP_GO_OBO, RAMP_PROBE_RESULTS, RAMP_POSTHOC_OUT (see README.md)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${RAMP_POSTHOC_OUT:-$HERE/output}"
mkdir -p "$OUT"
export RAMP_POSTHOC_OUT="$OUT"
run() { echo "== $1"; python3 "$HERE/$1" "$OUT/$2" > /dev/null; }
run diagnostics.py         diagnostics_out.json
run diagnostics2.py        diagnostics2_out.json
python3 "$HERE/dev_diag.py" > "$OUT/dev_diag_out.json"; echo "== dev_diag.py"
run verify_v1.py           verify_v1_out.json
run top1_definitions.py    top1_definitions_out.json
run paired_frequency.py    paired_frequency_out.json
run paired_frequency_grid.py paired_frequency_grid_out.json
run curve_data.py          curve_data.csv
echo "== make_figures.py"; python3 "$HERE/make_figures.py"
cp "$OUT/figures/fig_external_values.json" "$OUT/fig_external_values.json"
python3 "$HERE/check_reproduction.py" "$OUT"

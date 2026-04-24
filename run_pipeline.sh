#!/bin/bash
# ============================================================
# VisiumHD Pipeline Launcher
# ============================================================
# Usage:
#   ./run_pipeline.sh my_run.yml                  # recommended
#   ./run_pipeline.sh my_run.yml -resume          # resume a failed run
#   ./run_pipeline.sh --input_dir ... --probe_set ...  # CLI flags also work
# ============================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Check prerequisites ──────────────────────────────────────
ok=true

if ! command -v nextflow &> /dev/null; then
    echo "[ERROR] Nextflow not found. Install: curl -s https://get.nextflow.io | bash"
    ok=false
fi

if ! command -v conda &> /dev/null && ! command -v mamba &> /dev/null; then
    echo "[ERROR] Neither conda nor mamba found in PATH"
    ok=false
fi

if ! command -v spaceranger &> /dev/null; then
    echo "[WARNING] spaceranger not in PATH — set spaceranger_bin in your params file if installed elsewhere"
fi

$ok || exit 1

# ── Detect params file vs CLI flags ──────────────────────────
PARAMS_ARGS=()
EXTRA_ARGS=()

if [[ $# -ge 1 && -f "$1" && "$1" == *.yml ]]; then
    PARAMS_FILE="$1"
    shift
    PARAMS_ARGS=(-params-file "$PARAMS_FILE")
    echo "Using params file: $PARAMS_FILE"
fi

# Remaining args passed through (e.g. -resume, --input_dir, etc.)
EXTRA_ARGS=("$@")

# ── Run ──────────────────────────────────────────────────────
echo ""
echo "==========================================="
echo " VisiumHD Processing Pipeline"
echo " Pipeline : ${PIPELINE_DIR}"
echo " Work dir : ${PWD}/work"
echo "==========================================="
echo ""

nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile conda \
    -w "${PWD}/work" \
    "${PARAMS_ARGS[@]}" \
    "${EXTRA_ARGS[@]}"

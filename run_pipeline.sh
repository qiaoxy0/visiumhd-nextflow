#!/bin/bash
# run_pipeline.sh — VisiumHD Pipeline launcher
# Usage: run_pipeline.sh --input_dir /path/to/samples --spaceranger_ref /path/to/ref --probe_set /path/to/probes [options]
#
# This script can be called from any directory. The Nextflow work/ directory
# and .nextflow/ logs will be created in the caller's current directory,
# keeping the pipeline installation read-only and shareable.

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "==========================================="
echo " VisiumHD Processing Pipeline"
echo " Pipeline location: ${PIPELINE_DIR}"
echo " Working directory: ${PWD}"
echo "==========================================="

# Check prerequisites
if ! command -v nextflow &> /dev/null; then
    echo "[ERROR] Nextflow is not installed or not in PATH"
    echo "Install with: curl -s https://get.nextflow.io | bash"
    exit 1
fi

if ! command -v conda &> /dev/null && ! command -v mamba &> /dev/null; then
    echo "[WARNING] Neither conda nor mamba found in PATH"
    echo "The pipeline requires conda for environment management"
fi

# Run the pipeline from the caller's directory
nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile conda \
    -w "${PWD}/work" \
    "$@"

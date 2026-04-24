#!/bin/bash
# ============================================================
# VisiumHD Pipeline — First-Time Setup
# ============================================================
# Run this once after cloning to check prerequisites and
# pre-build conda environments so pipeline runs start fast.
#
# Usage:  ./setup.sh
# ============================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
PASS="\033[32m✓\033[0m"
FAIL="\033[31m✗\033[0m"
WARN="\033[33m!\033[0m"

echo ""
echo "==========================================="
echo " VisiumHD Pipeline — Setup"
echo "==========================================="
echo ""

# ── 1. Check Nextflow ────────────────────────────────────────
echo -n "Checking Nextflow... "
if command -v nextflow &> /dev/null; then
    echo -e "${PASS} $(nextflow -version 2>&1 | head -3 | tail -1 | xargs)"
else
    echo -e "${FAIL} Not found"
    echo "  Install: curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/"
fi

# ── 2. Check Conda ──────────────────────────────────────────
echo -n "Checking Conda... "
if command -v mamba &> /dev/null; then
    echo -e "${PASS} mamba $(mamba --version 2>&1 | head -1)"
elif command -v conda &> /dev/null; then
    echo -e "${PASS} $(conda --version 2>&1)"
else
    echo -e "${FAIL} Not found"
    echo "  Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
fi

# ── 3. Check SpaceRanger ────────────────────────────────────
echo -n "Checking SpaceRanger... "
if command -v spaceranger &> /dev/null; then
    echo -e "${PASS} $(spaceranger --version 2>&1 | head -1)"
else
    echo -e "${WARN} Not in PATH"
    echo "  Download from: https://www.10xgenomics.com/support/software/space-ranger/downloads"
    echo "  Then add to PATH or set spaceranger_bin in your params file"
fi

# ── 4. Check R (optional, conda will provide it) ────────────
echo -n "Checking R... "
if command -v Rscript &> /dev/null; then
    echo -e "${PASS} R $(Rscript --version 2>&1 | head -1)"
else
    echo -e "${WARN} Not in system PATH (conda env will provide R 4.3)"
fi

echo ""

# ── 5. Pre-build conda environments ─────────────────────────
echo "==========================================="
echo " Pre-building conda environments"
echo " (this may take 10-30 minutes on first run)"
echo "==========================================="
echo ""

CONDA_CMD="conda"
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
fi

CACHE_DIR="${PIPELINE_DIR}/.conda_cache"
mkdir -p "${CACHE_DIR}"

# Build Python environment (bin2cell + dependencies)
echo "Building Python environment (visiumhd)..."
if [ -d "${CACHE_DIR}/visiumhd" ] || ls "${CACHE_DIR}"/visiumhd-* 1>/dev/null 2>&1; then
    echo -e "  ${PASS} Already exists, skipping"
else
    $CONDA_CMD env create -f "${PIPELINE_DIR}/env/visiumhd.yml" -p "${CACHE_DIR}/visiumhd" -y
    echo -e "  ${PASS} Done"
fi

# Build R/Seurat environment
echo "Building R/Seurat environment (seurat)..."
if [ -d "${CACHE_DIR}/seurat" ] || ls "${CACHE_DIR}"/seurat-* 1>/dev/null 2>&1; then
    echo -e "  ${PASS} Already exists, skipping"
else
    $CONDA_CMD env create -f "${PIPELINE_DIR}/env/seurat.yml" -p "${CACHE_DIR}/seurat" -y
    echo -e "  ${PASS} Done"
fi

# ── 6. Set permissions ──────────────────────────────────────
echo ""
echo "Setting file permissions..."
chmod +x "${PIPELINE_DIR}/run_pipeline.sh"
chmod +x "${PIPELINE_DIR}/setup.sh"
chmod +x "${PIPELINE_DIR}/bin/"*
chmod -R 775 "${CACHE_DIR}" 2>/dev/null || true
echo -e "${PASS} Done"

# ── 7. Done ─────────────────────────────────────────────────
echo ""
echo "==========================================="
echo -e " ${PASS} Setup complete!"
echo "==========================================="
echo ""
echo " Next steps:"
echo "   1. Copy the params template:"
echo "        cp params.yml.template my_run.yml"
echo ""
echo "   2. Edit with your paths:"
echo "        nano my_run.yml"
echo ""
echo "   3. Run the pipeline:"
echo "        ./run_pipeline.sh my_run.yml"
echo ""
echo "   To resume a failed run:"
echo "        ./run_pipeline.sh my_run.yml -resume"
echo ""

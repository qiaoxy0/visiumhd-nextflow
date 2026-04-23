// SpaceRanger processing module

process SPACERANGER_COUNT {
    tag "$sample_id"
    publishDir "${params.output_dir}/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sample_dir)
    path(spaceranger_ref)
    path(probe_set)
    
    output:
    tuple val(sample_id), path("spaceranger_${sample_id}/outs"), path("original_image.tif"), emit: spaceranger_output
    
    script:
    """
    # Convert sample_dir to absolute path
    SAMPLE_DIR_ABS=\$(realpath "${sample_dir}")
    if [[ ! -d "\$SAMPLE_DIR_ABS" ]]; then
        echo "[ERROR] Directory not found: \$SAMPLE_DIR_ABS"
        exit 1
    fi

    # Check if SpaceRanger is available
    if ! command -v spaceranger &> /dev/null; then
        echo "ERROR: SpaceRanger not found!"
        echo "Download from: https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest"
        exit 1
    fi

    # Display SpaceRanger version
    spaceranger --version
    
    # Calculate memory for SpaceRanger (reserve 10% for system)
    # SpaceRanger expects memory in GB as integer
    TOTAL_MEM_GB=\$(echo "${task.memory}" | sed 's/\\.GB//' | sed 's/GB//' | sed 's/ //')
    SPACERANGER_MEM=\$(python3 -c "import math; print(int(math.floor(float('\$TOTAL_MEM_GB') * 0.9)))")
    
    echo "[INFO] Total allocated memory: ${task.memory}"
    echo "[INFO] Memory for SpaceRanger: \${SPACERANGER_MEM} GB"
    echo "[INFO] CPUs allocated: ${task.cpus}"

    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] Sample dir: ${sample_dir}"
    echo "[INFO] Probe set: ${probe_set}"

    # Find required input files
    FASTQ_DIR="\$(find "\$SAMPLE_DIR_ABS" -type f -name '*.fastq.gz' -print -quit | xargs -r dirname || true)"
    if [[ -z "\$FASTQ_DIR" ]]; then
        echo "[ERROR] No FASTQ files (*.fastq.gz) found under \$SAMPLE_DIR_ABS"
        exit 2
    fi
    echo "[DEBUG] FASTQ_DIR: \$FASTQ_DIR"

    # Loupe alignment JSON (if present)
    LOUPE_ALIGNMENT="\$(find "\$SAMPLE_DIR_ABS" -type f -name '*.json' -print -quit || true)"
    echo "[DEBUG] LOUPE_ALIGNMENT: \$LOUPE_ALIGNMENT"

    # Select hi-res image = largest tif; pick CytAssist by name if present, otherwise smallest tif
    HI_RES_IMAGE="\$(find "\$SAMPLE_DIR_ABS" -type f \\( -name '*.tif' -o -name '*.tiff' \\) -printf '%s\\t%p\\n' | sort -nr | head -n1 | cut -f2- || true)"
    CYTASSIST_IMAGE="\$(find "\$SAMPLE_DIR_ABS" -type f \\( -name '*[Cc]yt[Aa]ssist*.tif' -o -name '*[Cc]yt[Aa]ssist*.tiff' \\) -print -quit || true)"
    if [[ -z "\$CYTASSIST_IMAGE" ]]; then
        CYTASSIST_IMAGE="\$(find "\$SAMPLE_DIR_ABS" -type f \\( -name '*.tif' -o -name '*.tiff' \\) -printf '%s\\t%p\\n' | sort -n | head -n1 | cut -f2- || true)"
    fi
    echo "[DEBUG] HI_RES_IMAGE: \$HI_RES_IMAGE"
    echo "[DEBUG] CYTASSIST_IMAGE: \$CYTASSIST_IMAGE"

    # Copy hi-res to stable name for downstream
    if [[ -n "\$HI_RES_IMAGE" ]]; then
        cp "\$HI_RES_IMAGE" original_image.tif
    else
        echo "[WARNING] No hi-res image found, creating empty original_image.tif"
        touch original_image.tif
    fi

    # Parse slide/area from the loupe basename if possible
    SLIDE_ID=""
    AREA=""
    if [[ -n "\$LOUPE_ALIGNMENT" ]]; then
        BASENAME=\$(basename "\$LOUPE_ALIGNMENT" .json)
        # Split basename by '-' and extract slide ID (first two fields) and area (third field)
        IFS='-' read -r SLIDE_PREFIX SLIDE_CODE SLIDE_AREA _ <<< "\$BASENAME" || true
        SLIDE_ID="\${SLIDE_PREFIX}-\${SLIDE_CODE}"
        AREA="\$SLIDE_AREA"
        if [[ -z "\$SLIDE_ID" || -z "\$AREA" ]]; then
            echo "[WARNING] Could not parse SLIDE_ID or AREA from \$LOUPE_ALIGNMENT"
        fi
    fi
    echo "[DEBUG] SLIDE_ID: \$SLIDE_ID"
    echo "[DEBUG] AREA: \$AREA"

    echo "===================================="
    echo "Processing sample: ${sample_id}"
    echo "FASTQ directory: \$FASTQ_DIR"
    echo "Probe set: ${probe_set}"
    echo "Loupe alignment: \$LOUPE_ALIGNMENT"
    echo "Hi-res image: \$HI_RES_IMAGE"
    echo "CytAssist image: \$CYTASSIST_IMAGE"
    echo "Slide ID: \$SLIDE_ID"
    echo "Area: \$AREA"
    echo "Memory: \${SPACERANGER_MEM} GB"
    echo "===================================="
    
    # Build the spaceranger command and only include optional args when available
    CMD=(${params.spaceranger_bin} count --id=spaceranger_${sample_id} --transcriptome=${spaceranger_ref} --probe-set=${probe_set} --fastqs="\$FASTQ_DIR" --create-bam=false --localcores=${task.cpus} --localmem=\${SPACERANGER_MEM})

    [[ -n "\$HI_RES_IMAGE" ]]      && CMD+=(--image "\$HI_RES_IMAGE")
    [[ -n "\$CYTASSIST_IMAGE" ]]   && CMD+=(--cytaimage "\$CYTASSIST_IMAGE")
    [[ -n "\$LOUPE_ALIGNMENT" ]]   && CMD+=(--loupe-alignment "\$LOUPE_ALIGNMENT")
    [[ -n "\$SLIDE_ID" ]]          && CMD+=(--slide "\$SLIDE_ID")
    [[ -n "\$AREA" ]]              && CMD+=(--area "\$AREA")

    echo "[INFO] Running: \${CMD[*]}"
    "\${CMD[@]}"

    echo "[INFO] SpaceRanger finished for ${sample_id}"
    """
}
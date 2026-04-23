// bin2cell segmentation module

process BIN2CELL_SEGMENTATION {
    tag "$sample_id"
    publishDir "${params.output_dir}/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(spaceranger_outs), path(original_image)
    val mpp
    
    output:
    tuple val(sample_id), path("bin2cell_${sample_id}"), emit: segmented_data
    
    script:
    """
    echo "[BIN2CELL] Starting segmentation for sample: ${sample_id}"
    echo "[BIN2CELL] Process ID: \$\$"
    echo "[BIN2CELL] Timestamp: \$(date)"

    # Create output directory
    mkdir -p bin2cell_${sample_id}
    
    # Find required SpaceRanger outputs
    BIN_2UM_DIR="${spaceranger_outs}/binned_outputs/square_002um"
    SOURCE_IMAGE=${original_image}
    SPATIAL_DIR="${spaceranger_outs}/spatial"
    
    echo "Processing bin2cell segmentation for sample: ${sample_id}"
    echo "2um bin directory: \$BIN_2UM_DIR"
    echo "Source image: \$SOURCE_IMAGE"
    echo "Spatial directory: \$SPATIAL_DIR"
    
    # Run bin2cell segmentation using the script from bin/ directory
    python ${projectDir}/bin/b2c_segmentation.py \\
        --path \$BIN_2UM_DIR \\
        --source_image_path ${original_image} \\
        --spaceranger_image_path \$SPATIAL_DIR \\
        --mpp ${mpp} \\
        --out_dir bin2cell_${sample_id} 
    
    echo "[BIN2CELL] Segmentation completed for sample ${sample_id}"
    echo "[BIN2CELL] Output saved to: bin2cell_${sample_id}"
    echo "[BIN2CELL] Timestamp: \$(date)"
    """
}
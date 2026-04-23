process SEURAT_CLUSTERING {
    tag "$sample_id"
    publishDir "${params.output_dir}/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(input_data)
    val option  // "Cellseg" or "Default"
    
    output:
    tuple val(sample_id), path("seurat_${sample_id}"), emit: seurat_output
    
    script:
    """
    # Create output directory
    mkdir -p seurat_${sample_id}
    
    # Run Seurat analysis using the unified script
    # Arguments: sample_id, input_path, output_dir, prefix, option
    Rscript ${projectDir}/bin/seurat_sketch_cluster.R \\
        ${sample_id} \\
        ${input_data} \\
        seurat_${sample_id} \\
        ${sample_id} \\
        ${option}
    
    echo "Seurat analysis completed for sample ${sample_id} (${option})"
    """
}
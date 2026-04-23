
nextflow.enable.dsl=2

def helpMessage() {
    log.info """
    VisiumHD Processing Pipeline
    
    Usage:
        nextflow run main.nf --input_dir <path> --spaceranger_ref <path> --probe_set <path> [options]
    
    Required arguments:
        --input_dir           Path to directory containing sample folders
        --spaceranger_ref     Path to Space Ranger reference genome
        --probe_set           Path to probe set CSV file
    
    Optional arguments:
        --output_dir          Output directory (default: ./results)
        --spaceranger_bin     Path to SpaceRanger executable (default: spaceranger)
        --processing_mode     Processing mode: 'Default' or 'Cellseg' (default: Default)
        --mpp                 Microns per pixel for bin2cell segmentation (default: 0.5)
        --help                Show this help message
    
    Processing modes:
        Default:  FASTQ -> SpaceRanger -> Seurat sketch-clustering -> Marker identification
        Cellseg:  FASTQ -> SpaceRanger -> bin2cell segmentation -> Seurat sketch-clustering -> Marker identification

    Profiles:
        -profile conda        Use Conda environments (recommended)

    Examples:
        # Run with default mode
        nextflow run main.nf -profile conda \\
            --input_dir /path/to/data \\
            --spaceranger_ref /path/to/reference \\
            --probe_set /path/to/probe_set.csv

        # Run with cell segmentation
        nextflow run main.nf -profile conda \\
            --input_dir /path/to/data \\
            --spaceranger_ref /path/to/reference \\
            --probe_set /path/to/probe_set.csv \\
            --processing_mode Cellseg \\
            --mpp 0.5
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

// Check required parameters
if (!params.input_dir) {
    error "Please provide an input directory using --input_dir"
}

if (!params.spaceranger_ref) {
    error "Please provide a SpaceRanger reference using --spaceranger_ref"
}

if (!params.probe_set) {
    error "Please provide a probe set CSV file using --probe_set"
}

// Import modules
include { SPACERANGER_COUNT } from './modules/spaceranger'
include { BIN2CELL_SEGMENTATION } from './modules/bin2cell'
include { SEURAT_CLUSTERING } from './modules/seurat'

workflow {
    
    log.info """
    =========================================
    VisiumHD Processing Pipeline
    =========================================
    Input directory     : ${params.input_dir}
    SpaceRanger ref     : ${params.spaceranger_ref}
    SpaceRanger exec    : ${params.spaceranger_bin}
    Probe set           : ${params.probe_set}
    Output directory    : ${params.output_dir}
    Processing mode     : ${params.processing_mode}
    Microns per pixel   : ${params.mpp}
    =========================================
    Starting analysis...
    """
    // Check if input_dir is a single sample directory or a parent directory
    def inputDir = file(params.input_dir)
    def sample_dirs
    if (inputDir.isDirectory()) {
        // Check if input_dir contains *.fastq.gz files (indicating it's a sample directory)
        def fastqFiles = file("${params.input_dir}/*.fastq.gz")
        def hasFastq = fastqFiles instanceof List ? fastqFiles.any { it.exists() } : fastqFiles.exists()
        if (hasFastq) {
            // Treat input_dir as a single sample directory
            sample_dirs = Channel.of(tuple(inputDir.getName(), inputDir.toAbsolutePath()))
        } else {
            // Treat input_dir as a parent directory containing sample folders
            sample_dirs = Channel.fromPath("${params.input_dir}/*", type: 'dir')
                .filter { it.isDirectory() }
                .map { dir -> tuple(dir.getName(), dir.toAbsolutePath()) }
        }
    } else {
        error "Input directory does not exist: ${params.input_dir}"
    }

    // Debug: Print the sample directories
    sample_dirs.view { sample_id, dir -> "[DEBUG] Processing sample: ${sample_id}, Directory: ${dir}" }

    // Run SpaceRanger count
    log.info "Starting SpaceRanger analysis..."
    spaceranger_results = SPACERANGER_COUNT(
        sample_dirs,
        params.spaceranger_ref,
        params.probe_set
    )

    // Add status monitoring for SpaceRanger
    spaceranger_results.subscribe onNext: { sample_id, outs, image ->
        log.info "SpaceRanger completed for sample: ${sample_id}"
    }
    
    if (params.processing_mode == "Cellseg") {
        // Option 2: Include bin2cell segmentation
        bin2cell_input = spaceranger_results.map { sample_id, spaceranger_outs, original_image ->
            tuple(sample_id, spaceranger_outs, original_image)
        }
        
        segmentation_results = BIN2CELL_SEGMENTATION(
            bin2cell_input,
            params.mpp
        )

        // Add status monitoring for bin2cell
        segmentation_results.subscribe onNext: { sample_id, segmented_data ->
            log.info "bin2cell segmentation completed for sample: ${sample_id}"
        }
        
        // Seurat clustering on segmented data
        clustering_results = SEURAT_CLUSTERING(
            segmentation_results.segmented_data,
            "Cellseg"
        )
        
    } else {
        // Option 1: Direct Seurat clustering on SpaceRanger output
        seurat_input = spaceranger_results.map { sample_id, spaceranger_outs, original_image ->
            tuple(sample_id, spaceranger_outs)
        }
        
        clustering_results = SEURAT_CLUSTERING(
            seurat_input,
            "Default"
        )

        // Add status monitoring for Seurat
        clustering_results.subscribe onNext: { sample_id, seurat_output ->
            log.info "Seurat clustering completed for sample: ${sample_id} (Default mode)"
        }
    }

}

// Workflow completion handler
workflow.onComplete {
    if (workflow.success) {
        log.info """
        =========================================
        Pipeline execution completed successfully!
        =========================================
        Duration    : ${workflow.duration}
        Results dir : ${params.output_dir}
        Mode        : ${params.processing_mode}
        
        Output structure:
        ${params.output_dir}/
        ├── [sample_name]/
        │   ├── spaceranger_[sample]/     # SpaceRanger outputs
        │   ${params.processing_mode == "Cellseg" ? "│   ├── bin2cell_[sample]/        # Cell segmentation results" : ""}
        │   └── seurat_[sample]/          # Seurat analysis results
        └── pipeline_info/
            ├── timeline.html              # Execution timeline
            ├── report.html                # Execution report
            └── trace.txt                  # Execution trace
        =========================================
        """
    } else {
        log.error """
        =========================================
        Pipeline execution failed!
        =========================================
        Error message: ${workflow.errorMessage}
        Check the logs for more details.
        
        To resume the pipeline, run:
        nextflow run main.nf -resume [your parameters]
        =========================================
        """
    }
}

workflow.onError {
    log.error """
    =========================================
    Pipeline error occurred!
    =========================================
    Error: ${workflow.errorMessage}
    Check the .nextflow.log file for details
    =========================================
    """
}

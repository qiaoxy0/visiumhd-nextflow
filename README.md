# VisiumHD Pipeline

A Nextflow pipeline for processing 10x Genomics Visium HD spatial transcriptomics data, from FASTQ processing through cell segmentation and clustering analysis.

## For Lab Members (Quick Start)

```bash
# 1. Clone the repo
git clone https://github.com/<your-org>/VisiumHD.git
cd VisiumHD

# 2. Make sure spaceranger is in your PATH
spaceranger --version

# 3. Run from your own working directory
cd ~/my_analysis
/path/to/VisiumHD/run_pipeline.sh \
  --input_dir /path/to/samples \
  --output_dir ./results \
  --spaceranger_ref /path/to/refdata-gex-GRCm39-2024-A \
  --probe_set /path/to/probe_set.csv \
  --processing_mode Default
```

For long runs, use `screen` or `tmux`.

## Prerequisites

Before running the pipeline, install these on the host machine:

1. **[Nextflow](https://www.nextflow.io/docs/latest/install.html)** (>=21.04.0)
   ```bash
   curl -s https://get.nextflow.io | bash
   sudo mv nextflow /usr/local/bin/
   ```

2. **[Conda](https://docs.conda.io/en/latest/miniconda.html)** (Miniconda or Mamba)

3. **[SpaceRanger](https://www.10xgenomics.com/support/software/space-ranger/downloads)** — add to your `PATH`:
   ```bash
   export PATH=/path/to/spaceranger-x.x.x:$PATH
   spaceranger --version
   ```

4. **Reference data**: genome reference and probe set from 10x Genomics.

## Running the Pipeline

```bash
# Default mode (SpaceRanger segmentation → Seurat clustering)
nextflow run main.nf -profile conda \
  --input_dir /path/to/samples \
  --output_dir /path/to/results \
  --spaceranger_ref /path/to/refdata-gex-GRCm39-2024-A \
  --probe_set /path/to/probe_set.csv

# Cellseg mode (SpaceRanger → bin2cell segmentation → Seurat clustering)
nextflow run main.nf -profile conda \
  --input_dir /path/to/samples \
  --output_dir /path/to/results \
  --spaceranger_ref /path/to/refdata-gex-GRCm39-2024-A \
  --probe_set /path/to/probe_set.csv \
  --processing_mode Cellseg \
  --mpp 0.5
```

Or use the wrapper script from any directory:

```bash
/path/to/VisiumHD/run_pipeline.sh \
  --input_dir /path/to/samples \
  --spaceranger_ref /path/to/ref \
  --probe_set /path/to/probes.csv
```

## Input Directory Structure

The pipeline expects each sample in its own directory under `--input_dir`:

```
input_dir/
├── Sample1/
│   ├── Sample_1_S1_L001_R1_001.fastq.gz
│   ├── Sample_1_S1_L001_R2_001.fastq.gz
│   ├── *.json      # Loupe alignment file
│   ├── *.tif       # High-resolution image
│   └── *.tif       # CytAssist image
├── Sample2/
│   └── ...
```

### Required Files per Sample

1. **FASTQ files**: R1 and R2 paired-end sequencing files
2. **Loupe alignment file**: JSON file from 10x Loupe Browser after alignment
3. **High-resolution image**: TIFF microscopy image
4. **CytAssist image**: CytAssist image

## Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input_dir` | Directory containing sample folders | `/data/visiumhd_samples` |
| `--spaceranger_ref` | Path to SpaceRanger reference genome | `/refs/refdata-gex-GRCh38-2020-A` |
| `--probe_set` | Path to probe set CSV file | `probe_sets/Visium_Mouse_v2.1.0.csv` |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output_dir` | `./results` | Output directory for all results |
| `--processing_mode` | `Default` | Processing mode: `Default` or `Cellseg` |
| `--mpp` | `0.5` | Microns per pixel for bin2cell segmentation |
| `--max_memory` | `256.GB` | Maximum memory allocation |
| `--max_cpus` | `32` | Maximum CPU cores |
| `--help` | `false` | Show help message |

## Processing Modes

### Default Mode
Standard SpaceRanger analysis followed by Seurat clustering:
- SpaceRanger default cell segmentation → Seurat sketch clustering → Marker identification

### Cellseg Mode
Enhanced single-cell segmentation using bin2cell:
- SpaceRanger count → bin2cell segmentation → Seurat sketch clustering → Marker identification

## Output Directory Structure

```
results/
├── Sample1/
│   ├── spaceranger_Sample1/
│   ├── bin2cell_Sample1/       # Only in Cellseg mode
│   │   ├── cell_adata.h5ad
│   │   ├── segmentation_whole.png
│   │   ├── hist_bin.png
│   │   ├── total_count.png
│   │   └── ...
│   └── seurat_Sample1/
│       ├── Sample1_seurat_object.rds
│       ├── Sample1_umap_sketch.png
│       ├── Sample1_markers.csv
│       ├── Sample1_seurat_meta.csv
│       └── ...
├── Sample2/
│   └── ...
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## Deployment on a Shared Server

```bash
# Clone to a shared location
sudo git clone https://github.com/<your-org>/VisiumHD.git /opt/pipelines/VisiumHD

# Pre-build conda environments
cd /opt/pipelines/VisiumHD
nextflow run main.nf -profile conda --help

# Set permissions for lab group
sudo chown -R root:labgroup /opt/pipelines/VisiumHD
chmod -R 755 /opt/pipelines/VisiumHD
chmod -R 775 /opt/pipelines/VisiumHD/.conda_cache
```

Users then run from their own directories — `work/` and results stay in their space.

## Clean Cache After Run

```bash
rm -rf work/ .nextflow/

# List all previous runs
nextflow log

# Resume a failed run
nextflow run main.nf -resume [your parameters]
```

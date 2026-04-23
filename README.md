# VisiumHD Pipeline

A Nextflow pipeline for processing 10x Genomics Visium HD spatial transcriptomics data — from FASTQ processing through cell segmentation and clustering analysis.

## Prerequisites

You should install these first to set up the environment before running the analysis pipeline. Run `setup.sh` script to check after your installation:

| Tool | Install |
|------|---------|
| **Nextflow** (>=21.04) | `curl -s https://get.nextflow.io \| bash` |
| **Conda** or **Mamba** | [Miniconda](https://docs.conda.io/en/latest/miniconda.html) |
| **SpaceRanger** | [10x Downloads](https://www.10xgenomics.com/support/software/space-ranger/downloads) — follow their instructions add SpaceRanger to `PATH` |

Reference genome and probe set files are available from 10x Genomics.

## Quick Start (3 steps)

```bash
# Step 1: Clone and set up (one-time)
git clone https://github.com/qiaoxy0/visiumhd-nextflow.git
cd visiumhd-nextflow
./setup.sh

# Step 2: Configure your run
cp params.yml.template my_run.yml
nano my_run.yml          # fill in your paths

# Step 3: Run
./run_pipeline.sh my_run.yml
```

The `setup.sh` checks prerequisites and pre-builds conda environments. The `params.yml` file is the only thing you need to edit.


## Input Directory Structure

```
input_dir/
├── Sample1/
│   ├── *_R1_001.fastq.gz     # Read 1
│   ├── *_R2_001.fastq.gz     # Read 2
│   ├── *.json                 # Loupe alignment file
│   ├── HE.tif                 # High-resolution H&E image
│   └── CytAssist.tif          # CytAssist image
├── Sample2/
│   └── ...
```

The pipeline auto-detects whether `--input_dir` points to a single sample or a parent directory containing multiple samples.


## What to Put in `my_run.yml`

```yaml
# Required — fill these in
input_dir:       "/data/visiumhd/my_samples"
spaceranger_ref: "/path/to/refdata-gex-mm10-2020-A"
probe_set:       "/path/to/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv"

# Optional — change if needed
output_dir:       "./results"       # where results go
processing_mode:  "Default"         # "Default" or "Cellseg"
mpp:              0.5               # microns per pixel (Cellseg mode)
max_memory:       "256.GB"          # adjust to your server
max_cpus:         32
```

## Common Commands

```bash
# Run the pipeline
./run_pipeline.sh my_run.yml

# Resume a failed/interrupted run
./run_pipeline.sh my_run.yml -resume

# Run in background with screen
screen -S visiumhd
./run_pipeline.sh my_run.yml
# Detach: Ctrl+A, D
# Reattach: screen -r visiumhd

# Show help
./run_pipeline.sh --help

# Clean up after a run
rm -rf work/ .nextflow/
```

## Processing Modes

**Default** — SpaceRanger built-in segmentation → Seurat sketch clustering → marker genes

**Cellseg** — SpaceRanger → bin2cell nuclear segmentation (StarDist) → Seurat sketch clustering → marker genes. Use this for higher-resolution single-cell assignments.

## Output Structure

```
results/
├── Sample1/
│   ├── spaceranger_Sample1/          # SpaceRanger outputs
│   ├── bin2cell_Sample1/             # Cell segmentation (Cellseg mode only)
│   │   ├── cell_adata.h5ad
│   │   ├── segmentation_whole.png
│   │   └── ...
│   └── seurat_Sample1/              # Clustering results
│       ├── Sample1_seurat_object.rds
│       ├── Sample1_umap_sketch.png
│       ├── Sample1_markers.csv
│       ├── Sample1_seurat_meta.csv
│       └── ...
└── pipeline_info/
    ├── report.html                   # Execution report
    ├── timeline.html                 # Timeline visualization
    └── trace.txt                     # Resource usage
```

## Parameters Reference

### Required

| Parameter | Description |
|-----------|-------------|
| `input_dir` | Directory containing sample folders |
| `spaceranger_ref` | SpaceRanger reference genome path |
| `probe_set` | Probe set CSV file path |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `output_dir` | `./results` | Output directory |
| `processing_mode` | `Default` | `Default` or `Cellseg` |
| `mpp` | `0.5` | Microns per pixel (Cellseg mode) |
| `spaceranger_bin` | `spaceranger` | SpaceRanger executable path |
| `max_memory` | `256.GB` | Maximum memory per process |
| `max_cpus` | `32` | Maximum CPU cores per process |
| `max_time` | `48.h` | Maximum wall time per process |

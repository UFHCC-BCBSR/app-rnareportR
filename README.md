# RNA-seq Reporter Differential Expression Analysis Pipeline

A comprehensive Shiny-based workflow for configuring and executing differential expression analysis on RNA-seq data with automated report generation on HiPerGator.

---

## Overview

This pipeline provides

* **Interactive Shiny app** for parameter configuration
* **Three DE methods** DESeq2, edgeR GLM, or limma-voom
* **Flexible gene filtering** with edgeR's filterByExpr
* **Batch effect correction** (optional)
* **Automated HTML reports** with visualizations and enrichment analysis
* **SLURM integration** for HPC execution

---

## Directory Structure

```
rnareport/
├── app.R                           # Shiny app for parameter configuration
├── render-rnaseq-report.sbatch     # SLURM submission script
├── RNAseq_report.Rmd               # Main report template
├── R/
│   ├── analysis.R                  # Main analysis orchestration
│   ├── run-DE.R                    # DE analysis functions (DESeq2/edgeR/limma)
│   ├── filter-genes.R              # Gene filtering functions
│   ├── HRK_funcs.R                 # Plotting and enrichment functions
│   └── helper-funcs.R              # Utility functions
├── output/                         # Generated params files and reports
└── logs/                           # SLURM job logs
```

---

## Quick Start

### 1. Launch the Shiny App

```bash
module load R/4.5
rserver
```

Then access via SSH tunnel or HiPerGator web interface.

### 2. Configure Analysis Parameters

In the Shiny app

* Login with your HiPerGator group credentials
* Browse or upload
  * RSEM directory (with .genes.results files)
  * Sample metadata CSV
  * Contrasts file (txt with Group1-Group2 format)
* Set analysis parameters
  * DE Method limma-voom (default), DESeq2, or edgeR GLM
  * Filtering min count (default 10), min proportion (default 0.7)
  * Grouping variable Column name in metadata (e.g., "Condition")
  * Batch variable Optional batch correction variable
  * CPU cores Number of cores for enrichment analysis (default 4)
* Validate and Generate params.txt
* Submit Job directly from the app

### 3. Manual Job Submission (Alternative)

```bash
sbatch render-rnaseq-report.sbatch 
  --params-file output/my_project_params.txt 
  --title "My RNA-seq Analysis" 
  --output-dir /blue/your-group/results 
  --email your.email@ufl.edu
```

---

## Analysis Methods

### Differential Expression Options

| Method | Best For | Key Features |
|--------|----------|--------------|
| limma-voom | Complex designs, multiple comparisons | Fast, flexible, good power |
| DESeq2 | Small sample sizes (<6 per group) | Conservative, robust, shrinkage |
| edgeR GLM | Very small samples (<4 per group) | Exact test, memory efficient |

### Gene Filtering

Uses edgeR's filterByExpr with configurable parameters

* min_count Minimum count per sample (default 10)
* min_prop Minimum proportion of samples expressing gene (default 0.7)

Automatically adjusts thresholds based on sample size and group structure.

### Batch Correction

Optional batch correction using limma's duplicateCorrelation and removeBatchEffect

* Specify batch variable in Shiny app
* Applied before differential expression testing
* PCA plots show before/after correction

---

## Input Files

### 1. RSEM Directory

Directory containing STAR-RSEM output files

```
rsem_dir/
├── Sample1.genes.results
├── Sample2.genes.results
└── Sample3.genes.results
```

### 2. Sample Metadata (CSV)

```
SampleName,Condition,Batch
Sample1,Control,Batch1
Sample2,Control,Batch1
Sample3,Treatment,Batch2
Sample4,Treatment,Batch2
```

* First column Sample names (must match RSEM filenames)
* Additional columns Experimental variables

### 3. Contrasts File (TXT)

```
Treatment-Control
DrugA-Vehicle
DrugB-Vehicle
```

* Each line defines one comparison (format Group2-Group1)

---

## Output

### Generated Files

```
output_directory/
├── ProjectID_params.txt          # Parameters used
├── ProjectID.Report.html         # Main HTML report
├── voom_plot.png                 # Voom mean-variance plot (limma only)
├── SampleData.csv                # Processed metadata
└── contrast.txt                  # Contrasts tested
```

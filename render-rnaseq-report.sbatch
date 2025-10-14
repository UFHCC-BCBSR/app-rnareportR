#!/bin/bash
#SBATCH --job-name=rnaseq-report
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=01:00:00
#SBATCH --output=logs/rnaseq-report_%j.log
#SBATCH --error=logs/rnaseq-report_%j.error
#SBATCH --account=cancercenter-dept
#SBATCH --qos=cancercenter-dept

# Usage: sbatch render-rnaseq-report.sbatch --params-file my_params.txt --title "My Report" --output-dir /path/to/output --email user@ufl.edu
# NOTE: Run this script from the directory containing RNAseq_report.Rmd

# Parse command line arguments
PARAMS_FILE=""
REPORT_TITLE="RNA-seq Differential Expression Analysis"
OUTPUT_DIR=""
USER_EMAIL=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --params-file)
      PARAMS_FILE="$2"
      shift 2
      ;;
    --title)
      REPORT_TITLE="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --email)
      USER_EMAIL="$2"
      shift 2
      ;;
    *)
      echo "Unknown parameter: $1"
      exit 1
      ;;
  esac
done

# Set email notification if provided
if [[ -n "$USER_EMAIL" ]]; then
  echo "Setting up email notifications to: $USER_EMAIL"
  scontrol update JobId=$SLURM_JOB_ID MailUser=$USER_EMAIL MailType=END,FAIL
fi

# Validate required arguments
if [[ -z "$PARAMS_FILE" ]]; then
  echo "Error: --params-file is required"
  exit 1
fi

if [[ ! -f "$PARAMS_FILE" ]]; then
  echo "Error: Params file not found: $PARAMS_FILE"
  exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
  echo "Error: --output-dir is required"
  exit 1
fi

echo "=========================================="
echo "RNA-seq Report Generation"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Params file: $PARAMS_FILE"
echo "Report title: $REPORT_TITLE"
echo "Output directory: $OUTPUT_DIR"
if [[ -n "$USER_EMAIL" ]]; then
  echo "Email notifications: $USER_EMAIL"
fi
echo "Start time: $(date)"
echo "=========================================="

# Create logs directory if it doesn't exist
mkdir -p logs

# Initialize module system
source /etc/profile.d/modules.sh

# Load modules
module purge
module load R/4.5
module load pandoc/2.2.2.1

# Debug: check if they work
which Rscript
which pandoc
pandoc --version

# Set RSTUDIO_PANDOC
export RSTUDIO_PANDOC=$(which pandoc | xargs dirname)

# Extract sample_id from params file
# Replace lines 76-80 with:
# Extract sample_id from params file (this is the unique project identifier)
SAMPLE_ID=$(grep "^--sample_id" "$PARAMS_FILE" | sed -E "s/^--sample_id[[:space:]]+['\"]?(.+)['\"]?$/\1/" | sed 's/[^a-zA-Z0-9_-]/_/g')
if [[ -z "$SAMPLE_ID" ]]; then
  # Fallback to report_title if sample_id not found
  SAMPLE_ID=$(grep "^--report_title" "$PARAMS_FILE" | sed -E "s/^--report_title[[:space:]]+['\"]?(.+)['\"]?$/\1/" | sed 's/[^a-zA-Z0-9_-]/_/g')
  if [[ -z "$SAMPLE_ID" ]]; then
    SAMPLE_ID="RNAseq_Report"
  fi
fi
if [[ -z "$SAMPLE_ID" ]]; then
  SAMPLE_ID="RNAseq_Report"
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_ID}.Report.html"
echo "Output file: $OUTPUT_FILE"

# Render the report with parameters
Rscript -e "
.libPaths(c('/blue/cancercenter-dept/shared/R/library/4.5/', '/usr/local/lib/R/site-library', '/usr/lib/R/site-library', '/usr/lib/R/library'))
Sys.setenv(RSTUDIO_PANDOC='$RSTUDIO_PANDOC')
rmarkdown::render(
  'RNAseq_report.Rmd',
  params = list(
    params_file = '$PARAMS_FILE',
    report_title = '$REPORT_TITLE'
  ),
  output_file = '$OUTPUT_FILE'
)
"

if [[ $? -eq 0 ]]; then
  echo "=========================================="
  echo "Report completed successfully!"
  echo "Output: $OUTPUT_FILE"
  echo "End time: $(date)"
  echo "=========================================="
else
  echo "=========================================="
  echo "Report generation failed!"
  echo "Check error logs for details"
  echo "End time: $(date)"
  echo "=========================================="
  exit 1
fi

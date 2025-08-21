module load nextflow
module load apptainer

nextflow run /blue/cancercenter-dept/PIPELINES/REPORTING/rmarkdown_report_module/test_rmarkdown.nf \
  -profile singularity \
  --sample_id CX_lung_test_contrast \
  --params_file CX_lung_report_params.txt  \
  --rmarkdown_container /blue/cancercenter-dept/PIPELINES/REPORTING/rmarkdown_report_module/rmarkdown_report.sif \
  --rmarkdown_outdir results/reports \
  --log CX_test.log

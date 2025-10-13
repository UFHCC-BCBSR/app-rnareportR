process RENDER_RNASEQ_REPORT {
  tag "${params.sample_id}"                  // <-- must use params.
  publishDir "${params.rmarkdown_outdir}", mode: 'copy' // <-- must use params.

  input:
  path "RNAseq_report.Rmd"
  path "HRK_funcs.R"
  path "analysis.R"
  path "params_file"
  val rmarkdown_container
  val sample_id
  val rmarkdown_outdir

  output:
  path "*.html"

script:
  """
  REPORT_TITLE=\$(grep "^--report_title" ${params_file} | cut -d' ' -f2- | tr -d '"')

  echo "Rendering report for: \$REPORT_TITLE"

  apptainer exec ${rmarkdown_container} Rscript -e "rmarkdown::render(
    'RNAseq_report.Rmd',
    output_file = '${sample_id}.Report.html',
    params = list(
      params_file = '${params_file}',
      report_title = '\$REPORT_TITLE'
    )
  )" > ${sample_id}.render.log 2>&1
  """

}


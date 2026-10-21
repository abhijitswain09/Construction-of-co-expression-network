nextflow.enable.dsl=2
params.input = null
params.outdir = 'results'
process WGCNA {
  tag 'WGCNA analysis'
  conda "${projectDir}/environment.yml"
  input:
  path expression
  path script
  output:
  path '*.pdf', optional: true, emit: plots
  path '*.RData', optional: true, emit: rdata
  path '*.txt', optional: true, emit: tables
  script:
  """
  cp "${expression}" count_gene_all_rpkm2.csv
  Rscript "${script}"
  """
}
workflow {
  if (!params.input) error 'Provide --input with the tab-delimited expression matrix'
  WGCNA(file(params.input, checkIfExists: true), file("${projectDir}/WGCNA.R", checkIfExists: true))
}
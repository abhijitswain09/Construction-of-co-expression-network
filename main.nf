nextflow.enable.dsl=2
params.input = null
params.outdir = 'results'
process WGCNA {
  tag 'WGCNA analysis'
  conda "${projectDir}/environment.yml"
  input: path expression
  output: path '*.pdf', optional: true; path '*.RData', optional: true; path '*.txt', optional: true
  script:
  """
  cp "${expression}" count_gene_all_rpkm2.csv
  Rscript WGCNA.R
  """
}
workflow {
  if (!params.input) error 'Provide --input with the tab-delimited expression matrix'
  WGCNA(file(params.input, checkIfExists: true))
}
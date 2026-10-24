nextflow.enable.dsl=2
params.input = null
params.sample_columns = '8,14,15,5,6,7'
params.outdir = 'results'
process WGCNA {
    tag 'WGCNA analysis'
    conda "${projectDir}/environment.yml"
    cpus 4
    memory '16 GB'
    input:
    path expression
    path script
    output:
    path '*.pdf', optional: true, emit: plots
    path '*.RData', optional: true, emit: rdata
    path '*.txt', optional: true, emit: tables
    script:
    """
    Rscript "${script}" "${expression}" "${params.sample_columns}"
    """
    stub:
    """
    touch sample-clustering.pdf network.RData module-table.txt
    """
}
workflow {
    if (!params.input) error 'Provide --input with the tab-delimited expression matrix'
    WGCNA(file(params.input, checkIfExists: true), file("${projectDir}/WGCNA.R", checkIfExists: true))
}

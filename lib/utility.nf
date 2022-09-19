process BLAST_MAKEBLASTDB {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
        'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    input:
    tuple val(meta),path(fasta)
    val(args)
    
    output:
     tuple val(meta), path("$meta.id*"), emit: db

    script:
    """
    makeblastdb -in $fasta -title $meta.id -out $meta.id  $args
    """
}


process KLEBORATE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kleborate:3.1.3--pyhdfd78af_0' :
        'biocontainers/kleborate:3.1.3--pyhdfd78af_0' }"
        
    publishDir "${params.output_dir}/kleborate", mode: 'copy'

    input:
    tuple val(meta), path(fastas), val(module)

    output:
    tuple val(meta), path("*.tsv"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kleborate \\
    --threads ${task.cpus} \\
         $args \\
         -p ${module} \\
         -o . \\
         -a $fastas

    mv *_output.txt ${prefix}.${module}.kleborate.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${module}.kleborate.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """
}
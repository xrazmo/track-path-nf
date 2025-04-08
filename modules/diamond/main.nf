process DIAMOND_BLASTX {
    tag "$meta.id against $db_name" 
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.8--h43eeafb_0' :
        'biocontainers/diamond:2.1.8--h43eeafb_0' }"

    publishDir "${params.output_dir}/diamond/${db_name}/*.tsv", mode: 'copy'
    
    input:
    tuple val(meta) , path(fasta), path(db_fa)

    output:
    tuple val(meta), path('*.tsv')  , optional: true, emit: tsv
    tuple val(meta), path("*.log")  , emit: log
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def db_name = db_fa.getBaseName()
    
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    diamond makedb --in ${db_fa} -d ${db_name}
    
    diamond \\
        blastx \\
        --threads ${task.cpus} \\
        --db ${db_name}.dmnd \\
        --query ${fasta_name} \\
        --outfmt 6 qseqid sseqid pident slen qlen length mismatch gapopen qstart qend sstart send evalue bitscore stitle \\
        --out ${meta.id}.${db_name}.diamond.tsv \\
        --log \\
        ${args}

    mv diamond.log ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
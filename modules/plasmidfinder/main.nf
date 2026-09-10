process PLASMIDFINDER {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plasmidfinder:2.1.6--py310hdfd78af_1':
        'biocontainers/plasmidfinder:2.1.6--py310hdfd78af_1' }"

    publishDir "${params.output_dir}/plasmidfinder/${meta.id}", mode: 'copy',
        saveAs: { filename ->
            if (!filename.endsWith('.fsa')) {
                return filename
            } else {
                null
            }
        }

    input:
    tuple val(meta), path(seqs)
    path(db_path)

    output:
    tuple val(meta), path("*.json")                 , emit: json
    tuple val(meta), path("*.txt")                  , emit: txt
    tuple val(meta), path("*.tsv")                  , emit: tsv
    tuple val(meta), path("*-hit_in_genome_seq.fsa"), emit: genome_seq
    tuple val(meta), path("*-plasmid_seqs.fsa")     , emit: plasmid_seq
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.6'
    def is_compressed_fasta = seqs.getName().endsWith(".gz") ? true : false
    fasta_name = seqs.getName().replace(".gz", "")
    """
    if [ "$is_compressed_fasta" == "true" ]; then
        gzip -c -d $seqs > $fasta_name
    fi

    plasmidfinder.py \\
        $args \\
        -i $fasta_name \\
        -p $db_path \\
        -o ./ \\
        -x

    # Rename hard-coded outputs with prefix to avoid name collisions
    mv data.json ${prefix}.json
    mv results.txt ${prefix}.txt
    mv results_tab.tsv ${prefix}.tsv
    mv Hit_in_genome_seq.fsa ${prefix}-hit_in_genome_seq.fsa
    mv Plasmid_seqs.fsa ${prefix}-plasmid_seqs.fsa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: $VERSION
    END_VERSIONS
    """
}

process PLASMIDFINDER_UPDATE{
    
    tag "Updating PlasmidFinder"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plasmidfinder:2.1.6--py310hdfd78af_1':
        'biocontainers/plasmidfinder:2.1.6--py310hdfd78af_1' }"

    publishDir "${params.dataCacheDir}", mode: 'copy'

    output:
        path("plasmidfinder_db/"),emit: db_path

    """ 
   
    curl -L -o plasmidfinder_db.zip https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/master.zip
    unzip plasmidfinder_db.zip
    mv genomicepidemiology-plasmidfinder_db* plasmidfinder_db
    cd plasmidfinder_db

    PLASMID_DB=\$(pwd)
   
    # Install PlasmidFinder database with executable kma_index program
    python3 INSTALL.py kma_index
   
    """

}
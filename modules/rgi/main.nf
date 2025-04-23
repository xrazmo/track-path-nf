process RGI_MAIN {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_1':
        'biocontainers/rgi:6.0.3--pyha8f3691_1' }"

    publishDir "${params.output_dir}/rgi/${meta.id}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.txt') || filename.endsWith('.json')) {
                return filename
            } else {
                null
            }
        }

    input:
    tuple val(meta), path(fasta)
    path(card)
    path(wildcard)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.txt") , emit: tsv
    tuple val(meta), path("temp/") , emit: tmp
    env RGI_VERSION                , emit: tool_version
    env DB_VERSION                 , emit: db_version
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // This customizes the command: rgi load
    def args2 = task.ext.args2 ?: '' // This customizes the command: rgi main
    def prefix = task.ext.prefix ?: "${meta.id}"
    def load_wildcard = ""

    if (wildcard) {
        load_wildcard = """ \\
            --wildcard_annotation ${wildcard}/wildcard_database_v\$DB_VERSION.fasta \\
            --wildcard_annotation_all_models ${wildcard}/wildcard_database_v\$DB_VERSION\\_all.fasta \\
            --wildcard_index ${wildcard}/wildcard/index-for-model-sequences.txt \\
            --amr_kmers ${wildcard}/wildcard/all_amr_61mers.txt \\
            --kmer_database ${wildcard}/wildcard/61_kmer_db.json \\
            --kmer_size 61
        """
    }

    """
    DB_VERSION=\$(ls ${card}/card_database_*_all.fasta | sed "s/${card}\\/card_database_v\\([0-9].*[0-9]\\).*/\\1/")

    rgi \\
        load \\
        $args \\
        --card_json ${card}/card.json \\
        --debug --local \\
        --card_annotation ${card}/card_database_v\$DB_VERSION.fasta \\
        --card_annotation_all_models ${card}/card_database_v\$DB_VERSION\\_all.fasta \\
        $load_wildcard

    rgi \\
        main \\
        $args2 \\
        --num_threads $task.cpus \\
        --output_file ${prefix}.rgi \\
        --input_sequence $fasta

    mkdir temp/
    for FILE in *.xml *.fsa *.{nhr,nin,nsq} *.draft *.potentialGenes *{variant,rrna,protein,predictedGenes,overexpression,homolog}.json; do [[ -e \$FILE ]] && mv \$FILE temp/; done

    RGI_VERSION=\$(rgi main --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p temp
    touch test.json
    touch test.txt

    RGI_VERSION=\$(rgi main --version)
    DB_VERSION=stub_version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """
}

process RGI_UPDATE {

    tag "updating RGI v${version}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_1':
        'biocontainers/rgi:6.0.3--pyha8f3691_1' }"
    
    publishDir path: "${params.database_references_dir}", mode: 'copy', saveAs: { filename -> filename.equals("CARD") ? "CARD" : null }, overwrite: true
    publishDir path: "${params.dataCacheDir}", mode: 'copy', saveAs: { filename -> filename.equals("wildcard") ? "wildcard" : null }


    input:
        val(version)

    output:
        path("wildcard"), emit: wildcard
        path("CARD"), emit: card

    script:
    """
    # Work in the local process directory
    mkdir -p CARD
    cd CARD
    wget https://card.mcmaster.ca/latest/data
    tar -xvf data ./card.json
    rm ./data
    rgi card_annotation -i ./card.json > card_annotation.log 2>&1
    cd ..

    mkdir -p wildcard
    wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
    tar -xjf wildcard_data.tar.bz2 -C wildcard
    gunzip wildcard/*.gz

    rgi wildcard_annotation -i wildcard --card_json CARD/card.json -v ${version} > wildcard_annotation.log 2>&1

    rm ./wildcard_data.tar.bz2
    """
}
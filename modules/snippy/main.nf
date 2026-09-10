process SNIPPY_RUN {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_6' :
        'biocontainers/snippy:4.6.0--hdfd78af_2' }"

    
    publishDir "${params.output_dir}/snippy", mode: 'copy'

    input:
    tuple val(meta), path(reads), path(reference)

    output:
    tuple val(meta), path("${prefix}/${prefix}.tab")              , emit: tab
    tuple val(meta), path("${prefix}/${prefix}.csv")              , emit: csv
    tuple val(meta), path("${prefix}/${prefix}.html")             , emit: html
    tuple val(meta), path("${prefix}/${prefix}.vcf")              , emit: vcf
    tuple val(meta), path("${prefix}/${prefix}.bed")              , emit: bed
    tuple val(meta), path("${prefix}/${prefix}.gff")              , emit: gff
    tuple val(meta), path("${prefix}/${prefix}.bam")              , emit: bam
    tuple val(meta), path("${prefix}/${prefix}.bam.bai")          , emit: bai
    tuple val(meta), path("${prefix}/${prefix}.log")              , emit: log
    tuple val(meta), path("${prefix}/${prefix}.aligned.fa")       , emit: aligned_fa
    tuple val(meta), path("${prefix}/${prefix}.consensus.fa")     , emit: consensus_fa
    tuple val(meta), path("${prefix}/${prefix}.consensus.subs.fa"), emit: consensus_subs_fa
    tuple val(meta), path("${prefix}/${prefix}.raw.vcf")          , emit: raw_vcf
    tuple val(meta), path("${prefix}/${prefix}.filt.vcf")         , emit: filt_vcf
    tuple val(meta), path("${prefix}/${prefix}.vcf.gz")           , emit: vcf_gz
    tuple val(meta), path("${prefix}/${prefix}.vcf.gz.csi")       , emit: vcf_csi
    tuple val(meta), path("${prefix}/${prefix}.txt")              , emit: txt
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def read_inputs = meta.single_end ? "--se ${reads[0]}" : "--R1 ${reads[0]} --R2 ${reads[1]}"
    """
    snippy \\
        $args \\
        --cpus $task.cpus \\
        --ram $task.memory \\
        --outdir $prefix \\
        --reference $reference \\
        --prefix $prefix \\
        $read_inputs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy: \$(echo \$(snippy --version 2>&1) | sed 's/snippy //')
    END_VERSIONS
    """
}

process SNIPPY_CONTIGS_RUN {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_6' :
        'biocontainers/snippy:4.6.0--hdfd78af_2' }"

    
    publishDir "${params.output_dir}/snippy_contig", mode: 'copy'

    input:
    tuple val(meta), path(contigs), path(ref_fa), path(ref_gbk)

    output:
    tuple val(meta), path("${prefix}/${prefix}.tab")              , emit: tab
    tuple val(meta), path("${prefix}/${prefix}.csv")              , emit: csv
    tuple val(meta), path("${prefix}/${prefix}.html")             , emit: html
    tuple val(meta), path("${prefix}/${prefix}.vcf")              , emit: vcf
    tuple val(meta), path("${prefix}/${prefix}.snpEff.vcf")       , emit: snpEff_vcf
    // tuple val(meta), path("${prefix}/${prefix}.snpEff.tab")       , emit: snpEff_tab
    tuple val(meta), path("${prefix}/${prefix}.bed")              , emit: bed
    tuple val(meta), path("${prefix}/${prefix}.gff")              , emit: gff
    tuple val(meta), path("${prefix}/${prefix}.bam")              , emit: bam
    tuple val(meta), path("${prefix}/${prefix}.bam.bai")          , emit: bai
    tuple val(meta), path("${prefix}/${prefix}.log")              , emit: log
    tuple val(meta), path("${prefix}/${prefix}.aligned.fa")       , emit: aligned_fa
    tuple val(meta), path("${prefix}/${prefix}.consensus.fa")     , emit: consensus_fa
    tuple val(meta), path("${prefix}/${prefix}.consensus.subs.fa"), emit: consensus_subs_fa
    tuple val(meta), path("${prefix}/${prefix}.raw.vcf")          , emit: raw_vcf
    tuple val(meta), path("${prefix}/${prefix}.filt.vcf")         , emit: filt_vcf
    tuple val(meta), path("${prefix}/${prefix}.vcf.gz")           , emit: vcf_gz
    tuple val(meta), path("${prefix}/${prefix}.vcf.gz.csi")       , emit: vcf_csi
    tuple val(meta), path("${prefix}/${prefix}.txt")              , emit: txt
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = contigs.getExtension() == "gz" ? true : false
    def contigs_name = is_compressed ? contigs.getBaseName() : contigs

    """
    
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${contigs} > ${contigs_name}
    fi

    # Create snpEff directory structure and config
    mkdir -p ./.snpEff_cache/data/ref
    cp $ref_gbk ./.snpEff_cache/data/ref/genes.gbk

    # Create proper snpEff config with correct syntax
    echo -e "ref.genome : Snippy Reference\ndata_dir = ./data" > ./.snpEff_cache/snpEff.config

    # Build snpEff database
    snpEff build -genbank -config ./.snpEff_cache/snpEff.config -v ref

    # Extract features from reference GenBank
    perl ${baseDir}/bin/extract_features.pl --reference ${ref_gbk} --refdir .

    # Run snippy with proper escaping for Nextflow
    snippy \
        $args \
        --cpus $task.cpus \
        --ram $task.memory \
        --outdir $prefix \
        --reference $ref_fa \
        --prefix $prefix \
        --ctgs $contigs_name

    # Run snpEff annotation
    snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr \
        -config ./.snpEff_cache/snpEff.config ref \
        "${prefix}/${prefix}.filt.vcf" > "${prefix}/${prefix}.snpEff.vcf"

    # Remove version numbers from contig names
  #sed -E 's/^([^#])([^\t]+)\\.[0-9]+(\\t)/\\1\\2\\3/; s/^(##contig=<ID=)([^,]+)\\.[0-9]+(,)/\\1\\2\\3/' \\
   # "${prefix}.snpEff.withversion.vcf" > "${prefix}/${prefix}.snpEff.vcf"

    # Run snippy-vcf_to_tab
   
    #/usr/local/bin/snippy-vcf_to_tab --gff ./ref.gff --ref $ref_fa --vcf "${prefix}/${prefix}.snpEff.vcf" > "${prefix}/${prefix}.snpEff.tab"


    # Generate versions.yml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy: \$(snippy --version 2>&1 | sed 's/snippy //')
    END_VERSIONS
    """
}
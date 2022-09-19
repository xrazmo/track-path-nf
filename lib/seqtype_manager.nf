include { PRODIGAL } from "$baseDir/modules/nf-core/prodigal/main"
include { DIAMOND_MAKEDB } from "$baseDir/modules/nf-core/diamond/makedb/main"
include { DIAMOND_BLASTX } from "$baseDir/modules/nf-core/diamond/blastx/main"
 
 process CUS_BLASTN {
    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
        'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    input:
    tuple val(type),val(meta), path(fasta),val(db)
    val(args)
       

    output:
    tuple val(type), val(meta), path("*.${type}.blastn.txt"), emit: txt
   
    script:
    
    """
    blastn -num_threads $task.cpus -db ${db} -query $fasta -out ${meta.id}.${type}.blastn.txt ${args}
    """
}

process PROFILE_FINDER {
    tag "$meta.id"
    label "vshort"
    
    input:
    tuple val(type),val(meta), path(blastout),path(profile_tsv)
           

    output:
    tuple val(type),val(meta), path('*.profile.txt'), emit: profiles
   
    script:
    
    """
    python $baseDir/bin/profile_finder.py --profile $profile_tsv --blastout $blastout --type $type
    """
}

 workflow SEQUENCE_TYPING{
    take:
        input_gz
        data_dir
        annotate
    
    main:
        genomes_gz=Channel.fromPath("$input_gz/**/*gz")
                   .map{file ->[["id":file.baseName.replace('.','-v-'),"species":file.parent.simpleName],file] }

        PRODIGAL(genomes_gz , "gff")

        ch_orfs = PRODIGAL.out.nucleotide_fasta
        if(annotate){
            def diamond_cols = "qseqid sseqid pident qcovhsp scovhsp mismatch gaps evalue bitscore length qlen slen qstart qend sstart send stitle"
            dmnd_ch = Channel.fromPath("${data_dir}/dmnd/*.dmnd")
            ch_orfs.combine(dmnd_ch).combine(Channel.from('txt')).combine(Channel.from(tuple(diamond_cols)))
                    .map{it -> [[id:it[0].id+"__"+it[2].simpleName],it[1],it[2],it[3],it[4]]}
                    .multiMap {it -> fasta: tuple(it[0],it[1])
                                    db: it[2]
                                    ext: it[3]
                                    col: it[4]}
                    .set{diaParam}

            DIAMOND_BLASTX(diaParam.fasta,diaParam.db,diaParam.ext,diaParam.col)
        }

        ch_blasts_db = Channel.fromPath("${data_dir}/mlst/blastdb/*.nto")
            .map{it -> ['mlst',it.simpleName.split('-')[0],it.parent/it.simpleName]}
            .concat(
                    Channel.fromPath("${data_dir}/cgmlst/blastdb/*.nto")
                     .map{it -> ['cgmlst',it.simpleName.split('-')[0],it.parent/it.simpleName]}
                   )

        inputs_db = ch_blasts_db.combine(ch_orfs).filter(it->it[1]==it[3].species).map{it-> [it[0],it[3],it[4],it[2]]}


        CUS_BLASTN(inputs_db,"-max_target_seqs 1 -outfmt '6 qseqid sacc pident length slen qlen'")

        ch_profiles_db = Channel.fromPath("${data_dir}/mlst/**/*.profiles.tsv")
                .map{file->['mlst',file.parent.simpleName.split('-')[0],file]}
                .concat(
                            Channel.fromPath("${data_dir}/cgmlst/**/*.loci.tsv")
                            .map{file->['cgmlst',file.parent.simpleName.split('-')[0],file]}
                        )
            
        input_profile = ch_profiles_db.combine(CUS_BLASTN.out.txt)
                    .filter(it->(it[0]==it[3]) && (it[1]==it[4].species))
                    .map{it-> [it[0],it[4],it[5],it[2]]}
        
        PROFILE_FINDER(input_profile)

    emit:
       gene_annotations = PRODIGAL.out.gene_annotations
       nucleotide_fasta = PRODIGAL.out.nucleotide_fasta
       amino_acid_fasta = PRODIGAL.out.amino_acid_fasta
       diamond_txt = DIAMOND_BLASTX.out.txt
       seqtypes =PROFILE_FINDER.out  
}
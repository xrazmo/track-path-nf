include { PRODIGAL } from "$baseDir/modules/nf-core/prodigal/main"
include { DIAMOND_MAKEDB } from "$baseDir/modules/nf-core/diamond/makedb/main"
include { DIAMOND_BLASTX } from "$baseDir/modules/nf-core/diamond/blastx/main"
include { TRIMGALORE } from "$baseDir/modules/nf-core/trimgalore/main"
include { SPADES } from "$baseDir/modules/nf-core/spades/main"
include { KRAKEN2_KRAKEN2 } from "$baseDir/modules/local/kraken2"

 process BLASTN {
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
    tuple val(type),val(meta), path(blastout), path(profile_tsv)
           

    output:
    tuple val(type),val(meta), path('*.profile.txt'), emit: profiles
   
    script:
    
    """
    python $baseDir/bin/profile_finder.py --profile $profile_tsv --blastout $blastout --type $type
    """
}

process KRAKEN2_PARSER{
    tag "$meta.id"
    label "vshort"
    
    input:
    tuple val(meta), path(report)

    output:
    path("*.taxa.csv"), emit: taxa
   
    script:
    
    """
    python $baseDir/bin/parse_kraken2.py -i $report
    """
}
 workflow ANNOTATE_ISOLATES{
    take:
        input_ch
        data_dir
        annotate
    
    main:
       
        paired_ch = input_ch.filter{it.seq_rep=='paired'}
                                .map{it->[[id:it.sample.replace('.','-v-'),single_end:false],[it.file1,it.file2]]}
        genome_ch = input_ch.filter{it.seq_rep=='genome'}
                                .map{it->[[id:it.sample.replace('.','-v-'),single_end:true],[it.file1]]}
        
        taxa_csv = FIND_TAXA(paired_ch.concat(genome_ch),"$data_dir/kraken2/minikraken2_v1_8GB")
        taxa_ch = taxa_csv.splitCsv(header:true)
        TRIMGALORE(paired_ch)  

        SPADES(TRIMGALORE.out.reads.map{it -> [it[0],it[1],[],[]]},[])

        genome_reshaped_ch = genome_ch.map{it-> [it[0].id,it[1][0]]}.join(taxa_ch.map{it->[it.id,it]})
                  .map{it->[it[2],it[1]]}
        contigs_ch = SPADES.out.contigs.map{it->[it[0].id,it[1]]}
                                       .join(taxa_ch.map{it->[it.id,it]})
                                       .map{it->[it[2],it[1]]}

        PRODIGAL(genome_reshaped_ch.concat(contigs_ch) , "gff")

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

        inputs_ch = ch_blasts_db.combine(ch_orfs).filter(it->it[1]==it[3].tag).map{it-> [it[0],it[3],it[4],it[2]]}


        BLASTN(inputs_ch,"-max_target_seqs 1 -outfmt '6 qseqid sacc pident length slen qlen'")

        ch_profiles_db = Channel.fromPath("${data_dir}/mlst/**/*.profiles.tsv")
                .map{file->['mlst',file.parent.simpleName.split('-')[0],file]}
                .concat(
                            Channel.fromPath("${data_dir}/cgmlst/**/*.loci.tsv")
                            .map{file->['cgmlst',file.parent.simpleName.split('-')[0],file]}
                        )
            
        input_profile = ch_profiles_db.combine(BLASTN.out.txt)
                    .filter(it->(it[0]==it[3]) && (it[1]==it[4].tag))
                    .map{it-> [it[0],it[4],it[5],it[2]]}
        
        PROFILE_FINDER(input_profile)

    emit:
       taxa_csv = taxa_csv 
       contigs = SPADES.out.contigs
       gene_annotations = PRODIGAL.out.gene_annotations
       diamond_txt = DIAMOND_BLASTX.out.txt
       seqtypes =PROFILE_FINDER.out  
       nucleotide_fasta = PRODIGAL.out.nucleotide_fasta
       amino_acid_fasta = PRODIGAL.out.amino_acid_fasta
      
}

workflow FIND_TAXA{
    take:
        input_ch
        kraken2_db
    main:
                                            
        KRAKEN2_KRAKEN2(input_ch,kraken2_db,false,false)
        KRAKEN2_PARSER(KRAKEN2_KRAKEN2.out.report)

    emit:
        taxa=KRAKEN2_PARSER.out.taxa
}
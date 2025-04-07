#!/usr/bin/nextflow

// Enable DSL2 syntax
nextflow.enable.dsl=2

include {FASTQC} from "$baseDir/modules/fastqc/main"
include {TRIMMOMATIC} from "$baseDir/modules/trimmomatic/main"
include {MLST} from "$baseDir/modules/mlst/main"
include {SPADES} from "$baseDir/modules/spades/main"
include {QUAST} from "$baseDir/modules/quast/main"
include {PROKKA} from "$baseDir/modules/prokka/main"
include {SNIPPY_RUN} from "$baseDir/modules/snippy/main"
include {AMRFINDERPLUS_RUN} from "$baseDir/modules/amrfinderplus/run/main"
include {AMRFINDERPLUS_UPDATE} from "$baseDir/modules/amrfinderplus/update/main"

// Parameters with default values that can be overridden
params.input_dir = ""
params.output_dir = ""
params.reference_gbk = "$baseDir/assets/NC_002516v2.gbk"
params.reference_fa = "$baseDir/assets/NC_002516v2.fasta"
params.prodigal_tf="$baseDir/assets/Pseudomonas_aeruginosa.trn"


workflow {

    fq_ch = Channel
    .fromFilePairs("${params.input_dir}/*/*_{1,2}.{fq.gz,fastq.gz}")
    .concat(Channel.fromFilePairs("${params.input_dir}/*_{1,2}.{fq.gz,fastq.gz}"))
    .map{it->[[id:it[0],single_end:false],it[1]]}
    .unique { it[0] }


    AMRFINDERPLUS_UPDATE()

    FASTQC(fq_ch)

    TRIMMOMATIC(fq_ch)

    SPADES(TRIMMOMATIC.out.trimmed_reads.map{it -> [it[0],it[1],[],[]]},[])

    QUAST(SPADES.out.contigs.combine(Channel.of(file(params.reference_fa))))

    MLST(SPADES.out.contigs)

    PROKKA(SPADES.out.contigs,[],params.prodigal_tf)

    AMRFINDERPLUS_RUN(SPADES.out.contigs,AMRFINDERPLUS_UPDATE.out.db)

    SNIPPY_RUN(fq_ch,params.reference_gbk)
}
workflow PrepRef{
    ref_contigs_ch = Channel.of(file(params.reference_fa)).map{it->[[id:it.simpleName],it]}
    PROKKA(ref_contigs_ch,[],params.prodigal_tf)
}
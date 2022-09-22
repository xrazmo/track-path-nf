nextflow.enable.dsl=2

include {ANNOTATE_ISOLATES} from "$baseDir/modules/local/annotate_isolate"
include {CREATD_REF_DB} from "$baseDir/modules/local/db_manager"
include {SAVE_TO_DB} from "$baseDir/modules/local/db_manager"

def data_dir = "$baseDir/data"
def input_gz = "$baseDir/gz"
def ref_db = "$data_dir/db.sqlite3"
def contigs_dir ="$data_dir/contigs"

workflow{

    CREATD_REF_DB(ref_db,true)

    input_ch = Channel.fromPath("$input_gz/samples.csv").splitCsv(header:true);
    ANNOTATE_ISOLATES(input_ch,data_dir,true)

    ch_into_db =ANNOTATE_ISOLATES.out.taxa_csv.map{it -> ["samples",[id:it.simpleName],it]}
                .concat(ANNOTATE_ISOLATES.out.seqtypes)
                .concat(ANNOTATE_ISOLATES.out.gene_annotations.map{it->['orfs',it[0],it[1]]})
                .concat(ANNOTATE_ISOLATES.out.diamond_txt.map{it->['annotations',it[0],it[1]]})

    file(contigs_dir).mkdir()
    ANNOTATE_ISOLATES.out.contigs.map{it-> it[1]}.flatten().collectFile(storeDir:contigs_dir)

    
    SAVE_TO_DB(ch_into_db,ref_db)

 
}
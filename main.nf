nextflow.enable.dsl=2

include {SEQUENCE_TYPING} from "$baseDir/lib/seqtype_manager.nf"
include {REMOTE_SEQ_FETCHER} from "$baseDir/lib/biodb_manager.nf"

def data_dir = "$baseDir/data"
def input_gz = "$baseDir/gz"
def ref_db = "$data_dir/db.sqlite3"


process SAVE_TO_DB{
    tag "$meta.id"
    label "vshort"

    maxForks 1
    input:
    tuple val(type),val(meta), path(intsv)
    val(db_path)

    script:
    """
    python $baseDir/bin/save_to_db.py --db $db_path --in $intsv --type $type
    """
}

workflow{
    
    //    REMOTE_SEQ_FETCHER("$data_dir/selected_assembly_records_smpl-1000.tsv",input_gz)
    
    SEQUENCE_TYPING(input_gz,data_dir,true)

    ch_into_db =SEQUENCE_TYPING.out.seqtypes
               .concat(SEQUENCE_TYPING.out.gene_annotations.map{it->['orfs',it[0],it[1]]})
                .concat(SEQUENCE_TYPING.out.diamond_txt.map{it->['annotations',it[0],it[1]]})

    SAVE_TO_DB(ch_into_db,ref_db)

 
}
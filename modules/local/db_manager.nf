 
 
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

process CREATD_REF_DB{
    label "vshort"
    input:
        val db_path
        val replace

    script:
    """
    python $baseDir/bin/create_ref_db.py --db $db_path --replace $replace
    """
}
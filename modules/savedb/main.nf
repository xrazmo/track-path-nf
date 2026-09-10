process SAVE_TO_DB {
    tag "${db_name}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val ready
    val results_dir
    val db_name

    output:
    path "${db_name}", emit: db

    script:
    def add_seq = params.save_db_add_seq ? '--add_seq' : ''
    """
    python ${baseDir}/bin/save_to_db.py \\
        --input_dir ${results_dir} \\
        --output_dir . \\
        --db_name ${db_name} \\
        --cpus ${task.cpus} \\
        ${add_seq}
    """
}

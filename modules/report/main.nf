process GENERATE_REPORT {
    tag "report"
    label 'process_medium'

    conda "${baseDir}/modules/savedb/environment.yml"

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val ready
    val results_dir

    output:
    path "report", emit: report_dir

    script:
    """
    python ${baseDir}/bin/generate_report.py \\
        --input_dir ${results_dir} \\
        --output_dir report \\
        --ui_dir ${baseDir}/assets/ui
    """
}

process GENE_DIFF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${baseDir}/modules/genediff/genediff/environment.yml"

    publishDir "${params.output_dir}/genediff", mode: 'copy', 
        saveAs: { filename ->
            if (filename.endsWith('.json')) {
                return filename
            } else {
                null
            }
        }

    input:
    tuple val(meta), path(ref_fna),path(qry_ffa)

    output:
    tuple val(meta), path("*.json"), emit: json

    """
        python ${baseDir}/modules/genediff/genediff/scripts/gene_diff.py \\
        --reference ${ref_fna} \\
        --query ${qry_ffa} \\
        --output .
       
    """

}
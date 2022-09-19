
process DOWNLOAD_REF_DBs{
    input:
    val(out_dir)
    val(ignore_content)

    script:
    def ignore = ignore_content?'--ignore': ''
    """
    python $baseDir/bin/download_ref_db.py --out_dir $out_dir --ref_json $dbs_config_file  $ignore
    """
   
}

workflow DATABASE_INIT{
   
    DOWNLOAD_REF_DBs(Channel.fromPath(params.out_dir),false)
    
    ch_diamond_fa = Channel.fromPath("${params.out_dir}/dmnd/*.{fasta,fa,fna,fas}")
    
    ch_cgmlst_alle = Channel.fromPath("${params.out_dir}/cgmlst/**/*alleles.fna",type:'file')
    .map{file -> [["id":file.simpleName,"type":'cgmlst'],file]}.view()  

    ch_mlst_alle = Channel.fromPath("${params.out_dir}/mlst/**/*alleles.fna",type:'file')
                   .map{file -> [["id":file.simpleName,"type":'mlst'],file]}

    ch_db_input = ch_mlst_alle.concat(ch_cgmlst_alle)

    def mlst_blast = "${params.out_dir}/mlst/blastdb"
    def cgmlst_blast = "${params.out_dir}/cgmlst/blastdb"
    def dmnd_dir = "${params.out_dir}/dmnd"
    BLAST_MAKEBLASTDB(ch_db_input,'-dbtype nucl')
    
    BLAST_MAKEBLASTDB.out.db.filter{it[0].type=='mlst'}.map{it -> it[1]}.flatten().collectFile(storeDir:"${mlst_blast}")
    BLAST_MAKEBLASTDB.out.db.filter{it[0].type=='cgmlst'}.map{it -> it[1]}.flatten().collectFile(storeDir:"${cgmlst_blast}")

    DIAMOND_MAKEDB(ch_diamond_fa)
    DIAMOND_MAKEDB.out.db.collectFile(storeDir:"${dmnd_dir}")

} 
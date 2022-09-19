
process DOWNLOAD_NCBI_GENOME{
    tag "$meta.id"
    
    input:
        tuple val(meta), val(ftp)
        val(outDir)

      
    script:
  
    """
        mkdir -p $outDir/${meta.name}
        dirname=`echo $ftp | awk -F/ '{print \$NF}'`
        eval "wget --output-document=${meta.id}.gz $ftp/\${dirname}_genomic.fna.gz"
        mv ${meta.id}.gz $outDir//${meta.name}/${meta.id}.gz
    """
}

workflow REMOTE_SEQ_FETCHER{

    //TODO: if(contigs or reads)
    take:
        genomes_tsv
        outDir

    main:
        ch_genomes = Channel.fromPath(genomes_tsv).splitCsv(header:true,sep:'\t')
               .map{row->[["id":row.assembly_accession,"name":row.species_name],row.ftp_path]}
      
        DOWNLOAD_NCBI_GENOME(ch_genomes,outDir)
}
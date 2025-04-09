#!/usr/bin/nextflow

// Enable DSL2 syntax
nextflow.enable.dsl=2

include {FASTQC} from "$baseDir/modules/fastqc/main"
include {TRIMMOMATIC} from "$baseDir/modules/trimmomatic/main"
include {MLST} from "$baseDir/modules/mlst/main"
include {SPADES} from "$baseDir/modules/spades/main"
include {QUAST} from "$baseDir/modules/quast/main"
include {BUSCO_BUSCO as BUSCO} from "$baseDir/modules/busco/main"
include {PROKKA} from "$baseDir/modules/prokka/main"
include {SNIPPY_RUN} from "$baseDir/modules/snippy/main"
include {AMRFINDERPLUS_RUN} from "$baseDir/modules/amrfinderplus/run/main"
include {AMRFINDERPLUS_UPDATE} from "$baseDir/modules/amrfinderplus/update/main"
include {DIAMOND_BLASTX} from "$baseDir/modules/diamond/main"

include {RGI_UPDATE} from "$baseDir/modules/rgi/main"
include {RGI_MAIN} from "$baseDir/modules/rgi/main"
include {KLEBORATE} from "$baseDir/modules/kleborate/main"
include {PLASMIDFINDER} from "$baseDir/modules/plasmidfinder/main"


// Parameters with default values that can be overridden
params.reads_dir = ""              // Directory containing fastq files
params.contigs_dir = ""            // Directory containing contigs
params.output_dir = ""
params.assets = "$baseDir/assets"
params.reference_dir = "$baseDir/assets/references"
params.database_references_dir = "$baseDir/assets/databases"
params.species_config = "${params.reference_dir}/species_references.config"
params.db_config = "${params.database_references_dir}/database_references.config"


// Function to load species references from config file
def loadSpeciesConfig() {

    // Load species reference configuration from external file
    def speciesReferences = [:]
    def defaultReference = [:]

    def configFile = file(params.species_config)
    if (configFile.exists()) {
        def slurper = new groovy.json.JsonSlurper()
        def config = slurper.parseText(configFile.text)
        // Load the species configurations
        config.species.each { species ->
            speciesReferences[species.name] = [
                gbk: (species.gbk && species.gbk != "") ? "${params.reference_dir}/${species.gbk}":null,
                fasta: (species.fasta && species.fasta != "") ? "${params.reference_dir}/${species.fasta}":null,
                trn: (species.trn && species.trn != "") ? "${params.reference_dir}/${species.trn}":null,
                amrfindopt: (species.amrfindopt && species.amrfindopt != "") ? species.amrfindopt : null
            ]
        }
        
        // Set default if available
        if (config.default) {
            defaultReference = [
                gbk: null,
                fasta: null,
                trn: null,
                amrfindopt: null
            ]
        }
    } else {
        log.warn "Species configuration file not found: ${params.species_config}"
        
        // Fallback default if config file not found
        defaultReference = [
            gbk: null,
            fasta: null,
            trn: null,
            amrfindopt: null
        ]
        
        speciesReferences["unknown"] = defaultReference
    }
    return [speciesReferences: speciesReferences, defaultReference: defaultReference]
}

// Function to load database references for Diamond
def loadDatabaseConfig() {
    def databases = []
    def configFile = file(params.db_config)
    if (configFile.exists()) {
        def slurper = new groovy.json.JsonSlurper()
        def config = slurper.parseText(configFile.text)
        
        // Load database configurations
        config.databases.each { db ->
            databases << "${params.database_references_dir}/${db.path}"
        }
    } else {
        log.warn "Database configuration file not found: ${params.db_config}"
        log.warn "No Diamond databases will be used"
    }
    return databases
}

def loadCardVersion(){
   
    def configFile = file(params.db_config)
    if (configFile.exists()) {
        def slurper = new groovy.json.JsonSlurper()
        def config = slurper.parseText(configFile.text)
        return config.card_version
    } else {
        log.warn "Database configuration file not found: ${params.db_config}"
        log.warn "No wildCard option"
    }
    return "NA"
}

// Load configurations when pipeline starts
def references = loadSpeciesConfig()
def speciesReferences = references.speciesReferences
def defaultReference = references.defaultReference

def refDiamondFa = loadDatabaseConfig()
def cardversion = loadCardVersion()

workflow {
    // Create channels for different input types
    
    // 1. Input channel for fastq files (if directory provided)
    fq_ch = Channel.empty()
    if (params.reads_dir) {
        fq_ch = Channel
            .fromFilePairs("${params.reads_dir}/*/*_{1,2}.{fq.gz,fastq.gz}")
            .concat(Channel.fromFilePairs("${params.reads_dir}/*_{1,2}.{fq.gz,fastq.gz}"))
            .map{it->[[id:it[0],single_end:false,source:"reads"],it[1]]}
            .unique { it[0].id }
    }

    // 2. Input channel for pre-assembled contigs (if directory provided)
    contigs_direct_ch = Channel.empty()
    if (params.contigs_dir) {
        contigs_direct_ch = Channel
            .fromPath("${params.contigs_dir}/**/*.contigs.fa.gz")
            .map { file -> 
                def id = file.simpleName.toString().replace('.contigs', '')
                return [[id: id, source: "contigs"], file]
            }
            .concat(Channel
                .fromPath("${params.contigs_dir}/*.contigs.fa.gz")
                .map { file -> 
                    def id = file.simpleName.toString().replace('.contigs', '')
                    return [[id: id, source: "contigs"], file]
                }
            )
    }
   
    // Update AMRFinderPlus database (shared resource)
    AMRFINDERPLUS_UPDATE()

    // Process fastq reads if provided
    fq_ch.count().set { fq_count }
    fq_count.map { count ->
            if (count > 0) {
                log.info "Found ${count} fastq files, proceeding with assembly"
                // Clone the original channel for use in processes
                def fq_clone = fq_ch.buffer(size: count)
                
                FASTQC(fq_clone)
                TRIMMOMATIC(fq_clone)
                SPADES(TRIMMOMATIC.out.trimmed_reads.map{it -> [it[0],it[1],[],[]]},[])
                
                // Return the assembled contigs channel
                return SPADES.out.contigs
            } else {
                log.warn "No fastq files in input channel"
                return []
            }
        }
        .set { contigs_assembled_ch }

    // Combine direct contigs and assembled contigs
    contigs_ch = contigs_direct_ch
                                .concat(contigs_assembled_ch).filter{it != []}
                                
    // Exit if no input is provided
    contigs_ch
        .count()
        .map{count ->
            if (count==0) {
                error "No input data provided. Please specify at least one of --reads_dir or --contigs_dir"
            }
    }

    // Run MLST to identify species for all contigs
    MLST(contigs_ch)
    
    // Process MLST output to extract species information
    species_ch = MLST.out.tsv
        .map { meta, mlst_file ->
            def species = "unknown"
            mlst_file.withReader { reader ->
                def line = reader.readLine()
                if (line) {
                    def fields = line.split('\t')
                    if (fields.size() > 2) {
                        species = fields[1].trim()
                        // Normalize species name to match potential keys
                        if (species.contains("ecoli")) species = "Escherichia_coli"
                        else if (species.contains("klebsiella")) species = "Klebsiella_pneumoniae_complex"
                        else if (species.contains("koxytoca")) species = "Klebsiella_oxytoca_complex"
                        else if (species.contains("ecloacae")) species = "Enterobacter_cloacae"
                        else if (species.contains("paeruginosa")) species = "Pseudomonas_aeruginosa"
                    }
                }
            }
            return [meta, species]
        }

    // Log the identified species for each sample
    species_ch.map { meta, species ->
        log.info "Sample ${meta.id} identified as ${species}"
        return [meta, species]
    }

    // Combine assembly with species information and get reference files
    assembly_species_ch = contigs_ch
        .join(species_ch)
        .map { meta, contigs, species ->
            def updated_meta = meta + [species: species]
            def ref_genome = speciesReferences.containsKey(species) ? 
                speciesReferences[species] : defaultReference
            
            return [updated_meta , contigs, species, ref_genome]
        }

    // Run PROKKA with appropriate training file for all contigs
    prokka_ch = assembly_species_ch.map { meta, contigs, species, ref_genome ->
        def trn_file = (ref_genome.trn && ref_genome.trn != "") ? file(ref_genome.trn) : []
        [meta, contigs, trn_file, []]
    }
    
    PROKKA(prokka_ch)
    
    amrfinder_ch = assembly_species_ch.map { meta, contigs, species, ref_genome -> [meta, ref_genome.amrfindopt, contigs]
    }
    // Run AMRFinderPlus for all contigs
    AMRFINDERPLUS_RUN(amrfinder_ch, AMRFINDERPLUS_UPDATE.out.db)

    // Run Diamond against all reference databases (e.g., VFDB)
    diamond_db_channel = channel.fromList(refDiamondFa.collect { path -> file(path) })
    DIAMOND_BLASTX(PROKKA.out.ffn.combine(diamond_db_channel))

    
    // Run rgi and CARD annotations

    def wildcardPath = file("${params.dataCacheDir}/wildcard")
    def cardPath = file("${params.database_references_dir}/CARD")

    // Check if the directory exists and is non-empty
    if (!wildcardPath.exists() || wildcardPath.list().size() == 0) {
        log.warn "Directory ${params.dataCacheDir}/wildcard does not exist or is empty. Updating the CARD database..."
        
        // Run update process and pass its outputs to RGI_MAIN
        RGI_UPDATE(cardversion)
        RGI_MAIN(PROKKA.out.faa, RGI_UPDATE.out.card, RGI_UPDATE.out.wildcard)
    } else {
        log.warn "Directory ${params.dataCacheDir}/wildcard already exists and contains files. Skipping CARD update."
        
        // Use existing paths directly with RGI_MAIN
        RGI_MAIN(PROKKA.out.faa, cardPath, wildcardPath)
    }

    // Conditional Kleborate run for Klebsiella species
    kleborate_contigs = assembly_species_ch
        .filter { meta, contigs, species, ref_genome -> 
            species.contains("Klebsiella_pneumoniae")  || species.contains("Klebsiella_oxytoca")
        }
        .map { meta, contigs, species, ref_genome ->
            if (species.contains("Klebsiella_pneumoniae")) {
            [meta, contigs, 'kpsc']
            } else if (species.contains("Klebsiella_oxytoca")) {
                [meta, contigs, 'kosc']
            } 
        }

        if (kleborate_contigs) {
            KLEBORATE(kleborate_contigs)
        }

    // Run QUAST on all contigs
    quast_ch = assembly_species_ch.map { meta, contigs, species, ref_genome ->
        def fasta_file = (ref_genome.fasta && ref_genome.fasta != "") ? file(ref_genome.fasta) : []
        [meta, contigs, fasta_file]
     }.view()

    QUAST(quast_ch)

    //Run Busco on all contigs
    
    //Busco

    // Only run SNIPPY for samples that came from fastq reads
    // if (!fq_ch.isEmpty()) {
    //     // Get species information for reads samples
    //     reads_species_ch = species_ch.filter { meta, species ->
    //         meta.source == "reads"
    //     }

    //     // Prepare input for SNIPPY
    //     snippy_input = fq_ch.join(reads_species_ch).map { meta, reads, species ->
    //         def ref_genome = speciesReferences.containsKey(species) ? 
    //             speciesReferences[species] : defaultReference
    //         return [meta, reads, file(ref_genome.gbk)]
    //     }
        
    //     SNIPPY_RUN(snippy_input)
    // }
}

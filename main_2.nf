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


// Parameters with default values that can be overridden
params.reads_dir = ""              // Directory containing fastq files
params.contigs_dir = ""            // Directory containing pre-assembled contigs
params.output_dir = ""
params.reference_dir = "$baseDir/assets/references"
params.species_config = "$baseDir/conf/species_references.config"
params.db_config = "$baseDir/conf/database_references.config"

// Load species reference configuration from external file
def speciesReferences = [:]
def defaultReference = [:]

// Function to load species references from config file
def loadSpeciesConfig() {
    def configFile = file(params.species_config)
    if (configFile.exists()) {
        def slurper = new groovy.json.JsonSlurper()
        def config = slurper.parseText(configFile.text)
        
        // Load the species configurations
        config.species.each { species ->
            speciesReferences[species.name] = [
                gbk: species.gbk,
                fasta: species.fasta,
                trn: species.trn
            ]
        }
        
        // Set default if available
        if (config.default) {
            defaultReference = [
                gbk: config.default.gbk,
                fasta: config.default.fasta,
                trn: config.default.trn
            ]
        }
    } else {
        log.warn "Species configuration file not found: ${params.species_config}"
        log.warn "Using built-in defaults for Pseudomonas aeruginosa"
        
        // Fallback default if config file not found
        defaultReference = [
            gbk: "$params.reference_dir/pseudomonas/NC_002516v2.gbk",
            fasta: "$params.reference_dir/pseudomonas/NC_002516v2.fasta",
            trn: "$params.reference_dir/pseudomonas/Pseudomonas_aeruginosa.trn"
        ]
        
        speciesReferences["Pseudomonas aeruginosa"] = defaultReference
    }
}

// Function to load database references for Diamond
def loadDatabaseConfig() {
    def databases = [:]
    def configFile = file(params.db_config)
    if (configFile.exists()) {
        def slurper = new groovy.json.JsonSlurper()
        def config = slurper.parseText(configFile.text)
        
        // Load database configurations
        config.databases.each { db ->
            databases[db.name] = db.path
        }
    } else {
        log.warn "Database configuration file not found: ${params.db_config}"
        log.warn "No Diamond databases will be used"
    }
    return databases
}

// Load configurations when pipeline starts
loadSpeciesConfig()
def diamondDatabases = loadDatabaseConfig()

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
    if (!fq_ch.isEmpty()) {
        FASTQC(fq_ch)
        TRIMMOMATIC(fq_ch)
        SPADES(TRIMMOMATIC.out.trimmed_reads.map{it -> [it[0],it[1],[],[]]},[])
        
        // Create channel with assembled contigs from reads
        contigs_assembled_ch = SPADES.out.contigs
    } else {
        contigs_assembled_ch = Channel.empty()
    }

    // Combine direct contigs and assembled contigs
    contigs_ch = contigs_direct_ch.mix(contigs_assembled_ch)

    // Exit if no input is provided
    if (contigs_ch.isEmpty()) {
        error "No input data provided. Please specify at least one of --reads_dir or --contigs_dir"
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
                        if (species.contains("Escherichia")) species = "Escherichia coli"
                        else if (species.contains("Klebsiella")) species = "Klebsiella pneumoniae"
                        else if (species.contains("Enterobacter")) species = "Enterobacter cloacae"
                        else if (species.contains("Pseudomonas")) species = "Pseudomonas aeruginosa"
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
            def ref_genome = speciesReferences.containsKey(species) ? 
                speciesReferences[species] : defaultReference
            return [meta, contigs, species, ref_genome]
        }

    // Run PROKKA with appropriate training file for all contigs
    prokka_ch = assembly_species_ch.map { meta, contigs, species, ref_genome ->
        [meta, contigs, [], file(ref_genome.trn)]
    }
    
    PROKKA(prokka_ch)

    // Run AMRFinderPlus for all contigs
    AMRFINDERPLUS_RUN(contigs_ch, AMRFINDERPLUS_UPDATE.out.db)

    // Run Diamond against all configured databases using Prokka's protein FASTA
    diamondDatabases.each { db_name, db_path ->
        DIAMOND_RUN(PROKKA.out.proteins, db_name, db_path)
    }

    // Conditional Kleborate run for Klebsiella species
    klebsiella_contigs = assembly_species_ch
        .filter { meta, contigs, species, ref_genome -> 
            species.contains("Klebsiella")
        }
        .map { meta, contigs, species, ref_genome ->
            [meta, contigs]
        }
    
    if (!klebsiella_contigs.isEmpty()) {
        KLEBORATE(klebsiella_contigs)
    }

    // Run QUAST on all contigs
    QUAST(assembly_species_ch.map { meta, contigs, species, ref_genome ->
        [meta, contigs, file(ref_genome.fasta)]
    })

    //Run Busco on all contigs
    
    //Busco

    // Only run SNIPPY for samples that came from fastq reads
    if (!fq_ch.isEmpty()) {
        // Get species information for reads samples
        reads_species_ch = species_ch.filter { meta, species ->
            meta.source == "reads"
        }

        // Prepare input for SNIPPY
        snippy_input = fq_ch.join(reads_species_ch).map { meta, reads, species ->
            def ref_genome = speciesReferences.containsKey(species) ? 
                speciesReferences[species] : defaultReference
            return [meta, reads, file(ref_genome.gbk)]
        }
        
        SNIPPY_RUN(snippy_input)
    }
}

// Diamond process definition
process DIAMOND_RUN {
    tag "${meta.id} against ${db_name}"
    publishDir "${params.output_dir}/diamond/${db_name}", mode: 'copy'

    input:
    tuple val(meta), path(proteins)
    val db_name
    path db_path

    output:
    tuple val(meta), path("${meta.id}.${db_name}.diamond.tsv"), emit: results
    path "${meta.id}.${db_name}.diamond.html", emit: html

    script:
    """
    # Run Diamond blastp against the database
    diamond blastp \\
        --query ${proteins} \\
        --db ${db_path} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \\
        --threads ${task.cpus} \\
        --out ${meta.id}.${db_name}.diamond.tsv

    # Generate HTML report
    echo "<html><head><title>${meta.id} vs ${db_name}</title></head><body>" > ${meta.id}.${db_name}.diamond.html
    echo "<h1>Diamond results for ${meta.id} against ${db_name}</h1>" >> ${meta.id}.${db_name}.diamond.html
    echo "<table border='1'><tr><th>Query</th><th>Subject</th><th>Identity %</th><th>Length</th><th>E-value</th><th>Bit Score</th><th>Title</th></tr>" >> ${meta.id}.${db_name}.diamond.html
    
    # Use awk to format top hits into HTML table
    awk -F'\\t' 'BEGIN {OFS="</td><td>"} {print "<tr><td>"$1, $2, $3, $4, $11, $12, $13"</td></tr>"}' ${meta.id}.${db_name}.diamond.tsv >> ${meta.id}.${db_name}.diamond.html
    
    echo "</table></body></html>" >> ${meta.id}.${db_name}.diamond.html
    """
}

// Kleborate process definition
process KLEBORATE {
    tag "${meta.id}"
    publishDir "${params.output_dir}/kleborate", mode: 'copy'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.kleborate.txt"), emit: results
    path "${meta.id}.kleborate_summary.txt", optional: true, emit: summary

    script:
    """
    # Run Kleborate
    kleborate \\
        --outfile ${meta.id}.kleborate.txt \\
        --assemblies ${contigs} \\
        --all

    # Create summary
    head -n1 ${meta.id}.kleborate.txt > ${meta.id}.kleborate_summary.txt
    grep -v "^#" ${meta.id}.kleborate.txt >> ${meta.id}.kleborate_summary.txt
    """
}
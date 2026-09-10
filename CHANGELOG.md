# Changelog

## [Unreleased] — planned for v1.1.0

- Automated test coverage (e.g. `nf-test`) for the major pipeline branches: reads-only, contigs-only, and the Klebsiella-specific path — currently the pipeline has no regression protection beyond manual end-to-end runs.
- Validate the remaining curated species (currently only *E. coli* and *K. pneumoniae* have been run end-to-end; *Klebsiella oxytoca/variicola*, *Citrobacter freundii*, the *Enterobacter* species, *Pseudomonas aeruginosa*, *Proteus mirabilis*, *Providencia rettgeri*, and *Staphylococcus aureus* have not) to catch species-specific reference/config issues like the Staph aureus `.trn` casing bug already found and fixed.
- A standalone, on-demand cgMLST typing pipeline (chewBBACA allele calling + distance tree), reading a contigs directory directly — kept separate from `main.nf` by design, run whenever a user wants it rather than as part of every sample's processing.
- Migrate `main.nf`'s `include` statements off dynamic (`$baseDir`) string interpolation so the pipeline runs on current Nextflow (25.x+) instead of requiring the `23.10.x`–`24.x` pin.
- Give `GENE_DIFF`, `SAVE_TO_DB`, and `GENERATE_REPORT` a `container` directive so they work under container-only profiles without also needing `conda.enabled = true`.
- Wire reference proteins (`.faa`) into `species_references.config` for the reference-vs-reference BLASTx step (`main.nf:348`).

## [1.0.0] — 2026-09-10

Initial tagged release. Nextflow (DSL2) pipeline for bacterial WGS analysis from Illumina paired-end reads: assembly, species identification, AMR/virulence gene detection, plasmid and MLST typing, variant calling, and aggregation into a SQLite database and interactive HTML report.

Verified end-to-end on real public data (one *E. coli* and one *K. pneumoniae* ENA read set) under Singularity + Conda on a single workstation, including:

- Assembly (SPAdes) and QC (FastQC, Trimmomatic, QUAST, BUSCO)
- Taxonomic classification (Kraken2/Bracken) and species typing (MLST)
- Annotation (Prokka) and AMR/virulence detection (AMRFinderPlus, RGI/CARD, Diamond BLASTX vs VFDB)
- Species-specific typing (Kleborate for *Klebsiella*, PlasmidFinder for Enterobacterales)
- Variant calling (Snippy, from contigs and reads) and gene-level mutation comparison (GeneDiff)
- Aggregation into a per-run SQLite database and a self-contained interactive HTML report

Notable fixes made getting here (see git history for detail): asset path resolution, Trimmomatic adapter path and quality-encoding detection, config scoping bugs, a deprecated Biopython API in the `genediff` submodule, and a BUSCO/Nextflow path-staging conflict.

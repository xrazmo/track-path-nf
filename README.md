<img src="https://cdn.jsdelivr.net/gh/xrazmo/ui-util/track-path-ui/img/ico/trackpath.png" width="80" alt="TRACK-PATH logo">

# TRACK-PATH-NF

A Nextflow (DSL2) pipeline for analyzing bacterial whole-genome sequencing data from Illumina paired-end reads — assembly, species identification, AMR/virulence gene detection, plasmid and MLST typing, variant calling, and aggregation into a browsable SQLite database and an interactive HTML report.

**For research use only.** This pipeline is not validated or intended for clinical diagnosis, treatment decisions, or any other clinical use.

## Pipeline overview

```
reads (fastq) ──► FastQC ──► Trimmomatic ──► SPAdes ──┐
                                                       ├──► contigs
contigs (pre-assembled) ─────────────────────────────┘
                                                       │
reads ──► Kraken2 ──► Bracken (taxonomic classification)
                                                       │
contigs ──► MLST ──► species identification ──► per-species reference lookup
                                                       │
        ┌──────────────────────────────────────────────┼───────────────────────────────┐
        ▼                     ▼                ▼        ▼            ▼        ▼         ▼
     Prokka             AMRFinderPlus         RGI     Kleborate*   QUAST   BUSCO   PlasmidFinder**
   (annotation)          (AMR genes)       (CARD DB)  (Klebsiella)  (assembly QC)        │
        │                                                                                 │
        ├──► Diamond BLASTX (VFDB + reference proteins)                                   │
        ├──► GeneDiff (gene-level mutation diff vs. reference)                            │
        └──► Snippy (variant calling, from contigs and/or reads)                          │
                                                                                            │
                              (all of the above, for every sample) ───────────────────────┘
                                                       │
                                                       ▼
                                    SAVE_TO_DB (SQLite)   GENERATE_REPORT (interactive HTML)
```

\* Kleborate only runs for *Klebsiella* species.
\*\* PlasmidFinder only runs for species flagged `plasmidAct: true` in the species config (Enterobacterales).

Species-specific behavior (which reference genome/proteins/training file to use, which AMRFinderPlus organism flag to pass, whether PlasmidFinder applies) is driven entirely by [`assets/references/species_references.config`](assets/references/species_references.config) — species identity comes from MLST on the assembly, or can be overridden per-sample via `--species_csv`.

Currently curated species: *Escherichia coli*, *Klebsiella pneumoniae* complex, *Klebsiella oxytoca* complex, *Klebsiella variicola*, *Citrobacter freundii*, *Enterobacter cloacae/hormaechei/bugandensis/asburiae*, *Pseudomonas aeruginosa*, *Proteus mirabilis*, *Providencia rettgeri*, *Staphylococcus aureus*. Any other/unrecognized species falls back to a reference-free default (skips reference-dependent steps like Snippy/QUAST-vs-reference, still runs reference-free tools).

## Requirements

- [Nextflow](https://www.nextflow.io/) — pinned to the `23.10.x`–`24.x` range (see [Known issues](#known-issues) below)
- Java 17+ (required by Nextflow itself)
- One of: [Singularity](https://sylabs.io/singularity/) / [Apptainer](https://apptainer.org/), Docker, or Conda/[Mamba](https://mamba.readthedocs.io/) — most modules ship a `container` directive; a couple (`GENE_DIFF`, `SAVE_TO_DB`, `GENERATE_REPORT`) are conda-only and need `conda.enabled = true`
- ~450MB for the curated reference assets (already in `assets/references/`), plus a few GB of scratch space for auto-downloaded databases (Kraken2, CARD, PlasmidFinder, BUSCO lineages) on first run

## Installation

1. **Install Java 17+** (required by Nextflow):
   ```bash
   sudo apt-get install openjdk-17-jre-headless   # or any Java 17+ distribution
   ```

2. **Install Nextflow, pinned to a compatible version.** This pipeline's `include` statements use a dynamic-path syntax that newer Nextflow (25.x+) rejects, so pin the runtime with `NXF_VER` rather than installing whatever is latest:
   ```bash
   curl -s https://get.nextflow.io | bash   # downloads/bootstraps the `nextflow` launcher
   mv nextflow ~/.local/bin/                # or anywhere on your PATH
   export NXF_VER=24.10.5                   # add to your shell profile to make this permanent
   ```
   Verify with `nextflow -version` — it should report `24.10.5` (or another `23.10.x`–`24.x` release).

3. **Install a container/environment runtime.** Most modules ship a `container` directive (Singularity/Apptainer or Docker); a few conda-only modules (`GENE_DIFF`, `SAVE_TO_DB`, `GENERATE_REPORT`) need Conda/[Mamba](https://mamba.readthedocs.io/) as well:
   - Singularity/Apptainer: see [sylabs.io](https://sylabs.io/singularity/) / [apptainer.org](https://apptainer.org/)
   - Docker: see [docs.docker.com](https://docs.docker.com/get-docker/)
   - Conda/Mamba: see [Miniforge](https://github.com/conda-forge/miniforge)

4. **Clone the repository with its submodule** (the `genediff` gene-mutation-comparison tool lives in a separate repo):
   ```bash
   git clone --recurse-submodules https://github.com/xrazmo/track-path-nf.git
   cd track-path-nf
   ```
   If you already cloned without `--recurse-submodules`:
   ```bash
   git submodule update --init --recursive
   ```

5. **Verify the setup** with a dry run (`-preview` validates the pipeline DAG and config without actually executing anything, so a nonexistent `--reads_dir` is fine here):
   ```bash
   nextflow run main.nf -c main.config -c local_test.config -profile local \
     --reads_dir ./any_dir --output_dir ./results --run_assembly true \
     -preview
   ```
   A clean exit with a printed process list (no `ERROR`) means Nextflow, the config, and the submodule are all wired up correctly. The reference assets under `assets/references/` ship with the repo; databases like Kraken2/CARD/PlasmidFinder/BUSCO lineages download automatically into `dataCacheDir` on first real run.

## Quick start

```bash
./run.sh --reads_dir /path/to/fastq_dir --output_dir ./results
```

`run.sh` is an interactive wrapper: it validates inputs, generates a run `--ticket`, prints a confirmation summary, and invokes `nextflow run`. See `./run.sh --help` for all flags (`-q/--reads_dir`, `-g/--contigs_dir`, `-s/--species_csv`, `-o/--output`, `-w/--work`, `-c/--config`, `-p/--profile`, `-r/--resume`, `-f/--free_param`).

### Running directly with Nextflow

```bash
nextflow run main.nf -c main.config \
  --reads_dir /path/to/fastq_dir \
  --output_dir ./results \
  --run_assembly true
```

For local development (no SLURM/HPC, running Singularity+conda directly on a workstation), layer `local_test.config` on top, which overrides the HPC-oriented cache/asset paths in `main.config` and caps per-process resources to fit a single machine:

```bash
nextflow run main.nf -c main.config -c local_test.config -profile local \
  --reads_dir ./test_data/reads \
  --output_dir ./test_data/results \
  --run_assembly true \
  -w ./test_data/work
```

Launch from inside the working directory you pass via `-w` (or pass an absolute `-w` path) so Nextflow's own `.nextflow.log`/`.nextflow/` cache land there instead of cluttering the repo root.

### Input layout

- `--reads_dir`: a directory of paired-end fastq files, either flat (`SAMPLE_1.fastq.gz` / `SAMPLE_2.fastq.gz`, also accepts `_R1/_R2` and `.fq.gz`) or one subdirectory per sample.
- `--contigs_dir`: a directory of pre-assembled contigs (`SAMPLE.contigs.fa.gz`), for skipping assembly entirely. Can be combined with `--reads_dir`.
- `--species_csv`: optional CSV with `sample_id,species` columns to override MLST-based species identification (species name must match an entry in `species_references.config`).
- `--run_assembly true`: required to actually run FastQC/Trimmomatic/SPAdes/Kraken2/Bracken/Snippy-from-reads. If only `--contigs_dir` is given, leave this `false`.

### Key parameters

| Parameter | Default | Purpose |
|---|---|---|
| `--reads_dir` / `--contigs_dir` | `""` | Input data (at least one required) |
| `--output_dir` | `""` | Where results are published |
| `--run_assembly` | `false` | Run the read-based branch (QC, assembly, Kraken2/Bracken, Snippy-from-reads) |
| `--species_csv` | `""` | Manual sample→species override |
| `--ticket` | auto-generated | Run identifier, used to name the output database |
| `--db_name` | `<ticket>.trackpath_results.db` | SQLite output filename |
| `--save_db_add_seq` | `false` | Also store predicted ORF nucleotide/protein sequences in the database |
| `--skip_save_db` / `--skip_report` | `false` | Skip the corresponding aggregation step |

Asset/database locations (`assetsDir`, `dataCacheDir`, `species_config`, `db_config`) are set in `main.config`/`local_test.config`, not typically overridden on the command line.

## Outputs

Everything is published under `--output_dir`, one subdirectory per tool (`prokka/`, `mlst/`, `amrfinder/`, `rgi/`, `quast/`, `busco/`, `snippy/`, `snippy_contig/`, `plasmidfinder/`, `kleborate/`, `diamond/`, `genediff/`, `kraken2_bracken/`, `fastqc/`, ...), plus two aggregated outputs produced after every sample finishes:

- **`<ticket>.trackpath_results.db`** — a SQLite database with one table per tool (`species_typing`, `taxonomic_classification`, `genome_annotation`, `amr_genes`, `card_resistance_genes`, `virulence_factors`, `reference_protein_blast`, `assembly_quality`, `plasmid_typing`, `klebsiella_typing`, `variant_calls`, `variant_calls_summary`, `gene_mutation_diff`), deduplicated by content hash so re-running/resuming doesn't create duplicate rows.
- **`report/`** — a self-contained interactive HTML viewer (`report/index.html`, reading `report/js/data.js` and `report/img/`) summarizing QC, assembly quality, species/MLST calls, AMR/virulence findings, and plasmid typing per sample. Open `report/index.html` directly in a browser.

## Downstream (outside Nextflow)

Two SLURM batch scripts cover analyses that run after the main pipeline, across a cohort of samples rather than per-sample:

- **`run_cgmlst.sh`** — chewBBACA-based cgMLST allele calling and a distance-based Newick tree, from a directory of contigs.
- **`run_save_db.sh`** — a SLURM wrapper around `bin/save_to_db.py` for saving/merging results into a shared database outside of a single Nextflow run (e.g., backfilling historical runs).

## Repository layout

```
main.nf              Pipeline entrypoint (DSL2 workflow)
main.config           Default (HPC/SLURM-oriented) configuration
local_test.config      Local single-workstation override profile
modules/              One Nextflow module per tool (main.nf + optional environment.yml)
bin/                   Python/Perl helper scripts invoked by modules (report/DB generation, reference DB building, VCF parsing, ...)
assets/references/    Curated per-species reference genomes (fasta/gbk/faa/fna/trn) + species_references.config
assets/databases/      Bulk-downloaded databases (CARD, VFDB) + database_references.config — CARD/VFDB content is gitignored and regenerated via bin/download_ref_db.py
assets/ui/             Interactive HTML report template, copied into <output_dir>/report/ by GENERATE_REPORT
submodules/genediff/   Gene-level mutation comparison tool (git submodule, github.com/xrazmo/genediff)
run.sh                 Interactive pipeline launcher
run_cgmlst.sh / run_save_db.sh   SLURM scripts for downstream cgMLST tree-building / DB saving
```

Modules present but not currently wired into `main.nf`: `roary` (pan-genome analysis), `snpeff` (variant annotation), `utility` (GenBank protein extraction).

## Known issues

- **Nextflow version**: `main.nf`'s `include` statements use dynamic (`$baseDir`) string interpolation, which newer Nextflow (25.x+) rejects under its stricter DSL2 parsing. `main.config` pins `manifest.nextflowVersion` to `23.10.0`–`24.x`; use `NXF_VER=24.10.5 nextflow run ...` (or an equivalent version manager) if your default Nextflow is newer.
- **Conda-only modules under Singularity profiles**: `GENE_DIFF`, `SAVE_TO_DB`, and `GENERATE_REPORT` declare a `conda` environment but no `container`, so a Singularity-only profile must also set `conda.enabled = true` (see `local_test.config`) or they'll run against the bare host Python.
- No automated test suite yet, and only *E. coli* and *K. pneumoniae* have been validated end-to-end so far — see [CHANGELOG.md](CHANGELOG.md) for the full list of planned follow-up work.

## License

[MIT](LICENSE)

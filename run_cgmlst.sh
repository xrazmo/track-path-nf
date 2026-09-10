#!/bin/bash

#SBATCH --job-name=chewie_job
#SBATCH --output=chewie_%j.out.log
#SBATCH --error=chewie_%j.err.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# Function to display usage
usage() {
    echo "Usage: $0 --input_dir|-i <input_dir> --schema_dir|-s <schema_dir> --trn|-t <trn> --contig_dir|-c <contig_dir> --schema_name|-n <schema_name> --output_dir|-o <output_dir> [--skip_build]"
    exit 1
}
SLURM_CPUS_PER_TASK=8
TREE_MAKER_PY="~/source/workflows/track-path-nf/bin/cgmlst_tree.py"
# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir|-i)
            INPUT_DIR="$2"
            shift 2
            ;;
        --schema_dir|-s)
            SCHEMA_DIR="$2"
            shift 2
            ;;
        --trn|-t)
            TRN="$2"
            shift 2
            ;;
        --contig_dir|-c)
            CONTIG_DIR="$2"
            shift 2
            ;;
        --schema_name|-n)
            SCHEMA_NAME="$2"
            shift 2
            ;;
        --output_dir|-o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --skip_build)
            SKIP_BUILD=true
            shift
            ;;
        *)
            usage
            ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$INPUT_DIR" ] || [ -z "$SCHEMA_DIR" ] || [ -z "$TRN" ] || [ -z "$CONTIG_DIR" ] || [ -z "$SCHEMA_NAME" ] || [ -z "$OUTPUT_DIR" ]; then
    usage
fi
# Echo input parameters for review
echo "Input Parameters:"
echo "-----------------"
echo "Input Directory: $INPUT_DIR"
echo "Schema Directory: $SCHEMA_DIR"
echo "TRN File: $TRN"
echo "Contig Directory: $CONTIG_DIR"
echo "Schema Name: $SCHEMA_NAME"
echo "Output Directory: $OUTPUT_DIR"
echo "Skip PrepExternalSchema: ${SKIP_BUILD:-false}"
echo "-----------------"

# Load modules


# Activate conda environment
source activate chewie

# Step 1: PrepExternalSchema
if [ "${SKIP_BUILD:-false}" != "true" ]; then
echo "Running PrepExternalSchema with:"
echo "Input Directory: $INPUT_DIR"
echo "Output Directory: $SCHEMA_DIR/$SCHEMA_NAME"
echo "PTF File: $TRN"
echo "CPU cores: $SLURM_CPUS_PER_TASK"
chewBBACA.py PrepExternalSchema -g "$INPUT_DIR" -o "$SCHEMA_DIR/$SCHEMA_NAME" --ptf "$TRN" --cpu $SLURM_CPUS_PER_TASK
fi

# Step 2: Unzip contigs in contig_dir
echo "Unzipping contig files in $CONTIG_DIR..."
    for file in "$CONTIG_DIR"/*.gz; do
        if [ -f "$file" ]; then
            gunzip "$file"
        else
            echo "No .gz files found in $CONTIG_DIR"
            break
        fi
done

# Step 3: AlleleCall
chewBBACA.py AlleleCall -i "$CONTIG_DIR" -g "$SCHEMA_DIR/$SCHEMA_NAME" -o "$OUTPUT_DIR" --cpu $SLURM_CPUS_PER_TASK

# Step 4: cgmlst-dists
# Find results_alleles.tsv in any subdirectory of OUTPUT_DIR
CONTIGS_INFO=$(find "$OUTPUT_DIR" -type f -name "results_alleles.tsv" | head -n 1)
if [ -z "$CONTIGS_INFO" ]; then
    echo "Error: results_alleles.tsv not found in any subdirectory of $OUTPUT_DIR"
    exit 1
fi
echo "Found results_alleles.tsv: $CONTIGS_INFO"
# Get the directory of CONTIGS_INFO
CONTIGS_DIR=$(dirname "$CONTIGS_INFO")
# Run cgmlst-dists and output distances.tsv to the same directory as CONTIGS_INFO
cgmlst-dists "$CONTIGS_INFO" > "$CONTIGS_DIR/distances_cgmlstdist.tsv"
echo "Generated distances.tsv in: $CONTIGS_DIR/distances_cgmlstdist.tsv"


# Step 5: Make the cgMLST tree in newick formate
echo "Generating cgMLST tree from $CONTIGS_DIR/distances_cgmlstdist.tsv..."
eval "python $TREE_MAKER_PY -i \"$CONTIGS_DIR/distances_cgmlstdist.tsv\" -o \"$CONTIGS_DIR/cgmlst_tree.nwk\""
if [ $? -eq 0 ]; then
    echo "Newick tree saved as $CONTIGS_DIR/cgmlst_tree.nwk"
else
    echo "Error: Failed to generate cgMLST tree"
    exit 1
fi


# Deactivate conda environment
conda deactivate

echo "Job completed"
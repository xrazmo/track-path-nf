#!/bin/bash

#SBATCH --job-name=SaveInDB_TRACKPATH
#SBATCH --output=saveindb_%j.out
#SBATCH --error=saveindb_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=01:30:00

usage() {
    echo "Usage: $0 --input_dir|-i <input_dir> --output_dir|-o <output_dir> --db_name|-d <db_name> --add_seq|-a <add_seq>"
    exit 1
}
SLURM_CPUS_PER_TASK=10
SAVE_DB_PY="~/source/workflows/track-path-nf/bin/save_to_db.py"
# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir|-i)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output_dir|-o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --db_name|-d)
            DB_NAME="$2"
            shift 2
            ;;
        --add_seq|-a)
            ADD_SEQ="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done
# Check if all required arguments are provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    usage
fi

# Convert to absolute paths
INPUT_DIR="$(cd "$(dirname "$INPUT_DIR")" && pwd)/$(basename "$INPUT_DIR")"
OUTPUT_DIR="$(cd "$(dirname "$OUTPUT_DIR")" && pwd)/$(basename "$OUTPUT_DIR")"

# Echo input parameters for review
echo "Input Parameters:"
echo "-----------------"
echo "Input Directory: $INPUT_DIR"
echo "Output Directory: $OUTPUT_DIR"
if [ -n "$DB_NAME" ]; then
    echo "DB Name: $DB_NAME"
fi
if [ -n "$ADD_SEQ" ]; then
    echo "Add Seq: $ADD_SEQ"
fi
echo "-----------------"

source activate common_env

CMD="python $SAVE_DB_PY --input_dir \"$INPUT_DIR\" --output_dir \"$OUTPUT_DIR\""
if [ -n "$DB_NAME" ]; then
    CMD="$CMD --db_name \"$DB_NAME\""
fi
if [ -n "$ADD_SEQ" ]; then
    CMD="$CMD --add_seq \"$ADD_SEQ\""
fi
CMD="$CMD --cpus \"$SLURM_CPUS_PER_TASK\""
eval "$CMD"
 
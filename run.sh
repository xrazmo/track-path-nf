#!/bin/bash

# Nextflow Pipeline Execution Script

# Halt on any error
set -e

# Determine the absolute path of the script
SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

# Color codes for formatting
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to get absolute path
get_absolute_path() {
    local path="$1"
    # If path is relative, convert to absolute
    if [[ "$path" != /* ]]; then
        path="$(readlink -f "$path")"
    fi
    echo "$path"
}

# Function to display usage information
usage() {
    echo -e "${YELLOW}Usage:${NC} $0 [options]"
    echo "Options:"
    echo "  -i, --input         Input directory containing fastq.gz files (REQUIRED)"
    echo "  -o, --output        Output directory for results (default: ./results)"
    echo "  -w, --work          Nextflow working directory (default: ./work)"
    echo "  -c, --config        Custom Nextflow configuration file (default: nextflow.config)"
    echo "  -p, --profile       Execution profile (default: standard)"
    echo "  -r, --resume        Resume previous run"
    echo "  -h, --help          Display this help message"
    exit 1
}

# Default values with absolute paths
INPUT_DIR=""
OUTPUT_DIR="$(get_absolute_path "./results")"
WORK_DIR="$(get_absolute_path "./work")"
CONFIG_FILE="$(get_absolute_path "${SCRIPT_DIR}/nextflow.config")"
MAIN_NF_FILE="$(get_absolute_path "${SCRIPT_DIR}/main.nf")"
PROFILE="standard"
RESUME_FLAG=""

# Parse command-line arguments
ARGS=$(getopt -o i:o:w:c:p:rh --long input:,output:,work:,config:,profile:,resume,help -n "$0" -- "$@")

# Check for invalid arguments
if [ $? -ne 0 ]; then
    usage
fi

eval set -- "$ARGS"

# Extract options and their arguments
while true; do
    case "$1" in
        -i|--input)
            INPUT_DIR="$(get_absolute_path "$2")"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$(get_absolute_path "$2")"
            shift 2
            ;;
        -w|--work)
            WORK_DIR="$(get_absolute_path "$2")"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$(get_absolute_path "$2")"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -r|--resume)
            RESUME_FLAG="-resume"
            shift
            ;;
        -h|--help)
            usage
            ;;
        --)
            shift
            break
            ;;
        *)
            echo -e "${RED}Error: Invalid argument${NC}"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$INPUT_DIR" ]; then
    echo -e "${RED}Error: Input directory is required${NC}"
    usage
fi

# Ensure directories exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORK_DIR"

# Pre-run checks
echo -e "${GREEN}Running Nextflow Pipeline Validation${NC}"

# Check input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo -e "${RED}Error: Input directory does not exist: $INPUT_DIR${NC}"
    exit 1
fi

# Check main.nf file
if [ ! -f "$MAIN_NF_FILE" ]; then
    echo -e "${RED}Error: Main Nextflow script not found: $MAIN_NF_FILE${NC}"
    exit 1
fi

# Check config file
if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${RED}Error: Configuration file not found: $CONFIG_FILE${NC}"
    exit 1
fi

# Check Nextflow installation
if ! command -v nextflow &> /dev/null; then
    echo -e "${RED}Error: Nextflow is not installed${NC}"
    exit 1
fi

# Display Run Configuration
echo -e "${YELLOW}Pipeline Configuration:${NC}"
echo -e "Input Directory:       ${GREEN}$INPUT_DIR${NC}"
echo -e "Output Directory:      ${GREEN}$OUTPUT_DIR${NC}"
echo -e "Work Directory:        ${GREEN}$WORK_DIR${NC}"
echo -e "Main NF Script:        ${GREEN}$MAIN_NF_FILE${NC}"
echo -e "Configuration File:    ${GREEN}$CONFIG_FILE${NC}"
echo -e "Execution Profile:     ${GREEN}$PROFILE${NC}"
echo -e "Resume Flag:           ${GREEN}${RESUME_FLAG:-None}${NC}"

# Confirmation Prompt
read -p "Do you want to proceed with the pipeline? (y/n): " confirm
if [[ $confirm != [yY] && $confirm != [yY][eE][sS] ]]; then
    echo -e "${RED}Pipeline execution cancelled.${NC}"
    exit 0
fi

# Run Nextflow Pipeline
echo -e "${GREEN}Starting Nextflow Pipeline...${NC}"
cd "$WORK_DIR"
nextflow run "$MAIN_NF_FILE" \
    -profile "$PROFILE" \
    -c "$CONFIG_FILE" \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR" \
    -w "$WORK_DIR" \
    $RESUME_FLAG \
    -process.log

# Pipeline Completion Status
PIPELINE_EXIT_CODE=$?

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}Pipeline completed successfully!${NC}"
else
    echo -e "${RED}Pipeline failed with exit code $PIPELINE_EXIT_CODE${NC}"
fi

# Optional: Generate summary or perform post-run tasks
echo -e "${YELLOW}Pipeline Run Summary:${NC}"
echo "- Input Directory: $INPUT_DIR"
echo "- Output Directory: $OUTPUT_DIR"
echo "- Execution Profile: $PROFILE"
echo "- Completion Time: $(date)"

exit $PIPELINE_EXIT_CODE
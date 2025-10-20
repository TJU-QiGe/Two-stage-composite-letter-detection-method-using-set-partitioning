#!/bin/bash
#-------------------------------------------------------------------------------
# File: full_readout_pipeline.sh
# Author: [Qi Ge]
# Version: [2.0]
# Organization: [School of Microelectronics, Tianjin University]
# Date: [2025/10/15]
#
# Description:
# Complete DNA storage data readout pipeline consisting of 5 major steps:
# 1. Random downsampling of sequencing reads
# 2. Primer identification and removal  
# 3. Index-based read grouping
# 4. Composite letter detection and consensus calling
# 5. Reed-Solomon decoding and data recovery
#
# Usage: bash full_readout_pipeline.sh
#-------------------------------------------------------------------------------

# Set unlimited stack size to prevent overflow
ulimit -s unlimited

# Configuration Parameters
COVERAGE=20
INPUT_FASTQ="../Demo_FASTQ_eight_letters/Example-paired-end-assembled_100x.fastq"
OUTPUT_BASE="../Eight_letter_result"
FORWARD_PRIMER="TCACCATCCACTCTAAACAC"
REVERSE_PRIMER="ATGAGTGGAGGTGTAAAGTG"
BASE_TYPES="ATGCRYMK"
ALPHABET=1
THRESHOLD=2

echo "=================================================================================="
echo "DNA Storage Complete Readout Pipeline"
echo "Version: 2.0 | Author: Qi Ge | Organization: Tianjin University"
echo "=================================================================================="
echo "Configuration:"
echo "  Input Coverage: ${COVERAGE}x"
echo "  Input FASTQ: ${INPUT_FASTQ}"
echo "  Output Directory: ${OUTPUT_BASE}"
echo "=================================================================================="

#-------------------------------------------------------------------------------
# PROGRAM COMPILATION
#-------------------------------------------------------------------------------
echo "[COMPILE] Compiling all required programs..."

echo "[INFO] Compiling primer identification program..."
g++ -o primer_identification primer_identification.c edlib.cpp -lm
if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to compile primer_identification"
    exit 1
fi

echo "[INFO] Compiling read grouping program..."
gcc -o groupingByIndex groupingByIndex.c -lm
if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to compile groupingByIndex"
    exit 1
fi

echo "[INFO] Compiling letter detection core..."
gcc -o letter_detection_core letter_detection_core.c -lm
if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to compile letter_detection_core"
    exit 1
fi

echo "[INFO] Compiling data analysis program..."
gcc -o data_analysis data_analysis.c -lm
if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to compile data_analysis"
    exit 1
fi

echo "[INFO] Compiling Reed-Solomon decoder..."
g++ -std=c++11 RS_decoder.cpp -o RS_decoder -lm
if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to compile RS_decoder"
    exit 1
fi

echo "[INFO] Compiling data recovery program..."
gcc recovery_poem.c -o recovery_poem
if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to compile recovery_poem"
    exit 1
fi

echo "[SUCCESS] All programs compiled successfully"
echo ""

#-------------------------------------------------------------------------------
# STEP 1: RANDOM DOWNSAMPLING
#-------------------------------------------------------------------------------
echo "[STEP 1] Starting random downsampling of sequencing reads..."

# Setup directories
FASTQ_DIR="${OUTPUT_BASE}/${COVERAGE}/downsampling"
mkdir -p "${FASTQ_DIR}"

# Calculate required number of reads
NUM_READS=$((COVERAGE * 126))
OUTPUT_FASTQ="${FASTQ_DIR}/sub_sample.fastq"

# Perform random sampling with seqtk
echo "[INFO] Downsampling from 100x to ${COVERAGE}x coverage (${NUM_READS} reads)"
./seqtk sample -s100 "${INPUT_FASTQ}" "${NUM_READS}" > "${OUTPUT_FASTQ}"

echo "[SUCCESS] Step 1 completed - Downsampled reads saved to: ${OUTPUT_FASTQ}"
echo ""

#-------------------------------------------------------------------------------  
# STEP 2: PRIMER IDENTIFICATION AND REMOVAL
#-------------------------------------------------------------------------------
echo "[STEP 2] Starting primer identification and removal..."

# Setup output paths
OUTPUT_READS="${FASTQ_DIR}/read_without_primer.fasta"

# Run primer identification
echo "[INFO] Removing primers from reads..."
echo "  Forward primer: ${FORWARD_PRIMER}"
echo "  Reverse primer: ${REVERSE_PRIMER}"
./primer_identification "${OUTPUT_FASTQ}" "${OUTPUT_READS}" "${FORWARD_PRIMER}" "${REVERSE_PRIMER}"

echo "[SUCCESS] Step 2 completed - Primer-free reads saved to: ${OUTPUT_READS}"
echo ""

#-------------------------------------------------------------------------------
# STEP 3: INDEX-BASED READ GROUPING  
#-------------------------------------------------------------------------------
echo "[STEP 3] Starting index-based read grouping..."

# Setup input files and output directory
FRONT_INDEX="./configureFiles/index_8nt.txt"
BACK_INDEX="./configureFiles/index_8nt_interleaved.txt"
GROUP_DIR="${OUTPUT_BASE}/${COVERAGE}/group_reads"
mkdir -p "${GROUP_DIR}"

# Define output files
PAYLOAD_FILE="${GROUP_DIR}/payload.txt"
MATCH_INFO="${GROUP_DIR}/matched_information_record.txt"  
FREQ_FILE="${GROUP_DIR}/addrFrequency.txt"

# Run grouping algorithm
echo "[INFO] Grouping reads by double-end indices (threshold: ${THRESHOLD})"
./groupingByIndex "${OUTPUT_READS}" "${FRONT_INDEX}" "${BACK_INDEX}" "${FREQ_FILE}" "${MATCH_INFO}" "${PAYLOAD_FILE}" "${THRESHOLD}"

echo "[SUCCESS] Step 3 completed - Grouped reads saved to: ${PAYLOAD_FILE}"
echo ""

#-------------------------------------------------------------------------------
# STEP 4: COMPOSITE LETTER DETECTION
#-------------------------------------------------------------------------------
echo "[STEP 4] Starting composite letter detection..."

# Setup paths
RESULT_DIR="${OUTPUT_BASE}/${COVERAGE}/readout_result"
mkdir -p "${RESULT_DIR}"
CONSENSUS_FILE="${RESULT_DIR}/consensus.txt"
REFERENCE_FILE="./configureFiles/composite_payload_ref.txt"

# Run letter detection
echo "[INFO] Running two-stage letter detection algorithm"
echo "  Base types: ${BASE_TYPES}"
echo "  Alphabet complexity: ${ALPHABET}"
./letter_detection_core "${PAYLOAD_FILE}" "${CONSENSUS_FILE}" "${ALPHABET}"

# Run accuracy analysis
echo "[INFO] Analyzing detection accuracy..."
./data_analysis "${CONSENSUS_FILE}" "${PAYLOAD_FILE}" "${REFERENCE_FILE}" "${BASE_TYPES}" "${RESULT_DIR}" "${ALPHABET}"

echo "[SUCCESS] Step 4 completed - Consensus sequences saved to: ${CONSENSUS_FILE}"
echo ""

#-------------------------------------------------------------------------------
# STEP 5: REED-SOLOMON DECODING AND DATA RECOVERY
#-------------------------------------------------------------------------------
echo "[STEP 5] Starting Reed-Solomon decoding and data recovery..."

# Setup output files
DECODED_BITS="${RESULT_DIR}/decoded_information_bits.txt"
RECOVERED_POEM="${RESULT_DIR}/Tang_poem.txt"

# Run Reed-Solomon decoding
echo "[INFO] Performing Reed-Solomon error correction..."
./RS_decoder "${CONSENSUS_FILE}" "${DECODED_BITS}" "${ALPHABET}"

# Recover original data
echo "[INFO] Converting decoded bits to original text..."
./recovery_poem "${DECODED_BITS}" "${RECOVERED_POEM}"

echo "[SUCCESS] Step 5 completed - Recovered text saved to: ${RECOVERED_POEM}"
echo ""

#-------------------------------------------------------------------------------
# PIPELINE COMPLETION
#-------------------------------------------------------------------------------
echo "=================================================================================="
echo "Pipeline Execution Summary:"
echo "  Coverage Level: ${COVERAGE}x"
echo "  Total Steps: 5"
echo "  Final Output: ${RECOVERED_POEM}"
echo ""
echo "Generated Files:"
echo "  1. Downsampled reads: ${OUTPUT_FASTQ}"
echo "  2. Primer-free reads: ${OUTPUT_READS}"  
echo "  3. Grouped payload: ${PAYLOAD_FILE}"
echo "  4. Consensus sequence: ${CONSENSUS_FILE}"
echo "  5. Recovered text: ${RECOVERED_POEM}"
echo ""
echo "[SUCCESS] DNA storage readout pipeline completed successfully!"
echo "=================================================================================="

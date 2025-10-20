#!/bin/bash
# DNA Storage Encoding Pipeline - Build and Run Script
# Author: Qi Ge
# Institution: Tianjin University
# Date: 2025/10/15

echo "DNA Storage Encoding Pipeline"
echo "============================="
# Create output directory
output_path=../Encoding_result
if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
    echo "Created directory: $output_path"
fi

# Step 1: Data to Bitstream
echo ""
echo "Step 1: Converting data to bitstream..."
g++ -std=c++11 txt2bits.cpp -o txt2bits -lm
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile data_to_bitstream"
    exit 1
fi
Input_File="../User_Data/Tang_Poems.txt"
Output_File="$output_path/Bitstream_of_userData.txt"

if [ ! -f "$Input_File" ]; then
    echo "Error: Input file not found: $Input_File"
    exit 1
fi

./txt2bits $Input_File $Output_File 18900
echo "Bitstream conversion completed: $Output_File"

# Step 2: DNA Encoding
echo "Step 2: DNA encoding with error correction..."

g++ -o RS_1575_1890_encoder RS_1575_1890_encoder.cpp -std=c++11 -lm
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile DNA encoder"
    exit 1
fi

SOURCE_FILE="$output_path/Bitstream_of_userData.txt"
OUTPUT_DNA_PATH="$output_path/Encoded_strands"
scramble_sequence="./configureFiles/watermark_18900bit.txt"
Left_index_seq="./configureFiles/index_8nt.txt"
Right_index_seq="./configureFiles/index_8nt_interleaved.txt"


if [ ! -d "$OUTPUT_DNA_PATH" ]; then
    mkdir -p "$OUTPUT_DNA_PATH"
    echo "Created directory: $OUTPUT_DNA_PATH"
fi

if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source bitstream file not found: $SOURCE_FILE"
    exit 1
fi

./RS_1575_1890_encoder $SOURCE_FILE $scramble_sequence $OUTPUT_DNA_PATH $Left_index_seq $Right_index_seq

if [ $? -ne 0 ]; then
    echo "Error: DNA encoding failed"
    exit 1
fi

echo ""
echo "Pipeline completed successfully!"
echo "Input:  $Input_File"
echo "Output: $OUTPUT_DNA_PATH"
echo "End Processing time :" `date`

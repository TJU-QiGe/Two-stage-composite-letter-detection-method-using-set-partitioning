/******************************************************************************
 * File: txt2bits.cpp
 * Author: [Qi Ge]
 * Version: [1.0]
 * Organization: [School of Microelectronics, Tianjin University]
 * Date: [2025/02/23]
 *
 * Description:
 * This program converts input text files to bit sequences for RS encoding.
 * It reads binary data from a file, converts each byte to its 8-bit 
 * representation, and outputs the result as a sequence of '0' and '1' 
 * characters. If the input is smaller than expected, it pads with zeros.
 *
 * Usage:
 * gcc -o txt2bits txt2bits.cpp
 * ./txt2bits <input_file> <output_file> <expected_total_bits>
 * 
 * License:
 * This program is released under the MIT License. See LICENSE file for details.
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILENAME_LENGTH 500
#define MAX_BITSTREAM_LENGTH (6750 * 8)

/**
 * XOR operation for two character bits ('0' or '1')
 * @param a First bit character
 * @param b Second bit character
 * @return Result of XOR operation as character
 */
char xor_bits(char a, char b) {
    if (a == b) {
        return '0';
    }
    return '1';
}

/**
 * Convert text file to bit sequence for RS encoding
 * Usage: txt2bits <input_file> <output_file> <expected_total_bits>
 */
int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 4) {
        fprintf(stderr, "ERROR: Invalid number of arguments\n");
        fprintf(stderr, "USAGE: %s <input_file> <output_file> <expected_total_bits>\n", argv[0]);
        return 1;
    }

    int expected_total_bits = atoi(argv[3]);
    FILE* input_file = NULL;
    FILE* output_file = NULL;
    
    char bitstream_vector[MAX_BITSTREAM_LENGTH] = {0};
    char input_filename[MAX_FILENAME_LENGTH] = {0};
    char output_filename[MAX_FILENAME_LENGTH] = {0};

    // Prepare file names
    snprintf(input_filename, sizeof(input_filename), "%s", argv[1]);
    snprintf(output_filename, sizeof(output_filename), "%s", argv[2]);

    // Open input file
    if ((input_file = fopen(input_filename, "rb")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open input file '%s'\n", input_filename);
        return 1;
    }
    printf("INFO: Starting text-to-bit conversion process\n");
    
    // Open output file
    if ((output_file = fopen(output_filename, "wb")) == NULL) {
        fprintf(stderr, "ERROR: Cannot create output file '%s'\n", output_filename);
        fclose(input_file);
        return 1;
    }

    // Read file and convert to bits
    size_t bytes_read = 0;
    unsigned char current_byte;
    
    while (fread(&current_byte, sizeof(unsigned char), 1, input_file) == 1) {
        // Convert each byte to 8 bits (MSB first)
        for (int j = 7; j >= 0; j--) {
            char bit = (current_byte >> j) & 0x1;
            bitstream_vector[bytes_read * 8 + (7 - j)] = bit + '0';
        }
        bytes_read++;
    }
    
    printf("INFO: Processed %zu bytes from input file '%s'\n", bytes_read, input_filename);
    
    // Padding with zeros if necessary
    if (bytes_read * 8 < expected_total_bits) {
        size_t missing_bits = expected_total_bits - (bytes_read * 8);
        printf("INFO: Padding with %zu zero bits to meet expected length (%d bits)\n", 
               missing_bits, expected_total_bits);
        
        for (size_t i = 0; i < missing_bits; i++) {
            bitstream_vector[bytes_read * 8 + i] = '0';
        }
    } else if (bytes_read * 8 > expected_total_bits) {
        printf("WARNING: Input file exceeds expected length, truncating to %d bits\n", 
               expected_total_bits);
    }
    
    // Write bit sequence to output file
    for (int i = 0; i < expected_total_bits; i++) {
        fputc(bitstream_vector[i], output_file);
    }

    // Clean up
    fclose(input_file);
    fclose(output_file);
    
    printf("SUCCESS: Conversion completed successfully\n");
    printf("INFO: Output bit sequence written to '%s' (%d bits)\n", 
           output_filename, expected_total_bits);
    
    return 0;
}

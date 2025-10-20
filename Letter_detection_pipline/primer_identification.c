/******************************************************************************
 * File: primer_identification.c
 * Author: [Qi Ge]
 * Version: [1.0]
 * Organization: [School of Microelectronics, Tianjin University]
 * Date: [2025/09/23]
 *
 * Description:
 * This program implements a primer alignment algorithm based on double-ended
 * primers. It is designed to compare and align primer sequences efficiently,
 * ensuring accurate matching and identification.
 *
 * Usage:
 * g++ -o primer_identification primer_identification.c edlib.cpp -lm
 * 
 * License:
 * This program is released under the MIT License. See LICENSE file for details.
 *
 ******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "edlib.h"

#define MAX_SEQUENCE_LENGTH 500
#define MIN_EDIT_DISTANCE 4

// Change STANDARD_PAYLOAD_INDEX_LENGTH to global variable
int STANDARD_PAYLOAD_INDEX_LENGTH = 76;

void reverse(char* sequence, int length)
{
    int i, j, mid = (length - 1) / 2;
    char temp;
    for (i = 0; i <= mid; i++) 
	{
        j = length - 1 - i;
        temp = sequence[i];
        sequence[i] = sequence[j];
        sequence[j] = temp;
    }
}

void calculate_complement(const char *input, char *output, int length) 
{
    for (int i = 0; i < length; i++) 
	{
        switch (input[i]) 
		{
            case 'A':
                output[i] = 'T';
                break;
            case 'T':
                output[i] = 'A';
                break;
            case 'G':
                output[i] = 'C';
                break;
            case 'C':
                output[i] = 'G';
                break;
            default:
                output[i] = 'N';
                break;
        }
    }
    output[length] = '\0'; // Ensure null-terminated string
}

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        printf("Usage: %s <input_file> <output_file> <forward_primer> <reverse_primer>\n", argv[0]);
        printf("Arguments:\n");
        printf("  input_file     : Input FASTQ file path\n");
        printf("  output_file    : Output file for extracted sequences\n");
        printf("  forward_primer : Forward primer sequence\n");
        printf("  reverse_primer : Reverse primer sequence\n");
        return 1;
    }

    const char *input_file = argv[1];
    const char *output_file = argv[2];
    const char *forward_primer = argv[3];
    const char *reverse_primer = argv[4];
    
    printf("[INFO] Using default payload length: %d\n", STANDARD_PAYLOAD_INDEX_LENGTH);

    char current_sequence[500] = {0};
    char reverse_complement_sequence[500] = {0};
    char reverse_complement_forward_primer[50] = {0};
    char reverse_complement_reverse_primer[50] = {0};
    
    // Add variables for reading FASTQ four-line format
    char fastq_header[500] = {0};
    char fastq_plus_line[500] = {0};
    char fastq_quality[500] = {0};

    int forward_edit_distance = 0, forward_start_pos = 0, forward_end_pos = 0;
    int reverse_edit_distance = 0, reverse_start_pos = 0, reverse_end_pos = 0;
    int rc_reverse_edit_distance = 0, rc_reverse_start_pos = 0, rc_reverse_end_pos = 0;
    int rc_forward_edit_distance = 0, rc_forward_start_pos = 0, rc_forward_end_pos = 0;

    FILE *input_file_ptr = fopen(input_file, "rt");
    if (!input_file_ptr) 
    {
        fprintf(stderr, "[ERROR] Cannot open input file: %s\n", input_file);
        return 1;
    }

    FILE *output_file_ptr = fopen(output_file, "wt");
    if (!output_file_ptr) 
    {
        fprintf(stderr, "[ERROR] Cannot open output file: %s\n", output_file);
        return 1;
    }

    // Count total number of reads instead of total lines
    long total_reads_count = 0;
    printf("[INFO] Counting total reads in input file...\n");
    while (fgets(fastq_header, sizeof(fastq_header), input_file_ptr)) 
    {
        if (fastq_header[0] == '@') {
            // Skip sequence line, + line and quality line
            fgets(current_sequence, sizeof(current_sequence), input_file_ptr);
            fgets(fastq_plus_line, sizeof(fastq_plus_line), input_file_ptr);
            fgets(fastq_quality, sizeof(fastq_quality), input_file_ptr);
            total_reads_count++;
        }
    }
    rewind(input_file_ptr);
    printf("[INFO] Total reads found: %ld\n", total_reads_count);

    // Calculate reverse complements
    strcpy(reverse_complement_forward_primer, forward_primer);
    reverse(reverse_complement_forward_primer, strlen(reverse_complement_forward_primer));
    calculate_complement(reverse_complement_forward_primer, reverse_complement_forward_primer, strlen(reverse_complement_forward_primer));

	strcpy(reverse_complement_reverse_primer, reverse_primer);
    reverse(reverse_complement_reverse_primer, strlen(reverse_complement_reverse_primer));
    calculate_complement(reverse_complement_reverse_primer, reverse_complement_reverse_primer, strlen(reverse_complement_reverse_primer));
	
    int processed_reads_count = 0;
    int double_primer_count = 0;
    int invalid_reads_count = 0;
    
    printf("[INFO] Starting primer identification analysis...\n");
    
    // Modify main loop to correctly read FASTQ format
    while(fgets(fastq_header, sizeof(fastq_header), input_file_ptr))
    {
        // Check if it's a FASTQ header line
        if (fastq_header[0] != '@') {
            continue;
        }
        
        // Read sequence line
        if (!fgets(current_sequence, sizeof(current_sequence), input_file_ptr)) {
            break;
        }
        
        // Read + line
        if (!fgets(fastq_plus_line, sizeof(fastq_plus_line), input_file_ptr)) {
            break;
        }
        
        // Read quality line
        if (!fgets(fastq_quality, sizeof(fastq_quality), input_file_ptr)) {
            break;
        }

        int sequence_length = strlen(current_sequence) - 1; // Remove newline character

        processed_reads_count++;
        
        // Flag to mark if valid primer match is found
        int found_valid_match = 0;
        
        // Directly align forward and reverse primers on the original sequence
        EdlibAlignResult forward_primer_result = edlibAlign(forward_primer, strlen(forward_primer), current_sequence, sequence_length, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        EdlibAlignResult reverse_primer_result = edlibAlign(reverse_primer, strlen(reverse_primer), current_sequence, sequence_length, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        
        // Reverse alignment: align reversed complement primers
        EdlibAlignResult rc_reverse_primer_result = edlibAlign(reverse_complement_reverse_primer, strlen(reverse_complement_reverse_primer), current_sequence, sequence_length, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        EdlibAlignResult rc_forward_primer_result = edlibAlign(reverse_complement_forward_primer, strlen(reverse_complement_forward_primer), current_sequence, sequence_length, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

        forward_edit_distance = forward_primer_result.editDistance;
        forward_start_pos = forward_primer_result.startLocations[0];
        forward_end_pos = forward_primer_result.endLocations[0];

        reverse_edit_distance = reverse_primer_result.editDistance;
        reverse_start_pos = reverse_primer_result.startLocations[0];
        reverse_end_pos = reverse_primer_result.endLocations[0];

        rc_reverse_edit_distance = rc_reverse_primer_result.editDistance;
        rc_reverse_start_pos = rc_reverse_primer_result.startLocations[0];
        rc_reverse_end_pos = rc_reverse_primer_result.endLocations[0];

        rc_forward_edit_distance = rc_forward_primer_result.editDistance;
        rc_forward_start_pos = rc_forward_primer_result.startLocations[0];
        rc_forward_end_pos = rc_forward_primer_result.endLocations[0];

        if(forward_edit_distance + reverse_edit_distance <= rc_reverse_edit_distance + rc_forward_edit_distance) // Forward sequence
        {
            if(forward_edit_distance <= MIN_EDIT_DISTANCE && reverse_edit_distance <= MIN_EDIT_DISTANCE)
            {
                if(reverse_start_pos > forward_end_pos)
                {
                    int payload_length = reverse_start_pos - forward_end_pos - 1;
                    
                    double_primer_count++;
                    found_valid_match = 1;
                    fprintf(output_file_ptr,">forward_double_read_%d forward_ed_%d_pos_%d reverse_ed_%d_pos_%d\n",processed_reads_count, forward_edit_distance,forward_start_pos,reverse_edit_distance,reverse_start_pos);

                    // Write sequence fragment in one operation
                    fwrite(&current_sequence[forward_end_pos + 1], 1, payload_length, output_file_ptr);
                    fprintf(output_file_ptr,"\n");
                }
            }
        }
        else 
        {
            if(rc_reverse_edit_distance <= MIN_EDIT_DISTANCE && rc_forward_edit_distance <= MIN_EDIT_DISTANCE)
            {
                if(rc_forward_start_pos > rc_reverse_end_pos)
                {
                    int payload_length = rc_forward_start_pos - rc_reverse_end_pos - 1;
                    double_primer_count++;
                    found_valid_match = 1;
                    fprintf(output_file_ptr,">reverse_double_read_%d rc_reverse_ed_%d_pos_%d rc_forward_ed_%d_pos_%d\n",processed_reads_count, rc_reverse_edit_distance,rc_reverse_start_pos,rc_forward_edit_distance,rc_forward_start_pos);
                    
                    for(int i = rc_reverse_end_pos + 1; i < rc_forward_start_pos; i++)
                    {
                        reverse_complement_sequence[i - rc_reverse_end_pos - 1] = current_sequence[i]; // Reversing and complementing
                    }
                    reverse_complement_sequence[payload_length] = '\0';

                    reverse(reverse_complement_sequence, strlen(reverse_complement_sequence));
                    calculate_complement(reverse_complement_sequence, reverse_complement_sequence, strlen(reverse_complement_sequence));
                    
                    // Write reverse complement sequence in one operation
                    fwrite(reverse_complement_sequence, 1, strlen(reverse_complement_sequence), output_file_ptr);
                    fprintf(output_file_ptr,"\n");
                }
            }
        }
        
        // Count invalid reads if no valid match was found
        if (!found_valid_match) {
            invalid_reads_count++;
        }
        
        edlibFreeAlignResult(forward_primer_result);
        edlibFreeAlignResult(reverse_primer_result);
        edlibFreeAlignResult(rc_reverse_primer_result);
        edlibFreeAlignResult(rc_forward_primer_result);

        if (processed_reads_count % 10 == 0)
        { 
            double progress = (double)processed_reads_count / total_reads_count * 100;
            printf("\r[PROGRESS] Processing: %.2f%% (%d/%ld reads)", progress, processed_reads_count, total_reads_count);
            fflush(stdout); 
        }
    }
    
    printf("\n[RESULTS] Analysis completed successfully!\n");
    printf("=========================================\n");
    printf("Total reads processed: %d\n", processed_reads_count);
    printf("Double primer matches: %d\n", double_primer_count);
    printf("Invalid reads (no matches): %d\n", invalid_reads_count);
    printf("Valid percentage: %.4f%%\n", (double)double_primer_count/(double)processed_reads_count * 100);
    printf("=========================================\n");

    // Create statistics results file
    char stats_filename[600];
    snprintf(stats_filename, sizeof(stats_filename), "%s_statistics.txt", output_file);
    
    FILE *stats_file_ptr = fopen(stats_filename, "wt");
    if (stats_file_ptr) 
    {
        fprintf(stats_file_ptr, "========================================\n");
        fprintf(stats_file_ptr, "    PRIMER IDENTIFICATION STATISTICS    \n");
        fprintf(stats_file_ptr, "========================================\n");
        fprintf(stats_file_ptr, "Input file: %s\n", input_file);
        fprintf(stats_file_ptr, "Output file: %s\n", output_file);
        fprintf(stats_file_ptr, "Forward primer: %s\n", forward_primer);
        fprintf(stats_file_ptr, "Reverse primer: %s\n", reverse_primer);
        fprintf(stats_file_ptr, "Min edit distance threshold: %d\n", MIN_EDIT_DISTANCE);
        fprintf(stats_file_ptr, "Standard payload length: %d\n", STANDARD_PAYLOAD_INDEX_LENGTH);
        fprintf(stats_file_ptr, "\n========== PROCESSING RESULTS ==========\n");
        fprintf(stats_file_ptr, "Total reads processed: %d\n", processed_reads_count);
        fprintf(stats_file_ptr, "Double primer matches: %d\n", double_primer_count);
        fprintf(stats_file_ptr, "Invalid reads (no valid matches): %d\n", invalid_reads_count);
        fprintf(stats_file_ptr, "Valid percentage: %.4f%%\n", (double)double_primer_count/(double)processed_reads_count * 100);
        fprintf(stats_file_ptr, "Invalid percentage: %.4f%%\n", (double)invalid_reads_count/(double)processed_reads_count * 100);
        fprintf(stats_file_ptr, "========================================\n");
        
        fclose(stats_file_ptr);
        printf("[INFO] Statistics saved to: %s\n", stats_filename);
    }
    else 
    {
        fprintf(stderr, "[WARNING] Could not create statistics file: %s\n", stats_filename);
    }

    fclose(input_file_ptr);
    fclose(output_file_ptr);
    
    printf("[SUCCESS] Program completed successfully!\n");
    return 0;
}


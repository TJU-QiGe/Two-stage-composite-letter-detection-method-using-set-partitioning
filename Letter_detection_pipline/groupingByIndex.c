/******************************************************************************
 * File: groupingByIndex.c
 * Author: [Qi Ge]
 * Version: [1.0]
 * Organization: [School of Microelectronics, Tianjin University]
 * Date: [2025/02/23]
 *
 * Description:
 * This program groups reads based on the index sequence at each end of the read, 
 * and is used to subsequently infer composite letters within the groups
 *
 * Usage:
 * gcc -o groupingByIndex groupingByIndex.c -lm
 * 
 * License:
 * This program is released under the MIT License. See LICENSE file for details.
 *
 ******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define MAXLEN 500
#define MAX_SEQ 5000
#define MAX_ADDRESS 126
#define ADDRESS_LEN 16
#define FORWARD_ADDRESS_LEN 8
#define REVERSE_ADDRESS_LEN 8
#define SEQ_CHUNK 8
#define MIN_HAMMING 2
#define STANDARD_LENGTH 76

// Global data structure
typedef struct {
    char addresses_forward[MAX_ADDRESS + 1][FORWARD_ADDRESS_LEN + 3];
    char addresses_reverse[MAX_ADDRESS + 1][REVERSE_ADDRESS_LEN + 3];
    int address_count;
    int match_counts[MAX_ADDRESS];
    char ***grouped_reads;
    int *grouped_alloc_size;
    int grouped_counts[MAX_ADDRESS];
    int read_num;
    int read_with_standard_len;
    int match_num;
    int unmatch_num;
    int drop_num;
    
    // Front and back end matching statistics
    int front_distance_distribution[SEQ_CHUNK + 1];  // Front-end Hamming distance distribution (0-8)
    int back_distance_distribution[SEQ_CHUNK + 1];   // Back-end Hamming distance distribution (0-8)
    int front_back_match_count;                       // Number of sequences with same index for front and back
    int front_back_total_count;                       // Total number with matches in both ends
    
    // Statistics for different matching strategies
    int both_ends_match_count;                        // Number of both-ends matches
    int front_only_match_count;                       // Number of front-only matches
    int back_only_match_count;                        // Number of back-only matches

    // Sequence decoding result record file pointer
    FILE *decode_record_file;                         // Decoding result record file
} ProcessData;

// Function to compute Hamming distance
int hamming_distance(const char* seq1, const char* seq2, int length) 
{
    int distance = 0;
    for (int i = 0; i < length; i++) 
    {
        if (seq1[i] != seq2[i]) 
        {
            distance++;
        }
    }
    return distance;
}


// Load front-end address sequences
int load_front_addresses(const char *front_address_file, ProcessData *data) 
{
    FILE *f_address = fopen(front_address_file, "rt");
    if (!f_address) 
    {
        fprintf(stderr, "[ERROR] Cannot open front address file: %s\n", front_address_file);
        return 0;
    }

    data->address_count = 0;
    char temp_address[FORWARD_ADDRESS_LEN + 3];
    
    while (fgets(temp_address, sizeof(temp_address), f_address)) 
    {
        // Remove trailing newline character
        temp_address[strcspn(temp_address, "\n")] = '\0';
        
        // Copy to forward primer array
        strncpy(data->addresses_forward[data->address_count], temp_address, FORWARD_ADDRESS_LEN);
        data->addresses_forward[data->address_count][FORWARD_ADDRESS_LEN] = '\0';
        
        data->address_count++;
    }
    fclose(f_address);
    printf("[INFO] Loaded %d front-end addresses\n", data->address_count);
    return 1;
}

// Load back-end address sequences
int load_back_addresses(const char *back_address_file, ProcessData *data) 
{
    FILE *f_address = fopen(back_address_file, "rt");
    if (!f_address) 
    {
        fprintf(stderr, "[ERROR] Cannot open back address file: %s\n", back_address_file);
        return 0;
    }

    int back_count = 0;
    char temp_address[REVERSE_ADDRESS_LEN + 3];
    
    while (fgets(temp_address, sizeof(temp_address), f_address)) 
    {
        // Remove trailing newline character
        temp_address[strcspn(temp_address, "\n")] = '\0';
        
        // Copy to reverse primer array
        strncpy(data->addresses_reverse[back_count], temp_address, REVERSE_ADDRESS_LEN);
        data->addresses_reverse[back_count][REVERSE_ADDRESS_LEN] = '\0';
        
        back_count++;
    }
    
    // Check if front and back address counts match
    if (back_count != data->address_count) 
    {
        fprintf(stderr, "[WARNING] Front address count (%d) does not match back address count (%d)!\n", 
                data->address_count, back_count);
        // Use the smaller count
        if (back_count < data->address_count) {
            data->address_count = back_count;
        }
    }
    
    printf("[INFO] Loaded %d back-end addresses (total address pairs: %d)\n", back_count, data->address_count);
    fclose(f_address);
    return 1;
}

// Initialize processing data structure
int initialize_process_data(ProcessData *data) 
{
    memset(data->match_counts, 0, sizeof(data->match_counts));
    memset(data->grouped_counts, 0, sizeof(data->grouped_counts));
    memset(data->front_distance_distribution, 0, sizeof(data->front_distance_distribution));
    memset(data->back_distance_distribution, 0, sizeof(data->back_distance_distribution));
    data->read_num = data->match_num = data->unmatch_num = data->drop_num = 0;
    data->read_with_standard_len = 0;
    data->front_back_match_count = data->front_back_total_count = 0;
    data->both_ends_match_count = data->front_only_match_count = data->back_only_match_count = 0;

    data->grouped_reads = (char ***)malloc(MAX_ADDRESS * sizeof(char **));
    data->grouped_alloc_size = (int *)malloc(MAX_ADDRESS * sizeof(int));

    if (!data->grouped_reads || !data->grouped_alloc_size) 
    {
        fprintf(stderr, "[ERROR] Memory allocation failed!\n");
        return 0;
    }

    for (int i = 0; i < MAX_ADDRESS; i++) 
    {
        data->grouped_reads[i] = NULL;  
        data->grouped_alloc_size[i] = 0; 
    }
    
    printf("[INFO] Data structures initialized successfully\n");
    return 1;
}

// Expand grouped reads array
int expand_grouped_reads(ProcessData *data, int address_index) 
{
    int new_size = data->grouped_alloc_size[address_index] == 0 ? 10 : data->grouped_alloc_size[address_index] * 2;
    data->grouped_reads[address_index] = realloc(data->grouped_reads[address_index], new_size * sizeof(char *));
    
    if (!data->grouped_reads[address_index]) 
    {
        fprintf(stderr, "[ERROR] Memory allocation failed during array expansion!\n");
        return 0;
    }
    
    for (int k = data->grouped_alloc_size[address_index]; k < new_size; k++) 
    {
        data->grouped_reads[address_index][k] = NULL;
    }
    
    data->grouped_alloc_size[address_index] = new_size;
    return 1;
}

// Enhanced version: Process individual sequence with separate front and back matching
void process_sequence(const char *sequence, int sequence_length, ProcessData *data, int predefined_distance) 
{
    char front_seq[SEQ_CHUNK + 1] = {0};
    char back_seq[SEQ_CHUNK + 1] = {0};
    
    // Extract first 8 and last 8 bases
    strncpy(front_seq, sequence, SEQ_CHUNK);
    strncpy(back_seq, sequence + sequence_length - SEQ_CHUNK, SEQ_CHUNK);

    int front_min_distance = SEQ_CHUNK + 1;
    int back_min_distance = SEQ_CHUNK + 1;
    int front_best_index = -1;
    int back_best_index = -1;

    // Front-end matching: compare with first 8 bases of each address
    for (int i = 0; i < data->address_count; i++) 
    {
        int distance = hamming_distance(front_seq, data->addresses_forward[i], SEQ_CHUNK);
        if (distance < front_min_distance) 
        {
            front_min_distance = distance;
            front_best_index = i;
            if (distance == 0)
                break;
        }
    }

    // Back-end matching: compare with last 8 bases of each address
    for (int i = 0; i < data->address_count; i++) 
    {
        int distance = hamming_distance(back_seq, data->addresses_reverse[i], SEQ_CHUNK);
        if (distance < back_min_distance) 
        {
            back_min_distance = distance;
            back_best_index = i;
            if (distance == 0)
                break;
        }
    }

    data->front_distance_distribution[front_min_distance]++;
    data->back_distance_distribution[back_min_distance]++;
    
    // Check if front and back ends match within threshold
    int front_matched = (front_best_index >= 0 && front_min_distance <= predefined_distance);
    int back_matched = (back_best_index >= 0 && back_min_distance <= predefined_distance);
    
    // Record decoding result for each sequence to file
    if (data->decode_record_file) {
        fprintf(data->decode_record_file, "%d %d\n", front_best_index, back_best_index);
    }
    fflush(data->decode_record_file);

    // Use combined matching strategy (can be adjusted as needed)
    int best_address_index = -1;
    
    if (front_matched && back_matched) 
    {
        // Both ends matched
        data->both_ends_match_count++;
        
        // If front and back match the same address, use that address
        if (front_best_index == back_best_index) 
        {
            best_address_index = front_best_index;
            data->front_back_match_count++;
        }
    }

    // If match found, save sequence
    if (best_address_index >= 0)
    {
        data->match_num++;
        data->match_counts[best_address_index]++;
        
        // Check if array expansion is needed
        if (data->grouped_counts[best_address_index] >= data->grouped_alloc_size[best_address_index]) 
        {
            if (!expand_grouped_reads(data, best_address_index)) 
            {
                exit(1);
            }
        }

        // Save middle part after removing address sequences
        int content_length = sequence_length - 2 * SEQ_CHUNK;
        data->grouped_reads[best_address_index][data->grouped_counts[best_address_index]] = malloc((content_length + 1) * sizeof(char));
        if (!data->grouped_reads[best_address_index][data->grouped_counts[best_address_index]]) 
        {
            fprintf(stderr, "[ERROR] Memory allocation failed for sequence storage!\n");
            exit(1);
        }
        strncpy(data->grouped_reads[best_address_index][data->grouped_counts[best_address_index]], sequence + SEQ_CHUNK, content_length);
        data->grouped_reads[best_address_index][data->grouped_counts[best_address_index]][content_length] = '\0';
        data->grouped_counts[best_address_index]++;
    }
    else
    {
        data->unmatch_num++;
    }
}

// Process FASTA file
int process_fasta_file(const char *input_fasta, ProcessData *data, int predefined_distance) 
{
    FILE *f_input = fopen(input_fasta, "rt");
    if (!f_input) 
    {
        fprintf(stderr, "[ERROR] Cannot open input file: %s\n", input_fasta);
        return 0;
    }

    printf("[INFO] Starting to process FASTA file: %s\n", input_fasta);
    printf("[INFO] Using Hamming distance threshold: %d\n", predefined_distance);

    char line[MAXLEN];
    char head_line[MAXLEN];
    char sequence[MAXLEN] = {0};
    int sequence_index = 0;
    
    while (1) 
    {
        memset(head_line, 0, sizeof(char)*MAXLEN);
        memset(line, 0, sizeof(char)*MAXLEN);
        memset(sequence, 0, sizeof(char)*MAXLEN);

        fgets(head_line, sizeof(head_line), f_input);
        fgets(line, sizeof(line), f_input);
        if(feof(f_input))
            break;
        int len = strlen(line) - 1;
        if(len == STANDARD_LENGTH)
        {
            strncpy(sequence, line, len);
            sequence_index += len;
        }
        data->read_num++;
        
        // Progress indicator every 10 reads
        if (data->read_num % 10 == 0) {
            printf("\r[PROGRESS] Processed %d reads", data->read_num);
            fflush(stdout);
        }

        if (sequence_index > 0)
        {
            data->read_with_standard_len++;
            process_sequence(sequence, sequence_index, data, predefined_distance);
        }
        sequence_index = 0;
    }

    fclose(f_input);
    printf("\n[INFO] Finished processing FASTA file\n");
    return 1;
}

// Calculate dropout count
void calculate_dropout(ProcessData *data) 
{
    for(int i = 0; i < data->address_count; i++)
    {
        if(data->grouped_counts[i] == 0)
        {
            data->drop_num++;
        }
    }
}

// Write results to files
int write_results(ProcessData *data, const char *match_count_file, const char *unmatched_file, const char *output_data_file) 
{
    // Write match counts
    FILE *f_match_count = fopen(match_count_file, "wt");
    if (!f_match_count) 
    {
        fprintf(stderr, "[ERROR] Cannot open output file: %s\n", match_count_file);
        return 0;
    }

    for (int i = 0; i < data->address_count; i++) 
    {
        fprintf(f_match_count, "%d %d\n", i, data->match_counts[i]);
    }
    fclose(f_match_count);

    // Write grouped data
    FILE *f_grouped = fopen(output_data_file, "wt");
    if (!f_grouped) 
    {
        fprintf(stderr, "[ERROR] Cannot open output file: %s\n", output_data_file);
        return 0;
    }

    for (int i = 0; i < data->address_count; i++) 
    {
        fprintf(f_grouped, ">Address %d\n", i);
        for (int j = 0; j < data->grouped_counts[i]; j++) 
        {
            fprintf(f_grouped, "%s\n", data->grouped_reads[i][j]);
        }
    }
    fclose(f_grouped);

    // Write statistics
    FILE *f_unmatched = fopen(unmatched_file, "wt");
    if (!f_unmatched) 
    {
        fprintf(stderr, "[ERROR] Cannot open output file: %s\n", unmatched_file);
        return 0;
    }

    fprintf(f_unmatched, "total_reads %d\n", data->read_num);
    fprintf(f_unmatched, "standard_length_reads %d\n", data->read_with_standard_len);
    fprintf(f_unmatched, "matched_reads %d\n", data->match_num);
    fprintf(f_unmatched, "unmatched_reads %d\n", data->unmatch_num);
    fprintf(f_unmatched, "match_percentage %.4f\n", (double)data->match_num/(double)data->read_num);
    fprintf(f_unmatched, "dropout_addresses %d\n", data->drop_num);
    fclose(f_unmatched);

    return 1;
}

// Clear memory
void cleanup_memory(ProcessData *data) 
{
    for (int i = 0; i < MAX_ADDRESS; i++)
    {
        if (data->grouped_reads[i] != NULL) 
        {
            for (int j = 0; j < data->grouped_counts[i]; j++) 
            {
                if (data->grouped_reads[i][j] != NULL) 
                {
                    free(data->grouped_reads[i][j]);
                }
            }
            free(data->grouped_reads[i]);
        }
    }
    free(data->grouped_reads);
    free(data->grouped_alloc_size);
}

int main(int argc, char *argv[]) 
{
    if (argc != 8)
    {
        printf("Usage: %s <input_fasta> <front_address_file> <back_address_file> <match_count_file> <unmatched_file> <output_data_file> <hamming_distance_threshold>\n", argv[0]);
        printf("Arguments:\n");
        printf("  input_fasta            : Input FASTA file path\n");
        printf("  front_address_file     : Front-end address sequences file\n");
        printf("  back_address_file      : Back-end address sequences file\n");
        printf("  match_count_file       : Output file for match counts\n");
        printf("  unmatched_file         : Output file for statistics\n");
        printf("  output_data_file       : Output file for grouped data\n");
        printf("  hamming_distance_threshold : Maximum allowed Hamming distance\n");
        return 1;
    }
    
    const char *input_fasta = argv[1];
    const char *front_address_file = argv[2];
    const char *back_address_file = argv[3];
    const char *match_count_file = argv[4];
    const char *unmatched_file = argv[5];
    const char *output_data_file = argv[6];
    int predefined_distance = atoi(argv[7]);

    printf("[INFO] Starting grouping by index analysis\n");
    printf("[INFO] Input parameters:\n");
    printf("  - Input FASTA: %s\n", input_fasta);
    printf("  - Front address file: %s\n", front_address_file);
    printf("  - Back address file: %s\n", back_address_file);
    printf("  - Hamming distance threshold: %d\n", predefined_distance);

    ProcessData data;
    
    // Open decoding result record file
    char decode_record_filename[256];
    snprintf(decode_record_filename, sizeof(decode_record_filename), "%s.decode_record.txt", output_data_file);
    data.decode_record_file = fopen(decode_record_filename, "wt");
    if (!data.decode_record_file) {
        fprintf(stderr, "[ERROR] Cannot open decode record file: %s\n", decode_record_filename);
        return 1;
    }
    printf("[INFO] Decode record will be saved to: %s\n", decode_record_filename);

    // 1. Load front-end address sequences
    if (!load_front_addresses(front_address_file, &data)) 
    {
        return 1;
    }

    // 2. Load back-end address sequences
    if (!load_back_addresses(back_address_file, &data)) 
    {
        return 1;
    }

    // 3. Initialize processing data structure
    if (!initialize_process_data(&data))
    {
        return 1;
    }

    // 4. Process FASTA file
    if (!process_fasta_file(input_fasta, &data, predefined_distance)) 
    {
        cleanup_memory(&data);
        return 1;
    }  
    
    // Close decode record file
    fclose(data.decode_record_file);

    printf("[INFO] Finished processing FASTA file: %s\n", input_fasta);
    
    // 5. Calculate dropout count
    calculate_dropout(&data);

    // 6. Write results
    if (!write_results(&data, match_count_file, unmatched_file, output_data_file)) 
    {
        cleanup_memory(&data);
        return 1;
    }
    printf("[INFO] Results written to:\n");
    printf("  - Match counts: %s\n", match_count_file);
    printf("  - Statistics: %s\n", unmatched_file);
    printf("  - Grouped data: %s\n", output_data_file);
    
    // 7. Print final statistics
    printf("\n[SUMMARY] Processing completed successfully!\n");
    printf("=========================================\n");
    printf("Total reads processed: %d\n", data.read_num);
    printf("Standard length reads: %d\n", data.read_with_standard_len);
    printf("Successfully matched: %d\n", data.match_num);
    printf("Unmatched reads: %d\n", data.unmatch_num);
    printf("Match percentage: %.4f%%\n", (double)data.match_num/(double)data.read_num * 100.0);
    printf("Dropout addresses: %d\n", data.drop_num);
    printf("=========================================\n");

    // 8. Cleanup memory
    cleanup_memory(&data);

    printf("[SUCCESS] Program completed successfully!\n");
    return 0;
}

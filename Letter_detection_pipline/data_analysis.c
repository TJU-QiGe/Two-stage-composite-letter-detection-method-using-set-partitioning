/******************************************************************************
 * File: data_analysis.c
 * Author: [Qi Ge]
 * Version: [2.0]
 * Organization: [School of Microelectronics, Tianjin University]
 * Date: [2025/10/15]
 *
 * Description:
 * Data analysis program: Analyze the accuracy of consensus sequences, 
 * calculate error rates and entropy distributions for different base types
 * Input: Consensus sequence file, grouped reads file, reference sequence file
 * Output: Various statistical analysis result files
 * 
 * Usage:
 * gcc -o data_analysis data_analysis.c -lm
 * ./data_analysis <consensus_file> <grouped_reads_file> <reference_file> <basetype> <output_dir> <num_bases>
 * 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXLEN 500
#define MAX_LINE 5000
#define MAX_ADDRESS 126
#define STAND_LEN 60

double calculateEntropy(double pA, double pT, double pC, double pG) 
{
    double entropy = 0.0;
    if (pA > 0) entropy -= pA * log2(pA);
    if (pT > 0) entropy -= pT * log2(pT);
    if (pC > 0) entropy -= pC * log2(pC);
    if (pG > 0) entropy -= pG * log2(pG);
    return entropy;
}

void read_sequences(const char *filename, char sequences[][MAXLEN], int max_count)
{
    FILE *file = fopen(filename, "r");
    if (!file) 
    {
        fprintf(stderr, "[ERROR] Failed to open file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    printf("[INFO] Reading sequences from: %s\n", filename);
    
    char line[MAXLEN];
    int count = 0;
    while (fgets(line, MAXLEN, file) && count < max_count) 
    {
        size_t len = strlen(line);
        if (line[len - 1] == '\n') line[--len] = '\0';
        strncpy(sequences[count], line, STAND_LEN);
        sequences[count][STAND_LEN] = '\0';
        count++;
    }

    fclose(file);
    printf("[INFO] Successfully loaded %d sequences\n", count);
}

void read_grouped_reads(const char *filename, char ***grouped_reads, int *grouped_counts) 
{
    FILE *file = fopen(filename, "r");
    if (!file) 
    {
        fprintf(stderr, "[ERROR] Failed to open grouped reads file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    printf("[INFO] Reading grouped reads from: %s\n", filename);

    char line[MAXLEN];
    int current_address = -1;

    while (fgets(line, MAXLEN, file)) 
    {
        size_t len = strlen(line);
        if (line[len - 1] == '\n') line[--len] = '\0';

        if (line[0] == '>') 
        {
            current_address++;
            grouped_counts[current_address] = 0;
        } 
        else 
        {
            if (current_address >= 0 && current_address < MAX_ADDRESS) 
            {
                int read_index = grouped_counts[current_address];
                grouped_reads[current_address][read_index] = strdup(line);
                grouped_counts[current_address]++;
            }
        }
    }

    fclose(file);
    printf("[INFO] Successfully loaded grouped reads from %d addresses\n", current_address + 1);
}

void analyze_entropy_distribution(char ***grouped_reads, int *grouped_counts, 
                                const char *reference_sequences, const char *basetype, 
                                const char *output_dir)
{
    printf("[INFO] Analyzing entropy distribution for base types: %s\n", basetype);
    
    int typenum = strlen(basetype);
    char reference[MAX_ADDRESS][MAXLEN];
    
    // Read reference sequences
    FILE *ref_file = fopen(reference_sequences, "r");
    if (!ref_file) 
    {
        fprintf(stderr, "[ERROR] Failed to open reference file: %s\n", reference_sequences);
        exit(EXIT_FAILURE);
    }
    
    int ref_count = 0;
    char line[MAXLEN];
    while (fgets(line, MAXLEN, ref_file) && ref_count < MAX_ADDRESS) 
    {
        size_t len = strlen(line);
        if (line[len - 1] == '\n') line[--len] = '\0';
        strncpy(reference[ref_count], line, STAND_LEN);
        reference[ref_count][STAND_LEN] = '\0';
        ref_count++;
    }
    fclose(ref_file);

    // Open entropy output files
    FILE *fw_entropy[10] = {NULL};
    for (int i = 0; i < typenum; i++)
    {
        char filename[500];
        sprintf(filename, "%s/Entropy-%c.txt", output_dir, basetype[i]);
        fw_entropy[i] = fopen(filename, "wt");
        if (!fw_entropy[i])
        {
            fprintf(stderr, "[ERROR] Failed to create entropy file: %s\n", filename);
            exit(EXIT_FAILURE);
        }
        printf("[INFO] Created entropy analysis file: %s\n", filename);
    }

    int processed_groups = 0;
    
    // Analyze each group
    for (int j = 0; j < MAX_ADDRESS; j++)
    {
        int n = grouped_counts[j];
        if (n < 4) continue;

        int numA[STAND_LEN] = {0}, numT[STAND_LEN] = {0}, numG[STAND_LEN] = {0}, numC[STAND_LEN] = {0};
        int actually_count = 0;

        // Count bases at each position
        for (int i = 0; i < n; i++)
        {
            if (strlen(grouped_reads[j][i]) != STAND_LEN) continue;
            actually_count++;
            for (int ii = 0; ii < STAND_LEN; ii++)
            {
                switch (grouped_reads[j][i][ii]) 
                {
                    case 'A': numA[ii]++; break;
                    case 'T': numT[ii]++; break;
                    case 'G': numG[ii]++; break;
                    case 'C': numC[ii]++; break;
                }
            }
        }

        if (actually_count < 4) continue;
        processed_groups++;

        // Calculate entropy for each position
        for (int i = 0; i < STAND_LEN; i++)
        {
            double total = numA[i] + numT[i] + numC[i] + numG[i];
            double pA = (double)numA[i] / total;
            double pT = (double)numT[i] / total;
            double pC = (double)numC[i] / total;
            double pG = (double)numG[i] / total;
            double entropy = calculateEntropy(pA, pT, pC, pG);

            // Write entropy to corresponding base type file
            for (int ii = 0; ii < typenum; ii++)
            {
                if (reference[j][i] == basetype[ii])
                {
                    fprintf(fw_entropy[ii], "%lf\n", entropy);
                    break;
                }
            }
        }
    }

    // Close files
    for (int i = 0; i < typenum; i++) 
    {
        fclose(fw_entropy[i]);
    }
    
    printf("[INFO] Entropy analysis completed for %d groups\n", processed_groups);
}

void analyze_consensus_accuracy(const char *consensus_file, const char *reference_file, 
                              const char *basetype, const char *output_dir)
{
    printf("[INFO] Analyzing consensus sequence accuracy\n");
    
    char consensus[MAX_ADDRESS][MAXLEN];
    char reference[MAX_ADDRESS][MAXLEN];
    int typenum = strlen(basetype);

    // Read consensus and reference sequences
    read_sequences(consensus_file, consensus, MAX_ADDRESS);
    read_sequences(reference_file, reference, MAX_ADDRESS);

    // Open output files
    char error_file[500], base_error_file[500], total_error_file[500];
    sprintf(error_file, "%s/Errors_per_strand.txt", output_dir);
    sprintf(base_error_file, "%s/Different_letter_errors.txt", output_dir);
    sprintf(total_error_file, "%s/Total_errors.txt", output_dir);

    FILE *fw_error = fopen(error_file, "wt");
    FILE *fw_base_error = fopen(base_error_file, "wt");
    FILE *fw_total_error = fopen(total_error_file, "wt");

    if (!fw_error || !fw_base_error || !fw_total_error) 
    {
        fprintf(stderr, "[ERROR] Failed to create analysis output files\n");
        exit(EXIT_FAILURE);
    }
    
    printf("[INFO] Created analysis files:\n");
    printf("  - %s\n", error_file);
    printf("  - %s\n", base_error_file);
    printf("  - %s\n", total_error_file);

    int consensus_errors[MAX_ADDRESS] = {0};
    int base_errors[10] = {0}; // Different base type errors
    int erasure_count = 0;
    int total_errors = 0;

    // Analyze each sequence
    for (int i = 0; i < MAX_ADDRESS; i++)
    {
        for (int j = 0; j < STAND_LEN; j++)
        {
            if (consensus[i][j] == 'E')
            {
                if (j == 0) erasure_count++; // Count erasure only once per sequence
                continue;
            }

            if (reference[i][j] != consensus[i][j])
            {
                consensus_errors[i]++;
                total_errors++;

                // Count errors for different base types
                for (int k = 0; k < typenum; k++)
                {
                    if (reference[i][j] == basetype[k])
                    {
                        base_errors[k]++;
                        break;
                    }
                }
            }
        }
    }

    // Write results
    for (int i = 0; i < MAX_ADDRESS; i++)
    {
        fprintf(fw_error, "%d\n", consensus_errors[i]);
    }

    int valid_sequences = MAX_ADDRESS - erasure_count;
    for (int i = 0; i < typenum; i++)
    {
        fprintf(fw_base_error, "%c %.6f\n", basetype[i], 
                (double)base_errors[i] / (STAND_LEN * valid_sequences));
    }

    fprintf(fw_total_error, "erasure %d error %d\n", erasure_count * STAND_LEN, total_errors);
    
    double error_rate = (double)total_errors / (STAND_LEN * valid_sequences);
    printf("[RESULT] Letter detection error rate: %.6f\n", error_rate);
    printf("[RESULT] Total erasures: %d sequences\n", erasure_count);
    printf("[RESULT] Total errors: %d positions\n", total_errors);

    fclose(fw_error);
    fclose(fw_base_error);
    fclose(fw_total_error);
}

int main(int argc, char *argv[])
{
    // Record start time
    clock_t start_time = clock();
    time_t start_wall_time = time(NULL);
    
    printf("================================================================================\n");
    printf("Composite Letter Consensus Sequence Analysis Tool\n");
    printf("================================================================================\n");
    
    if (argc != 7) 
    {
        fprintf(stderr, "[ERROR] Invalid number of arguments\n");
        fprintf(stderr, "Usage: %s <consensus_file> <grouped_reads_file> <reference_file> <basetype> <output_dir> <num_bases>\n", argv[0]);
        fprintf(stderr, "Parameters:\n");
        fprintf(stderr, "  consensus_file      : Input consensus sequence file\n");
        fprintf(stderr, "  grouped_reads_file  : Input grouped reads file\n");
        fprintf(stderr, "  reference_file      : Reference sequence file for comparison\n");
        fprintf(stderr, "  basetype           : Base types to analyze (e.g., 'ATCG')\n");
        fprintf(stderr, "  output_dir         : Directory for output analysis files\n");
        fprintf(stderr, "  num_bases          : Number of bases in composite letters (1-3)\n");
        return EXIT_FAILURE;
    }

    const char *consensus_file = argv[1];
    const char *grouped_reads_file = argv[2];
    const char *reference_file = argv[3];
    const char *basetype = argv[4];
    const char *output_dir = argv[5];
    int num_bases = atoi(argv[6]);

    printf("[INFO] Consensus file: %s\n", consensus_file);
    printf("[INFO] Grouped reads file: %s\n", grouped_reads_file);
    printf("[INFO] Reference file: %s\n", reference_file);
    printf("[INFO] Base types: %s\n", basetype);
    printf("[INFO] Output directory: %s\n", output_dir);
    printf("[INFO] Composite base complexity: %d\n", num_bases);

    // Allocate memory for grouped reads
    char ***grouped_reads = (char ***)malloc(MAX_ADDRESS * sizeof(char **));
    int *grouped_counts = (int *)malloc(MAX_ADDRESS * sizeof(int));
    
    if (!grouped_reads || !grouped_counts) {
        fprintf(stderr, "[ERROR] Memory allocation failed\n");
        return EXIT_FAILURE;
    }
    
    for (int i = 0; i < MAX_ADDRESS; i++)
    {
        grouped_reads[i] = (char **)malloc(MAX_LINE * sizeof(char *));
        if (!grouped_reads[i]) {
            fprintf(stderr, "[ERROR] Memory allocation failed for address %d\n", i);
            return EXIT_FAILURE;
        }
        grouped_counts[i] = 0;
    }

    // Read grouped reads
    read_grouped_reads(grouped_reads_file, grouped_reads, grouped_counts);

    // Perform analyses
    printf("[INFO] Starting entropy distribution analysis...\n");
    analyze_entropy_distribution(grouped_reads, grouped_counts, reference_file, basetype, output_dir);
    
    printf("[INFO] Starting consensus accuracy analysis...\n");
    analyze_consensus_accuracy(consensus_file, reference_file, basetype, output_dir);

    // Calculate execution time
    clock_t end_time = clock();
    time_t end_wall_time = time(NULL);
    double cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    double wall_time = difftime(end_wall_time, start_wall_time);

    printf("================================================================================\n");
    printf("Analysis Summary:\n");
    printf("  CPU Time: %.3f seconds\n", cpu_time);
    printf("  Wall Time: %.0f seconds\n", wall_time);
    printf("  Analyzed Addresses: %d\n", MAX_ADDRESS);
    printf("  Sequence Length: %d bp\n", STAND_LEN);
    printf("  Base Types: %s\n", basetype);
    printf("[SUCCESS] All analyses completed successfully\n");
    printf("[OUTPUT] Results saved to directory: %s\n", output_dir);
    printf("================================================================================\n");

    // Free memory
    for (int i = 0; i < MAX_ADDRESS; i++) 
    {
        for (int j = 0; j < grouped_counts[i]; j++) 
        {
            free(grouped_reads[i][j]);
        }
        free(grouped_reads[i]);
    }
    free(grouped_reads);
    free(grouped_counts);

    printf("[INFO] Program completed successfully\n");
    return EXIT_SUCCESS;
}

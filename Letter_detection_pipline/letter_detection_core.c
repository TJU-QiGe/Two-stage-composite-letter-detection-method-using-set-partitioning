/******************************************************************************
 * File: letter_detection_core.c
 * Author: [Qi Ge]
 * Version: [2.0]
 * Organization: [School of Microelectronics, Tianjin University]
 * Date: [2025/10/15]
 *
 * Description:
 * Core letter detection program: Composite letter detection based on information entropy 
 * and maximum likelihood estimation
 * Input: Grouped reads file
 * Output: Consensus sequence file
 * 
 * Usage:
 * gcc -o letter_detection_core letter_detection_core.c -lm
 * ./letter_detection_core <grouped_reads_file> <output_consensus_file> <num_bases>
 * 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define MAXLEN 500
#define MAX_LINE 5000
#define MAX_ADDRESS 126
#define STAND_LEN 60

// Define the frequency vector of a composite letter
typedef struct 
{
    char symbol;
    double sigma[4]; // A, C, G, T
} DegenerateBase;

// Define concatenated base list (Mixture of 2 bases)
DegenerateBase degenerate_bases_one[] = 
{
    {'R', {1, 0, 1, 0}},  // A, G
    {'Y', {0, 1, 0, 1}},  // C, T
    {'M', {1, 1, 0, 0}},  // A, C
    {'K', {0, 0, 1, 1}},  // G, T
};

DegenerateBase degenerate_bases_two[] = 
{
    {'R', {1, 0, 1, 0}},  // A, G
    {'Y', {0, 1, 0, 1}},  // C, T
    {'M', {1, 1, 0, 0}},  // A, C
    {'K', {0, 0, 1, 1}},  // G, T
    {'S', {0, 1, 1, 0}},  // G, C
    {'W', {1, 0, 0, 1}},  // A, T
};

// Define concatenated base list (Mixture of 3 bases)
DegenerateBase degenerate_bases_three[] = 
{
    {'H', {1, 1, 0, 1}},  // A, C, T
    {'B', {0, 1, 1, 1}},  // G, C, T
    {'V', {1, 1, 1, 0}},  // A, G, C
    {'D', {1, 0, 1, 1}},  // A, G, T
};

#define NUM_DEGENERATE_BASES_ONE (sizeof(degenerate_bases_one) / sizeof(degenerate_bases_one[0]))
#define NUM_DEGENERATE_BASES_TWO (sizeof(degenerate_bases_two) / sizeof(degenerate_bases_two[0]))
#define NUM_DEGENERATE_BASES_THREE (sizeof(degenerate_bases_three) / sizeof(degenerate_bases_three[0]))

char find_most_likely_letter(int count_A, int count_C, int count_G, int count_T, int num_bases) 
{
    double epsilon = 1e-6;
    char best_base = '\0';
    double max_log_likelihood = -DBL_MAX;
    DegenerateBase *bases;
    int num_degenerate_bases;

    if(num_bases == 1)
    {
        bases = degenerate_bases_one;
        num_degenerate_bases = NUM_DEGENERATE_BASES_ONE;
    }
    else if (num_bases == 2) 
    {
        bases = degenerate_bases_two;
        num_degenerate_bases = NUM_DEGENERATE_BASES_TWO;
    } 
    else 
    {
        bases = degenerate_bases_three;
        num_degenerate_bases = NUM_DEGENERATE_BASES_THREE;
    }

    for (int i = 0; i < num_degenerate_bases; i++)
    {
        double p_A = (bases[i].sigma[0] + epsilon) / (bases[i].sigma[0] + bases[i].sigma[1] + bases[i].sigma[2] + bases[i].sigma[3] + 4 * epsilon);
        double p_C = (bases[i].sigma[1] + epsilon) / (bases[i].sigma[0] + bases[i].sigma[1] + bases[i].sigma[2] + bases[i].sigma[3] + 4 * epsilon);
        double p_G = (bases[i].sigma[2] + epsilon) / (bases[i].sigma[0] + bases[i].sigma[1] + bases[i].sigma[2] + bases[i].sigma[3] + 4 * epsilon);
        double p_T = (bases[i].sigma[3] + epsilon) / (bases[i].sigma[0] + bases[i].sigma[1] + bases[i].sigma[2] + bases[i].sigma[3] + 4 * epsilon);

        double log_likelihood = count_A * log(p_A) + count_C * log(p_C) + count_G * log(p_G) + count_T * log(p_T);

        if (log_likelihood > max_log_likelihood) 
        {
            max_log_likelihood = log_likelihood;
            best_base = bases[i].symbol;
        }
    }

    return best_base;
}

double calculateEntropy(double pA, double pT, double pC, double pG) 
{
    double entropy = 0.0;
    if (pA > 0) entropy -= pA * log2(pA);
    if (pT > 0) entropy -= pT * log2(pT);
    if (pC > 0) entropy -= pC * log2(pC);
    if (pG > 0) entropy -= pG * log2(pG);
    return entropy;
}

void swap(double *a, double *b, char *c, char *d) 
{
    double temp = *a;
    *a = *b;
    *b = temp;

    char temp2 = *c;
    *c = *d;
    *d = temp2;
}

int partition(double arr[], char bases[], int low, int high) 
{
    double pivot = arr[high];
    int i = (low - 1);
  
    for (int j = low; j <= high - 1; j++) 
    {
        if (arr[j] > pivot)
        {
            i++;
            swap(&arr[i], &arr[j], &bases[i], &bases[j]);
        }
    }
    swap(&arr[i + 1], &arr[high], &bases[i + 1], &bases[high]);
    return (i + 1);
}

void quickSort(double arr[], char bases[], int low, int high)
{
    if (low < high) 
    {
        int pi = partition(arr, bases, low, high);
        quickSort(arr, bases, low, pi - 1);
        quickSort(arr, bases, pi + 1, high);
    }
}

char determine_consensus_base(int *numA, int *numT, int *numG, int *numC, int position, int n, int num_bases)
{
    double pA = (double)numA[position] / n;
    double pT = (double)numT[position] / n;
    double pG = (double)numG[position] / n;
    double pC = (double)numC[position] / n;
    double entropy = calculateEntropy(pA, pT, pC, pG);
    
    if (num_bases == 1)
    {
        if (entropy < 0.5) // ATGC
        {
            char bases[4] = {'A', 'T', 'C', 'G'};
            double probabilities[4] = {pA, pT, pC, pG};
            quickSort(probabilities, bases, 0, 3);
            return bases[0];
        }
        else // RYMK
        {
            return find_most_likely_letter(numA[position], numC[position], numG[position], numT[position], num_bases);
        }
    }
    else if (num_bases == 2)
    {
        if (entropy < 0.50) // ATGC
        {
            return (pA >= pT) ? 'A' : 'T';
        }
        else // RYMKSW
        {
            return find_most_likely_letter(numA[position], numC[position], numG[position], numT[position], num_bases);
        }
    }
    else if (num_bases == 3)
    {
        if (entropy < 1) // ATGC
        {
            char bases[4] = {'A', 'T', 'C', 'G'};
            double probabilities[4] = {pA, pT, pC, pG};
            quickSort(probabilities, bases, 0, 3);
            return bases[0];
        }
        else // BDHV
        {
            return find_most_likely_letter(numA[position], numC[position], numG[position], numT[position], num_bases);
        }
    }
    
    return 'N'; // Default case
}

void process_grouped_reads(char ***grouped_reads, int *grouped_counts, const char *output_file, int num_bases) 
{
    FILE *fw_consensus = fopen(output_file, "wt");
    if (!fw_consensus) 
    {
        fprintf(stderr, "Error opening output file: %s\n", output_file);
        exit(EXIT_FAILURE);
    }

    for (int j = 0; j < MAX_ADDRESS; j++)
    {
        int numA[STAND_LEN] = {0}, numT[STAND_LEN] = {0}, numG[STAND_LEN] = {0}, numC[STAND_LEN] = {0};
        char consensusBase[STAND_LEN + 1] = {0};
        int n = grouped_counts[j];
        int actually_count = 0;

        if (n < 4)
        {
            // Mark as erasure
            for (int i = 0; i < STAND_LEN; i++)
            {
                consensusBase[i] = 'E';
            }
        }
        else
        {
            // Count bases at each position
            for (int i = 0; i < n; i++)
            {
                if (strlen(grouped_reads[j][i]) != STAND_LEN) 
                {
                    continue;
                }
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

            if (actually_count < 4)
            {
                for (int i = 0; i < STAND_LEN; i++)
                {
                    consensusBase[i] = 'E';
                }
            }
            else
            {
                // Determine consensus base for each position
                for (int i = 0; i < STAND_LEN; i++)
                {
                    consensusBase[i] = determine_consensus_base(numA, numT, numG, numC, i, actually_count, num_bases);
                }
            }
        }

        // Write consensus sequence
        for (int i = 0; i < STAND_LEN; i++)
        {
            fprintf(fw_consensus, "%c", consensusBase[i]);
        }
        fprintf(fw_consensus, "\n");
    }

    fclose(fw_consensus);
}

void read_grouped_reads(const char *filename, char ***grouped_reads, int *grouped_counts) 
{
    FILE *file = fopen(filename, "r");
    if (!file) 
    {
        fprintf(stderr, "[ERROR] Failed to open input file: %s\n", filename);
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

int main(int argc, char *argv[])
{
    // Record start time
    clock_t start_time = clock();
    time_t start_wall_time = time(NULL);
    
    printf("================================================================================\n");
    printf("Letter Detection via Set Partitioning\n");
    printf("================================================================================\n");
    
    if (argc != 4) 
    {
        fprintf(stderr, "[ERROR] Invalid number of arguments\n");
        fprintf(stderr, "Usage: %s <grouped_reads_file> <output_consensus_file> <num_bases>\n", argv[0]);
        fprintf(stderr, "Parameters:\n");
        fprintf(stderr, "  grouped_reads_file    : Input file containing grouped DNA reads\n");
        fprintf(stderr, "  output_consensus_file : Output file for consensus sequences\n");
        fprintf(stderr, "  num_bases            : Number of bases in composite letters (1-3)\n");
        return EXIT_FAILURE;
    }

    const char *grouped_reads_file = argv[1];
    const char *output_consensus_file = argv[2];
    int num_bases = atoi(argv[3]);
    
    printf("[INFO] Input file: %s\n", grouped_reads_file);
    printf("[INFO] Output file: %s\n", output_consensus_file);
    printf("[INFO] Composite base complexity: %d\n", num_bases);
    printf("[INFO] Starting consensus sequence generation...\n");

    // Allocate memory
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

    // Read grouped reads and process
    read_grouped_reads(grouped_reads_file, grouped_reads, grouped_counts);
    
    printf("[INFO] Processing sequences and generating consensus...\n");
    process_grouped_reads(grouped_reads, grouped_counts, output_consensus_file, num_bases);

    // Calculate execution time
    clock_t end_time = clock();
    time_t end_wall_time = time(NULL);
    double cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    double wall_time = difftime(end_wall_time, start_wall_time);

    printf("[SUCCESS] Consensus sequences successfully generated\n");
    printf("[OUTPUT] Results saved to: %s\n", output_consensus_file);
    printf("================================================================================\n");
    printf("Execution Summary:\n");
    printf("  CPU Time: %.3f seconds\n", cpu_time);
    printf("  Wall Time: %.0f seconds\n", wall_time);
    printf("  Processed Addresses: %d\n", MAX_ADDRESS);
    printf("  Sequence Length: %d bp\n", STAND_LEN);
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

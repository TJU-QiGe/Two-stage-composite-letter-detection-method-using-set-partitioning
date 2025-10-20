/******************************************************************************
 * File: recovery_poem.c
 * Author: [Qi Ge]
 * Version: [2.0]
 * Organization: [School of Microelectronics, Tianjin University]
 * Date: [2025/10/15]
 *
 * Description:
 * Data recovery program: Recovers original text data from encoded binary files
 * using XOR descrambling with watermark data. Converts binary representation
 * back to character format and handles GBK to UTF-8 encoding conversion.
 * Input: Encoded binary file, watermark file
 * Output: Recovered original text file
 * 
 * Usage:
 * gcc -o recovery_poem recovery_poem.c
 * ./recovery_poem <encoded_bit_file> <output_recovered_file>
 * 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <locale.h>
#include <iconv.h>
#include <time.h>

#define WATERMARK_SIZE 18900
#define CHAR_COUNT 2355
#define BITS_PER_CHAR 8
#define TOTAL_BITS (CHAR_COUNT * BITS_PER_CHAR)

/**
 * Performs XOR operation on two binary characters
 * @param a First binary character ('0' or '1')
 * @param b Second binary character ('0' or '1')  
 * @return Result of XOR operation ('0' or '1')
 */
char dexor(char a, char b) 
{
    if(a=='0' && b == '0') return '0';
    if(a=='0' && b == '1') return '1';
    if(a=='1' && b == '0') return '1';
    return '0';
}

/**
 * Converts GBK encoded text to UTF-8 and prints to stdout
 * @param input Input string in GBK encoding
 * @param in_len Length of input string
 */
void convert_gbk_to_utf8_and_print(const char *input, size_t in_len) 
{
    setlocale(LC_ALL, "");

    iconv_t cd = iconv_open("UTF-8", "GBK");
    if (cd == (iconv_t)-1) {
        fprintf(stderr, "[ERROR] Failed to initialize character encoding conversion\n");
        return;
    }

    size_t out_len = in_len * 2;
    char *outbuf = malloc(out_len);
    if (!outbuf) {
        fprintf(stderr, "[ERROR] Memory allocation failed for encoding conversion\n");
        iconv_close(cd);
        return;
    }
    
    char *outptr = outbuf;
    const char *inptr = input;
    memset(outbuf, 0, out_len);

    if (iconv(cd, (char **)&inptr, &in_len, &outptr, &out_len) == (size_t)-1) {
        fprintf(stderr, "[ERROR] Character encoding conversion failed\n");
    } else {
        printf("[RECOVERED TEXT]\n");
        fwrite(outbuf, 1, strlen(outbuf), stdout);
        printf("\n");
    }

    free(outbuf);
    iconv_close(cd);
}

int main(int argc, char *argv[]) 
{
    // Record start time
    clock_t start_time = clock();
    time_t start_wall_time = time(NULL);
    
    printf("================================================================================\n");
    printf("DNA Storage Data Recovery\n");
    printf("================================================================================\n");
    
    if (argc != 3) {
        fprintf(stderr, "[ERROR] Invalid number of arguments\n");
        fprintf(stderr, "Usage: %s <encoded_bit_file> <output_recovered_file>\n", argv[0]);
        fprintf(stderr, "Parameters:\n");
        fprintf(stderr, "  encoded_bit_file      : Input file containing encoded binary data\n");
        fprintf(stderr, "  output_recovered_file : Output file for recovered text data\n");
        return EXIT_FAILURE;
    }

    setlocale(LC_ALL, "");

    printf("[INFO] Input encoded file: %s\n", argv[1]);
    printf("[INFO] Output recovered file: %s\n", argv[2]);
    printf("[INFO] Expected data size: %d characters (%d bits)\n", CHAR_COUNT, TOTAL_BITS);

    // Allocate memory for decoded results
    char *dec_results = (char *)malloc(TOTAL_BITS * sizeof(char));
    if (dec_results == NULL) {
        fprintf(stderr, "[ERROR] Memory allocation failed for decoded results\n");
        return EXIT_FAILURE;
    }

    // Read watermark file
    char buff_scramble[WATERMARK_SIZE] = {0};
    FILE *fp_water = fopen("./configureFiles/watermark_18900bit.txt", "r");
    if (fp_water == NULL) {
        fprintf(stderr, "[ERROR] Failed to open watermark file: watermark_18900bit.txt\n");
        free(dec_results);
        return EXIT_FAILURE;
    }
    
    int len_water = fread(buff_scramble, sizeof(char), WATERMARK_SIZE, fp_water);
    fclose(fp_water);
    
    if (len_water != WATERMARK_SIZE) {
        fprintf(stderr, "[WARNING] Watermark file size mismatch: expected %d, got %d\n", 
                WATERMARK_SIZE, len_water);
    } else {
        printf("[INFO] Successfully loaded watermark data (%d bits)\n", len_water);
    }

    // Read encoded input file
    FILE *file = fopen(argv[1], "r");
    if (file == NULL) {
        fprintf(stderr, "[ERROR] Failed to open input file: %s\n", argv[1]);
        free(dec_results);
        return EXIT_FAILURE;
    }

    size_t result = fread(dec_results, sizeof(char), TOTAL_BITS, file);
    if (result != TOTAL_BITS) {
        fprintf(stderr, "[ERROR] Input file size mismatch: expected %zu bits, got %zu bits\n", 
                (size_t)TOTAL_BITS, result);
        free(dec_results);
        fclose(file);
        return EXIT_FAILURE;
    }
    fclose(file);
    printf("[INFO] Successfully read %zu bits from input file\n", result);

    // Perform XOR descrambling
    for(int i = 0; i < TOTAL_BITS; i++) {
        dec_results[i] = dexor(dec_results[i], buff_scramble[i]);
    }

    // Convert binary to characters
    char audio[CHAR_COUNT] = {0};
    int t1 = 0;
    for (int i = 0; i < CHAR_COUNT; i++) {
        for (int j = BITS_PER_CHAR * i; j < BITS_PER_CHAR * (i + 1); j++) {
            audio[t1] += (dec_results[j] - '0') * (1 << (7 - (j % BITS_PER_CHAR)));
        }
        t1++;
    }

    // Save recovered file
    printf("[INFO] Saving recovered data to output file...\n");
    FILE *fpoutput1 = fopen(argv[2], "wb");
    if (fpoutput1 == NULL) {
        fprintf(stderr, "[ERROR] Failed to create output file: %s\n", argv[2]);
        free(dec_results);
        return EXIT_FAILURE;
    }
    
    fwrite(audio, sizeof(char), CHAR_COUNT, fpoutput1);
    fclose(fpoutput1);

    // Display recovered content
    convert_gbk_to_utf8_and_print(audio, sizeof(audio));

    // Calculate execution time
    clock_t end_time = clock();
    time_t end_wall_time = time(NULL);
    double cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    double wall_time = difftime(end_wall_time, start_wall_time);

    printf("================================================================================\n");
    printf("Recovery Summary:\n");
    printf("  CPU Time: %.3f seconds\n", cpu_time);
    printf("  Wall Time: %.0f seconds\n", wall_time);
    printf("  Processed Data: %d characters (%d bits)\n", CHAR_COUNT, TOTAL_BITS);
    printf("  Watermark Size: %d bits\n", WATERMARK_SIZE);
    printf("[SUCCESS] Data recovery completed successfully\n");
    printf("[OUTPUT] Recovered data saved to: %s\n", argv[2]);
    printf("================================================================================\n");

    // Cleanup
    free(dec_results);
    printf("[INFO] Program completed successfully\n");

    return EXIT_SUCCESS;
}
#ifndef DECODE_RS_H
#define DECODE_RS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <float.h>
#include <limits.h>

// 常量定义
#define MAX_RANDOM LONG_MAX
#define ROWS 126
#define COLS 60
#define informLen 18900
#define rsCodeNum 1
#define rsCodeLen 12
#define check_length 315
#define rs_k 1575
#define rs_n 1890
#define payloadLen 60
#define NUM_GROUPS_PER_ROW 15
#define BASES_PER_GROUP 4
#define TOTAL_BASES_PER_ROW (NUM_GROUPS_PER_ROW * BASES_PER_GROUP)
#define TOTAL_BITS_PER_ROW 180
#define GROUP_SIZE 15

extern int ri, t, m, n, length, k, t2, d, parity_len, init_zero;
extern int* p;
extern int numerr, errpos[5000], errval[5000], decerror;
extern int biterror, error;
extern int numera, era[5000], eraval[5000];
extern int* g, * rdata, * index_of, * alpha_to, * recd, * b;

// Memory management functions
void init_rs_parameters();
void cleanup_rs_memory();
int** allocate_2d_int_array(int rows, int cols);
void free_2d_int_array(int** array, int rows);
char** allocate_2d_char_array(int rows, int cols);
void free_2d_char_array(char** array, int rows);

// Command line parsing function
int parse_command_line(int args, char** argv, char* noisy_composite_dna, char* output_bit_stream);

// DNA sequence processing functions
int process_dna_sequences(const char* filename, int rule, int bitSequence[rs_n][rsCodeLen], int* ERROR);
void convert_bits_to_rs_codewords(int bitSequence[rs_n][rsCodeLen], int** rsCodeWithErr, int* ERROR);

// Reed-Solomon decoding functions
void setup_erasure_positions(int* ERROR, int* era, int* numera);
int perform_rs_decoding(int** rsCodeWithErr, int** codewordRS, FILE* fw_bit_output, int* ERROR);

// Utility functions
void mapCharacterToBits(char c, int rule, char* bits);
char dexor(char a, char b);
void Dec2Bin(int num10, int* output);

// Reed-Solomon core functions
void read_p(void);
void generate_gf(void);
void gen_poly(void);
void encode_rs(void);
void decode_rs(void);

#endif

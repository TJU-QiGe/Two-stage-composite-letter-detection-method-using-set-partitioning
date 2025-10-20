/******************************************************************************
 * File: RS_1575_1890_encoder.cpp
 * Author: [Qi Ge]
 * Version: [1.0]
 * Organization: [School of Microelectronics, Tianjin University]
 * Date: [2025/02/23]
 *
 * Description:
 * Reed-Solomon (1575,1890) encoder for DNA storage system with composite letter
 * encoding. This program performs RS encoding on input data with scrambling
 * and converts the result to composite letter strands (A,T,G,C,R,Y,M,K).
 *
 * Usage:
 * g++ -o RS_encoder RS_1575_1890_encoder.cpp -lm
 * ./RS_encoder <input_data_file> <scramble_sequence_file> <output_directory> <left_index_file> <right_index_file>
 * 
 * License:
 * This program is released under the MIT License. See LICENSE file for details.
 *
 ******************************************************************************/

#define _CRT_SECURE_NO_WARNINGS
#define _cplusplus

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <float.h>
#include <limits.h>
#include <sys/stat.h>

using namespace std;
//code rate 1575/1890=0.833
//
#define N_address 20  // Fixed number of bits
#define MAX_RANDOM LONG_MAX    // Maximum value of random() 
#define MAXLEN 20000000
#define informLen 18900
#define rsCodeNum 1
#define rsCodeLen 12
#define rs_k 1575
#define rs_n 1890

#define payloadLen 60
#define NUM_GROUPS_PER_ROW 15
#define BASES_PER_GROUP 4
#define TOTAL_BASES_PER_ROW (NUM_GROUPS_PER_ROW * BASES_PER_GROUP)
#define GROUP_SIZE 15
#define check_length 315 
char buff[informLen] = { 0 };
char buff_scramble[informLen] = {0};
char buff_scramble_char[informLen] = { 0 };
int Pre_RS_array[informLen] = {0}; // 15*27*16200


int ri;
int t;
int m, n, length, k, t2, d, parity_check;
int init_zero;
int* p;
int numerr, errpos[600], errval[600], decerror;
int biterror, error;
char filename[40], name2[40];
int numera;
int era[600], eraval[600];
int* g, * rdata, * index_of, * alpha_to, * recd, * b;


void read_p(void);
void generate_gf(void);
void gen_poly(void);
void encode_rs(void);

void dec2bin(int num10, int* output) // Convert decimal to binary
{
	int arr[N_address] = { 0 }; // Initialize array to 0 (pad with zeros when insufficient bits)
	int i;


	for (i = N_address - 1; i >= 0; i--)  // Assign values to array from back to front
	{
		arr[i] = num10 % 2;
		num10 = num10 / 2;
	}


	for (i = 0; i <= N_address - 1; i++)
	{

		output[i] = arr[i];

	}

}
char dexor(char a, char b) // XOR operation
{
	//char output;
	if (a == '0' && b == '0')
	{
		// output='0';
		return '0';
	}
	if (a == '0' && b == '1')
	{
		// output='1';
		return '1';
	}
	if (a == '1' && b == '0')
	{
		// output='1';
		return '1';
	}
	if (a == '1' && b == '1')
	{
		// output='0';
		return '0';
	}
	return '0';
}

/**
 * Check if directory exists
 * @param path Directory path to check
 * @return true if exists, false otherwise
 */
bool directory_exists(const char* path) {
    struct stat info;
    return (stat(path, &info) == 0 && (info.st_mode & S_IFDIR));
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 6) {
        fprintf(stderr, "ERROR: Invalid number of arguments\n");
        fprintf(stderr, "USAGE: %s <input_data_file> <scramble_sequence_file> <output_directory> <left_index_file> <right_index_file>\n", argv[0]);
        fprintf(stderr, "  input_data_file: Text file containing data to be encoded (18900 bits)\n");
        fprintf(stderr, "  scramble_sequence_file: Text file containing scrambling sequence (18900 bits)\n");
        fprintf(stderr, "  output_directory: Directory to store output files\n");
        fprintf(stderr, "  left_index_file: File containing left-end index sequences (126 x 8nt)\n");
        fprintf(stderr, "  right_index_file: File containing right-end index sequences (126 x 8nt)\n");
        return 1;
    }

    // Parse command line arguments
    char input_data_file[500];
    char scramble_file[500];
    char output_dir[500];
    char left_index_file[500];
    char right_index_file[500];
    
    snprintf(input_data_file, sizeof(input_data_file), "%s", argv[1]);
    snprintf(scramble_file, sizeof(scramble_file), "%s", argv[2]);
    snprintf(output_dir, sizeof(output_dir), "%s", argv[3]);
    snprintf(left_index_file, sizeof(left_index_file), "%s", argv[4]);
    snprintf(right_index_file, sizeof(right_index_file), "%s", argv[5]);

    // Check if output directory exists
    if (!directory_exists(output_dir)) {
        fprintf(stderr, "ERROR: Output directory '%s' does not exist\n", output_dir);
        return 1;
    }

    printf("INFO: Starting RS encoding process\n");
    printf("INFO: Input data file: %s\n", input_data_file);
    printf("INFO: Scramble sequence file: %s\n", scramble_file);
    printf("INFO: Output directory: %s\n", output_dir);
    printf("INFO: Left index file: %s\n", left_index_file);
    printf("INFO: Right index file: %s\n", right_index_file);

    int i, j, k_1 = 0;

    /***************************** RS Encoding ********************************/
    // Load pseudo-random sequence
    FILE* fp_water;
    fp_water = fopen(scramble_file, "r");
    if (fp_water == NULL) {
        fprintf(stderr, "ERROR: Cannot open scramble sequence file '%s'\n", scramble_file);
        return 1;
    }
    int len1 = 0;
    len1 = fread(buff_scramble, sizeof(char), informLen, fp_water);
    fclose(fp_water);
    printf("INFO: Loaded %d bits from scramble sequence file\n", len1);

    // Load stored data
    FILE* fp_scramble;
    fp_scramble = fopen(input_data_file, "r");
    if (fp_scramble == NULL) {
        fprintf(stderr, "ERROR: Cannot open input data file '%s'\n", input_data_file);
        return 1;
    }

    len1 = fread(buff, sizeof(char), informLen, fp_scramble);
    printf("INFO: Loaded %d bits from input data file\n", len1);

    for (i = 0; i < informLen; i++) {
        buff[i] = dexor(buff[i], buff_scramble[i]);
    }
    // Use memcpy function for data copying
    memcpy(buff_scramble_char, buff, informLen);
    fclose(fp_scramble);
    printf("INFO: Data scrambling completed\n");

    // Prepare output file paths
    char oligopool_file[600];
    char original_bits_file[600];
    char rs_codeword_file[600];
    char rs_bits_file[600];
    char payload_file[600];
    
    snprintf(oligopool_file, sizeof(oligopool_file), "%s/126composite_strands_116letters.txt", output_dir);
    snprintf(original_bits_file, sizeof(original_bits_file), "%s/original_bits_18900.txt", output_dir);
    snprintf(rs_codeword_file, sizeof(rs_codeword_file), "%s/rs_encoded_codeword.txt", output_dir);
    snprintf(rs_bits_file, sizeof(rs_bits_file), "%s/rscode_bits.txt", output_dir);
    snprintf(payload_file, sizeof(payload_file), "%s/composite_payload_ref.txt", output_dir);

    FILE* fw_oligopool = fopen(oligopool_file, "wt");
    if (fw_oligopool == NULL) {
        fprintf(stderr, "ERROR: Cannot create output file '%s'\n", oligopool_file);
        return 1;
    }


	int** Pre_RS_matrix;
	Pre_RS_matrix = (int**)malloc(sizeof(int*) * rs_k); // 1575 rows
	for (int l0 = 0; l0 < rs_k; l0++)
	{
		Pre_RS_matrix[l0] = (int*)malloc(sizeof(int) * (rsCodeNum * rsCodeLen)); // 12 columns
	}

	int** Pre_RS_Symbol_matrix;
	Pre_RS_Symbol_matrix = (int**)malloc(sizeof(int*) * rs_k); // 1575 rows
	for (int r = 0; r < rs_k; r++)
	{
		Pre_RS_Symbol_matrix[r] = (int*)malloc(sizeof(int) * rsCodeNum); // 1 column
	}


	int** Multi_RS_Codes; // Codeword values
	Multi_RS_Codes = (int**)malloc(sizeof(int*) * rsCodeNum); // 1 row
	for (int i1 = 0; i1 < rsCodeNum; i1++)
	{
		Multi_RS_Codes[i1] = (int*)malloc(sizeof(int) * rs_n); // 1890 columns
	}

	int** rsEncodedBits;
	rsEncodedBits = (int**)malloc(sizeof(int*) * rs_n);
	for (int l1 = 0; l1 < rs_n; l1++)
	{
		rsEncodedBits[l1] = (int*)malloc(sizeof(int) * rsCodeLen * rsCodeNum);
		memset(rsEncodedBits[l1], 0, sizeof(int) * rsCodeLen * rsCodeNum); // 12 bits
	}
	for (int x = 0; x < rs_n; x++)
	{
		for (int y = 0; y < rsCodeLen * rsCodeNum; y++)
		{
			rsEncodedBits[x][y] = 0;
		}
	} // Initialize RS code binary array

	int** base15_symbols = (int**)malloc(rs_n * sizeof(int*));
	if (base15_symbols == NULL)
	{
		fprintf(stderr, "ERROR: Memory allocation failed\n");
		return 1;
	}
	for (int i = 0; i < rs_n; i++)
	{
		base15_symbols[i] = (int*)malloc((TOTAL_BASES_PER_ROW + 1) * sizeof(int));
		if (base15_symbols[i] == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation failed\n");
			return 1;
		}
		// Initialize array
		memset(base15_symbols[i], 0, (TOTAL_BASES_PER_ROW + 1) * sizeof(int));
	}

	FILE* outfile = fopen(original_bits_file, "w");
    if (outfile == NULL) {
        fprintf(stderr, "ERROR: Cannot create output file '%s'\n", original_bits_file);
        fclose(fw_oligopool);
        return 1;
    }

    FILE* outfile2 = fopen(rs_codeword_file, "w");
    if (outfile2 == NULL) {
        fprintf(stderr, "ERROR: Cannot create output file '%s'\n", rs_codeword_file);
        fclose(fw_oligopool);
        fclose(outfile);
        return 1;
    }
	
	g = (int*)malloc(600 * sizeof(int)); // Parity bits
	rdata = (int*)malloc(4095 * sizeof(int));
	index_of = (int*)malloc(4095 * sizeof(int));
	alpha_to = (int*)malloc(4095 * sizeof(int));
	recd = (int*)malloc(4095 * sizeof(int));
	b = (int*)malloc(600 * sizeof(int));
	p = (int*)malloc(600 * sizeof(int));

	m = rsCodeLen; // 12
	length = 4095; // 4095
	parity_check = rs_n - rs_k; // Parity bits (315)
	init_zero = 0;
	k = length - parity_check; // 3780
	t = parity_check / 2;
	t2 = 2 * t;

	for (int group = 0; group < rsCodeNum * rsCodeLen * rs_k; group++) // 1*12*1575
	{
		Pre_RS_array[group] = buff_scramble_char[group] - '0';
	}

	read_p();        // Read m
	generate_gf();   // Construct the Galois Field GF(2^m)
	gen_poly();      // Compute the generator polynomial of RS code
	//=======================================================================//


	for (i = 0; i < rs_k; i++) // Convert rows to matrix
	{
		for (j = 0; j < (rsCodeLen * rsCodeNum); j++)
		{

			// Assign values to Pre_RS_matrix based on Pre_RS_array values
			Pre_RS_matrix[i][j] = (Pre_RS_array[i * (rsCodeLen * rsCodeNum) + j] == 0) ? 0 : 1;
			fprintf(outfile, "%d", Pre_RS_matrix[i][j]);

		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);

	for (int jj = 0; jj < rs_k; jj++) // Convert bits to symbols
	{
		for (int aa = 0; aa < rsCodeNum; aa++) // Convert bits to symbols
		{
			unsigned long int sum12 = 0;
			for (int bb = 0; bb < 12; bb++)
			{
				sum12 += (int)(pow(2, 11 - bb) * Pre_RS_matrix[jj][aa * 12 + bb]);  // Convert every 12 bits to one decimal number
			}
			Pre_RS_Symbol_matrix[jj][aa] = sum12;

		}
	}
	

	for (int cc = 0; cc < rsCodeNum; cc++)
	{

		for (int dd1 = 0; dd1 < 2205; dd1++)
		{
			rdata[dd1] = 0;  // Pad first (4095-1890=2205) bits of RS code with 0
		}
		for (int dd2 = 0; dd2 < rs_k; dd2++)
		{
			rdata[dd2 + 2205] = Pre_RS_Symbol_matrix[dd2][cc];  // Remaining rs_k information bits
		}
		

		encode_rs();    // Encode data - RS encoding 
		
		for (ri = 0; ri < rs_k; ri++)
		{
			recd[ri] = rdata[ri + 2205];
		}
		for (ri = 0; ri < parity_check; ri++) // 315 Check bits
		{
			recd[ri + rs_k] = b[ri]; // Parity bits at the end
			// printf("%d ", b[ri]);
		}
		// printf("\n");

		for (int ss = 0; ss < rs_n; ss++)
		{

			Multi_RS_Codes[cc][ss] = recd[ss];
			fprintf(outfile2, "%d ", Multi_RS_Codes[cc][ss]);

		}
		fprintf(outfile2, "\n");
		
		// RS code verification
		int CodeCheck[4096] = { 0 };
		for (int i1 = 0; i1 < 3780; i1++)
		{
			CodeCheck[i1] = *(rdata + i1);
			//printf("%d ", code[i1]);
		}
		for (int i1 = 0; i1 < parity_check; i1++)
		{
			CodeCheck[i1 + 3780] = *(b + i1);
			//printf("%d ", code[i1+251]);
		}
		
		int syn_error = 0;
		int s[5500] = { 0 };
		for (int mm = 1; mm <= t2; mm++)
		{
			s[mm] = 0;
			for (int kk = 0; kk < length; kk++)
				if (CodeCheck[kk] != 0)
					s[mm] ^= alpha_to[(index_of[CodeCheck[kk]] + (mm + init_zero - 1) * kk) % n];
			// convert syndrome from vector form to log form
			if (s[mm] != 0)
				syn_error = 1;         // set flag if non-zero syndrome => error
			//
			// Note:    If the code is used only for ERROR DETECTION, then
			//          exit program here indicating the presence of errors.
			//
			s[mm] = index_of[s[mm]];
		}
		if (syn_error == 0) {
			printf("INFO: RS code num %d, encoding successful, check result %4d\n", cc, syn_error);
		} else {
			printf("INFO: RS code num %d, check result %4d\n", cc, syn_error);
		}
		//
	}

	fclose(outfile2);

	// Call conversion function
	for (int xx = 0; xx < rsCodeNum; xx++)
	{

		for (int yy = 0; yy < rs_n; yy++)
		{
			// First initialize corresponding bits to 0
			for (int pos = 0; pos < 12; pos++)
			{
				rsEncodedBits[yy][12 * xx + pos] = 0;
			}

			int value_rs = Multi_RS_Codes[xx][yy];

			for (int xxx = 0; xxx < 12; xxx++)
			{
				rsEncodedBits[yy][12 * xx + 11 - xxx] = (value_rs >> xxx) & 0x01;

			}

			//int temp_rs = 0;
			//int pos_rs = 0;
			//do
			//{
			//	temp_rs = value_rs % 2;
			//	rsEncodedBits[yy][12 * xx + 11 - pos_rs] = temp_rs;
			//	value_rs = (value_rs - temp_rs) / 2;
			//	pos_rs++;
			//} while (value_rs != 0); // Convert RS code back to binary, write in 12 columns each, get 1890*12 bits
		}

	}
	FILE* rs_bits = fopen(rs_bits_file, "w");
    if (rs_bits == NULL) {
        fprintf(stderr, "ERROR: Cannot create output file '%s'\n", rs_bits_file);
        return 1;
    }
    printf("INFO: RS bits output file '%s' opened for writing\n", rs_bits_file);
	
	for (int yy = 0; yy < rs_n; yy++)
	{

		for (int xxx = 0; xxx < rsCodeNum*rsCodeLen; xxx++)
		{
			fprintf(rs_bits,"%d", rsEncodedBits[yy][xxx]);

		}
		
	}
	fclose(rs_bits);

	/*****************************Composite Base Encoding********************************/
	// Create a 126x180 array for the composite base encoding
	int** compositeBase = (int**)malloc(sizeof(int*) * 126);
	for (int i = 0; i < 126; i++) 
	{
		compositeBase[i] = (int*)malloc(sizeof(int) * 180);
		memset(compositeBase[i], 0, sizeof(int) * 180); // Initialize to zeros
	}

	// Reshape the data from rsEncodedBits (1890x12) to compositeBase (126x180)
	int idx = 0;
	for (int row = 0; row < 126; row++) 
	{
		for (int col = 0; col < 180; col++) 
		{
			int srcRow = idx / 12;
			int srcCol = idx % 12;
			if (srcRow < rs_n) 
			{
				compositeBase[row][col] = rsEncodedBits[srcRow][srcCol];
			}
			idx++;
		}
	}

	printf("INFO: Composite base projection completed\n");
	/*****************************Composite Base Encoding********************************/

    FILE* payloadfile = fopen(payload_file, "w");
    if (payloadfile == NULL) {
        fprintf(stderr, "ERROR: Cannot create payload file '%s'\n", payload_file);
        return 1;
    }
    printf("INFO: Payload file '%s' opened for writing\n", payload_file);

    FILE* left_indexfile;
    FILE* right_indexfile;
    
    if ((left_indexfile = fopen(left_index_file, "rt")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open left index file '%s'\n", left_index_file);
        fclose(payloadfile);
        fclose(fw_oligopool);
        return 1;
    }

    if ((right_indexfile = fopen(right_index_file, "rt")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open right index file '%s'\n", right_index_file);
        fclose(left_indexfile);
        fclose(payloadfile);
        fclose(fw_oligopool);
        return 1;
    }

    // Primer sequences (can be made configurable in future versions)
    char primer1[21] = "AGCCTTGTGTCCATCAATCC";
    char primer2[21] = "TGCGCTATGGTTTGGCTAAT";

	char sequence[126][60] = {0};
	char left_indextemp[10] = {0};
	char right_indextemp[10] = {0};

	for(int item=0;item<126;item++)
	{
		int code[180] = { 0 };
		for ( j = 0; j < 180; j++)
		{
			code[j] = compositeBase[item][j];
		}
		// for (j = 0; j < 180; j++)
		// {
		// 	printf("%d", code[j]);
		// }
		// printf("\n");
		
		for ( j = 0; j < 60; j++)
		{
			
			if(  *(code +j*3)  ==0   &&     *(code +j*3 +1)  ==0  &&     *(code +j*3 +2)  ==0   )
			{
				  sequence[item][j] = 'A';
			}

			if(  *(code +j*3)  ==0   &&     *(code +j*3 +1 )  ==0  &&     *(code +j*3 +2)  ==1   )
			{
				  sequence[item][j] = 'T';
			}
			
			if(  *(code +j*3)  ==0   &&     *(code +j*3 +1 )  ==1  &&     *(code +j*3 +2 )  ==0   )
			{
				  sequence[item][j] = 'G';
			}			   
		
			if(  *(code +j*3)  ==0   &&     *(code +j*3 +1)  ==1  &&     *(code +j*3 +2)  ==1   )
			{
				  sequence[item][j] = 'C';
			}					   


			if(  *(code +j*3)  ==1   &&     *(code +j*3 +1 )  ==0  &&     *(code +j*3 +2 )  ==0   )
			{
				  sequence[item][j] = 'R';
			}

			if(  *(code +j*3)  ==1   &&     *(code +j*3 +1)  ==0  &&     *(code +j*3 +2 )  ==1   )
			{
				  sequence[item][j] = 'Y';
			}
			
			if(  *(code +j*3)  ==1   &&     *(code +j*3 +1 )  ==1  &&     *(code +j*3 +2 )  ==0   )
			{
				  sequence[item][j] = 'M';
			}			   
		
			if(  *(code +j*3)  ==1   &&     *(code +j*3 +1 )  ==1  &&     *(code +j*3 +2 )  ==1   )
			{
				  sequence[item][j] = 'K';
			}
		}
		int aaa_left = fscanf(left_indexfile, "%s\n", left_indextemp);
		int aaa_right = fscanf(right_indexfile, "%s\n", right_indextemp);
		
		// printf("%d Left: %s, Right: %s\n", item, left_indextemp, right_indextemp);
		
		for(j=0;j<20;j++)
        {
            fputc(primer1[j],fw_oligopool);
        }
        for(j=0;j<8;j++)
        {
            fputc(left_indextemp[j],fw_oligopool);
        }

        for(j=0;j<60;j++)
        {
            fputc(sequence[item][j],fw_oligopool);
        }		
        for (j=0;j<60;j++)
        {
            fputc(sequence[item][j],payloadfile);
        }
        fputc('\n',payloadfile);
        
        for(j=0;j<8;j++)
        {
            fputc(right_indextemp[j],fw_oligopool);
        }

        for(j=0;j<20;j++)
        {
            fputc(primer2[j],fw_oligopool);
        }
        fputc('\n',fw_oligopool);
	}
	fclose(payloadfile);
	fclose(left_indexfile);
	fclose(right_indexfile);
	fclose(fw_oligopool);

    printf("SUCCESS: RS encoding completed successfully\n");
    printf("INFO: All output files saved to directory '%s'\n", output_dir);
    printf("INFO: Generated files:\n");
    printf("  - %s (composite DNA sequences)\n", oligopool_file);
    printf("  - %s (original bits)\n", original_bits_file);
    printf("  - %s (RS codewords)\n", rs_codeword_file);
    printf("  - %s (RS bits)\n", rs_bits_file);
    printf("  - %s (payload sequences)\n", payload_file);

	for (int i = 0; i < 126; i++) {
		free(compositeBase[i]);
	}
	free(compositeBase);


	for (int i = 0; i < rs_n; i++)
	{
		free(base15_symbols[i]);
	}
	free(base15_symbols);


	for (int ll0 = 0; ll0 < rs_k; ll0++)
		free(Pre_RS_matrix[ll0]);
	free(Pre_RS_matrix);

	for (int rr2 = 0; rr2 < rs_k; rr2++)
		free(Pre_RS_Symbol_matrix[rr2]);
	free(Pre_RS_Symbol_matrix);

	for (int ii1 = 0; ii1 < rsCodeNum; ii1++)
		free(Multi_RS_Codes[ii1]);
	free(Multi_RS_Codes);


	for (int ll1 = 0; ll1 < rs_n; ll1++)
		free(rsEncodedBits[ll1]);
	free(rsEncodedBits);


	free(b);
	free(g);
	free(p);
	free(alpha_to);


	return 0;
	
}




void read_p()
//      Read m, the degree of a primitive polynomial p(x) used to compute the
//      Galois field GF(2**m). Get precomputed coefficients p[] of p(x). Read
//      the code length.
{
	int i, ninf = 0;
	/*
	  printf("\nSimulation of RS codes \n");
	  printf("Copyright 2002 (c)Robert Morelos-Zaragoza. All rights reserved.\n\n");*/

	for (i = 1; i < m; i++)
		p[i] = 0;
	p[0] = p[m] = 1;

	if (m == 2)             p[1] = 1;
	else if (m == 3)        p[1] = 1;
	else if (m == 4)        p[3] = 1;
	// else if (m == 4)        p[1] = 1;  // Commented out to match example p. 68
	else if (m == 5)        p[2] = 1;
	else if (m == 6)        p[1] = 1;
	else if (m == 7)        p[1] = 1;
	else if (m == 8)        p[4] = p[5] = p[6] = 1;
	else if (m == 9)        p[4] = 1;
	else if (m == 10)       p[3] = 1;
	else if (m == 11)       p[2] = 1;
	else if (m == 12)       p[3] = p[4] = p[7] = 1;
	else if (m == 13)       p[1] = p[3] = p[4] = 1;
	else if (m == 14)       p[1] = p[11] = p[12] = 1;
	else if (m == 15)       p[1] = 1;
	else if (m == 16)       p[2] = p[3] = p[5] = 1;
	else if (m == 17)       p[3] = 1;
	else if (m == 18)       p[7] = 1;
	else if (m == 19)       p[1] = p[5] = p[6] = 1;
	else if (m == 20)       p[3] = 1;
	/*  printf("Primitive polynomial of GF(2^%d), (LSB first)   p(x) = ",m);*/

	n = 1;
	for (i = 0; i <= m; i++)
	{
		n *= 2;
		//     printf("%1d", p[i]);
	}
	//  printf("\n");

	n = n / 2 - 1;
}

void generate_gf()
// generate GF(2^m) from the irreducible polynomial p(X) in p[0]..p[m]
//
// lookup tables:  log->vector form           alpha_to[] contains j=alpha**i;
//                 vector form -> log form    index_of[j=alpha**i] = i
// alpha=2 is the primitive element of GF(2^m)
{

	register int i, mask;
	mask = 1;
	alpha_to[m] = 0;
	for (i = 0; i < m; i++)
	{
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;
		if (p[i] != 0)
			alpha_to[m] ^= mask;
		mask <<= 1;
	}

	index_of[alpha_to[m]] = m;
	mask >>= 1;
	for (i = m + 1; i < n; i++)
	{
		if (alpha_to[i - 1] >= mask)
			alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
		else alpha_to[i] = alpha_to[i - 1] << 1;
		index_of[alpha_to[i]] = i;
	}
	index_of[0] = -1;
	//#define PRINT_GF

#ifdef PRINT_GF/*
printf("Table of GF(%d):\n",n+1);
printf("----------------------\n");
printf("   i\tvector \tlog\n");
printf("----------------------\n");
for (i=0; i<=n; i++)
printf("%4d\t%4x\t%4d\n", i, alpha_to[i], index_of[i]);*/
#endif
}





void gen_poly()
// Compute the generator polynomial of the t-error correcting, length
// n=(2^m -1) Reed-Solomon code from the product of (X+alpha^i), for
// i = init_zero, init_zero + 1, ..., init_zero+length-k-1
{
	register int i, j;


	g[0] = alpha_to[init_zero];  //  <--- vector form of alpha^init_zero
	g[1] = 1;     // g(x) = (X+alpha^init_zero)
	for (i = 2; i <= length - k; i++)
	{
		g[i] = 1;
		for (j = i - 1; j > 0; j--)
			if (g[j] != 0)
				g[j] = g[j - 1] ^ alpha_to[(index_of[g[j]] + i + init_zero - 1) % n];
			else
				g[j] = g[j - 1];
		g[0] = alpha_to[(index_of[g[0]] + i + init_zero - 1) % n];
	}

	// convert g[] to log form for quicker encoding 
	for (i = 0; i <= length - k; i++)  g[i] = index_of[g[i]];


	//#define PRINT_POLY
#ifdef PRINT_POLY
//printf("Generator polynomial (independent term first):\ng(x) = ");
for (i=0; i<=length-k; i++) printf("%5d", g[i]);
printf("\n");
#endif
}





void encode_rs()
// Compute the 2t parity symbols in b[0]..b[2*t-1]
// data[] is input and b[] is output in polynomial form.
// Encoding is done by using a feedback shift register with connections
// specified by the elements of g[].

{
	register int i, j;
	int feedback;
	
	for (i = 0; i < length - k; i++)
		b[i] = 0;
	
	for (i = k - 1; i >= 0; i--)
	{
		feedback = index_of[rdata[i] ^ b[length - k - 1]];
		if (feedback != -1)
		{
			for (j = length - k - 1; j > 0; j--)
				if (g[j] != -1)
					b[j] = b[j - 1] ^ alpha_to[(g[j] + feedback) % n];
				else
					b[j] = b[j - 1];
			b[0] = alpha_to[(g[0] + feedback) % n];
		}
		else
		{
			for (j = length - k - 1; j > 0; j--)
				b[j] = b[j - 1];
			b[0] = 0;
		}
	}
}





static double r[98];
static int iff,ix1,ix2,ix3;

double ran1(int& idum)
{
    double rm1,rm2,t;
	int m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3,j;
	 
    m1 = 259200, ia1 = 7141, ic1 = 54773, rm1 = 0.0000038580247,
    m2 = 134456, ia2 = 8121, ic2 = 28411, rm2 = 0.0000074373773,
    m3 = 243000, ia3= 4561, ic3 = 51349;
	
    if ((idum < 0) || (iff == 0))
	{
        iff = 1;
        ix1 =  (ic1 - idum) % m1;
        ix1 = (ia1 * ix1 + ic1) % m1;
        ix2 = ix1 % m2;
        ix1 = (ia1 * ix1 + ic1) % m1;
        ix3 = ix1 % m3;
        for (j = 1; j<=97; j++)
		{
            ix1 = (ia1 * ix1 + ic1)% m1;
            ix2 = (ia2 * ix2 + ic2)% m2;
            r[j] = (double(ix1) + double(ix2) * rm2) * rm1;
        }
        idum = 1;
	}
      ix1 = (ia1 * ix1 + ic1) % m1;
      ix2 = (ia2 * ix2 + ic2) % m2;
      ix3 = (ia3 * ix3 + ic3) % m3;
      j = 1 + int((97 * ix3) / m3);
      if ((j > 97) || (j < 1)) 
	  {
		  cout<< "abnormal exit"<<endl;
	      return 1;
	  }
      t = r[j];
      r[j] = (float(ix1) + float(ix2) * rm2) * rm1;
	  return t;
}





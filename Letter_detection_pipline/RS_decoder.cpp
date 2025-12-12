#define _CRT_SECURE_NO_WARNINGS
#include "RS_decoder.h"

using namespace std;

// Global variable definitions
int ri, t, m, n, length, k, t2, d, parity_len, init_zero;
int* p;
int numerr, errpos[5000], errval[5000], decerror;
int biterror, error;
int numera = 0;
int era[5000] = {}, eraval[5000] = {};
int* g, * rdata, * index_of, * alpha_to, * recd, * b;

// Utility functions
char dexor(char a, char b)
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
void Dec2Bin(int num10, int* output)
{
    int arr[8] = { 0 };
    int i;


    for (i = 8 - 1; i >= 0; i--)  
    {
        arr[i] = num10 % 2;
        num10 = num10 / 2;
    }


    for (i = 0; i <= 8 - 1; i++)
    {

        output[i] = arr[i];

    }

}
void mapCharacterToBits(char c, int rule, char* bits) 
{
    if (rule == 1) 
    {
        // Rule 1 mapping
        switch (c) 
        {
            case 'A': strcpy(bits, "000"); break;
            case 'T': strcpy(bits, "001"); break;
            case 'G': strcpy(bits, "010"); break;
            case 'C': strcpy(bits, "011"); break;
            case 'R': strcpy(bits, "100"); break;
            case 'Y': strcpy(bits, "101"); break;
            case 'M': strcpy(bits, "110"); break;
            case 'K': strcpy(bits, "111"); break;
            case 'E': strcpy(bits, "EEE"); break;
            default: strcpy(bits, "000"); break;  // Default case
        }
    } 
    else if (rule == 2) 
    {
        // Rule 2 mapping
        switch (c) 
        {
            case 'A': strcpy(bits, "000"); break;
            case 'T': strcpy(bits, "001"); break;
            case 'S': strcpy(bits, "010"); break;
            case 'W': strcpy(bits, "011"); break;
            case 'R': strcpy(bits, "100"); break;
            case 'Y': strcpy(bits, "101"); break;
            case 'M': strcpy(bits, "110"); break;
            case 'K': strcpy(bits, "111"); break;
            case 'E': strcpy(bits, "EEE"); break;
            default: strcpy(bits, "000"); break;  // Default case
        }
    } 
    else if (rule == 3) 
    {
        // Rule 3 mapping
        switch (c) 
        {
            case 'A': strcpy(bits, "000"); break;
            case 'T': strcpy(bits, "001"); break;
            case 'G': strcpy(bits, "010"); break;
            case 'C': strcpy(bits, "011"); break;
            case 'H': strcpy(bits, "100"); break;
            case 'B': strcpy(bits, "101"); break;
            case 'V': strcpy(bits, "110"); break;
            case 'D': strcpy(bits, "111"); break;
            case 'E': strcpy(bits, "EEE"); break;
            default: strcpy(bits, "000"); break;  // Default case
        }
    }
}

// Memory management module
void init_rs_parameters()
{
    g = (int*)malloc(600 * sizeof(int));
    rdata = (int*)malloc(4096 * sizeof(int));
    index_of = (int*)malloc(4096 * sizeof(int));
    alpha_to = (int*)malloc(4096 * sizeof(int));
    recd = (int*)malloc(4096 * sizeof(int));
    b = (int*)malloc(600 * sizeof(int));
    p = (int*)malloc(600 * sizeof(int));

    m = rsCodeLen;
    length = 4095;
    parity_len = rs_n - rs_k;
    init_zero = 0;
    k = length - parity_len;
    t = parity_len / 2;
    t2 = 2 * t;

    read_p();
    generate_gf();
    gen_poly();
}

void cleanup_rs_memory()
{
    free(g);
    free(rdata);
    free(index_of);
    free(alpha_to);
    free(recd);
    free(b);
    free(p);
}

int** allocate_2d_int_array(int rows, int cols)
{
    int** array = (int**)malloc(sizeof(int*) * rows);
    for (int i = 0; i < rows; i++) {
        array[i] = (int*)malloc(sizeof(int) * cols);
        memset(array[i], 0, sizeof(int) * cols);
    }
    return array;
}

void free_2d_int_array(int** array, int rows)
{
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}

char** allocate_2d_char_array(int rows, int cols)
{
    char** array = (char**)malloc(sizeof(char*) * rows);
    for (int i = 0; i < rows; i++) {
        array[i] = (char*)malloc(sizeof(char) * cols);
    }
    return array;
}

void free_2d_char_array(char** array, int rows)
{
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}

// Command line parameter parsing module
int parse_command_line(int args, char** argv, char* noisy_composite_dna, 
                      char* output_bit_stream)
{
    if (args < 4) {
        printf("Usage: %s <noisy_composite_dna> <output_bit_stream> <rule>\n", argv[0]);
        printf("  noisy_composite_dna: Input DNA sequence file\n");
        printf("  output_bit_stream: Output decoded bit stream file\n");
        printf("  rule: Mapping rule (1, 2, or 3)\n");
        return -1;
    }

    sscanf(argv[1], "%s", noisy_composite_dna);
    sscanf(argv[2], "%s", output_bit_stream);
    
    int rule = atoi(argv[3]);
    if (rule != 1 && rule != 2 && rule != 3) {
        printf("ERROR: Invalid mapping rule %d. Must be 1, 2, or 3.\n", rule);
        return -1;
    }
    
    printf("INFO: Using mapping rule %d\n", rule);
    return rule;
}

// DNA sequence processing module
int process_dna_sequences(const char* filename, int rule, int bitSequence[rs_n][rsCodeLen], int* ERROR)
{
    FILE* fpQuery = fopen(filename, "rt");
    if (!fpQuery) {
        printf("ERROR: Cannot open DNA sequence file: %s\n", filename);
        return -1;
    }

    printf("INFO: Processing DNA sequences from file: %s\n", filename);

    char readWithPrimer[62] = {0};
    char readWithoutPrimer[62] = {0};
    int reads_num = 0;

    while (reads_num < 126) 
	{
        fscanf(fpQuery, "%s", readWithPrimer);
        
        for (int group = 0; group < 60; group++) {
            readWithoutPrimer[group] = readWithPrimer[group];
        }
        readWithoutPrimer[60] = '\0';
        
        int tempBits[180] = {0};
        
        if (readWithoutPrimer[0] == 'E') 
		{
            for (int row = 0; row < 15; row++) 
			{
                int destRow = reads_num * 15 + row;
                ERROR[destRow] = 1;
                for (int j = 0; j < 12; j++) 
				{
                    bitSequence[destRow][j] = 0;
                }
            }
        } 
		else
		{
            char bits[4] = {0};
            int bit_index = 0;
            for (int j = 0; j < 60; j++) {
                mapCharacterToBits(readWithoutPrimer[j], rule, bits);
                for (int k = 0; k < 3; k++) {
                    tempBits[bit_index] = bits[k] - '0';
                    bit_index++;
                }
            }
            
            for (int row = 0; row < 15; row++) {
                int destRow = reads_num * 15 + row;
                for (int j = 0; j < 12; j++) {
                    bitSequence[destRow][j] = tempBits[row * 12 + j];
                }
            }
        }
        reads_num++;
    }
    
    fclose(fpQuery);
    printf("INFO: Successfully processed %d DNA sequences and converted to RS codewords\n", reads_num);
    return 0;
}

// Bit sequence conversion module
void convert_bits_to_rs_codewords(int bitSequence[rs_n][rsCodeLen], int** rsCodeWithErr, int* ERROR)
{
    for (int i = 0; i < rs_n; i++)
	{
        if (ERROR[i] == 1) 
		{
            for (int j = 0; j < rsCodeNum; j++) 
			{
                rsCodeWithErr[j][i] = 0;
            }
        } 
		else 
		{
            int sym_sum12 = 0;
            for (int g2 = 0; g2 < rsCodeNum; g2++) 
			{
                sym_sum12 = 0;
                for (int gg = 0; gg < 12; gg++) {
                    sym_sum12 += (int)(pow(2, 11 - gg) * bitSequence[i][12 * g2 + gg]);
                }
                rsCodeWithErr[g2][i] = sym_sum12;
            }
        }
    }
}

// Erasure position setup module
void setup_erasure_positions(int* ERROR, int* era, int* numera)
{
    *numera = 0;
    for (int i = 0; i < rs_n; i++) 
    {
        if (ERROR[i] == 1) {
            era[*numera] = i + 2205;
            (*numera)++;
        }
    }
    printf("INFO: Found %d erasure positions\n", *numera);
}

// Reed-Solomon decoding module
int perform_rs_decoding(int** rsCodeWithErr, int** codewordRS, 
                       FILE* fw_bit_output, int* ERROR)
{
    int erasure_num;
    int sequence0[12] = {0};
    int value_rs, pos_rs;

    printf("INFO: Starting Reed-Solomon decoding process\n");

    for (int i = 0; i < rsCodeNum; i++) 
    {
        erasure_num = 0;
        
        // Prepare decoding data
        for (int j = 0; j < 2205; j++) {
            recd[j] = 0;
        }
        for (int j = 0; j < 1575; j++) {
            recd[j + 2205] = rsCodeWithErr[i][j];
        }
        for (int j = 3780; j < 4095; j++) {
            recd[j] = rsCodeWithErr[i][j - 2205];
        }

        // Calculate erasure count
        for (int j = 0; j < rs_n; j++) {
            if(ERROR[j] == 1) {
                erasure_num++;
            }
        }
        
        printf("INFO: Codeword %d - Erasures detected: %d\n", i + 1, erasure_num);

        if (numera >= 315)
        {
            printf("WARNING: Erasure count (%d) exceeds correction capability (315). Skipping decoding.\n", numera);
        }
        else
        {
            printf("INFO: Erasure count (%d) within correction capability. Proceeding with decoding.\n", numera);
            decode_rs();

            for (int j = 0; j < rs_k; j++) {
                codewordRS[j][i] = recd[j + 2205];
            }
            
            printf("INFO: Codeword %d decoding completed successfully\n", i + 1);
            
            // Output decoding results
            for (int jj = 0; jj < rs_k; jj++) 
            {
				memset(sequence0, 0, sizeof(sequence0));
                for (int aa = 0; aa < rsCodeNum; aa++) 
                {
                    value_rs = codewordRS[jj][aa];
                    pos_rs = 11;
                    
                    while (value_rs > 0 && pos_rs >= 0) 
                    {
                        sequence0[12 * aa + pos_rs] = value_rs % 2;
                        value_rs = value_rs / 2;
                        pos_rs--;
                    }
                }

                for (int i = 0; i < 12; i++) 
                {
                    fprintf(fw_bit_output, "%d", sequence0[i]);
                }
                fprintf(fw_bit_output, "\n");
            }
        }
    }
    
    printf("INFO: Reed-Solomon decoding process completed\n");
    return 0;
}

// Main function
int main(int args, char** argv)
{
    printf("Reed-Solomon Decoder for DNA Storage\n");
    printf("====================================\n");

    // Filename variables
    char noisy_composite_dna[200] = {0};
    char output_bit_stream[200] = {0};
    
    // Declare FILE pointer variable
    FILE* fw_bit_output = NULL;
    
    // Parse command line arguments
    int rule = parse_command_line(args, argv, noisy_composite_dna, output_bit_stream);
    if (rule < 0) {
        return 1;
    }

    // Initialize RS parameters and memory
    printf("INFO: Initializing Reed-Solomon parameters\n");
    init_rs_parameters();

    // Allocate memory
    printf("INFO: Allocating memory for data structures\n");
    int** rsCodeWithErr = allocate_2d_int_array(rsCodeNum, rs_n);
    int** codewordRS = allocate_2d_int_array(rs_k, rsCodeNum);
    int** seq_after_bch = allocate_2d_int_array(rs_n, rsCodeLen * rsCodeNum);
    int** sequence11 = allocate_2d_int_array(rs_k, rsCodeLen * rsCodeNum);
    char** original = allocate_2d_char_array(rs_k, rsCodeLen * rsCodeNum);
    int** rsEncodedBits = allocate_2d_int_array(rs_n, TOTAL_BITS_PER_ROW);

    // Error flag array
    int ERROR[4000] = {0};
    for (int x = 0; x < 600; x++) {
        era[x] = 0;
    }
    numera = 0;

    // Process DNA sequences
    int bitSequence[rs_n][rsCodeLen] = {0};
    if (process_dna_sequences(noisy_composite_dna, rule, bitSequence, ERROR) < 0) {
        printf("ERROR: Failed to process DNA sequences\n");
        return 1;
    }

    // Convert bit sequences to RS codewords
    printf("INFO: Converting bit sequences to RS codewords\n");
    convert_bits_to_rs_codewords(bitSequence, rsCodeWithErr, ERROR);

    // Setup erasure positions
    setup_erasure_positions(ERROR, era, &numera);

    // Open output file
    printf("INFO: Opening output file: %s\n", output_bit_stream);
    fw_bit_output = fopen(output_bit_stream, "wt");
    
    if (!fw_bit_output) {
        printf("ERROR: Cannot open output file: %s\n", output_bit_stream);
        return 1;
    }

    // Perform RS decoding
    perform_rs_decoding(rsCodeWithErr, codewordRS, fw_bit_output, ERROR);

    // Close file
    fclose(fw_bit_output);
    printf("INFO: Output written to: %s\n", output_bit_stream);

    // Free memory
    printf("INFO: Cleaning up allocated memory\n");
    free_2d_int_array(rsCodeWithErr, rsCodeNum);
    free_2d_int_array(codewordRS, rs_k);
    free_2d_int_array(seq_after_bch, rs_n);
    free_2d_int_array(sequence11, rs_k);
    free_2d_char_array(original, rs_k);
    free_2d_int_array(rsEncodedBits, rs_n);
    
    cleanup_rs_memory();

    printf("INFO: Reed-Solomon decoding completed successfully\n");
    return 0;
}

* NOTE:    
*          The authors are not responsible for any malfunctioning of
*          this program, nor for any damage caused by it. Please include the
*          original program along with these comments in any redistribution.
*
* Portions of this program are from a Reed-Solomon encoder/decoder
* in C, written by Simon Rockliff (simon@augean.ua.oz.au) on 21/9/89.
*
* COPYRIGHT NOTICE: This computer program is free for non-commercial purposes.
* You may implement this program for any non-commercial application. You may 
* also implement this program for commercial purposes, provided that you
* obtain my written permission. Any modification of this program is covered
* by this copyright.
*
* Copyright (c) 1995.  Robert Morelos-Zaragoza and Hari Thirumoorthy.
*                      All rights reserved.
// Reed-Solomon related function implementations
void read_p(void)
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
// Generate GF(2^m) from the irreducible polynomial p(X) in p[0]..p[m]
//
// Lookup tables:  log->vector form           alpha_to[] contains j=alpha**i;
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

	// Convert g[] to log form for quicker encoding 
	for (i = 0; i <= length - k; i++)  g[i] = index_of[g[i]];


	//#define PRINT_POLY
#ifdef PRINT_POLY
	//printf("Generator polynomial (independent term first):\ng(x) = ");
	/*for (i=0; i<=length-k; i++) printf("%5d", g[i]);
	printf("\n");*/
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

void decode_rs()
{
	register int i, j, u, q;
	int count = 0, syn_error = 0;
	int degphi, ell, temp = 0;


	int* d;
	d = (int*)malloc((600+2) * sizeof(int));
	int* l;
	l = (int*)malloc((600+2) * sizeof(int));
	int* u_lu;
	u_lu = (int*)malloc((600+2) * sizeof(int));
	int* s, * forney, * err, * omega, * phi, * phiprime;
	s = (int*)malloc((600+2) * sizeof(int));
	forney = (int*)malloc((600+2) * sizeof(int));
	err = (int*)malloc((4100+2) * sizeof(int));
	omega = (int*)malloc((600+2) * sizeof(int));
	phi = (int*)malloc((4100+2) * sizeof(int));
	phiprime = (int*)malloc((4100+2) * sizeof(int));
	int* tau, * root, * loc, * z, * reg, * aux;
	tau = (int*)malloc((600+2) * sizeof(int));
	root = (int*)malloc((600+2) * sizeof(int));
	loc = (int*)malloc((600+2) * sizeof(int));
	z = (int*)malloc((600+2) * sizeof(int));
	reg = (int*)malloc((600+2) * sizeof(int));
	aux = (int*)malloc((600+2) * sizeof(int));
	/*
	int **elp = new int*[10536];
		for (i = 0; i < 10536; i++)
			elp[i] = new int[65538];

	*/


	int** elp;
	elp = (int**)malloc(sizeof(int*) * (600+2) );
	for (i = 0; i < (600+2) ; i++)
		elp[i] = (int*)malloc(sizeof(int) * (600+2) );

	// Compute the syndromes

 //#ifdef PRINT_SYNDROME
 //   printf("\ns =         0 ");
 //#endif
	for (i = 1; i <= t2; i++)
	{
		s[i] = 0;
		for (j = 0; j < length; j++)
			if (recd[j] != 0)
				s[i] ^= alpha_to[(index_of[recd[j]] + (i + init_zero - 1) * j) % n];
		// convert syndrome from vector form to log form  */
		if (s[i] != 0)
			syn_error = 1;         // set flag if non-zero syndrome => error
		  //
		  // Note:    If the code is used only for ERROR DETECTION, then
		  //          exit program here indicating the presence of errors.
		  //
		s[i] = index_of[s[i]];
		//#define PRINT_SYNDROME
#ifdef PRINT_SYNDROME
//   printf("%4d ", s[i]);
#endif
	}
	//printf("\nsrr =         0 ");
    printf("INFO: Syndrome check result: %s\n", syn_error ? "Errors detected" : "No errors detected");
	if (syn_error)       // If syndromes are nonzero then try to correct
	{

		s[0] = 0;  // S(x) = 1 + s_1x + ...

		 // TO HANDLE ERASURES

		if (numera)
			// If erasures are present, compute the erasure locator polynomial, tau(x)
		{
			for (i = 0; i <= t2; i++)
			{
				tau[i] = 0; aux[i] = 0;
			}

			aux[1] = alpha_to[era[0]];
			aux[0] = 1;       // (X + era[0]) 

			if (numera > 1)
				for (i = 0; i < numera; i++)
				{
					p[1] = era[i];
					p[0] = 0;
					for (j = 0; j < 2; j++)
						for (ell = 0; ell <= i; ell++)
							// Line below added 8/17/2003
							if ((p[j] != -1) && (aux[ell] != 0))
								tau[j + ell] ^= alpha_to[(p[j] + index_of[aux[ell]]) % n]; 
								// printf("%4d  %4d\n", ell, tau[ell]);
					if (i != (numera - 1))
						for (ell = 0; ell <= (i + 1); ell++)
						{
							aux[ell] = tau[ell];
							tau[ell] = 0;
						}
				}

			else {
				tau[0] = aux[0]; tau[1] = aux[1];
			}


			// Put in index (log) form
			for (i = 0; i <= numera; i++)
				tau[i] = index_of[tau[i]]; /* tau in log form */


#ifdef PRINT_SYNDROME
//printf("\ntau =    ");
//for (i=0; i<=numera; i++)
 // printf("%4d ", tau[i]);
//printf("\nforney = ");
#endif


		// Compute FORNEY modified syndrome:
		//            forney(x) = [ s(x) tau(x) + 1 ] mod x^{t2}

			for (i = 0; i <= n - k; i++) forney[i] = 0;
			for (i = 0; i <= n - k; i++)
				for (j = 0; j <= numera; j++)
					if (i + j <= (n - k)) // mod x^{n-k+1}
						if ((s[i] != -1) && (tau[j] != -1))
							forney[i + j] ^= alpha_to[(s[i] + tau[j]) % n];

			forney[0] ^= 1;

			for (i = 0; i <= n - k; i++)
				forney[i] = index_of[forney[i]];

#ifdef PRINT_SYNDROME
			/*for (i=0; i<=n-k; i++)
			  printf("%4d ", forney[i]);*/
#endif

		}

		else // No erasures
		{
			tau[0] = 0;
			for (i = 1; i <= n - k; i++) forney[i] = s[i];
		}

#ifdef PRINT_SYNDROME
		//printf("\n");
#endif

  // --------------------------------------------------------------
  //    THE BERLEKAMP-MASSEY ALGORITHM FOR ERRORS AND ERASURES
  // --------------------------------------------------------------

	  // Initialize table entries
		d[0] = 0;                // log form
		d[1] = forney[numera + 1]; // log form
		elp[0][0] = 0;           // log form 
		elp[1][0] = 1;           // vector form 
		for (i = 1; i < t2; i++)
		{
			elp[0][i] = -1;   // log form
			elp[1][i] = 0;    // vector form 
		}
		l[0] = 0;
		l[1] = 0;
		u_lu[0] = -1;
		u_lu[1] = 0;
		u = 0;

		if (numera < t2) {  // If errors can be corrected
			do
			{
				u++;

				if (d[u] == -1)
				{
#ifdef PRINT_SYNDROME
					//printf("d[%d] is zero\n",u);
#endif


					l[u + 1] = l[u];

					for (i = 0; i <= l[u]; i++)
					{
						//	printf("elp[u][i]=%d,index_of[elp[u][i]]=%d",elp[u][i],index_of[elp[u][i]]);


						elp[u + 1][i] = elp[u][i];

						elp[u][i] = index_of[elp[u][i]];

					}

				}
				else
					// search for words with greatest u_lu[q] for which d[q]!=0
				{
					q = u - 1;
					while ((d[q] == -1) && (q > 0))    q--;
					// have found first non-zero d[q] 
					if (q > 0)
					{
						j = q;
						do
						{
							j--;
							if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
								q = j;
						} while (j > 0);
					}

#ifdef PRINT_SYNDROME
					//printf("u = %4d, q = %4d, d[q] = %4d d[u] = %4d\n", u, q, d[q],d[u]);
#endif

			// have now found q such that d[u]!=0 and u_lu[q] is maximum 
			// store degree of new elp polynomial 
					if (l[u] > l[q] + u - q)  l[u + 1] = l[u];
					else  l[u + 1] = l[q] + u - q;

#ifdef PRINT_SYNDROME
					//printf("l[q] = %4d, l[u] = %4d\n", l[q], l[u]);
#endif

			// compute new elp(x) 
			// for (i=0; i<t2-numera; i++)    elp[u+1][i] = 0; 
					for (i = 0; i < t2; i++)    elp[u + 1][i] = 0;
					for (i = 0; i <= l[q]; i++)
						if (elp[q][i] != -1)
							elp[u + 1][i + u - q] = alpha_to[(d[u] + n - d[q] + elp[q][i]) % n];
					for (i = 0; i <= l[u]; i++)
					{
						elp[u + 1][i] ^= elp[u][i];
						elp[u][i] = index_of[elp[u][i]];
					}

#ifdef PRINT_SYNDROME
					//printf("l[u+1] = %4d, elp[u+1] = ", l[u+1]);
					/*for (i=0;  i<=l[u+1]; i++) printf("%4d ",index_of[elp[u+1][i]]);
					printf("\n");*/
#endif

				}


				u_lu[u + 1] = u - l[u + 1];
				// compute (u+1)th discrepancy 
				if (u < (t2 - numera)) // no discrepancy computed on last iteration 
				// --- if ( u < (l[u+1]+t-1-(numera/2)) ) 
				{
					// if (s[u+1]!=-1)
					if (forney[numera + u + 1] != -1)
						d[u + 1] = alpha_to[forney[numera + u + 1]];
					else
						d[u + 1] = 0;
#ifdef PRINT_SYNDROME
					//printf("discrepancy for u = %d: d[u+1] = %4d\n", u, index_of[d[u+1]]);
#endif

					for (i = 1; i <= l[u + 1]; i++)
						// if ((s[u+1-i]!=-1) && (elp[u+1][i]!=0))
						//   d[u+1] ^= alpha_to[(s[numera+u+1-i]
						if ((forney[numera + u + 1 - i] != -1) && (elp[u + 1][i] != 0)) {
							d[u + 1] ^= alpha_to[(forney[numera + u + 1 - i]
								+ index_of[elp[u + 1][i]]) % n];

#ifdef PRINT_SYNDROME
							//printf("i=%d, forney[%d] = %4d, d[u+1] = %4d\n",i,numera+u+1-i,
											 //  forney[numera+u+1-i],index_of[d[u+1]]);
#endif
						}
					d[u + 1] = index_of[d[u + 1]];     // put d[u+1] into index form 

#ifdef PRINT_SYNDROME
//printf("d[u+1] = %4d\n", d[u+1]);
#endif
				}
				//	printf("u=%d,(t2-numera)=%d,l[u+1]=%d,(t2-numera)/2)=%d\n",u,t2-numera,l[u+1],(t2-numera)/2);
			} while ((u < (t2 - numera)) && (l[u + 1] <= ((t2 - numera) / 2)));


		}

		// else // case of 2t erasures
		// {
		//   elp[1][0] = 0;
		//   count = 0;
		// }

		u++;
		//     printf("%d %d %d %d",l[u],t-numera/2,t,numera);
		if (l[u] <= t - numera / 2)         // can correct errors
		{

			// put elp into index form 
			for (i = 0; i <= l[u]; i++)   elp[u][i] = index_of[elp[u][i]];
			/*
			printf("\nBM algorithm, after %d iterations:\nsigma = ", (u-1));
			for (i=0; i<=l[u]; i++)   printf("%4d ", elp[u][i]);
			printf("\n");*/

			// find roots of the error location polynomial 
			for (i = 1; i <= l[u]; i++)
				reg[i] = elp[u][i];
			count = 0;
			for (i = 1; i <= n; i++)
			{
				q = 1;
				for (j = 1; j <= l[u]; j++)
					if (reg[j] != -1)
					{
						reg[j] = (reg[j] + j) % n;
						q ^= alpha_to[reg[j]];
					}
				if (!q)        // store root and error location number indices 
				{
					root[count] = i;
					loc[count] = n - i;

#ifdef PRINT_SYNDROME
					//printf("loc[%4d] = %4d\n", count, loc[count]);
#endif

					count++;
				}
			}

			if (count == l[u])    // no. roots = degree of elp hence <= t errors 
			{

				// Compute the errata evaluator polynomial, omega(x)

				forney[0] = 0;  // as a log, to construct 1+T(x)
				for (i = 0; i <= t2; i++)
					omega[i] = 0;
				for (i = 0; i <= t2; i++)
				{
					for (j = 0; j <= l[u]; j++)\
					{
						if (i + j <= t2) // mod x^{t2}
							if ((forney[i] != -1) && (elp[u][j] != -1))
								omega[i + j] ^= alpha_to[(forney[i] + elp[u][j]) % n];
					}
				}
				for (i = 0; i <= t2; i++)
					omega[i] = index_of[omega[i]];

#ifdef PRINT_SYNDROME
				/*printf("\nomega =    ");
				for (i=0; i<=t2; i++)   printf("%4d ", omega[i]);
				printf("\n");*/
#endif

				// Compute the errata locator polynomial, phi(x)

				degphi = numera + l[u];
				for (i = 0; i <= degphi; i++) phi[i] = 0;
				for (i = 0; i <= numera; i++)
					for (j = 0; j <= l[u]; j++)
						if ((tau[i] != -1) && (elp[u][j] != -1))
							phi[i + j] ^= alpha_to[(tau[i] + elp[u][j]) % n];
				for (i = 0; i <= degphi; i++)
					phi[i] = index_of[phi[i]];


#ifdef PRINT_SYNDROME
				/*printf("phi =      ");
				for (i=0; i<=degphi; i++)   printf("%4d ", phi[i]);
				printf("\n");*/
#endif

				// Compute the "derivative" of phi(x): phiprime

				for (i = 0; i <= degphi; i++) phiprime[i] = -1; // as a log
				for (i = 0; i <= degphi; i++)
					if (i % 2)  // Odd powers of phi(x) give terms in phiprime(x)
						phiprime[i - 1] = phi[i];

#ifdef PRINT_SYNDROME
				/*printf("phiprime = ");
				for (i=0; i<=degphi; i++)   printf("%4d ", phiprime[i]);
				printf("\n\n");*/
#endif

				//printf("**** count = %d\n", count);

				if (numera)
					// Add erasure positions to error locations
					for (i = 0; i < numera; i++) {
						loc[count + i] = era[i];
						root[count + i] = (n - era[i]) % n;
					}

				// evaluate errors at locations given by errata locations, loc[i]
				//                             for (i=0; i<l[u]; i++)    
				for (i = 0; i < degphi; i++)
				{
					// compute numerator of error term  
					err[loc[i]] = 0;
					for (j = 0; j <= t2; j++)
						if ((omega[j] != -1) && (root[i] != -1))
							err[loc[i]] ^= alpha_to[(omega[j] + j * root[i]) % n];

					// -------  The term loc[i]^{2-init_zero}
					if ((err[loc[i]] != 0) && (loc[i] != -1))
						err[loc[i]] = alpha_to[(index_of[err[loc[i]]]
							+ loc[i] * (2 - init_zero + n)) % n];
					if (err[loc[i]] != 0)
					{
						err[loc[i]] = index_of[err[loc[i]]];
						// compute denominator of error term 
						q = 0;
						for (j = 0; j <= degphi; j++)
							if ((phiprime[j] != -1) && (root[i] != -1))
								q ^= alpha_to[(phiprime[j] + j * root[i]) % n];

						// Division by q
						err[loc[i]] = alpha_to[(err[loc[i]] - index_of[q] + n) % n];

#ifdef PRINT_SYNDROME
						//printf("errata[%4d] = %4d (%4d) \n",loc[i],index_of[err[loc[i]]],err[loc[i]]);
#endif

						recd[loc[i]] ^= err[loc[i]];
					}
				}

				//    printf("\nCor =");
		  //         for (i=0; i<length; i++) {
		  //             printf("%4d ", index_of[recd[i]]);
		  //            }
		  //          printf("\n     ");
			 /*       for (i=0; i<length; i++) {
					   printf("%4d ", recd[i]);
					   }
					printf("\n");*/


			}

			else
				// no. roots != degree of elp => >t errors and cannot solve
				;

		}

		else         // elp has degree has degree >t hence cannot solve 
			;

	}

	else       // no non-zero syndromes => no errors: output received codeword 
		;

	for (i = 0; i < (parity_len+2) ; i++)
		free(elp[i]);
	free(elp);
	free(d);
	free(l);
	free(u_lu);
	free(s);
	free(forney);
	free(err);
	free(omega);
	free(phi);
	free(phiprime);
	free(tau);
	free(root);
	free(loc);
	free(z);
	free(reg);
	free(aux);
	/*
	   for (i=0;i<10536;i++)
	   {
		delete []elp[i];
	   }
	  delete []elp;
	*/
}


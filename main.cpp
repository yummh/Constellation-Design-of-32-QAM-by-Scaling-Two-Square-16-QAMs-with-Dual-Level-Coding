#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "include/mod2sparse.h"
#include "include/rcode.h"
#include "include/enc.h"
#include "include/dec.h"
#include "include/check.h"

#include "include/proto_say.h"
#include "include/encoding_say.h"
#include "include/DMOFDM_AWGN_1.h"
#include "include/random_say.h"

int main(int argc, char* argv[]) {
	int myid = 0;
	int Nfft = 1024;

	int i;
	int seed;

	double scaling = 0;
	double EbN0;
	double EsN0;
	int NCW = 10000;
	int MaxIteration;
	int channel;
	int decoder;
	int SEED;

	long err_cnt_1st = 0;
	long err_cnt_2nd = 0;

	long bit_err_1st = 0;
	long bit_err_2nd = 0;
	long bit_err_cnt_1st = 0;
	long bit_err_cnt_2nd = 0;
	long block_err_1st = 0;
	long block_err_2nd = 0;
	long check_1st = 0;
	long check_2nd = 0;
	long undetected_err_1st = 0;
	long undetected_err_2nd = 0;
	long cw_block_err = 0;

	double BER_1st, BER_2nd;
	double FER_1st, FER_2nd;
	double BER, FER;

	double R_1st;
	double R_2nd;
	double R;

	int B, E;
	int BG_1st, BG_2nd;
	int B_1st, B_2nd;
	int E_1st, E_2nd;
	int M_1st, M_2nd;
	int N_1st, N_2nd;
	int M_BASE_1st, M_BASE_2nd;         // BG row size
	int N_BASE_1st, N_BASE_2nd;         // BG column size
	int Z_factor_1st, Z_factor_2nd;
	int bps = 5;                        // bit per symbol: 32-QAM
	char* BG_type_1st = 0;              // 1: BG1, 2: BG2
	char* BG_type_2nd = 0;

	int Kb_1st, Kb_2nd, BofKb_1st, BofKb_2nd, iLS_1st = 0, iLS_2nd = 0;

	int cw;

	char FileName_pchk_1st[255];
	char FileName_proto_1st[255];
	char FileName_pchk_2nd[255];
	char FileName_proto_2nd[255];
	int** proto_matrix_1st;
	int** proto_matrix_2nd;
	char* codeword_inter_1st;
	char* codeword_inter_2nd;
	char* codeword_inter;

	int k1, k2, k3, k4;      // 1st layer
	int kk1, kk2, kk3, kk4;  // 2nd layer

	char* message_1st;
	char* codeword_1st;
	double* received_LR_1st;
	char* bit_stream_trans_1st;
	char* chks_1st;
	long* iter_1st;
	long* iter_bit_1st;
	long* iter_frame_1st;
	long* iter_check_1st;
	double* iter_ber_1st; //total
	double* iter_fer_1st; //total
	int* rand_arr_1st;
	long* arr_err_1st;

	char* message_2nd;
	char* codeword_2nd;
	double* received_LR_2nd;
	char* bit_stream_trans_2nd;
	char* chks_2nd;
	long* iter_2nd;
	long* iter_bit_2nd;
	long* iter_frame_2nd;
	long* iter_check_2nd;
	double* iter_ber_2nd; //total
	double* iter_fer_2nd; //total
	int* rand_arr_2nd;
	long* arr_err_2nd;


	// sets of LDPC lifting size Z
	int lifting_size[2][51] = {
		{2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,28,30,32,36,40,44,48,52,56,60,64,72,80,88,96,104,112,120,128,144,160,176,192,208,224,240,256,288,320,352,384},
		{0,1,0,2,1,3,0,4, 2, 5, 1, 6, 3, 7, 0, 4, 2, 5, 1, 6, 3, 7, 0, 4, 2, 5, 1, 6, 3, 7, 0, 4, 2, 5, 1,  6,  3,  7,  0,  4,  2,  5,  1,  6,  3,  7,  0,  4,  2,  5,  1} };

	if (argc == 7)
	{
		scaling = atof(argv[1]);
		EsN0 = atof(argv[2]);
		MaxIteration = atoi(argv[3]);

		//channel = 0;                      // 0: AWGN, 1: Rayleigh fading
		//decoder = 0;                      //

		B = atoi(argv[4]);
		B_1st = atoi(argv[5]);            // 1st layer B (message bit)
		E_1st = Nfft;                     // 1st layer E = 1024 bits
		B_2nd = B - B_1st;                // 2st layer B
		E_2nd = Nfft * (bps - 1);         // 2st layer E = 4096 bits

		SEED = atoi(argv[6]);

		// BG1 or BG2
		R_1st = (double)B_1st / (double)E_1st;
		R_2nd = (double)B_2nd / (double)E_2nd;
		R = R_1st + R_2nd * 4;

		if (B_1st <= 292) BG_1st = 2;
		else if (B_1st <= 3824 && R_1st <= 0.67) BG_1st = 2;
		else if (R_1st <= 0.25) BG_1st = 2;
		else BG_1st = 1;

		if (B_2nd <= 292) BG_2nd = 2;
		else if (B_2nd <= 3824 && R_2nd <= 0.67) BG_2nd = 2;
		else if (R_2nd <= 0.25) BG_2nd = 2;
		else BG_2nd = 1;

		// determine Kb
		if (BG_1st == 1) Kb_1st = 22;
		else {
			if (B_1st > 640) Kb_1st = 10;
			else if (B_1st > 560) Kb_1st = 9;
			else if (B_1st > 192) Kb_1st = 8;
			else Kb_1st = 6;
		}

		if (BG_2nd == 1) Kb_2nd = 22;
		else {
			if (B_2nd > 640) Kb_2nd = 10;
			else if (B_2nd > 560) Kb_2nd = 9;
			else if (B_2nd > 192) Kb_2nd = 8;
			else Kb_2nd = 6;
		}

		double BofKb_1st = (double)B_1st / (double)Kb_1st;
		double BofKb_2nd = (double)B_2nd / (double)Kb_2nd;

		for (int i = 0; i < 51; i++) {
			if (lifting_size[0][i] >= BofKb_1st) {
				Z_factor_1st = lifting_size[0][i];
				iLS_1st = lifting_size[1][i];
				break;
			}
		}
		for (int i = 0; i < 51; i++) {
			if (lifting_size[0][i] >= BofKb_2nd) {
				Z_factor_2nd = lifting_size[0][i];
				iLS_2nd = lifting_size[1][i];
				break;
			}
		}
		if (Z_factor_1st == -1 || Z_factor_2nd == -1) {
			printf("Z_factor = -1");
			return 0;
		}

		sprintf(FileName_pchk_1st, "proto/BG%d_original_iLs%d_%d.pchk", BG_1st, iLS_1st, Z_factor_1st);			// pchk file name
		sprintf(FileName_proto_1st, "proto/BG%d_original_iLs%d_%d.proto", BG_1st, iLS_1st, Z_factor_1st);	    // proto file name
		sprintf(FileName_pchk_2nd, "proto/BG%d_original_iLs%d_%d.pchk", BG_2nd, iLS_2nd, Z_factor_2nd);			// pchk file name
		sprintf(FileName_proto_2nd, "proto/BG%d_original_iLs%d_%d.proto", BG_2nd, iLS_2nd, Z_factor_2nd);		// proto file name

		if (BG_1st == 1)
		{
			BG_1st = 22;
			M_BASE_1st = 48;
			N_BASE_1st = 68;
			BG_type_1st = "BG1";
		}
		else
		{
			BG_1st = 10;
			M_BASE_1st = 42;
			N_BASE_1st = 52;
			BG_type_1st = "BG2";
		}

		if (BG_2nd == 1)
		{
			BG_2nd = 22;
			M_BASE_2nd = 48;
			N_BASE_2nd = 68;
			BG_type_2nd = "BG1";
		}
		else
		{
			BG_2nd = 10;
			M_BASE_2nd = 42;
			N_BASE_2nd = 52;
			BG_type_2nd = "BG2";
		}

		EbN0 = EsN0 - 10 * log10((double)((double)B / (double)Nfft));
	}
	else   //print error
	{
		if (myid == 0)
		{
			printf("Not enough input arguments. Check again.\n");
			printf("[scaling_value] [Es/No] [max_iter] [B] [B_1st] [SEED]\n");
		}
		exit(1);
	}


	// H matrix
	read_pchk(FileName_pchk_1st);          // read the H matrix (PCM)
	mod2sparse* H_1st = H;                 // H matrix ¡æ sparse matrix
	M_1st = mod2sparse_rows(H_1st);   
	N_1st = mod2sparse_cols(H_1st);     

	read_pchk(FileName_pchk_2nd);
	mod2sparse* H_2nd = H;
	M_2nd = mod2sparse_rows(H_2nd);
	N_2nd = mod2sparse_cols(H_2nd);


	// proto parameter load
	proto_matrix_1st = (int**)malloc(sizeof(int*) * M_BASE_1st);
	for (int i = 0; i < M_BASE_1st; i++)
		proto_matrix_1st[i] = (int*)malloc(sizeof(int) * N_BASE_1st);
	load_protograph(FileName_proto_1st, M_BASE_1st, N_BASE_1st, proto_matrix_1st);

	proto_matrix_2nd = (int**)malloc(sizeof(int*) * M_BASE_2nd);
	for (int i = 0; i < M_BASE_2nd; i++)
		proto_matrix_2nd[i] = (int*)malloc(sizeof(int) * N_BASE_2nd);
	load_protograph(FileName_proto_2nd, M_BASE_2nd, N_BASE_2nd, proto_matrix_2nd);


	// section setting
	k1 = 2 * Z_factor_1st;
	k2 = B_1st;
	k3 = BG_1st * Z_factor_1st;
	k4 = E_1st + k1 + (k3 - k2);

	kk1 = 2 * Z_factor_2nd;
	kk2 = B_2nd;
	kk3 = BG_2nd * Z_factor_2nd;
	kk4 = E_2nd + kk1 + (kk3 - kk2);


	// memory allocation
	message_1st = (char*)calloc(N_1st - M_1st, sizeof(char));
	codeword_1st = (char*)calloc(N_1st, sizeof(char));
	received_LR_1st = (double*)calloc(N_1st, sizeof(double));
	bit_stream_trans_1st = (char*)calloc(N_1st, sizeof(char));
	chks_1st = (char*)calloc(M_1st, sizeof(char));
	iter_1st = (long*)calloc(MaxIteration + 1, sizeof(long));
	iter_bit_1st = (long*)calloc(MaxIteration + 1, sizeof(long));
	iter_frame_1st = (long*)calloc(MaxIteration + 1, sizeof(long));
	iter_check_1st = (long*)calloc(MaxIteration + 1, sizeof(long));
	iter_ber_1st = (double*)calloc(MaxIteration + 1, sizeof(double));
	iter_fer_1st = (double*)calloc(MaxIteration + 1, sizeof(double));
	rand_arr_1st = (int*)calloc(E_1st, sizeof(int));
	arr_err_1st = (long*)calloc(N_1st, sizeof(long));

	message_2nd = (char*)calloc(N_2nd - M_2nd, sizeof(char));
	codeword_2nd = (char*)calloc(N_2nd, sizeof(char));
	received_LR_2nd = (double*)calloc(N_2nd, sizeof(double));
	bit_stream_trans_2nd = (char*)calloc(N_2nd, sizeof(char));
	chks_2nd = (char*)calloc(M_2nd, sizeof(char));
	iter_2nd = (long*)calloc(MaxIteration + 1, sizeof(long));
	iter_bit_2nd = (long*)calloc(MaxIteration + 1, sizeof(long));
	iter_frame_2nd = (long*)calloc(MaxIteration + 1, sizeof(long));
	iter_check_2nd = (long*)calloc(MaxIteration + 1, sizeof(long));
	iter_ber_2nd = (double*)calloc(MaxIteration + 1, sizeof(double));
	iter_fer_2nd = (double*)calloc(MaxIteration + 1, sizeof(double));
	rand_arr_2nd = (int*)calloc(E_2nd, sizeof(int));
	arr_err_2nd = (long*)calloc(N_2nd, sizeof(long));

	char* est_1st = (char*)calloc(N_1st, sizeof(char));


	// load random file
	seed = myid + SEED;     
	srand(seed);
	shuffle(rand_arr_1st, E_1st);
	shuffle(rand_arr_2nd, E_2nd);

	for (cw = 1; cw < NCW + 1; cw++) {
		//////////////////// Message Generator ////////////////////
		for (int i = 0; i < N_1st - M_1st; i++) {
			if (i < B_1st)
				message_1st[i] = rand() % 2;          // message index < B ¡æ random binary 
			else
				message_1st[i] = (char)0;             // message index > B ¡æ 0 (shortening)
		}
		for (int i = 0; i < N_2nd - M_2nd; i++) {
			if (i < B_2nd)
				message_2nd[i] = rand() % 2;
			else
				message_2nd[i] = (char)0;
		}


		//////////////////// Encoding ////////////////////
		encoding(message_1st, codeword_1st, proto_matrix_1st, M_1st, N_1st, Z_factor_1st, BG_type_1st);
		encoding(message_2nd, codeword_2nd, proto_matrix_2nd, M_2nd, N_2nd, Z_factor_2nd, BG_type_2nd);

		mod2sparse_mulvec(H_1st, codeword_1st, chks_1st);
		mod2sparse_mulvec(H_2nd, codeword_2nd, chks_2nd);

		for (i = 0; i < M_1st; i++)  
		{
			if (chks_1st[i] == 1) {  
				fprintf(stderr, "Output 1st block is not a codeword!  (Fails check 1)\n");
				exit(1);
			}
		}
		for (i = 0; i < M_2nd; i++)
		{
			if (chks_2nd[i] == 1) {
				fprintf(stderr, "Output 2nd block is not a codeword!  (Fails check 1)\n");
				exit(1);
			}
		}


		//////////////////// Interleaving ////////////////////
		codeword_inter_2nd = (char*)calloc(N_2nd, sizeof(char));

		for (int i = 0; i < N_2nd; i++) codeword_inter_2nd[i] = codeword_2nd[i];

		section_division_inter(codeword_2nd, codeword_inter_2nd, rand_arr_2nd, kk1, kk2, kk3, kk4);


		//////////////////// OFDM ////////////////////
		char* codeword_msg_1st = (char*)malloc(sizeof(char) * E_1st);
		char* codeword_msg_2nd = (char*)malloc(sizeof(char) * E_2nd);
		double* LR_msg_1st = (double*)malloc(sizeof(double) * E_1st);
		double* LR_msg_2nd = (double*)malloc(sizeof(double) * E_2nd);
		complex* received_symbol = (complex*)malloc(sizeof(complex) * Nfft);
		complex* channel = (complex*)malloc(sizeof(complex) * Nfft);


		// Save interleaving codeword copy of section E range
		int k, kk;
		for (int i = k1, k = 0; i < k2; i++, k++) codeword_msg_1st[k] = codeword_1st[i];
		for (int i = k3, k = k2 - k1; i < k4; i++, k++) codeword_msg_1st[k] = codeword_1st[i];
		for (int i = kk1, kk = 0; i < kk2; i++, kk++) codeword_msg_2nd[kk] = codeword_inter_2nd[i];
		for (int i = kk3, kk = kk2 - kk1; i < kk4; i++, kk++) codeword_msg_2nd[kk] = codeword_inter_2nd[i];


		// Calculate LLR1
		OFDM_1(LR_msg_1st, received_symbol, codeword_msg_1st, codeword_msg_2nd, EsN0, scaling, channel);

		for (int i = k1, k = 0; i < k2; i++, k++) received_LR_1st[i] = LR_msg_1st[k];        // message bit (E¡æN)
		for (int i = k3, k = k2 - k1; i < k4; i++, k++) received_LR_1st[i] = LR_msg_1st[k];  // message + parity bit 

		iter_1st[prprp_decode(H_1st, received_LR_1st, bit_stream_trans_1st, chks_1st, &err_cnt_1st, iter_bit_1st, iter_frame_1st, iter_check_1st, codeword_1st, MaxIteration, BG_1st, B_1st, E_1st, Z_factor_1st)]++;

		for (int i = 0; i < N_1st; i++)
			est_1st[i] = bit_stream_trans_1st[i];

		char* est_bit_1st = (char*)calloc(E_1st, sizeof(char));

		for (int i = k1, k = 0; i < k2; i++, k++) est_bit_1st[k] = est_1st[i];        // message bit
		for (int i = k3, k = k2 - k1; i < k4; i++, k++) est_bit_1st[k] = est_1st[i];  // message + parity bit


		// Calculate LLR2
		OFDM_2(LR_msg_2nd, received_symbol, est_bit_1st, EsN0, scaling, channel);

		for (int i = kk1, kk = 0; i < kk2; i++, kk++) received_LR_2nd[i] = LR_msg_2nd[kk];           // message bit
		for (int i = kk3, kk = kk2 - kk1; i < kk4; i++, kk++) received_LR_2nd[i] = LR_msg_2nd[kk];   // message + parity bit

		section_division_deinter(received_LR_2nd, rand_arr_2nd, kk1, kk2, kk3, kk4);

		iter_2nd[prprp_decode(H_2nd, received_LR_2nd, bit_stream_trans_2nd, chks_2nd, &err_cnt_2nd, iter_bit_2nd, iter_frame_2nd, iter_check_2nd, codeword_2nd, MaxIteration, BG_2nd, B_2nd, E_2nd, Z_factor_2nd)]++;


		// bit error count
		bit_err_1st += err_cnt_1st;
		bit_err_2nd += err_cnt_2nd;

		int check_1st = 0;
		int check_2nd = 0;

		if (err_cnt_1st)
		{
			block_err_1st++;																	// Accumulate the number of block errors in cw1
			check_1st++;																	    //			  the number of block errors in cw1
			if (check(H_1st, bit_stream_trans_1st, chks_1st) == 0)   undetected_err_1st++;      // count undetected error seperatly
		}

		if (err_cnt_2nd)
		{
			block_err_2nd++;
			check_2nd++;
			if (check(H_2nd, bit_stream_trans_2nd, chks_2nd) == 0)   undetected_err_2nd++; 
		}

		if (check_1st + check_2nd) {
			cw_block_err++;
		}

		if (cw_block_err >= 100) {
			cw += 1;
			break;
		}
	}

	cw -= 1;
	BER_1st = ((double)bit_err_1st) / ((double)((double)cw * (double)B_1st));
	BER_2nd = ((double)bit_err_2nd) / ((double)((double)cw * (double)B_2nd));
	BER = ((double)bit_err_1st + (double)bit_err_2nd) / ((double)((double)cw * ((double)(B_1st)+(double)(B_2nd))));   // calc BER & FER after decoding all frame
	FER = ((double)cw_block_err) / ((double)cw);
	FER_1st = ((double)block_err_1st) / ((double)cw);
	FER_2nd = ((double)block_err_2nd) / ((double)cw);


	printf("Es/N0 : %.2f \t R1 : %.3f \t B1 : %d \t BER1 : %.8e \t BLER1 : %.8e \t bit_error1 : %5d   /   block_error1 : %3d   /   undetected_error1 : %ld\n", EsN0, R_1st, B_1st, BER_1st, FER_1st, bit_err_1st, block_err_1st, undetected_err_1st);
	printf("Es/N0 : %.2f \t R2 : %.3f \t B2 : %d \t BER2 : %.8e \t BLER2 : %.8e \t bit_error2 : %5d   /   block_error2 : %3d   /   undetected_error2 : %ld\n", EsN0, R_2nd, B_2nd, BER_2nd, FER_2nd, bit_err_2nd, block_err_2nd, undetected_err_2nd);
	printf("Es/N0 : %.2f \t R  : %.3f \t B  : %d \t BER  : %.8e \t BLER  : %.8e \t codeword : %5d\n\n", EsN0, R, B, BER, FER, cw);
}
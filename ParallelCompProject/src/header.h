/*
 * header.h
 *
 *  Created on: Aug 17, 2020
 *      Author: linuxu
 */

#ifndef HEADER_H_
#define HEADER_H_

#define MAX_SEQ1 3000
#define MAX_SEQ2 2000

const char* conservativeGroup[] = {"NDEQ",	"NEQK",	"STA",
								   "MILV",	"QHRK",	"NHQK",
								   "FYW",	"HY",	"MILF"};

const int conservativeCount = 9;

const char* semiConservativeGroup[] = {"SAG",		"ATV",		"CSA",
									   "SGND",		"STPA",		"STNK",
									   "NEQHRK",	"NDEQHK",	"SNDEQK",
										 	 	 	"HFY",		"FVLIM"};

const int semiConservativeCount = 9;

enum conserve_type { SEMI_CONSERVATIVE = 0, CONSERVATIVE = 1 };

typedef struct {
	float weights[4];
	char* seq1;
	int numOfSeq2;
	char** seq2;
} DataSeq;

int readInputFile(DataSeq* data);
int getLine(char* line, FILE* file);
void saveResultsOutputFile(int* results_n, int* results_k, int numOfResults);

void broadcastGroundData(DataSeq* data, int my_rank);
int sendSeq2(DataSeq* data, int numOfProcs);
void recieveSeq2(DataSeq* data, MPI_Status* status);

void produceMaxAlignScore(char* seq1, char* seq2, float* weights, int* n, int* k);
float directScoreCalc(char* seq1, char* seq2, float* weights, int n, int k);
int checkConservativation(char c1, char c2, const char* group[], int groupSize);

void freeData(DataSeq* data);

#endif /* HEADER_H_ */

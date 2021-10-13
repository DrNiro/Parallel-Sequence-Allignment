/*
 ============================================================================
 Name        : ParallelCompProject.c
 Author      : Nir Druk
 Description : Final project for Parallel Computations course - String compare with score.
 ============================================================================
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include <omp.h>
#include "header.h"

#define ROOT 0
#define TAG 0
#define INFINITY 999999

int main(int argc, char* argv[]) {

	int  my_rank; 			/* rank of process */
	int  p;       			/* number of processes */
	MPI_Status status ;   	/* return status for receive */

	double startTime, finishTime;
	int n, k;
	int* result_n;			// array to hold offset results.
	int* result_k;			// array to hold mutate-index results.

	DataSeq data; 			// struct to hold relevant information for each process.
	int leftoverForRoot;	// number of jobs left for root process after distribution (will be calculated later on).

	/* start up MPI */
	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	startTime = MPI_Wtime();

//	master process. here its job is to read from file and distribute work to all processes (including himself?).
	if(my_rank == ROOT)
	{
		if(!readInputFile(&data))
		{
			printf("Couldn't read file.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if(data.numOfSeq2 < p)
		{
			printf("Invalid number of processes: %d. Must be lower or equal to the amount of 'Seq2's: %d.\n", p, data.numOfSeq2);
			MPI_Abort(MPI_COMM_WORLD, 2);
		}
	}

//	broadcast common data (weights and seq1)
	broadcastGroundData(&data, my_rank);

	if(my_rank == ROOT) // master send data and start process leftover data
	{
		result_n = (int*)malloc(sizeof(int)*data.numOfSeq2); // prepare result offset array.
		result_k = (int*)malloc(sizeof(int)*data.numOfSeq2); // prepare result mutate-index array.
		leftoverForRoot = sendSeq2(&data, p); // send seq2s to sub-processes.

		printf("total number of seq2s: %d\n", data.numOfSeq2);
		int i, seqIdx;
		for(i = 0; i < leftoverForRoot; i++) // calculate data in master and save in arrays.
		{
			seqIdx = data.numOfSeq2-leftoverForRoot+i;
			produceMaxAlignScore(data.seq1, data.seq2[seqIdx], data.weights, &n, &k);
			result_n[seqIdx] = n;
			result_k[seqIdx] = k;
		}
	}
	else // sub-process receive sent data and sending results to master.
	{
		recieveSeq2(&data, &status);
		int i;
		for(i = 0; i < data.numOfSeq2; i++)
		{
			produceMaxAlignScore(data.seq1, data.seq2[i], data.weights, &n, &k);
			MPI_Send(&n, 1, MPI_INT, ROOT, TAG, MPI_COMM_WORLD);
			MPI_Send(&k, 1, MPI_INT, ROOT, TAG, MPI_COMM_WORLD);
		}
	}

//	measure and show time for each process.
	finishTime = MPI_Wtime();
	printf("Process #%d finish time: %f\n", my_rank, finishTime - startTime);

//	receive results from sub-processes into master.
//	assemble results in order of the input file, and save in result array.
//	save result to output file.
	if(my_rank == ROOT) {
		int i, j;
		int idx = 0;
		for(i = 1; i < p; i++) {
			for(j = 0; j < data.numOfSeq2/p; j++) {
				MPI_Recv(&result_n[idx], 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&result_k[idx], 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
				idx++;
			}
		}

//		measure and show total time after receiving all results.
		finishTime = MPI_Wtime();
		printf("Finished all, Total time: %f\n", finishTime - startTime);

//		show all results.
		printf("\n\n");
		for(i = 0; i < data.numOfSeq2; i++) {
			printf("n: %d, k: %d\n", result_n[i], result_k[i]);
		}

		saveResultsOutputFile(result_n, result_k, data.numOfSeq2);
	}

//	free allocations.
	freeData(&data);
	if(my_rank == ROOT) {
		free(result_n);
		free(result_k);
	}

	/* shut down MPI */
	MPI_Finalize();

	return 0;
}

// process data from input file into data structure.
int readInputFile(DataSeq* data)
{
	int i;
	int lineLen;
//	char* file_name = "//home//linuxu//eclipse-workspace//ParallelCompProject//input.txt";
//	char* file_name = "input.txt";

	FILE* file;
	file = fopen("input.txt", "r");
	if (file == NULL)
	{
		printf("Could not open file\n");
		return 0;
	}

	//read weights and skip '\n'
	fscanf(file, "%f %f %f %f %*[\n]", &data->weights[0], &data->weights[1], &data->weights[2], &data->weights[3]);

	//read seq1
	char seq1_holder[MAX_SEQ1];
	lineLen = getLine(seq1_holder, file);
	data->seq1 = (char*)malloc(lineLen);
	strcpy(data->seq1, seq1_holder);

	//read numOfSeq2 and skip '\n'
	fscanf(file, "%d %*[\n]", &data->numOfSeq2);

	//read all seq2
	data->seq2 = (char**)malloc(sizeof(char*)*(data->numOfSeq2));
	char seq2_holder[data->numOfSeq2][MAX_SEQ2];
	for(i = 0; i < data->numOfSeq2; i++)
	{
		lineLen = getLine(seq2_holder[i], file);
		data->seq2[i] = (char*)malloc(lineLen + 1);
		strcpy(data->seq2[i], seq2_holder[i]);
		data->seq2[i][lineLen] = '\0';
	}

	fclose(file);
	printf("read file successfully\n");
	return 1;
}

void saveResultsOutputFile(int* results_n, int* results_k, int numOfResults) {
	FILE * outputFile;
	int i;
	/* open the file for writing*/
	outputFile = fopen ("output.txt","w");

	for(i = 0; i < numOfResults; i++)
	{
	   fprintf(outputFile, "n: %d, k: %d\n",results_n[i], results_k[i]);
	}

	fclose(outputFile);
}

// get a line from file and return its length.
int getLine(char* line, FILE* file)
{
	char c;
	int i;
	for(i = 0; i < MAX_SEQ1; i++)
	{
		c = fgetc(file);
		if(c == EOF || c == '\n')
		{
			line[i] = '\0';
			break;
		} else {
			line[i] = c;
		}
	}
	return i-1;
}

// broadcast data that all threads needs.
void broadcastGroundData(DataSeq* data, int my_rank)
{
	int s1_len;

	if(my_rank == ROOT)
	{
		s1_len = strlen(data->seq1);
	}

	MPI_Bcast(&s1_len, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	if(my_rank != ROOT)
	{
		data->seq1 = (char*)malloc(s1_len+1);
	}

	MPI_Bcast(data->weights, 4, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(data->seq1, s1_len+1, MPI_CHAR, ROOT, MPI_COMM_WORLD);
}

// call from root process only. distribute all seq2 between processes.
// return the amount of tasks for master process.
int sendSeq2(DataSeq* data, int numOfProcs)
{
	int i, j, sendCount, leftover;
	int seq_len;

	sendCount = data->numOfSeq2 / numOfProcs; // calculate evenly how many seq2s each process receives.
	leftover = data->numOfSeq2 % numOfProcs; // calculate how many additional seq2s master process need to receive.

//	send all tasks seq2s.
	for(i = 1; i < numOfProcs; i++)
	{
		MPI_Send(&sendCount, 1, MPI_INT, i, TAG, MPI_COMM_WORLD);
		for(j = i; j < sendCount + i; j++)
		{
			seq_len = strlen(data->seq2[j-1]);
			MPI_Send(&seq_len, 1, MPI_INT, i, TAG, MPI_COMM_WORLD);
			MPI_Send(data->seq2[(j - 1)], seq_len+1, MPI_CHAR, i, TAG, MPI_COMM_WORLD);
		}
	}
	return sendCount + leftover;
}

// receive seq2
void recieveSeq2(DataSeq* data, MPI_Status* status)
{
	int i, seq2_len, sentCount;

	MPI_Recv(&sentCount, 1, MPI_INT, ROOT, TAG, MPI_COMM_WORLD, status);
	data->seq2 = (char**)malloc(sizeof(char*)*sentCount);
	data->numOfSeq2 = sentCount;

	for(i = 0; i < sentCount; i++)
	{
		MPI_Recv(&seq2_len, 1, MPI_INT, ROOT, TAG, MPI_COMM_WORLD, status);
		data->seq2[i] = (char*)malloc(seq2_len + 1);
		MPI_Recv(data->seq2[i], seq2_len + 1, MPI_CHAR, ROOT, TAG, MPI_COMM_WORLD, status);
	}
}

// calculate the score between 2 sequences and saves the offset and mutation-index in n and k respectively.
void produceMaxAlignScore(char* seq1, char* seq2, float* weights, int* n, int* k) {
	int i, j, offsetAmount, mutesAmount;
	float maxScore = -INFINITY; // initial maxScore to a large negative value.
	float score;
	mutesAmount = strlen(seq2);
	offsetAmount = strlen(seq1) - mutesAmount;
// 	using openMP parallel for here to parallelize the most heavy task 'directScoreCalc'.
//	this task include in it many for loops, but only the 2 outer for loops are parallelized.
//	there is no other 'parallel for' parallelization because any use of it causes a negative impact on the finish time.
#pragma omp parallel for collapse(2) // parallelize the heavy task 'directScoreCalc' within 2 nested for-loops.
	for(i = 0; i < offsetAmount; i++) {
		for(j = 1; j < mutesAmount; j++) {
			score = directScoreCalc(seq1, seq2, weights, i, j);

			#pragma omp critical // This block need to be accessed by only 1 thread at a time to prevent race-condition.
			{
				if(score > maxScore) // if condition satisfied - swap n and k, set a new maxScore
				{
					maxScore = score;
					*n = i;
					*k = j;
				}
			}
		}
	}
}

// calculate the score by taking in account the directions from the assignment in the fastest method.
float directScoreCalc(char* seq1, char* seq2, float* weights, int n, int k) {
	int astriks = 0, colons = 0, dots = 0, spaces = 0;
	int i;
	int seenK = 0;
	char c1, c2;

	int seq2len = strlen(seq2);
	//	validate that seq2 doesn't pass seq1 after offset.
	if(seq2len + n > strlen(seq1))
	{
		printf("offset not legal (seq2 passed seq1)\n");
		return weights[0]*seq2len + 1; //return max score possible + 1 to indicate error.
	}

	for(i = 0; i < seq2len; i++)
	{
		if(i == k) {
			seenK = 1;
			spaces++;
		}
		c1 = seq1[i+n+seenK];
		c2 = seq2[i];

		if(c1 == c2)
			astriks++;
		else if(checkConservativation(c1, c2, conservativeGroup, conservativeCount))
			colons++;
		else if(checkConservativation(c1, c2, semiConservativeGroup, semiConservativeCount))
			dots++;
		else
			spaces++;
	}

	return weights[0]*astriks - weights[1]*colons - weights[2]*dots - weights[3]*spaces;
}

// function checking if 2 characters are in the same conservative group.
int checkConservativation(char c1, char c2, const char* group[], int groupSize) {
	int i;
	for(i = 0; i < groupSize; i++) {
		if(strchr(group[i], c1) && strchr(group[i], c2)) {
			return 1;
		}
	}
	return 0;
}

void freeData(DataSeq* data) {
	int i;
	free(data->seq1);
	for(i = 0; i < data->numOfSeq2; i++) {
		free(data->seq2[i]);
	}
	free(data->seq2);
}

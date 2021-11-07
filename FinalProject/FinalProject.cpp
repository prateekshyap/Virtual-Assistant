// FinalProject.cpp : Defines the entry point for the console application.
//
using namespace std;

#include "stdafx.h"
#include <iostream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "commonvar.h"
#include "lbg.h"
#include "forwardbackward.h"
#include "viterbi.h"
#include "baumwelch.h"
#include "functions.h"

extern long double *** xi;

extern long double ** alpha, ** beta, ** delta, ** gamma, ** codebook;

extern long double ** A, ** AComplement, ** B, ** BComplement;

extern long double * Pi, * PiComplement;

extern int ** psi;

extern int * O, * qStar, * qStarComplement;

extern FILE * AComplementFile, * BComplementFile, * PiComplementFile;

extern char * resultWord;

FILE * AFile = NULL, * BFile = NULL, * PiFile = NULL;

extern int T, N, M, R, duration, p;

int D = 0, F = 0;
int ** fullO = NULL;

void trainBeginningModel(int);
void recognize(char *, char **, char **, char **, int, int, int);

int _tmain(int argc, _TCHAR* argv[])
{
	int d = 0, r = 0, range = 0, i = 0;
	char * temp = NULL;
	int trainingStatus = -1; //set a variable to know if the model has to be trained or not
	int buildingStatus = -1; //set a variable to know if codebook needs to be built or not
	FILE * file = fopen("data/info.txt","r"); //open info.txt
	fscanf(file,"%d",&trainingStatus); //read from the file
	fclose(file);
	
	//read the duration
	file = fopen("data/duration.txt","r");
	fscanf(file,"%d",&duration);
	fclose(file);

	//read N, M, range of amplitudes
	file = fopen("data/N.txt","r");
	fscanf(file,"%d",&N);
	fclose(file);
	file = fopen("data/M.txt","r");
	fscanf(file,"%d",&M);
	fclose(file);
	file = fopen("data/range.txt","r");
	fscanf(file,"%d",&range);
	fclose(file);
	
	//initialize
	A = new long double *[N]; //A Matrix
	B = new long double *[N]; //B Matrix
	Pi = new long double[N]; //Pi Array

	for (i = 0; i < N; ++i)
		A[i] = new long double[N];
	for (i = 0; i < N; ++i)
		B[i] = new long double[M];
	
	alpha = new long double *[N];
	beta = new long double *[N];

	for (i = 0; i < N; ++i)
		alpha[i] = new long double[T];
	for (i = 0; i < N; ++i)
		beta[i] = new long double[T];


	if (trainingStatus == 0) //if file says 0
	{
		printf("Training the initial Model\n\n");
		file = fopen("data/build.txt","r");
		fscanf(file,"%d",&buildingStatus);
		fclose(file);
		trainBeginningModel(buildingStatus); //call the model training function
		file = fopen("data/info.txt","w");
		fprintf(file,"1"); //replace 0 with 1 to avoid training in further executions
		fclose(file);
		printf("\nTraining completed\n\n");
	}

	//fill the file details
	file = fopen("data/D.txt","r");
	fscanf(file,"%d",&D);
	char * folder = "HMM/";
	char ** digits = new char * [D];
	for (d = 0; d < D; ++d)
		digits[d] = (char *)malloc(sizeof(char *));
	for (d = 0; d < D; ++d)
		fscanf(file,"%s",digits[d]);
	char ** files = NULL, ** dataFiles = NULL;
	fclose(file);
	dataFiles = new char * [R];
	for (r = 0; r < R; ++r)
		dataFiles[r] = (char *)malloc(sizeof(char *));
	for (r = 0; r < R; ++r)
		sprintf(dataFiles[r],"%d",r);
	files = new char * [4];
	for (r = 0; r < 3; ++r)
		files[r] = (char *)malloc(sizeof(char *));
	files[0] = "A";
	files[1] = "B";
	files[2] = "Pi";

	//switch case will be added here for all kind of operations (training, testing etc)
	
	//live training

	
	
	//add image display feature or youtube link feature here, don't modify the existing code

	
	
	//testing
	
	char command[200];
	sprintf(command,"Recording_Module.exe %d data/o.wav data/o.txt",duration);
	std :: system(command);
	recognize(folder, digits, files, dataFiles, D, R, range);
	printf("You spoke %s",resultWord);
	

	printf("Breakpoint\n");
	return 0;
}

void trainBeginningModel(int buildingStatus)
{
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	bool isUpdated = true;
	int i = 0, j = 0, k = 0, r = 0, c = 0, d = 0, f = 0, t = 0, m = 0, modelCount = 0;
	long double ** refA = NULL, ** refB = NULL; 
	long double * refPi = NULL;
	FILE * file = NULL;
	char temp[100];
	int range = 0, progressCount = 0, percentageSize = 25, modelLoop = 4, updationLimit = 0;
	long double value = 0;

	


	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Reading file information
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//Read the updated information from info.txt
	file = fopen("data/F.txt","r");
	fscanf(file,"%d",&R); F = R;
	fclose(file);

	file = fopen("data/D.txt","r");
	fscanf(file,"%d",&D);
	
	char * folder = "HMM/";
	char ** digits = new char * [D];
	for (d = 0; d < D; ++d)
		digits[d] = (char *)malloc(sizeof(char *));
	for (d = 0; d < D; ++d)
		fscanf(file,"%s",digits[d]);
	char ** files = NULL, ** dataFiles = NULL;
	fclose(file);

	file = fopen("data/duration.txt","r");
	fscanf(file,"%d",&duration);
	fclose(file);

	file = fopen("data/range.txt","r");
	fscanf(file,"%d",&range);
	fclose(file);

	file = fopen("data/N.txt","r");
	fscanf(file,"%d",&N);
	fclose(file);

	file = fopen("data/M.txt","r");
	fscanf(file,"%d",&M);
	fclose(file);
	

	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Building the universe and the codebook
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	dataFiles = new char * [R];
	for (r = 0; r < R; ++r)
		dataFiles[r] = (char *)malloc(sizeof(char *));
	for (r = 0; r < R; ++r)
		sprintf(dataFiles[r],"%d",r);
	files = new char * [4];
	for (r = 0; r < 4; ++r)
		files[r] = (char *)malloc(sizeof(char *));
	files[0] = "A";
	files[1] = "B";
	files[2] = "O";
	files[3] = "Pi";
	resetModel(folder, digits, dataFiles, D, R);
	if (buildingStatus == 0)
	{
		printf("Building the Universe\n");
		buildUniverse(folder, digits, dataFiles, D, R, range);
		printf("\n\nUniverse building completed\n\nBuilding the Codebook\n");
		buildCodebook();
		printf("\n\nCodebook building completed\n");
	}
	
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Define variables for HMM
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

	//define everything

	delta = new long double *[N];
	psi = new int *[N];
	qStar = new int[T];
	qStarComplement = new int[T];

	for (i = 0; i < N; ++i)
		delta[i] = new long double[T];
	for (i = 0; i < N; ++i)
		psi[i] = new int[T];

	xi = new long double ** [N];
	for (i = 0; i < N; ++i)
		xi[i] = new long double * [N];
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			xi[i][j] = new long double[T-1];

	gamma = new long double * [N];
	for (i = 0; i < N; ++i)
		gamma[i] = new long double[T];

	PiComplement = new long double [N];

	AComplement = new long double * [N];
	for (i = 0; i < N; ++i)
		AComplement[i] = new long double[N];

	BComplement = new long double * [N];
	for (i = 0; i < N; ++i)
		BComplement[i] = new long double[M];

	O = new int [N]; //Observation Sequence
	fullO = new int *[R];

	for (i = 0; i < R; ++i)
		fullO[i] = new int[T];

	refA = new long double *[N]; //A Matrix
	refB = new long double *[N]; //B Matrix
	refPi = new long double[N]; //Pi Array

	for (i = 0; i < N; ++i)
		refA[i] = new long double[N];
	for (i = 0; i < N; ++i)
		refB[i] = new long double[M];
	
	char * buffer = new char[1024];
	char skipLine[1024];
	
	char fileName[100], AComplementFileName[100], BComplementFileName[100], PiComplementFileName[100];




	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Build the initial HMM
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//run the hmm building process
	printf("\nBuilding the model\n");
	while (modelCount != modelLoop)
	{
		++modelCount;
		//printf("%d\n",modelCount);
		//converge the model for each observatio sequence
		for (d = 0; d < D; ++d)
		{
			generateObservationSequences(folder, digits, dataFiles, d, R, range); //generate the full observation sequence
			for (f = 0; f < 4; ++f)
			{
				strcpy(fileName,folder);
				strcat(fileName,digits[d]);
				strcat(fileName,"/");
				strcat(fileName,files[f]);
				strcat(fileName,".txt");
			
				file = fopen(fileName,"r");

				if (f == 0)
				{
					for (i = 0; i < N; ++i)
					{
						for (j = 0; j < N; ++j)
						{
							fscanf(file,"%lf",&A[i][j]);
							refA[i][j] = A[i][j];
						}
					}
				}
				else if (f == 1)
				{
					for (i = 0; i < N; ++i)
					{
						for (j = 0; j < M; ++j)
						{
							fscanf(file,"%lf",&B[i][j]);
							refB[i][j] = B[i][j];
						}
					}
				}
				else if (f == 2)
				{
					for (i = 0; i < R; ++i)
					{
						//fscanf(file,"%s %d %s",&skipLine,&k,&skipLine);
						for (j = 0; j < T; ++j)
							fscanf(file,"%d",&fullO[i][j]);
					}				
				}
				else if (f == 3)
				{
					for (i = 0; i < N; ++i)
					{
						fscanf(file,"%lf",&Pi[i]);
						refPi[i] = Pi[i];				
					}
				}
			
				fclose(file);
			}
			for (k = 0; k < R; ++k)
			{
				sprintf(temp,"%d",k);
			
				strcpy(AComplementFileName,folder);
				strcat(AComplementFileName,digits[d]);
				strcat(AComplementFileName,"/");
				strcat(AComplementFileName,"AComplement");
				strcat(AComplementFileName,temp);
				strcat(AComplementFileName,".txt");

				strcpy(BComplementFileName,folder);
				strcat(BComplementFileName,digits[d]);
				strcat(BComplementFileName,"/");
				strcat(BComplementFileName,"BComplement");
				strcat(BComplementFileName,temp);
				strcat(BComplementFileName,".txt");

				strcpy(PiComplementFileName,folder);
				strcat(PiComplementFileName,digits[d]);
				strcat(PiComplementFileName,"/");
				strcat(PiComplementFileName,"PiComplement");
				strcat(PiComplementFileName,temp);
				strcat(PiComplementFileName,".txt");

				AComplementFile = fopen(AComplementFileName,"w+");
				BComplementFile = fopen(BComplementFileName,"w+");
				PiComplementFile = fopen(PiComplementFileName,"w+");

				//reset A
				for (i = 0; i < N; ++i)
					for (j = 0; j < N; ++j)
						A[i][j] = refA[i][j];
			
				//reset B
				for (i = 0; i < N; ++i)
					for (m = 0; m < M; ++m)
						B[i][m] = refB[i][m];
			
				//reset Pi
				for (j = 0; j < N; ++j)
					Pi[j] = refPi[j];

				//set new O
				O = fullO[k];
				
				/*
				//reset A
				for (i = 0; i < N; ++i)
				{
					for (j = 0; j < N; ++j)
						printf("%g ",A[i][j]);
					printf("\n");
				}
			
				//reset B
				for (i = 0; i < N; ++i)
				{
					for (m = 0; m < M; ++m)
						printf("%g ",B[i][j]);
					printf("\n");
				}
			
				//reset Pi
				for (j = 0; j < N; ++j)
					printf("%g ",Pi[j]);
				printf("\n");

				for (j = 0; j < T; ++j)
					printf("%d ",O[j]);
				printf("\n");
				*/


				isUpdated = true;

				//update new pStar
				runViterbi(0);
				updationLimit = 0;

				//run till the model is being updated
				while(isUpdated && updationLimit++ <= 200)
				{
					runForwardBackward();
					runBaumWelch();
					runViterbi(1);
					isUpdated = compareAndUpdateModel();
				}
			
				//print the new models to respective files
				for (i = 0; i < N; ++i)
				{
					for (j = 0; j < N; ++j)
						fprintf(AComplementFile,"%g ",A[i][j]);
					fprintf(AComplementFile,"\n");
				}
			
				for (i = 0; i < N; ++i)
				{
					for (m = 0; m < M; ++m)
						fprintf(BComplementFile,"%g ",B[i][m]);
					fprintf(BComplementFile,"\n");
				}
			
				for (j = 0; j < N; ++j)
					fprintf(PiComplementFile,"%g ",Pi[j]);
				fprintf(PiComplementFile,"\n");

				fclose(AComplementFile);
				fclose(BComplementFile);
				fclose(PiComplementFile);
			}
		}
		//printf("\nModel converged\n");
		//replace the initial model with the new average model
		for (d = 0; d < D; ++d)
		{
			strcpy(fileName,folder);
			strcat(fileName,digits[d]);
			strcat(fileName,"/");
			strcat(fileName,"A");
			strcat(fileName,".txt");
			AFile = fopen(fileName,"w");

			strcpy(fileName,folder);
			strcat(fileName,digits[d]);
			strcat(fileName,"/");
			strcat(fileName,"B");
			strcat(fileName,".txt");
			BFile = fopen(fileName,"w");

			strcpy(fileName,folder);
			strcat(fileName,digits[d]);
			strcat(fileName,"/");
			strcat(fileName,"Pi");
			strcat(fileName,".txt");
			PiFile = fopen(fileName,"w");

			for (i = 0; i < N; ++i)
				for (j = 0; j < N; ++j)
					AComplement[i][j] = 0;

			for (i = 0; i < N; ++i)
				for (j = 0; j < M; ++j)
					BComplement[i][j] = 0;

			for (i = 0; i < N; ++i)
				PiComplement[i] = 0;

			for (k = 0; k < R; ++k)
			{
				sprintf(temp,"%d",k);
			
				strcpy(AComplementFileName,folder);
				strcat(AComplementFileName,digits[d]);
				strcat(AComplementFileName,"/");
				strcat(AComplementFileName,"AComplement");
				strcat(AComplementFileName,temp);
				strcat(AComplementFileName,".txt");

				strcpy(BComplementFileName,folder);
				strcat(BComplementFileName,digits[d]);
				strcat(BComplementFileName,"/");
				strcat(BComplementFileName,"BComplement");
				strcat(BComplementFileName,temp);
				strcat(BComplementFileName,".txt");

				strcpy(PiComplementFileName,folder);
				strcat(PiComplementFileName,digits[d]);
				strcat(PiComplementFileName,"/");
				strcat(PiComplementFileName,"PiComplement");
				strcat(PiComplementFileName,temp);
				strcat(PiComplementFileName,".txt");

				AComplementFile = fopen(AComplementFileName,"r");
				BComplementFile = fopen(BComplementFileName,"r");
				PiComplementFile = fopen(PiComplementFileName,"r");

				for (i = 0; i < N; ++i)
				{
					for (j = 0; j < N; ++j)
					{
						fscanf(AComplementFile,"%lf",&value);
						AComplement[i][j] += value;
					}
				}

				for (i = 0; i < N; ++i)
				{
					for (j = 0; j < M; ++j)
					{
						fscanf(BComplementFile,"%lf",&value);
						BComplement[i][j] += value;
					}
				}

				for (i = 0; i < N; ++i)
				{
					fscanf(PiComplementFile,"%lf",&value);
					PiComplement[i] += value;
				}

				fclose(AComplementFile);
				fclose(BComplementFile);
				fclose(PiComplementFile);
			}

			for (i = 0; i < N; ++i)
				for (j = 0; j < N; ++j)
					AComplement[i][j] /= R;

			for (i = 0; i < N; ++i)
				for (j = 0; j < M; ++j)
					BComplement[i][j] /= R;

			for (i = 0; i < N; ++i)
				PiComplement[i] /= R;

			ensureStochastic();

			for (i = 0; i < N; ++i)
			{
				for (j = 0; j < N; ++j)
					fprintf(AFile,"%g ",AComplement[i][j]);
				fprintf(AFile,"\n");
			}

			for (i = 0; i < N; ++i)
			{
				for (j = 0; j < M; ++j)
					fprintf(BFile,"%g ",BComplement[i][j]);
				fprintf(BFile,"\n");
			}

			for (i = 0; i < N; ++i)
				fprintf(PiFile,"%g ",PiComplement[i]);
			fprintf(PiFile,"\n");

			fclose(AFile);
			fclose(BFile);
			fclose(PiFile);
		}

		//printing the progress bar
		progressCount += percentageSize;
		printf("\r"); //move to the beginning of the line
		/*print the progress*/
		i = 0;
		printf("[");
		for (i = 1; i <= progressCount; ++i)
			printf("|");
		for (; i <= 100; ++i)
			printf(" ");
		printf("]%d",(progressCount)); //print the percentage
	}
	printf("\n\nModel building completed\n\n");
}



/*

////////////////////////////////////////////////////
Working main code for viterbi
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

char * folders[] = {"Samples/Digit ","Samples/Digit ","Samples/Digit "};
	char * digits[] = {"1","7","8"};
	char * files[] = {"A_","B_","obs_seq_","Pi_"};

	char * buffer = new char[1024];
	char skipLine[1024];
	A = new long double *[N]; //A Matrix
	B = new long double *[N]; //B Matrix
	O = new int *[N]; //Observation Sequence
	Pi = new long double[N]; //Pi Array
	int i = 0, j = 0, k = 0, r = 0, c = 0, d = 0, f = 0, t = 0;
	char fileName[100], outputFileName[100];
	FILE * file = NULL;
	
	for (i = 0; i < N; ++i)
		A[i] = new long double[N];
	for (i = 0; i < N; ++i)
		B[i] = new long double[M];
	for (i = 0; i < N; ++i)
		O[i] = new int[T];

	for (d = 0; d < D; ++d)
	{
		strcpy(outputFileName,folders[d]);
		strcat(outputFileName,digits[d]);
		strcat(outputFileName,"/");
		strcat(outputFileName,"Output.txt");

		dataOutputFile = fopen(outputFileName,"w+");
		printf("-----------------------------------------------------------\n");
		printf("                         Digit %s\n",digits[d]);
		printf("-----------------------------------------------------------\n");
		fprintf(dataOutputFile,"-----------------------------------------------------------\n");
		fprintf(dataOutputFile,"                         Digit %s\n",digits[d]);
		fprintf(dataOutputFile,"-----------------------------------------------------------\n");

		for (f = 0; f < F; ++f)
		{
			strcpy(fileName,folders[d]);
			strcat(fileName,digits[d]);
			strcat(fileName,"/");
			strcat(fileName,files[f]);
			strcat(fileName,digits[d]);
			strcat(fileName,".txt");
			
			file = fopen(fileName,"r");

			if (f == 0)
			{
				printf("Reading %s for A matrix\n",fileName);
				fprintf(dataOutputFile,"Reading %s for A matrix\n",fileName);
				for (i = 0; i < N; ++i)
					for (j = 0; j < N; ++j)
						fscanf(file,"%lf",&A[i][j]);

				/*
				//remove this comment to print on console
				for (i = 0; i < N; ++i)
				{
					for (j = 0; j < N; ++j)
						printf("%.20lf\t", A[i][j]);
					printf("\n");
				}
				
			}
			else if (f == 1)
			{
				printf("Reading %s for B matrix\n",fileName);
				fprintf(dataOutputFile,"Reading %s for B matrix\n",fileName);
				for (i = 0; i < N; ++i)
					for (j = 0; j < M; ++j)
						fscanf(file,"%lf",&B[i][j]);
				/*
				//remove this comment to print on console
				for (i = 0; i < N; ++i)
				{
					for (j = 0; j < M; ++j)
						printf("%.20lf\t", B[i][j]);
					printf("\n");
				}
				
			}
			else if (f == 2)
			{
				printf("Reading %s for observation sequences\n",fileName);
				fprintf(dataOutputFile,"Reading %s for observation sequences\n",fileName);
				for (i = 0; i < 5; ++i)
				{
					fscanf(file,"%s %d %s",&skipLine,&k,&skipLine);
					for (j = 0; j < T; ++j)
						fscanf(file,"%d",&O[i][j]);
				}

				/*
				//remove this comment to print on console
				for (i = 0; i < N; ++i)
				{
					for (j = 0; j < T; ++j)
						printf("%d ",O[i][j]);
					printf("\n");
				}
				
			}
			else if (f == 3)
			{
				printf("Reading %s for Pi\n",fileName);
				fprintf(dataOutputFile,"Reading %s for Pi\n",fileName);
				for (i = 0; i < N; ++i)
					fscanf(file,"%lf",&Pi[i]);

				/*
				//remove this comment to print on console
				for (i = 0; i < N; ++i)
					printf("%.20lf\t",Pi[i]);
				printf("\n");
				
			}
			
			fclose(file);
		}
		for (i = 0; i < N; ++i)
			runViterbi(i);
		printf("Required outputs are printed in %s\n",outputFileName);
		fclose(dataOutputFile);
	}

	printf("Breakpoint\n");
	return 0;


////////////////////////////////////////////////////
Working main code for three problems to verify data
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


char * folders[] = {"Samples/Digit ","Samples/Digit ","Samples/Digit "};
	char * digits[] = {"1","7","8"};
	char * files[] = {"A_","B_","obs_seq_","Pi_"};

	char * buffer = new char[1024];
	char skipLine[1024];
	A = new long double *[N]; //A Matrix
	B = new long double *[N]; //B Matrix
	O = new int [N]; //Observation Sequence
	fullO = new int *[N];
	Pi = new long double[N]; //Pi Array
	int i = 0, j = 0, k = 0, r = 0, c = 0, d = 0, f = 0, t = 0;
	char fileName[100], dataOutputFileName[100], modelOutputFileName[100];
	FILE * file = NULL;
	
	for (i = 0; i < N; ++i)
		A[i] = new long double[N];
	for (i = 0; i < N; ++i)
		B[i] = new long double[M];
	for (i = 0; i < N; ++i)
		fullO[i] = new int[T];

	for (d = 0; d < D; ++d)
	{
		strcpy(dataOutputFileName,folders[d]);
		strcat(dataOutputFileName,digits[d]);
		strcat(dataOutputFileName,"/");
		strcat(dataOutputFileName,"DataOutput.txt");

		strcpy(modelOutputFileName,folders[d]);
		strcat(modelOutputFileName,digits[d]);
		strcat(modelOutputFileName,"/");
		strcat(modelOutputFileName,"ModelOutput.txt");

		dataOutputFile = fopen(dataOutputFileName,"w+");
		modelOutputFile = fopen(modelOutputFileName,"w+");
		printf("-----------------------------------------------------------\n");
		printf("                         Digit %s\n",digits[d]);
		printf("-----------------------------------------------------------\n");
		fprintf(dataOutputFile,"-----------------------------------------------------------\n");
		fprintf(dataOutputFile,"                         Digit %s\n",digits[d]);
		fprintf(dataOutputFile,"-----------------------------------------------------------\n");
		fprintf(modelOutputFile,"-----------------------------------------------------------\n");
		fprintf(modelOutputFile,"                         Digit %s\n",digits[d]);
		fprintf(modelOutputFile,"-----------------------------------------------------------\n");

		for (f = 0; f < F; ++f)
		{
			strcpy(fileName,folders[d]);
			strcat(fileName,digits[d]);
			strcat(fileName,"/");
			strcat(fileName,files[f]);
			strcat(fileName,digits[d]);
			strcat(fileName,".txt");
			
			file = fopen(fileName,"r");

			if (f == 0)
			{
				printf("Reading %s for A matrix\n",fileName);
				fprintf(dataOutputFile,"Reading %s for A matrix\n",fileName);
				for (i = 0; i < N; ++i)
					for (j = 0; j < N; ++j)
						fscanf(file,"%lf",&A[i][j]);	
			}
			else if (f == 1)
			{
				printf("Reading %s for B matrix\n",fileName);
				fprintf(dataOutputFile,"Reading %s for B matrix\n",fileName);
				for (i = 0; i < N; ++i)
					for (j = 0; j < M; ++j)
						fscanf(file,"%lf",&B[i][j]);
			}
			else if (f == 2)
			{
				printf("Reading %s for observation sequences\n",fileName);
				fprintf(dataOutputFile,"Reading %s for observation sequences\n",fileName);
				for (i = 0; i < 5; ++i)
				{
					fscanf(file,"%s %d %s",&skipLine,&k,&skipLine);
					for (j = 0; j < T; ++j)
						fscanf(file,"%d",&fullO[i][j]);
				}				
			}
			else if (f == 3)
			{
				printf("Reading %s for Pi\n",fileName);
				fprintf(dataOutputFile,"Reading %s for Pi\n",fileName);
				for (i = 0; i < N; ++i)
					fscanf(file,"%lf",&Pi[i]);				
			}
			
			fclose(file);
		}
		for (i = 0; i < N; ++i)
		{
			O = fullO[i];
			runForwardBackward();
			runViterbi();
			runBaumWelch();
		}
		printf("Required outputs are printed in %s\n",dataOutputFileName);
		printf("Other information is printed in %s\n",modelOutputFileName);
		fclose(dataOutputFile); fclose(modelOutputFile);
	}

	printf("Breakpoint\n");
	return 0;


////////////////////////////////////////////////////////////////////////
Working main code for single convergence for each observation sequence
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

bool isUpdated = true;
	int i = 0, j = 0, k = 0, r = 0, c = 0, d = 0, f = 0, t = 0, m = 0;
	long double ** refA = NULL, ** refB = NULL; 
	long double * refPi = NULL;
	FILE * file = NULL;

	/*Read the updated information from info.txt
	file = fopen("Info.txt","r");
	fscanf(file,"%d",&R);

	/*define everything
	alpha = new long double *[N];
	beta = new long double *[N];

	for (i = 0; i < N; ++i)
		alpha[i] = new long double[T];
	for (i = 0; i < N; ++i)
		beta[i] = new long double[T];

	delta = new long double *[N];
	psi = new int *[N];
	qStar = new int[T];
	qStarComplement = new int[T];

	for (i = 0; i < N; ++i)
		delta[i] = new long double[T];
	for (i = 0; i < N; ++i)
		psi[i] = new int[T];

	xi = new long double ** [N];
	for (i = 0; i < N; ++i)
		xi[i] = new long double * [N];
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			xi[i][j] = new long double[T-1];

	gamma = new long double * [N];
	for (i = 0; i < N; ++i)
		gamma[i] = new long double[T];

	PiComplement = new long double [N];

	AComplement = new long double * [N];
	for (i = 0; i < N; ++i)
		AComplement[i] = new long double[N];

	BComplement = new long double * [N];
	for (i = 0; i < N; ++i)
		BComplement[i] = new long double[M];

	A = new long double *[N]; //A Matrix
	B = new long double *[N]; //B Matrix
	O = new int [N]; //Observation Sequence
	fullO = new int *[R];
	Pi = new long double[N]; //Pi Array

	for (i = 0; i < N; ++i)
		A[i] = new long double[N];
	for (i = 0; i < N; ++i)
		B[i] = new long double[M];
	for (i = 0; i < R; ++i)
		fullO[i] = new int[T];

	refA = new long double *[N]; //A Matrix
	refB = new long double *[N]; //B Matrix
	refPi = new long double[N]; //Pi Array

	for (i = 0; i < N; ++i)
		refA[i] = new long double[N];
	for (i = 0; i < N; ++i)
		refB[i] = new long double[M];
	
	char * folders[] = {"ConvergenceCheck/Digit ","ConvergenceCheck/Digit ","ConvergenceCheck/Digit "};
	char * digits[] = {"1","7","8"};
	char * files[] = {"A_","B_","obs_seq_","Pi_"};

	char * buffer = new char[1024];
	char skipLine[1024];
	
	char fileName[100], AComplementFileName[100], BComplementFileName[100], PiComplementFileName[100], temp[100];

	for (d = 0; d < D; ++d)
	{
		printf("-----------------------------------------------------------\n");
		printf("                         Digit %s\n",digits[d]);
		printf("-----------------------------------------------------------\n");
		
		for (f = 0; f < F; ++f)
		{
			strcpy(fileName,folders[d]);
			strcat(fileName,digits[d]);
			strcat(fileName,"/");
			strcat(fileName,files[f]);
			strcat(fileName,digits[d]);
			strcat(fileName,".txt");
			
			file = fopen(fileName,"r");

			if (f == 0)
			{
				printf("Reading %s for A matrix\n",fileName);
				for (i = 0; i < N; ++i)
				{
					for (j = 0; j < N; ++j)
					{
						fscanf(file,"%lf",&A[i][j]);
						refA[i][j] = A[i][j];
					}
				}
			}
			else if (f == 1)
			{
				printf("Reading %s for B matrix\n",fileName);
				for (i = 0; i < N; ++i)
				{
					for (j = 0; j < M; ++j)
					{
						fscanf(file,"%lf",&B[i][j]);
						refB[i][j] = B[i][j];
					}
				}
			}
			else if (f == 2)
			{
				printf("Reading %s for observation sequences\n",fileName);
				for (i = 0; i < 5; ++i)
				{
					fscanf(file,"%s %d %s",&skipLine,&k,&skipLine);
					for (j = 0; j < T; ++j)
						fscanf(file,"%d",&fullO[i][j]);
				}				
			}
			else if (f == 3)
			{
				printf("Reading %s for Pi\n",fileName);
				for (i = 0; i < N; ++i)
				{
					fscanf(file,"%lf",&Pi[i]);
					refPi[i] = Pi[i];				
				}
			}
			
			fclose(file);
		}
		for (k = 0; k < N; ++k)
		{
			sprintf(temp,"%d",k);
			
			strcpy(AComplementFileName,folders[d]);
			strcat(AComplementFileName,digits[d]);
			strcat(AComplementFileName,"/");
			strcat(AComplementFileName,digits[d]);
			strcat(AComplementFileName,"AComplement");
			strcat(AComplementFileName,temp);
			strcat(AComplementFileName,".txt");

			strcpy(BComplementFileName,folders[d]);
			strcat(BComplementFileName,digits[d]);
			strcat(BComplementFileName,"/");
			strcat(BComplementFileName,digits[d]);
			strcat(BComplementFileName,"BComplement");
			strcat(BComplementFileName,temp);
			strcat(BComplementFileName,".txt");

			strcpy(PiComplementFileName,folders[d]);
			strcat(PiComplementFileName,digits[d]);
			strcat(PiComplementFileName,"/");
			strcat(PiComplementFileName,digits[d]);
			strcat(PiComplementFileName,"PiComplement");
			strcat(PiComplementFileName,temp);
			strcat(PiComplementFileName,".txt");

			AComplementFile = fopen(AComplementFileName,"w+");
			BComplementFile = fopen(BComplementFileName,"w+");
			PiComplementFile = fopen(PiComplementFileName,"w+");

			//reset A
			for (i = 0; i < N; ++i)
				for (j = 0; j < N; ++j)
					A[i][j] = refA[i][j];
			
			//reset B
			for (i = 0; i < N; ++i)
				for (m = 0; m < M; ++m)
					B[i][m] = refB[i][m];
			
			//reset Pi
			for (j = 0; j < N; ++j)
				Pi[j] = refPi[j];

			//set new O
			O = fullO[k];

			isUpdated = true;

			//update new pStar
			runViterbi(0);

			//run till the model is being updated
			while(isUpdated)
			{
				runForwardBackward();
				runBaumWelch();
				runViterbi(1);
				isUpdated = compareAndUpdateModel();
			}
			
			//print the new models to respective files
			for (i = 0; i < N; ++i)
			{
				for (j = 0; j < N; ++j)
					fprintf(AComplementFile,"%g ",A[i][j]);
				fprintf(AComplementFile,"\n");
			}
			
			for (i = 0; i < N; ++i)
			{
				for (m = 0; m < M; ++m)
					fprintf(BComplementFile,"%g ",B[i][m]);
				fprintf(BComplementFile,"\n");
			}
			
			for (j = 0; j < N; ++j)
				fprintf(PiComplementFile,"%g ",Pi[j]);
			fprintf(PiComplementFile,"\n");

			fclose(AComplementFile);
			fclose(BComplementFile);
			fclose(PiComplementFile);
		}
	}

	printf("Breakpoint\n");
	return 0;
*/
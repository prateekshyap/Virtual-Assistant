// FinalProject.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "commonvar.h"
#include "forwardbackward.h"
#include "viterbi.h"
#include "baumwelch.h"
#include "functions.h"

extern long double *** xi;

extern long double ** alpha, ** beta, ** delta, ** gamma;

extern long double ** A, ** AComplement, ** B, ** BComplement;

extern long double * Pi, * PiComplement;

extern int ** psi;

extern int * O, * qStar, * qStarComplement;

extern FILE * AComplementFile, * BComplementFile, * PiComplementFile;

extern int T, N, M;

int D = 3, F = 4;
int ** fullO = NULL;

int _tmain(int argc, _TCHAR* argv[])
{
	bool isUpdated = true;
	int i = 0, j = 0, k = 0, r = 0, c = 0, d = 0, f = 0, t = 0, m = 0;
	long double ** refA = NULL, ** refB = NULL; 
	long double * refPi = NULL;

	/*define everything*/
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
	fullO = new int *[N];
	Pi = new long double[N]; //Pi Array

	for (i = 0; i < N; ++i)
		A[i] = new long double[N];
	for (i = 0; i < N; ++i)
		B[i] = new long double[M];
	for (i = 0; i < N; ++i)
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
	FILE * file = NULL;

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
*/
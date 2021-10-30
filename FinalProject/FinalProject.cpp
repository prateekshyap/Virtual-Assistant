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
//#include "viterbi.h"
//#include "baumwelch.h"

//for viterbi
int D = 3, F = 4;
extern FILE * viterbiOutputFile;

//for baum welch

int _tmain(int argc, _TCHAR* argv[])
{

	/*Find Ci*/


	/*LBG*/


	/*Forward and Backward Process*/


	/*Viterbi*/
	//viterbi algorithm will be called here


	/*Baum Welch*/


	
}



/*

////////////////////////////////////////
Working main code for viterbi-
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

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

		viterbiOutputFile = fopen(outputFileName,"w+");
		printf("-----------------------------------------------------------\n");
		printf("                         Digit %s\n",digits[d]);
		printf("-----------------------------------------------------------\n");
		fprintf(viterbiOutputFile,"-----------------------------------------------------------\n");
		fprintf(viterbiOutputFile,"                         Digit %s\n",digits[d]);
		fprintf(viterbiOutputFile,"-----------------------------------------------------------\n");

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
				fprintf(viterbiOutputFile,"Reading %s for A matrix\n",fileName);
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
				fprintf(viterbiOutputFile,"Reading %s for B matrix\n",fileName);
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
				fprintf(viterbiOutputFile,"Reading %s for observation sequences\n",fileName);
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
				fprintf(viterbiOutputFile,"Reading %s for Pi\n",fileName);
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
		fclose(viterbiOutputFile);
	}

	printf("Breakpoint\n");
	return 0;


*/
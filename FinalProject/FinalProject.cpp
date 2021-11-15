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
#include <direct.h>


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
int trainFlag = 1;
char trainUserString[128];

void trainBeginningModel(int);
void generateTestReport(char *, char **, char **, int);
int preProcess();
void postProcess();
void trainNewModel();
void searchWord(char **);

int _tmain(int argc, _TCHAR* argv[])
{
	int d = 0, r = 0, range = 0, i = 0, testChoice = 1, userChoice = -1, currTrainCount = 0;
	char * temp = NULL;
	int trainingStatus = -1; //set a variable to know if the model has to be trained or not
	int buildingStatus = -1; //set a variable to know if codebook needs to be built or not
	int readyStatus = -1;
	char command[200];
	char * testFileName = NULL, * userString = NULL;
	char trainFileName[256], filePath[128];

	range = preProcess();

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

	FILE * file = fopen("data/info.txt","r"); //open info.txt
	fscanf(file,"%d",&trainingStatus); //read from the file
	fclose(file);

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

	//operations
	while (true)
	{
		printf("============================================================\n");
		printf("                   Tutorial for Kids\n");
		printf("============================================================\n");
		printf("1 - Training mode\n");
		printf("2 - Testing mode\n");
		printf("3 - Quit\n");
		printf("Enter the corresponding number-\n");

		scanf("%d",&userChoice);

		switch (userChoice)
		{
		case 1:
			searchWord(digits);
			if (trainFlag == 0) //existing word
			{
				sprintf(trainFileName,"HMM/%s/count.txt",trainUserString);
				file = fopen(trainFileName,"r");
				fscanf(file,"%d",&currTrainCount);
				fclose(file);
				sprintf(trainFileName,"HMM/%s/%d.txt",trainUserString,currTrainCount);
				sprintf(command,"Recording_Module.exe %d data/o.wav %s",duration,trainFileName);
				std :: system(command);
				trainModel(trainUserString,range,trainFileName);
			}
			else //new word
			{
				printf("--------------Warning----------------\n");
				printf("It will take a few minutes\n");
				printf("Please record %d times\n",R);
				sprintf(filePath,"HMM/%s",trainUserString);
				mkdir(filePath);
				for (i = 0; i < R; ++i)
				{
					printf("Enter 1 when you're ready to record\n");
					scanf("%d",&readyStatus);
					sprintf(command,"Recording_Module.exe %d data/o.wav %s/%d.txt",duration,filePath,i);
					std :: system(command);
					printf("%d files recorded\n",(i+1));
				}

				//update the set of words in file as well as in the array
				file = fopen("data/D.txt","w");
				++D;
				fprintf(file,"%d\n",D);
				for (i = 0; i < D-1; ++i)
					fprintf(file,"%s\n",digits[i]);
				fprintf(file,"%s\n",trainUserString);
				fclose(file);
				file = fopen("data/D.txt","r");
				digits = new char * [D];
				for (d = 0; d < D; ++d)
					digits[d] = (char *)malloc(sizeof(char *));
				fscanf(file,"%d",&D);
				for (d = 0; d < D; ++d)
					fscanf(file,"%s",digits[d]);
				fclose(file);

				//write bakis model in the folder
				writeGeneralModel(filePath,R);

				//train the model from the beginning along with the live trained existing words
				range = preProcess();
				trainNewModel();
			}
			break;

		case 2:
			testFileName = "data/o.txt";
			sprintf(command,"Recording_Module.exe %d data/o.wav data/o.txt",duration);
			std :: system(command);
			recognize(folder, digits, files, dataFiles, D, R, range, testFileName);
			printf("You spoke %s\n",resultWord);
			/* Put image display feature here */
			break;

		case 3:
			printf("Thank you for using the app.");
			return 0;

		default:
			printf("Invalid input! Please try again!");
		}
	}
	
	return 0;
}

void searchWord(char ** digits)
{
	int i = 0;
	printf("Enter the word for which you want to train-\n");
	scanf("%s",&trainUserString);
	trainFlag = 1;
	for (i = 0; i < D; ++i)
	{
		if (strcmp(trainUserString,digits[i]) == 0)
		{
			trainFlag = 0;
			break;
		}
	}
}

int preProcess()
{	
	int range = 0;

	//read the duration
	FILE * file = fopen("data/duration.txt","r");
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
	file = fopen("data/F.txt","r");
	fscanf(file,"%d",&R);
	fclose(file);

	return range;
}

void postProcess()
{
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
	int range = 0, progressCount = 0, percentageSize = 50, modelLoop = 2, updationLimit = 0;
	long double value = 0;

	


	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Reading file information
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//Read the updated information from info.txt
	file = fopen("data/F.txt","r");
	fscanf(file,"%d",&R); F = R;
	fclose(file);
	
	//reading the existing words
	file = fopen("data/D.txt","r");
	fscanf(file,"%d",&D);

	char * folder = "HMM/";
	char ** digits = new char * [D];
	for (d = 0; d < D; ++d)
		digits[d] = (char *)malloc(sizeof(char *));
	for (d = 0; d < D; ++d)
		fscanf(file,"%s",digits[d]);
	fclose(file);

	//resetting the D file for the first 15 words
	D = 17;
	file = fopen("data/D.txt","w");
	fprintf(file,"%d\n",D);
	for (d = 0; d < D; ++d)
		fprintf(file,"%s\n",digits[d]);
	fclose(file);

	//updating the digits array with the old 15 words
	file = fopen("data/D.txt","r");
	fscanf(file,"%d",&D);
	digits = new char * [D];
	for (d = 0; d < D; ++d)
		digits[d] = (char *)malloc(sizeof(char *));
	for (d = 0; d < D; ++d)
		fscanf(file,"%s",digits[d]);
	fclose(file);

	char ** files = NULL, ** dataFiles = NULL;

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
	resetCount(folder, digits, dataFiles, D, R);
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
	define();

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
	printf("\nBuilding the model\n\n");
	while (modelCount != modelLoop) //for a few times
	{
		++modelCount;
		//printf("%d\n",modelCount);
		//converge the model for each observation sequence
		for (d = 0; d < D; ++d) //for each digit
		{
			generateObservationSequences(folder, digits, dataFiles, d, R, range); //generate the full observation sequence
			//read the model
			for (f = 0; f < 4; ++f) //for each file
			{
				sprintf(fileName,"%s%s/%s.txt",folder,digits[d],files[f]);
			
				file = fopen(fileName,"r"); //open the file and read it

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
				sprintf(AComplementFileName,"%s%s/AComplement%d.txt",folder,digits[d],k);
				sprintf(BComplementFileName,"%s%s/BComplement%d.txt",folder,digits[d],k);
				sprintf(PiComplementFileName,"%s%s/PiComplement%d.txt",folder,digits[d],k);

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
			sprintf(fileName,"%s%s/A.txt",folder,digits[d]);
			AFile = fopen(fileName,"w");

			sprintf(fileName,"%s%s/B.txt",folder,digits[d]);
			BFile = fopen(fileName,"w");

			sprintf(fileName,"%s%s/Pi.txt",folder,digits[d]);
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
				sprintf(AComplementFileName,"%s%s/AComplement%d.txt",folder,digits[d],k);
				sprintf(BComplementFileName,"%s%s/BComplement%d.txt",folder,digits[d],k);
				sprintf(PiComplementFileName,"%s%s/PiComplement%d.txt",folder,digits[d],k);

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

void trainNewModel()
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
	int range = 0, progressCount = 0, percentageSize = 50, modelLoop = 2, updationLimit = 0;
	long double value = 0;

	


	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Reading file information
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//Read the updated information from info.txt
	file = fopen("data/F.txt","r");
	fscanf(file,"%d",&R);
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
	printf("Building the Universe\n");
	buildUniverse(folder, digits, dataFiles, D, R, range);
	printf("\n\nUniverse building completed\n\nBuilding the Codebook\n");
	buildCodebook();
	printf("\n\nCodebook building completed\n");
	
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Define variables for HMM
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//define everything
	define();

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
	printf("\nBuilding the model\n\n");
	while (modelCount != modelLoop) //for a few times
	{
		++modelCount;
		//printf("%d\n",modelCount);
		//converge the model for each observation sequence
		for (d = 0; d < D; ++d) //for each digit
		{
			sprintf(fileName,"%s%s/count.txt",folder,digits[d]);

			file = fopen(fileName,"r");

			fscanf(file,"%d",&F);

			fclose(file);

			generateObservationSequences(folder, digits, dataFiles, d, F, range); //generate the full observation sequence
			//read the model
			for (f = 0; f < 4; ++f) //for each file
			{
				sprintf(fileName,"%s%s/%s.txt",folder,digits[d],files[f]);
			
				file = fopen(fileName,"r"); //open the file and read it

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
					fullO = new int *[F];
					for (i = 0; i < F; ++i)
						fullO[i] = new int[T];
					for (i = 0; i < F; ++i)
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

			for (k = 0; k < F; ++k)
			{
				sprintf(AComplementFileName,"%s%s/AComplement%d.txt",folder,digits[d],k);
				sprintf(BComplementFileName,"%s%s/BComplement%d.txt",folder,digits[d],k);
				sprintf(PiComplementFileName,"%s%s/PiComplement%d.txt",folder,digits[d],k);

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
			sprintf(fileName,"%s%s/count.txt",folder,digits[d]);

			file = fopen(fileName,"r");

			fscanf(file,"%d",&F);

			fclose(file);

			sprintf(fileName,"%s%s/A.txt",folder,digits[d]);
			AFile = fopen(fileName,"w");

			sprintf(fileName,"%s%s/B.txt",folder,digits[d]);
			BFile = fopen(fileName,"w");

			sprintf(fileName,"%s%s/Pi.txt",folder,digits[d]);
			PiFile = fopen(fileName,"w");

			for (i = 0; i < N; ++i)
				for (j = 0; j < N; ++j)
					AComplement[i][j] = 0;

			for (i = 0; i < N; ++i)
				for (j = 0; j < M; ++j)
					BComplement[i][j] = 0;

			for (i = 0; i < N; ++i)
				PiComplement[i] = 0;

			for (k = 0; k < F; ++k)
			{
				sprintf(AComplementFileName,"%s%s/AComplement%d.txt",folder,digits[d],k);
				sprintf(BComplementFileName,"%s%s/BComplement%d.txt",folder,digits[d],k);
				sprintf(PiComplementFileName,"%s%s/PiComplement%d.txt",folder,digits[d],k);

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
					AComplement[i][j] /= F;

			for (i = 0; i < N; ++i)
				for (j = 0; j < M; ++j)
					BComplement[i][j] /= F;

			for (i = 0; i < N; ++i)
				PiComplement[i] /= F;

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
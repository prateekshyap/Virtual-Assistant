extern FILE * AComplementFile, * BComplementFile, * PiComplementFile;
extern long double pStar, pStarComplement, pOfOGivenLambda;
extern long double ** AComplement, ** BComplement, ** A, ** B, ** codebook;
extern long double * PiComplement, * Pi;
extern int * O;
extern int N, M, T, duration, p, universeSize;
extern char * resultWord;

int noOfFrames = 100, samplingRate = 16000;

double * getRValues(double * x, int sampleSize, int p)
{
	double * R = new double[13];
	int i = 0, j = 0;
	for (i = 0; i <= p; ++i)
	{
		R[i] = 0;
		for (j = 0; j < sampleSize-i; ++j)
			R[i] += x[j]*x[i+j];
	}
	return R;
}

double * getAValues(double * R, int sampleSize, int p)
{
	//following Durbin's Algorithm
	int i = 0, j = 0;
	double * AA = new double[p];
	double * E = new double[p+1];
	double * K = new double[p+1];
	double ** tempA = new double*[p+1];

	//initializing everything to 0
	for (i = 0; i < p+1; ++i)
		tempA[i] = new double[p+1];
	for (i = 0; i < p+1; ++i)
	{
		E[i] = 0;
		K[i] = 0;
	}
	for (i = 0; i < p+1; ++i)
		for (j = 0; j < p+1; ++j)
			tempA[i][j] = 0;

	//Algorithm steps
	E[0] = R[0];
	for (i = 1; i <= p; ++i)
	{
		K[i] = 0;
		for (j = 1; j < i; ++j)
			K[i] += tempA[i-1][j]*R[i-j];
		K[i] = R[i]-K[i];
		K[i] /= E[i-1];
		tempA[i][i] = K[i];
		for (j = 1; j < i; ++j)
			tempA[i][j] = tempA[i-1][j]-(K[i]*tempA[i-1][i-j]);
		E[i] = (1-(K[i]*K[i]))*E[i-1];
	}
	for (i = 1; i <= p; ++i)
		AA[i-1] = tempA[p][i];
	return AA;
}

double * getCValues(double * AA, int sampleSize, int p, double r0)
{
	int i = 0, j = 0;
	double * C = new double[p+1];
	C[0] = log(r0*r0);
	for (i = 1; i < p+1; ++i)
	{
		C[i] = 0;
		double m = i;
		for (j = 1; j < i; ++j)
		{
			double k = j;
			C[i] += (k/m)*C[j]*AA[i-j-1];
		}
		C[i] += AA[i-1];
	}
	return C;
}

void buildUniverse(char * folder, char * digits[], char * files[], int D, int R, int range)
{
	int d = 0, f = 0, r = 0, i = 0, j = 0, k = 0, q = p, index = 0, prog = 0;
	int xCount = duration*samplingRate, dcCount = 0, shift = 0, temp = 0, sampleSize = xCount/noOfFrames;
	char * buffer = new char[1024];
	double * x = new double[xCount], * RR = NULL, * AA = NULL, * C = NULL, * trainingData = new double[sampleSize];
	double DCShift = 0, maxData = 0, minData = 0, normalizationFactor = 0, copy = 0;
	char fileName[100];
	FILE * universeFile = fopen("data/universe.txt","w+"), * file = NULL;

	/////////////////////////////////////////////////////////////////
	//Implementing a progress bar
	int nop = D*R*noOfFrames;
	int percentageBreak = D*R;
	int progressCount = 0;
	/////////////////////////////////////////////////////////////////
	printf("\n");
	for (d = 0; d < D; ++d)
	{		
		for (f = 0; f < R; ++f)
		{
			sprintf(fileName,"%s%s/%s.txt",folder,digits[d],files[f]);

			file = fopen(fileName,"r");
			dcCount = 0;
			for (i = 0; i < 1000; ++i) //for all x values
			{
				fgets(buffer, 1024, file); //read one line from the file
				x[i] = atof(buffer); //convert it into float
				DCShift += x[i]; ++dcCount; //add it to dc shift
			}
			DCShift /= dcCount; //find out dcshift
			for (i = 0; i < 1000; ++i)
			{
				x[i] -= DCShift; //subtract dcshift from already stored values
				if (x[i] > maxData) maxData = x[i]; //update maxData
				if (x[i] < minData) minData = x[i]; //update minData
			}
			for (i = 1000; i < xCount; ++i) //for rest values
			{
				fgets(buffer, 1024, file); //read one line from the file
				x[i] = atof(buffer); //convert it into float
				x[i] -= DCShift; //subtract dcshift
				if (x[i] > maxData) maxData = x[i]; //update maxData
				if (x[i] < minData) minData = x[i]; //update minData
			}
			fclose(file);

			//find out normalization factor
			minData = abs(minData); //get the absolute value of the minimum
			normalizationFactor = range/((maxData+minData)/2); //get the normalization factor
			//truncate the normalizationFactor to one decimal digit
			copy = normalizationFactor; //store the normalizationFactor in copy
			while (copy < 1) //till copy is less than 1
			{
				++shift; //increment shift
				copy *= 10; //keep multiplying copy with 10
			}
			temp = (int)copy; //store the integer part of copy in temp
			copy -= temp; //subtract temp from copy
			while (shift-- > 0) //reverse the shift process
				copy /= 10;
			/*the above process helps in getting all the digits except the first decimal digit*/
			normalizationFactor -= copy; //subtract copy from normalizationFactor
			/*this will ensure that only one digit remains after decimal point*/

			//multiply normalization factor
			for (i = 0; i < xCount; ++i)
				x[i] *= normalizationFactor;

			index = 0;

			for (r = 0; r < noOfFrames; ++r) //indicates frame number
			{
				for (i = 0; i < sampleSize; ++i)
					trainingData[i] = x[index++];

				RR = getRValues(trainingData,sampleSize,p);
				AA = getAValues(RR,sampleSize,p);
				C = getCValues(AA,sampleSize,p,RR[0]);
				
				//raised sine window on C and adding it to the universe file
				for (j = 1; j <= p; ++j)
				{
					C[j] *= 1+((q/2)*sin((3.14*j)/q));
					fprintf(universeFile,"%f ",C[j]);
				}
				++universeSize;
				fprintf(universeFile,"\n");
				
				//printing the progress bar
				++progressCount;
				if (progressCount%percentageBreak == 0)
				{
					printf("\r"); //move to the beginning of the line
					/*print the progress*/
					prog = 0;
					printf("[");
					for (prog = percentageBreak; prog <= progressCount; prog += percentageBreak)
						printf("|");
					for (; prog <= nop; prog += percentageBreak)
						printf(" ");
					printf("]%d",(progressCount/percentageBreak)); //print the percentage
				}
			}
		}
	}
	
	fclose(universeFile);
	delete(x);
	delete(buffer);
	delete(RR);
	delete(AA);
	delete(C);
	delete(trainingData);
}

bool compareAndUpdateModel()
{
	int i = 0, j = 0, t = 0, m = 0;
	if (pStarComplement > pStar)
	//if (pStarComplement < 1e-30)
	//if (abs(pStarComplement-pStar) > 1e-2)
	{
		//update A
		for (i = 0; i < N; ++i)
			for (j = 0; j < N; ++j)
				A[i][j] = AComplement[i][j];
		
		//update B
		for (i = 0; i < N; ++i)
			for (m = 0; m < M; ++m)
				B[i][m] = BComplement[i][m];
		
		//update Pi
		for (j = 0; j < N; ++j)
			Pi[j] = PiComplement[j];
	
		//update pStar
		pStar = pStarComplement;
		//printf("%g\n",pStarComplement);

		/*//reinitialize complement matrices
		PiComplement = new long double [N];

		AComplement = new long double * [N];
		for (i = 0; i < N; ++i)
			AComplement[i] = new long double[N];

		BComplement = new long double * [N];
		for (i = 0; i < N; ++i)
			BComplement[i] = new long double[M];
			*/
		return true;
	}
	else
	{
		//printf("PStar updated- %g\n",pStar);
		return false;
	}
}

void generateObservationSequences(char * folder, char * digits[], char * files[], int digitIndex, int F, int range)
{
	int f = 0, r = 0, i = 0, j = 0, k = 0, q = p, index = 0, codeBookIndex = 0;
	int xCount = duration*samplingRate, dcCount = 0, shift = 0, temp = 0, sampleSize = xCount/noOfFrames, minIndex = -1;
	char * buffer = new char[1024];
	double * x = new double[xCount], * RR = NULL, * AA = NULL, * C = NULL, * trainingData = new double[sampleSize];
	double DCShift = 0, maxData = 0, minData = 0, normalizationFactor = 0, copy = 0;
	int ** Ob = new int*[F];
	long double * distances = new long double[M]; //stores the tokhura distances
	long double w[] = {1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0}; //given Tokhura weights
	long double min = 1e30, temp1 = 0, temp2 = 0;
	char fileName[100];
	FILE * file = NULL;
	for (i = 0; i < F; ++i)
		Ob[i] = new int[noOfFrames];	
	for (f = 0; f < F; ++f)
	{
		sprintf(fileName,"%s%s/%d.txt",folder,digits[digitIndex],f);

		file = fopen(fileName,"r");
		dcCount = 0;
		for (i = 0; i < 1000; ++i) //for all x values
		{
			fgets(buffer, 1024, file); //read one line from the file
			x[i] = atof(buffer); //convert it into float
			DCShift += x[i]; ++dcCount; //add it to dc shift
		}
		DCShift /= dcCount; //find out dcshift
		for (i = 0; i < 1000; ++i)
		{
			x[i] -= DCShift; //subtract dcshift from already stored values
			if (x[i] > maxData) maxData = x[i]; //update maxData
			if (x[i] < minData) minData = x[i]; //update minData
		}
		for (i = 1000; i < xCount; ++i) //for rest values
		{
			fgets(buffer, 1024, file); //read one line from the file
			x[i] = atof(buffer); //convert it into float
			x[i] -= DCShift; //subtract dcshift
			if (x[i] > maxData) maxData = x[i]; //update maxData
			if (x[i] < minData) minData = x[i]; //update minData
		}
		fclose(file);

		//find out normalization factor
		minData = abs(minData); //get the absolute value of the minimum
		normalizationFactor = range/((maxData+minData)/2); //get the normalization factor

		//truncate the normalizationFactor to one decimal digit
		copy = normalizationFactor; //store the normalizationFactor in copy
		while (copy < 1) //till copy is less than 1
		{
			++shift; //increment shift
			copy *= 10; //keep multiplying copy with 10
		}
		temp = (int)copy; //store the integer part of copy in temp
		copy -= temp; //subtract temp from copy
		while (shift-- > 0) //reverse the shift process
			copy /= 10;
		/*the above process helps in getting all the digits except the first decimal digit*/
		normalizationFactor -= copy; //subtract copy from normalizationFactor
		/*this will ensure that only one digit remains after decimal point*/

		//multiply normalization factor
		for (i = 0; i < xCount; ++i)
			x[i] *= normalizationFactor;

		index = 0;
		for (r = 0; r < noOfFrames; ++r) //indicates frame number
		{
			for (i = 0; i < sampleSize; ++i)
				trainingData[i] = x[index++];

			RR = getRValues(trainingData,sampleSize,p);
			AA = getAValues(RR,sampleSize,p);
			C = getCValues(AA,sampleSize,p,RR[0]);
			
			//raised sine window on C and adding it to the universe file
			for (j = 1; j <= p; ++j)
				C[j] *= 1+((q/2)*sin((3.14*j)/q));

			//Tokhura Distance from codebook
			file = fopen("data/codebook.txt","r");
			for (j = 0; j < M; ++j) //for each vector in the codebook
			{
				distances[j] = 0;
				for (k = 0; k < p; ++k) //for each p value
				{
					temp1 = (long double)C[k+1];
					fscanf(file,"%lf",&temp2);
					distances[j] += w[k]*(temp1-temp2)*(temp1-temp2); //find out the distance
				}
			}
			fclose(file);
			//Minimum index
			min = 1e30; //reset the minimum
			for (j = 0; j < M; ++j)
			{
				if (min > distances[j])
				{
					min = distances[j];
					minIndex = j;
				}
			}

			Ob[f][r] = minIndex+1; //+1 is done to maintain the 1-32 range format and not 0-31 since this subtraction is being done at later steps
		}
	}

	sprintf(fileName,"%s%s/O.txt",folder,digits[digitIndex]);
	file = fopen(fileName,"w+");
	for (i = 0; i < F; ++i)
	{
		for (j = 0; j < noOfFrames; ++j)
		{
			fprintf(file,"%d ",Ob[i][j]);
		}
		fprintf(file,"\n");
	}
	fclose(file);
	delete(x);
	delete(buffer);
	delete(Ob);
	delete(RR);
	delete(AA);
	delete(C);
	delete(trainingData);
}

void resetCount(char * folder, char * digits[], char * files[], int D, int R)
{
	int d = 0;
	char fileName[100];
	FILE * countFile = NULL;
	for (d = 0; d < D; ++d)
	{
		sprintf(fileName,"%s%s/count.txt",folder,digits[d]);
		countFile = fopen(fileName,"w");

		fprintf(countFile,"%d\n",R);

		fclose(countFile);
	}
}

void resetModel(char * folder, char * digits[], char * files[], int D, int R)
{
	int d = 0, k = 0, i = 0, j = 0;
	FILE * AFile = NULL, * BFile = NULL, * PiFile = NULL;
	char fileName[100];
	long double * row = new long double[N];
	long double tempB = 0.03125, one = 1, zero = 0;
	row[0] = 0.8;
	row[1] = 0.2;
	row[2] = row[3] = row[4] = 0;
	for (d = 0; d < D; ++d)
	{
		sprintf(fileName,"%s%s/A.txt",folder,digits[d]);
		AFile = fopen(fileName,"w");

		sprintf(fileName,"%s%s/B.txt",folder,digits[d]);
		BFile = fopen(fileName,"w");

		sprintf(fileName,"%s%s/Pi.txt",folder,digits[d]);
		PiFile = fopen(fileName,"w");

		tempB = 0.03125, one = 1, zero = 0;
		row[0] = 0.8;
		row[1] = 0.2;
		row[2] = row[3] = row[4] = 0;

		for (i = 0; i < N-1; ++i)
		{
			for (j = 0; j < N; ++j)
				fprintf(AFile,"%g ",row[j]);
			fprintf(AFile,"\n");
			for (j = N-2; j >= 0; --j)
				row[j+1] = row[j];
			row[0] = 0;
		}
		for (i = 0; i < N-1; ++i)
			fprintf(AFile,"%g ",zero);
		fprintf(AFile,"%g \n",one);

		for (i = 0; i < N; ++i)
		{
			for (j = 0; j < M; ++j)
				fprintf(BFile,"%g ",tempB);
			fprintf(BFile,"\n");
		}

		fprintf(PiFile,"%g ",one);
		for (i = 1; i < N; ++i)
			fprintf(PiFile,"%g ",zero);
		fprintf(PiFile,"\n");

		fclose(AFile);
		fclose(BFile);
		fclose(PiFile);
	}
}

void recognize(char * folder, char * digits[], char * files[], char * dataFiles[], int D, int R, int range, char * testFileName)
{
	int d = 0, f = 0, r = 0, i = 0, j = 0, k = 0, q = p, index = 0;
	int xCount = duration*samplingRate, dcCount = 0, shift = 0, temp = 0, sampleSize = xCount/noOfFrames, minIndex = -1;
	char * buffer = new char[1024];
	double * x = new double[xCount], * RR = NULL, * AA = NULL, * C = NULL, * trainingData = new double[sampleSize];
	double DCShift = 0, maxData = 0, minData = 0, normalizationFactor = 0, copy = 0;
	int * Ob = new int[noOfFrames];
	long double * distances = new long double[M]; //stores the tokhura distances
	long double w[] = {1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0}; //given Tokhura weights
	long double min = 1e30, temp1 = 0, temp2 = 0, maxProb = 0;
	char fileName[200], * resultantWord = NULL;
	FILE * file = NULL, * AFile = NULL, * BFile = NULL, * PiFile = NULL;	

	printf("%s\n",testFileName);
	file = fopen(testFileName,"r");
	dcCount = 0;
	for (i = 0; i < 1000; ++i) //for all x values
	{
		fgets(buffer, 1024, file); //read one line from the file
		x[i] = atof(buffer); //convert it into float
		DCShift += x[i]; ++dcCount; //add it to dc shift
	}
	DCShift /= dcCount; //find out dcshift
	for (i = 0; i < 1000; ++i)
	{
		x[i] -= DCShift; //subtract dcshift from already stored values
		if (x[i] > maxData) maxData = x[i]; //update maxData
		if (x[i] < minData) minData = x[i]; //update minData
	}
	for (i = 1000; i < xCount; ++i) //for rest values
	{
		fgets(buffer, 1024, file); //read one line from the file
		x[i] = atof(buffer); //convert it into float
		x[i] -= DCShift; //subtract dcshift
		if (x[i] > maxData) maxData = x[i]; //update maxData
		if (x[i] < minData) minData = x[i]; //update minData
	}
	fclose(file);

	//find out normalization factor
	minData = abs(minData); //get the absolute value of the minimum
	normalizationFactor = range/((maxData+minData)/2); //get the normalization factor

	//truncate the normalizationFactor to one decimal digit
	copy = normalizationFactor; //store the normalizationFactor in copy
	while (copy < 1) //till copy is less than 1
	{
		++shift; //increment shift
		copy *= 10; //keep multiplying copy with 10
	}
	temp = (int)copy; //store the integer part of copy in temp
	copy -= temp; //subtract temp from copy
	while (shift-- > 0) //reverse the shift process
		copy /= 10;
	/*the above process helps in getting all the digits except the first decimal digit*/
	normalizationFactor -= copy; //subtract copy from normalizationFactor
	/*this will ensure that only one digit remains after decimal point*/

	//multiply normalization factor
	for (i = 0; i < xCount; ++i)
		x[i] *= normalizationFactor;
	
	index = 0;
	
	for (r = 0; r < noOfFrames; ++r) //indicates frame number
	{
		for (i = 0; i < sampleSize; ++i) //copy the next frame to training data
			trainingData[i] = x[index++];

		//find out Ci
		RR = getRValues(trainingData,sampleSize,p);
		AA = getAValues(RR,sampleSize,p);
		C = getCValues(AA,sampleSize,p,RR[0]);
			
		//raised sine window on C and adding it to the universe file
		for (j = 1; j <= p; ++j)
			C[j] *= 1+((q/2)*sin((3.14*j)/q));

		//Tokhura Distance from codebook
		file = fopen("data/codebook.txt","r");
		for (j = 0; j < M; ++j) //for each vector in the codebook
		{
			distances[j] = 0;
			for (k = 0; k < p; ++k) //for each p value
			{
				temp1 = (long double)C[k+1];
				fscanf(file,"%lf",&temp2);
				distances[j] += w[k]*(temp1-temp2)*(temp1-temp2); //find out the distance
			}
		}
		
		fclose(file);
		//Minimum index
		min = 1e30; //reset the minimum
		for (j = 0; j < M; ++j)
		{
			if (min > distances[j])
			{
				min = distances[j];
				minIndex = j;
			}
		}

		Ob[r] = minIndex+1; //+1 is done to maintain the 1-32 range format and not 0-31 since this subtraction is being done at later steps
	}

	O = Ob;
	for (d = 0; d < D; ++d)
	{
		strcpy(fileName,folder);
		strcat(fileName,digits[d]);
		strcat(fileName,"/A.txt");
		AFile = fopen(fileName,"r");
		strcpy(fileName,folder);
		strcat(fileName,digits[d]);
		strcat(fileName,"/B.txt");
		BFile = fopen(fileName,"r");
		strcpy(fileName,folder);
		strcat(fileName,digits[d]);
		strcat(fileName,"/Pi.txt");
		PiFile = fopen(fileName,"r");

		for (i = 0; i < N; ++i)
			for (j = 0; j < N; ++j)
				fscanf(AFile,"%lf",&A[i][j]);

		for (i = 0; i < N; ++i)
			for (j = 0; j < M; ++j)
				fscanf(BFile,"%lf",&B[i][j]);

		for (i = 0; i < N; ++i)
			fscanf(PiFile,"%lf",&Pi[i]);

		runForwardBackward();

		printf("d = %s, %g\n",digits[d],pOfOGivenLambda);
		if (pOfOGivenLambda > maxProb)
		{
			maxProb = pOfOGivenLambda;
			resultantWord = digits[d];
		}

		fclose(AFile);
		fclose(BFile);
		fclose(PiFile);
	}

	delete(x);
	delete(buffer);
	delete(Ob);
	delete(RR);
	delete(AA);
	delete(C);
	delete(trainingData);
	resultWord = resultantWord;
}

void trainModel(char * word, int range, char * trainFileName)
{
	bool isUpdated = true;
	int xCount = duration*samplingRate, dcCount = 0, shift = 0, temp = 0, q = p, sampleSize = xCount/noOfFrames, minIndex = -1, index = 0, codebookIndex = 0, updationLimit = 0, currCount = -1;
	char * buffer = new char[1024];
	double * x = new double[xCount], * RR = NULL, * AA = NULL, * C = NULL, * trainingData = new double[sampleSize];
	double DCShift = 0, maxData = 0, minData = 0, normalizationFactor = 0, copy = 0;
	O = new int[noOfFrames];
	long double * distances = new long double[M]; //stores the tokhura distances
	long double w[] = {1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0}; //given Tokhura weights
	long double min = 1e30, temp1 = 0, temp2 = 0, tempValue = 0;
	FILE * file = NULL;
	
	int i = 0, j = 0, r = 0, k = 0, m = 0;
	char filePath[128];
	sprintf(filePath,"HMM/%s/",word);
	char fileName[256];

	//define variables
	define();

	//set A B and Pi to Bakis model
	long double * row = new long double[N];
	long double tempB = 0.03125, one = 1, zero = 0;
	row[0] = 0.8;
	row[1] = 0.2;
	row[2] = row[3] = row[4] = 0;

	for (i = 0; i < N-1; ++i)
	{
		for (j = 0; j < N; ++j)
			A[i][j] = row[j];
		for (j = N-2; j >= 0; --j)
			row[j+1] = row[j];
		row[0] = 0;
	}
	for (i = 0; i < N-1; ++i)
		A[N-1][i] = zero;
	A[N-1][N-1] = one;

	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < M; ++j)
			B[i][j] = tempB;
	}

	Pi[0] = 1;
	for (i = 1; i < N; ++i)
		Pi[i] = 0;

	//generate observation sequence for the spoken word
	file = fopen(trainFileName,"r");
	dcCount = 0;
	for (i = 0; i < 1000; ++i) //for all x values
	{
		fgets(buffer, 1024, file); //read one line from the file
		x[i] = atof(buffer); //convert it into float
		DCShift += x[i]; ++dcCount; //add it to dc shift
	}
	DCShift /= dcCount; //find out dcshift
	for (i = 0; i < 1000; ++i)
	{
		x[i] -= DCShift; //subtract dcshift from already stored values
		if (x[i] > maxData) maxData = x[i]; //update maxData
		if (x[i] < minData) minData = x[i]; //update minData
	}
	for (i = 1000; i < xCount; ++i) //for rest values
	{
		fgets(buffer, 1024, file); //read one line from the file
		x[i] = atof(buffer); //convert it into float
		x[i] -= DCShift; //subtract dcshift
		if (x[i] > maxData) maxData = x[i]; //update maxData
		if (x[i] < minData) minData = x[i]; //update minData
	}
	fclose(file);

	//find out normalization factor
	minData = abs(minData); //get the absolute value of the minimum
	normalizationFactor = range/((maxData+minData)/2); //get the normalization factor

	//truncate the normalizationFactor to one decimal digit
	copy = normalizationFactor; //store the normalizationFactor in copy
	while (copy < 1) //till copy is less than 1
	{
		++shift; //increment shift
		copy *= 10; //keep multiplying copy with 10
	}
	temp = (int)copy; //store the integer part of copy in temp
	copy -= temp; //subtract temp from copy
	while (shift-- > 0) //reverse the shift process
		copy /= 10;
	/*the above process helps in getting all the digits except the first decimal digit*/
	normalizationFactor -= copy; //subtract copy from normalizationFactor
	/*this will ensure that only one digit remains after decimal point*/

	//multiply normalization factor
	for (i = 0; i < xCount; ++i)
		x[i] *= normalizationFactor;

	index = 0;
	for (r = 0; r < noOfFrames; ++r) //indicates frame number
	{
		for (i = 0; i < sampleSize; ++i)
			trainingData[i] = x[index++];

		RR = getRValues(trainingData,sampleSize,p);
		AA = getAValues(RR,sampleSize,p);
		C = getCValues(AA,sampleSize,p,RR[0]);
			
		//raised sine window on C and adding it to the universe file
		for (j = 1; j <= p; ++j)
			C[j] *= 1+((q/2)*sin((3.14*j)/q));

		//Tokhura Distance from codebook
		file = fopen("data/codebook.txt","r");
		for (j = 0; j < M; ++j) //for each vector in the codebook
		{
			distances[j] = 0;
			for (k = 0; k < p; ++k) //for each p value
			{
				temp1 = (long double)C[k+1];
				fscanf(file,"%lf",&temp2);
				distances[j] += w[k]*(temp1-temp2)*(temp1-temp2); //find out the distance
			}
		}
		fclose(file);
		//Minimum index
		min = 1e30; //reset the minimum
		for (j = 0; j < M; ++j)
		{
			if (min > distances[j])
			{
				min = distances[j];
				minIndex = j;
			}
		}

		O[r] = minIndex+1; //+1 is done to maintain the 1-32 range format and not 0-31 since this subtraction is being done at later steps
	}

	//train the model
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
	strcpy(fileName,filePath);
	strcat(fileName,"count.txt");
	file = fopen(fileName,"r");
	fscanf(file,"%d",&currCount);
	fclose(file);

	file = fopen(fileName,"w");
	fprintf(file,"%d\n",(currCount+1));
	fclose(file);

	strcpy(fileName,filePath);
	strcat(fileName,"A.txt");

	file = fopen(fileName,"r");

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			fscanf(file,"%lf",&AComplement[i][j]);

	fclose(file);

	file = fopen(fileName,"w");

	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < N; ++j)
		{
			tempValue = ((currCount*AComplement[i][j])+A[i][j])/(currCount+1);
			fprintf(file,"%g ",tempValue);
		}
		fprintf(file,"\n");
	}
	fclose(file);
			
	strcpy(fileName,filePath);
	strcat(fileName,"B.txt");
	
	file = fopen(fileName,"r");

	for (i = 0; i < N; ++i)
		for (j = 0; j < M; ++j)
			fscanf(file,"%lf",&BComplement[i][j]);

	fclose(file);

	file = fopen(fileName,"w");

	for (i = 0; i < N; ++i)
	{
		for (m = 0; m < M; ++m)
		{
			tempValue = ((currCount*BComplement[i][m])+B[i][m])/(currCount+1);
			fprintf(file,"%g ",tempValue);
		}
		fprintf(file,"\n");
	}
	fclose(file);
			
	strcpy(fileName,filePath);
	strcat(fileName,"Pi.txt");

	file = fopen(fileName,"r");

	for (i = 0; i < N; ++i)
		fscanf(file,"%lf",&PiComplement[i]);

	fclose(file);

	file = fopen(fileName,"w");

	for (j = 0; j < N; ++j)
	{
		tempValue = ((currCount*PiComplement[j])+Pi[j])/(currCount+1);
		fprintf(file,"%g ",tempValue);
	}
	fprintf(file,"\n");
	fclose(file);
}

void writeGeneralModel(char * filePath, int R)
{
	char fileName[256];
	FILE * file = NULL;
	int d = 0, k = 0, i = 0, j = 0;
	long double * row = new long double[N];
	long double tempB = 0.03125, one = 1, zero = 0;
	row[0] = 0.8;
	row[1] = 0.2;
	row[2] = row[3] = row[4] = 0;

	sprintf(fileName,"%s/A.txt",filePath);
	file = fopen(fileName,"w");

	for (i = 0; i < N-1; ++i)
	{
		for (j = 0; j < N; ++j)
			fprintf(file,"%g ",row[j]);
		fprintf(file,"\n");
		for (j = N-2; j >= 0; --j)
			row[j+1] = row[j];
		row[0] = 0;
	}
	for (i = 0; i < N-1; ++i)
		fprintf(file,"%g ",zero);
	fprintf(file,"%g \n",one);

	fclose(file);

	sprintf(fileName,"%s/B.txt",filePath);
	file = fopen(fileName,"w");

	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < M; ++j)
			fprintf(file,"%g ",tempB);
		fprintf(file,"\n");
	}

	fclose(file);

	sprintf(fileName,"%s/Pi.txt",filePath);
	file = fopen(fileName,"w");

	fprintf(file,"%g ",one);
	for (i = 1; i < N; ++i)
		fprintf(file,"%g ",zero);
	fprintf(file,"\n");

	fclose(file);

	sprintf(fileName,"%s/count.txt",filePath);
	file = fopen(fileName,"w");
	fprintf(file,"%d\n",R);
	fclose(file);
}
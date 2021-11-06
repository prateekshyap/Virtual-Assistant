extern int M, universeSize;
extern long double ** codebook;
int K = M;

void KMeans(char * fileName, int currK, long double delta, long double w[])
{
	int minIndex = 9999999; //minIndex stores the index of the code vector with minimum distortion
	/*previous distortion stores the distortion of iteration m-1
	current distortion stores the distortion of iteration m
	min stores the distortion value of the code vector with minimum distortion
	totalDistortion stores the total distortion of the universe*/
	long double previousDistortion = 0.0, currentDistortion = 1000.0, min = 1000000000.0, totalDistortion = 0.0, temp = 0;
	long long int iterationCount = 0; //total number of iterations in kmeans
	/*bucketNumber stores values from 0 to N-1 included i.e. the respective bucket numbers of all the vectors in the universe
	bucketSize stores the number of x vectors of each bucket*/
	int * bucketNumber = new int[universeSize], * bucketSize = new int[K]; //should lie between 0 to K-1 included
	long double * distances = new long double[K]; //stores the tokhura distances
	long double ** bucketSum = new long double*[K]; //stores the sum of buckets
	long double * currentUniverseVector = new long double[p]; // holds the current universe vector
	int i = 0, j = 0, k = 0, m = 0; //loop variables
	FILE * file = NULL;

	//initialize everything to 0
	for (i = 0; i < currK; ++i)
		bucketSum[i] = new long double[p];
	for (i = 0; i < currK; ++i)
		bucketSize[i] = 0;
	for (i = 0; i < currK; ++i)
		for (j = 0; j < p; ++j)
			bucketSum[i][j] = 0;
	
	i = 0;
	j = 0;
	while (abs(currentDistortion-previousDistortion) >= delta) //till the difference is more than the threshold
	{
	++iterationCount;
	totalDistortion = 0;
	for (i = 0; i < currK; ++i)
		bucketSize[i] = 0;
	for (i = 0; i < currK; ++i)
		for (j = 0; j < p; ++j)
			bucketSum[i][j] = 0;
	file = fopen(fileName,"r");
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*classification*/
	for (i = 0; i < universeSize; ++i) //for each vector in the universe
	{
		//read the current vector from the file
		for (j = 0; j < p; ++j) //for each value in current row of the file
		{
			fscanf(file,"%lf",&temp); //read the next value
			currentUniverseVector[j] = temp; //add it to the corresponding index of the current universe vector
		}

		min = 1000000000.0; //reset the minimum

		//find out the tokhura distances from each vector in the codebook
		for (j = 0; j < currK; ++j) //for each vector in the codebook
		{
			distances[j] = 0;
			for (k = 0; k < p; ++k) //for each p value
			{
				distances[j] += w[k]*(currentUniverseVector[k]-codebook[j][k])*(currentUniverseVector[k]-codebook[j][k]); //find out the distance
			}
		}

		//find out the minimum
		for (j = 0; j < currK; ++j)
		{
			if (min > distances[j])
			{
				min = distances[j];
				minIndex = j;
			}
		}

		//store the minIndex in the bucketNumber array
		bucketNumber[i] = minIndex;

		//increase the corresponding bucket size
		++bucketSize[minIndex];

		/*//update the maximum bucket size
		if (bucketSize[minIndex] > M)
			M = bucketSize[minIndex];*/

		//add the current vector to the bucket sum
		for (j = 0; j < p; ++j)
			bucketSum[minIndex][j] += currentUniverseVector[j];

		//add the current distortion to the distortion sum
		totalDistortion += min;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*code vector updation*/

	//taking the average and updating the codebook
	for (i = 0; i < currK; ++i)
		for (j = 0; j < p; ++j)
			codebook[i][j] = bucketSum[i][j]/bucketSize[i];

	previousDistortion = currentDistortion; //store the previous distortion
	currentDistortion = totalDistortion / universeSize; //find out the current distortion
	fclose(file);
	}
	delete(bucketNumber);
	delete(bucketSize);
	delete(distances);
	delete(bucketSum);
	delete(currentUniverseVector);
}

void LBG(char * fileName, long double * initialCodebook, long double delta, long double epsilon, long double w[])
{
	long double ** tempCodebook = NULL;
	codebook = NULL;
	int codebookSize = 1;
	int i = 0, j = 0, progressCount = 0, printGap = 0, percentageSize = 0;

	codebook = new long double*[1];
	for (i = 0; i < 1; ++i)
		codebook[i] = new long double[p];

	//store the initial codebook in codebook
	for (i = 0; i < 1; ++i)
		for (j = 0; j < p; ++j)
			codebook[i][j] = initialCodebook[j];
	
	//find out the progress printing details
	while(codebookSize != K)
	{
		codebookSize *= 2;
		++printGap;
	}
	percentageSize = 100/printGap;

	codebookSize = 1;
	
	printf("\n");
	while(codebookSize != K) //till the desired size is not reached
	{
		tempCodebook = new long double*[2*codebookSize]; //create a double sized codebook
		for (j = 0; j < 2*codebookSize; ++j)
			tempCodebook[j] = new long double[p];

		//split
		for (i = 0; i < codebookSize; ++i)
		{
			for (j = 0; j < p; ++j)
			{
				tempCodebook[2*i][j] = codebook[i][j]*(1+epsilon);
				tempCodebook[(2*i)+1][j] = codebook[i][j]*(1-epsilon);
			}
		}

		codebook = tempCodebook; //store the new codebook
		codebookSize *= 2; //update the new codebook size
		KMeans(fileName, codebookSize, delta, w); //run k means for optimal codebook

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
	//printing the progress bar if it is still incomplete
	printf("\r"); //move to the beginning of the line
	/*print the progress*/
	i = 0;
	printf("[");
	for (i = 1; i <= 100; ++i)
		printf("|");
	printf("]%d",100); //print the percentage
	//delete(tempCodebook);
}

long double * getInitialCodebook(char * fileName, int rowSize)
{
	FILE * file = NULL;
	char * buffer = new char[rowSize];
	int i = 0, j = 0;
	long double temp = 0;

	long double * initialCodebook = new long double[p];

	//store 0 in the entire codebook
	for (j = 0; j < p; ++j)
		initialCodebook[j] = 0;

	//add all the vectors from the universe file
	file = fopen(fileName,"r");

	for (i = 0; i < universeSize; ++i)
	{
		for (j = 0; j < p; ++j)
		{
			fscanf(file,"%lf",&temp); //read the next value
			initialCodebook[j] += temp; //add it to the corresponding index of the initial codebook
		}
	}

	fclose(file);

	//take the average
	for (j = 0; j < p; ++j)
		initialCodebook[j] /= universeSize;

	delete(buffer);
	return initialCodebook;
}

void buildCodebook()
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*variable declaration*/
	char * universeFileName = "data/universe.txt"; //universe file name
	/*universeSize is the total number of vectors given
	p stores the value of p
	rowSize stores the size of one row of the .csv file
	K is the codebook size
	*/
	int p = 12, rowSize = 1024;
	K = M;
	long double epsilon = 0.03, delta = 0.00001; // or 0.0001
	long double w[] = {1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0}; //given Tokhura weights

	long double * initialCodebook = NULL;
	int i = 0, j = 0; //loop variables

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*execution*/
	initialCodebook = getInitialCodebook(universeFileName, rowSize); //generate initial codebook	
	LBG(universeFileName, initialCodebook, delta, epsilon, w); //run LBG
	FILE * codebookFile = fopen("data/codebook.txt","w+");
	for (i = 0; i < M; ++i)
	{
		for (j = 0; j < p; ++j)
			fprintf(codebookFile,"%lf ",codebook[i][j]);
		fprintf(codebookFile,"\n");
	}
	fclose(codebookFile);
	//delete(codebook);
	//delete(universeFileName);
	//delete(initialCodebook);
}
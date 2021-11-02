extern long double *** xi;
extern long double ** A, ** B, ** alpha, ** beta, ** gamma;
extern long double ** AComplement, ** BComplement;
extern long double * PiComplement;
extern long double floorB;
extern int * O;
extern int N, M, T;
extern FILE * dataOutputFile, * modelOutputFile;

long double ** expectedTransitionsFromSiToSj = NULL, ** expectedObservationOfKthSymbolAtSj = NULL;

long double * pOfOGivenLambdaInDenom = NULL, * expectedTransitionsFromSi = NULL;

long double gammaDenom = 0;

void ensureStochastic()
{
	long double rowMax = 0, rowSum = 0, normalizer = 0;
	int maxIndex = 0;
	int i = 0, j = 0;
	for (i = 0; i < N; ++i)
	{
		rowMax = rowSum = 0;
		maxIndex = -1;
		for (j = 0; j < N; ++j)
		{
			if (AComplement[i][j] > rowMax)
			{
				rowMax = AComplement[i][j];
				maxIndex = j;
			}
			rowSum += AComplement[i][j];
		}
		if (rowSum != 1)
		{
			normalizer = abs(rowSum-1);
			AComplement[i][maxIndex] = rowSum > 1 ? AComplement[i][maxIndex]-normalizer : AComplement[i][maxIndex]+normalizer;
		}
	}

	for (i = 0; i < N; ++i)
	{
		rowMax = rowSum = 0;
		maxIndex = -1;
		for (j = 0; j < M; ++j)
		{
			if (BComplement[i][j] > rowMax)
			{
				rowMax = BComplement[i][j];
				maxIndex = j;
			}
			rowSum += BComplement[i][j];
		}
		if (rowSum != 1)
		{
			normalizer = abs(rowSum-1);
			BComplement[i][maxIndex] = rowSum > 1 ? BComplement[i][maxIndex]-normalizer : BComplement[i][maxIndex]+normalizer;
		}
	}
}

void runBaumWelch()
{
	int i = 0, j = 0, k = 0, t = 0, m = 0;

	expectedTransitionsFromSi = new long double[N];
	
	expectedTransitionsFromSiToSj = new long double * [N];
	for (i = 0; i < N; ++i)
		expectedTransitionsFromSiToSj[i] = new long double[N];

	expectedObservationOfKthSymbolAtSj = new long double * [N];
	for (i = 0; i < N; ++i)
		expectedObservationOfKthSymbolAtSj[i] = new long double[M];

	pOfOGivenLambdaInDenom = new long double[T-1];
	for (i = 0; i < T; ++i)
		pOfOGivenLambdaInDenom[i] = 0;

	//find out p of O given lambda for each t
	for (t = T-2; t >= 0; --t)
		for (i = 0; i < N; ++i)
			for (j = 0; j < N; ++j)
				pOfOGivenLambdaInDenom[t] += (alpha[i][t]*A[i][j]*B[j][O[t+1]-1]*beta[j][t+1]);

	//find out xi
	for (t = T-2; t >= 0; --t)
		for (i = 0; i < N; ++i)
			for (j = 0; j < N; ++j)
				xi[i][j][t] = (alpha[i][t]*A[i][j]*B[j][O[t+1]-1]*beta[j][t+1])/pOfOGivenLambdaInDenom[t];


	//find out gamma
	/*
	//using xi
	for (i = 0; i < N; ++i)
	{
		for (t = 0; t < T; ++t)
		{
			gamma[i][t] = 0;
			for (j = 0; j < N; ++j)
				gamma[i][t] += xi[i][j][t];
		}
	}
	*/
	
	//using alpha and beta
	for (t = 0; t < T; ++t)
	{
		gammaDenom = 0;
		for (j = 0; j < N; ++j)
			gammaDenom += alpha[j][t]*beta[j][t];
		for (i = 0; i < N; ++i)
			gamma[i][t] = (alpha[i][t]*beta[i][t])/gammaDenom;
	}
	

	//find out the expectations
	for (i = 0; i < N; ++i)
	{
		expectedTransitionsFromSi[i] = 0;
		for (t = 0; t < T; ++t)
			expectedTransitionsFromSi[i] += gamma[i][t];
	}

	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < N; ++j)
		{
			expectedTransitionsFromSiToSj[i][j] = 0;
			for (t = 0; t < T-1; ++t)
				expectedTransitionsFromSiToSj[i][j] += xi[i][j][t];
		}
	}

	for (j = 0; j < N; ++j)
	{
		for (k = 1; k <= M; ++k)
		{
			expectedObservationOfKthSymbolAtSj[j][k-1] = 0;
			for (t = 0; t < T; ++t)
			{
				if (O[t] == k)
					expectedObservationOfKthSymbolAtSj[j][k-1] += gamma[j][t];
			}
		}
	}

	//find out new Pi, A and B
	for (i = 0; i < N; ++i)
		PiComplement[i] = gamma[i][0];

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			AComplement[i][j] = expectedTransitionsFromSiToSj[i][j]/expectedTransitionsFromSi[i];

	for (j = 0; j < N; ++j)
		for (k = 0; k < M; ++k)
		{
			BComplement[j][k] = expectedObservationOfKthSymbolAtSj[j][k]/expectedTransitionsFromSi[j];
			if (BComplement[j][k] == 0)
				BComplement[j][k] = floorB;
		}

	ensureStochastic(); //for A and B both

	/*compare lambda and lambdaComplement here*/

}


/*
void ensureStochastic()
{
	long double rowMax = 0, rowSum = 0, normalizer = 0;
	int maxIndex = 0;
	int i = 0, j = 0;
	for (i = 0; i < N; ++i)
	{
		rowMax = rowSum = 0;
		maxIndex = -1;
		for (j = 0; j < N; ++j)
		{
			if (AComplement[i][j] > rowMax)
			{
				rowMax = AComplement[i][j];
				maxIndex = j;
			}
			rowSum += AComplement[i][j];
		}
		if (rowSum != 1)
		{
			normalizer = abs(rowSum-1);
			AComplement[i][maxIndex] = rowSum > 1 ? AComplement[i][maxIndex]-normalizer : AComplement[i][maxIndex]+normalizer;
		}
	}

	for (i = 0; i < N; ++i)
	{
		rowMax = rowSum = 0;
		maxIndex = -1;
		for (j = 0; j < M; ++j)
		{
			if (BComplement[i][j] > rowMax)
			{
				rowMax = BComplement[i][j];
				maxIndex = j;
			}
			rowSum += BComplement[i][j];
		}
		if (rowSum != 1)
		{
			normalizer = abs(rowSum-1);
			BComplement[i][maxIndex] = rowSum > 1 ? BComplement[i][maxIndex]-normalizer : BComplement[i][maxIndex]+normalizer;
		}
	}
}

void runBaumWelch()
{
	int i = 0, j = 0, k = 0, t = 0, m = 0;

	expectedTransitionsFromSi = new long double[N];
	
	expectedTransitionsFromSiToSj = new long double * [N];
	for (i = 0; i < N; ++i)
		expectedTransitionsFromSiToSj[i] = new long double[N];

	expectedObservationOfKthSymbolAtSj = new long double * [N];
	for (i = 0; i < N; ++i)
		expectedObservationOfKthSymbolAtSj[i] = new long double[M];

	pOfOGivenLambdaInDenom = new long double[T-1];
	for (i = 0; i < T; ++i)
		pOfOGivenLambdaInDenom[i] = 0;

	/*xi = new long double ** [N];
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

	//find out p of O given lambda for each t
	for (t = T-2; t >= 0; --t)
		for (i = 0; i < N; ++i)
			for (j = 0; j < N; ++j)
				pOfOGivenLambdaInDenom[t] += (alpha[i][t]*A[i][j]*B[j][O[t+1]-1]*beta[j][t+1]);

	//find out xi
	for (t = T-2; t >= 0; --t)
		for (i = 0; i < N; ++i)
			for (j = 0; j < N; ++j)
				xi[i][j][t] = (alpha[i][t]*A[i][j]*B[j][O[t+1]-1]*beta[j][t+1])/pOfOGivenLambdaInDenom[t];


	//find out gamma
	/*
	//using xi
	for (i = 0; i < N; ++i)
	{
		for (t = 0; t < T; ++t)
		{
			gamma[i][t] = 0;
			for (j = 0; j < N; ++j)
				gamma[i][t] += xi[i][j][t];
		}
	}
	
	
	//using alpha and beta
	for (t = 0; t < T; ++t)
	{
		gammaDenom = 0;
		for (j = 0; j < N; ++j)
			gammaDenom += alpha[j][t]*beta[j][t];
		for (i = 0; i < N; ++i)
			gamma[i][t] = (alpha[i][t]*beta[i][t])/gammaDenom;
	}
	

	//find out the expectations
	for (i = 0; i < N; ++i)
	{
		expectedTransitionsFromSi[i] = 0;
		for (t = 0; t < T; ++t)
			expectedTransitionsFromSi[i] += gamma[i][t];
	}

	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < N; ++j)
		{
			expectedTransitionsFromSiToSj[i][j] = 0;
			for (t = 0; t < T-1; ++t)
				expectedTransitionsFromSiToSj[i][j] += xi[i][j][t];
		}
	}

	for (j = 0; j < N; ++j)
	{
		for (k = 1; k <= M; ++k)
		{
			expectedObservationOfKthSymbolAtSj[j][k-1] = 0;
			for (t = 0; t < T; ++t)
			{
				if (O[t] == k)
					expectedObservationOfKthSymbolAtSj[j][k-1] += gamma[j][t];
			}
		}
	}

	//find out new Pi, A and B
	for (i = 0; i < N; ++i)
		PiComplement[i] = gamma[i][0];

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			AComplement[i][j] = expectedTransitionsFromSiToSj[i][j]/expectedTransitionsFromSi[i];

	for (j = 0; j < N; ++j)
		for (k = 0; k < M; ++k)
		{
			BComplement[j][k] = expectedObservationOfKthSymbolAtSj[j][k]/expectedTransitionsFromSi[j];
			if (BComplement[j][k] == 0)
				BComplement[j][k] = floorB;
		}

	ensureStochastic(); //for A and B both

	/*compare lambda and lambdaComplement here

	//print into file
	fprintf(modelOutputFile,"\n\nxi matrix-\n");
	for (i = 0; i < N; ++i)
	{
		fprintf(modelOutputFile,"\ni = %d\n",i);
		for (t = 0; t < T-1; ++t)
		{
			for (j = 0; j < N; ++j)
				fprintf(modelOutputFile,"%g ",xi[i][j][t]);
			fprintf(modelOutputFile,"\n");
		}
	}
	fprintf(modelOutputFile,"\n\ngamma matrix-\n");
	for (t = 0; t < T; ++t)
	{
		for (i = 0; i < N; ++i)
			fprintf(modelOutputFile,"%g ",gamma[i][t]);
		fprintf(modelOutputFile,"\n");
	}
	fprintf(modelOutputFile,"\n\nNew model generated-\n");
	fprintf(modelOutputFile,"\nA matrix-\n");
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < N; ++j)
			fprintf(modelOutputFile,"%g\t\t",AComplement[i][j]);
		fprintf(modelOutputFile,"\n");
	}
	fprintf(modelOutputFile,"\nB matrix-\n");
	for (m = 0; m < M; ++m)
	{
		for (i = 0; i < N; ++i)
			fprintf(modelOutputFile,"%g\t",BComplement[i][m]);
		fprintf(modelOutputFile,"\n");
	}
	fprintf(modelOutputFile,"\nPi array-\n");
	for (i = 0; i < N; ++i)
	{
		fprintf(modelOutputFile,"%g\t",PiComplement[i]);
	}
	fprintf(modelOutputFile,"\n");
}



*/
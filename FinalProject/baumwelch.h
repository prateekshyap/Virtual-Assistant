extern long double *** xi;
extern long double ** gamma;
extern long double ** AComplement, ** BComplement;

extern long double * PiComplement;

extern int N, M, T;

long double ** expectedTransitionsFromSiToSj = NULL, ** expectedObservationOfKthSymbolAtSj = NULL;

long double * pOfOGivenLambda = NULL, * expectedTransitionsFromSi = NULL;

void runBaumWelch()
{
	int i = 0, j = 0, k = 0, t = 0;

	expectedTransitionsFromSi = new long double[N];
	
	expectedTransitionsFromSiToSj = new long double * [N];
	for (i = 0; i < N; ++i)
		expectedTransitionsFromSiToSj[i] = new long double[N];

	expectedObservationOfKthSymbolAtSj = new long double * [N];
	for (i = 0; i < N; ++i)
		expectedObservationOfKthSymbolAtSj[i] = new long double[T];

	pOfOGivenLambda = new long double[T-1];
	for (i = 0; i < T; ++i)
		pOfOGivenLambda[i] = 0;

	xi = new long double ** [N];
	for (i = 0; i < N; ++i)
		xi[i] = new long double * [N];
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			xi[i][j] = new long double[T-1];

	gamma = new long double * [N];
	for (i = 0; i < N; ++i)
		gamma[i] = new long double[T-1];

	PiComplement = new long double [N];

	AComplement = new long double * [N];
	for (i = 0; i < N; ++i)
		AComplement[i] = new long double[N];

	BComplement = new long double * [N];
	for (i = 0; i < N; ++i)
		BComplement = new long double[T];


	//find out p of O given lambda for each t
	for (t = T-2; t >= 0; --t)
		for (j = 0; j < N; ++j)
			for (i = 0; i < N; ++i)
				pOfOGivenLambda[t] += (alpha[i][t]*A[i][j]*B[j][t+1]*beta[j][t+1]);

	//find out xi
	for (t = T-2; t >= 0; --t)
		for (j = 0; j < N; ++j)
			for (i = 0; i < N; ++i)
				xi[i][j][t] = (alpha[i][t]*A[i][j]*B[j][O[t+1]-1]*beta[j][t+1])/pOfOGivenLambda[t];


	//find out gamma
	for (i = 0; i < N; ++i)
	{
		for (t = 0; t < T-1; ++t)
		{
			gamma[i][t] = 0;
			for (j = 0; j < N; ++j)
				gamma[i][t] += xi[i][j][t];
		}
	}

	//find out the expectations
	for (i = 0; i < N; ++i)
	{
		expectedTransitionsFromSi[i] = 0;
		for (t = 0; t < T-2; ++t)
		{
			expectedTransitionsFromSi[i] += gamma[i][t];
		}
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

	/*Not sure about this loop, will update later*/
	for (j = 0; j < N; ++j)
	{
		for (k = 1; k <= T; ++k)
		{
			expectedObservationOfKthSymbolAtSj[j][k-1] = 0;
			for (t = 0; t < T; ++t)
				if (O[t] == k)
					expectedObservationOfKthSymbolAtSj[j][k-1] += gamma[j][t];
		}
	}

	//find out new Pi, A and B
	for (i = 0; i < N; ++i)
		PiComplement[i] = gamma[0][i];

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			AComplement[i][j] = expectedTransitionsFromSiToSj[i][j]/expectedTransitionsFromSi[i];

	for (j = 0; j < N; ++j)
		for (k = 0; k < T; ++k)
			BComplement[j][k] = expectedObservationOfKthSymbolAtSj[j][k]/expectedTransitionsFromSi[j];

	/*compare lambda and lambdaComplement here*/
}
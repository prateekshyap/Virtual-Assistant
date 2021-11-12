extern long double ** A, ** B, ** alpha, ** beta;
extern long double * Pi;
extern long double pOfOGivenLambda;
extern int * O;
extern int N, M, T;
extern FILE * dataOutputFile, * modelOutputFile;

void runForwardBackward()
{
	int i = 0, j = 0, t = 0;

	//initialization
	for (i = 0; i < N; ++i)
		alpha[i][0] = Pi[i]*B[i][O[0]-1];

	//induction
	for (t = 0; t < T-1; ++t)
	{
		for (j = 0; j < N; ++j)
		{
			alpha[j][t+1] = 0;
			for (i = 0; i < N; ++i)
				alpha[j][t+1] += alpha[i][t]*A[i][j];
			alpha[j][t+1] *= B[j][O[t+1]-1];
		}
	}

	//finding out probability
	pOfOGivenLambda = 0;
	for (i = 0; i < N; ++i)
		pOfOGivenLambda += alpha[i][T-1];

	//initialization
	for (i = 0; i < N; ++i)
		beta[i][T-1] = 1;

	//induction
	for (t = T-2; t >= 0; --t)
	{
		for (i = 0; i < N; ++i)
		{
			beta[i][t] = 0;
			for (j = 0; j < N; ++j)
				beta[i][t] += A[i][j]*B[j][O[t+1]-1]*beta[j][t+1];
		}
	}
}
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



/*


int i = 0, j = 0, t = 0;
	alpha = new long double *[N];
	beta = new long double *[N];

	for (i = 0; i < N; ++i)
		alpha[i] = new long double[T];
	for (i = 0; i < N; ++i)
		beta[i] = new long double[T];

	fprintf(dataOutputFile,"\nObservation Sequence:\n");
	for (i = 0; i < T; ++i)
		fprintf(dataOutputFile,"%d ",O[i]);
	fprintf(dataOutputFile,"\n");

	fprintf(modelOutputFile,"\nObservation Sequence:\n");
	for (i = 0; i < T; ++i)
		fprintf(modelOutputFile,"%d ",O[i]);
	fprintf(modelOutputFile,"\n");

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
			{
				alpha[j][t+1] += alpha[i][t]*A[i][j];
			}
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
			{
				beta[i][t] += A[i][j]*B[j][O[t+1]-1]*beta[j][t+1];
			}
		}
	}

	//print into file
	fprintf(modelOutputFile,"\n\nalpha matrix-\n");
	for (t = 0; t < T; ++t)
	{
		for (i = 0; i < N; ++i)
			fprintf(modelOutputFile,"%g\t",alpha[i][t]);
		fprintf(modelOutputFile,"\n");
	}
	fprintf(modelOutputFile,"\n\n\nbeta matrix-\n");
	for (t = 0; t < T; ++t)
	{
		for (i = 0; i < N; ++i)
			fprintf(modelOutputFile,"%g\t",beta[i][t]);
		fprintf(modelOutputFile,"\n");
	}
	fprintf(dataOutputFile,"Probability of Observation Sequence Given the model Lambda-\n");
	fprintf(dataOutputFile,"%g\n",pOfOGivenLambda);

*/
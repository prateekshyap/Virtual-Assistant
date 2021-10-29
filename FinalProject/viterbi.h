int N = 5, M = 32, D = 3, F = 4, T = 85;

long double ** A = NULL, ** B = NULL, ** delta = NULL;
long double * Pi = NULL;
long double pStar = 0, temp = 0;
int ** O = NULL, ** psi = NULL;
int * qStar = NULL;
FILE * outputFile = NULL;

void runViterbi(int row)
{
	int i = 0, j = 0, t = 0;
	delta = new long double *[N];
	psi = new int *[N];
	pStar = temp = 0;
	qStar = new int[T];

	for (i = 0; i < N; ++i)
		delta[i] = new long double[T];
	for (i = 0; i < N; ++i)
		psi[i] = new int[T];

	fprintf(outputFile,"\nObservation Sequence:\n");
	for (i = 0; i < T; ++i)
		fprintf(outputFile,"%d ",O[row][i]);

	//initialization
	for (i = 0; i < N; ++i)
	{
		delta[i][0] = Pi[i]*B[i][O[row][0]-1];
		psi[i][0] = 0;
	}

	//recursion
	for (t = 1; t < T; ++t)
	{
		for (j = 0; j < N; ++j)
		{
			delta[j][t] = 0;
			psi[j][t] = -1;
			for (i = 0; i < N; ++i)
			{
				temp = delta[i][t-1]*A[i][j];
				if (temp > delta[j][t])
				{
					delta[j][t] = temp;
					psi[j][t] = i;
				}
			}
			delta[j][t] *= B[j][O[row][t]-1];
		}
	}

	/*
	//remove this comment to print delta
	for (i = 0; i < T; ++i)
	{
		for (j = 0; j < N; ++j)
			printf("%g ",delta[j][i]);
		printf("\n");
	}
	*/

	/*
	//remove this comment to print psi
	for (i = 0; i < T; ++i)
	{
		for (j = 0; j < N; ++j)
			printf("%d ",psi[j][i]);
		printf("\n");
	}
	*/

	//termination
	qStar[T-1] = -1;
	for (i = 0; i < N; ++i)
	{
		if (delta[i][T-1] > pStar)
		{
			pStar = delta[i][T-1];
			qStar[T-1] = i;
		}
	}

	//backtracking
	for (t = T-2; t >= 0; --t)
		qStar[t] = psi[qStar[t+1]][t+1];

	/*
	//modify this part after getting values till this point
	fprintf(outputFile,"\nValue of P_*:%g\n",pStar);
	fprintf(outputFile,"q_*:\n");
	for (t = 0; t < T; ++t)
		fprintf(outputFile,"%d ",(qStar[t]+1));
	fprintf(outputFile,"\n");
	*/
}
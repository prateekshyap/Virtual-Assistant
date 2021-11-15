extern long double ** A, ** B, ** AComplement, **BComplement, ** delta;
extern long double * Pi, * PiComplement;
extern long double pStar, pStarComplement;
extern int * O, ** psi;
extern int * qStar, * qStarComplement;
extern int N, M, T;
extern FILE * dataOutputFile, * modelOutputFile;

long double temp = 0;

void runViterbi(int status)
{
	int i = 0, j = 0, t = 0;
	
	if (status == 0)
	{
		//initialization
		for (i = 0; i < N; ++i)
		{
			delta[i][0] = Pi[i]*B[i][O[0]-1];
			psi[i][0] = 0;
		}

		//recursion
		for (t = 1; t < T; ++t)
		{
			for (j = 0; j < N; ++j)
			{
				delta[j][t] = (long double)INT_MIN; //suspicious
				psi[j][t] = -1;
				for (i = 0; i < N; ++i)
				{
					temp = delta[i][t-1]*A[i][j];
					if (temp >= delta[j][t]) //suspicious
					{
						delta[j][t] = temp;
						psi[j][t] = i;
					}
				}
				delta[j][t] *= B[j][O[t]-1];
			}
		}

		//termination
		qStar[T-1] = -1;
		pStar = (long double)INT_MIN;
		for (i = 0; i < N; ++i)
		{
			if (delta[i][T-1] >= pStar) //suspicious
			{
				pStar = delta[i][T-1];
				qStar[T-1] = i;
			}
		}

		//backtracking
		for (t = T-2; t >= 0; --t)
			qStar[t] = psi[qStar[t+1]][t+1];

		//printf("pStar = %g\n",pStar);
	}
	else if (status == 1)
	{
		//initialization
		for (i = 0; i < N; ++i)
		{
			delta[i][0] = PiComplement[i]*BComplement[i][O[0]-1];
			psi[i][0] = 0;
		}

		//recursion
		for (t = 1; t < T; ++t)
		{
			for (j = 0; j < N; ++j)
			{
				delta[j][t] = (long double)INT_MIN; //suspicious
				psi[j][t] = -1;
				for (i = 0; i < N; ++i)
				{
					temp = delta[i][t-1]*AComplement[i][j];
					if (temp >= delta[j][t]) //suspicious
					{
						delta[j][t] = temp;
						psi[j][t] = i;
					}
				}
				delta[j][t] *= BComplement[j][O[t]-1];
			}
		}

		//termination
		qStarComplement[T-1] = -1;
		pStarComplement = (long double)INT_MIN;
		for (i = 0; i < N; ++i)
		{
			if (delta[i][T-1] >= pStarComplement) //suspicious
			{
				pStarComplement = delta[i][T-1];
				qStarComplement[T-1] = i;
			}
		}

		//backtracking
		for (t = T-2; t >= 0; --t)
		{
			//printf("%d ",qStarComplement[t+1]);
			qStarComplement[t] = psi[qStarComplement[t+1]][t+1];
		}
	}
}
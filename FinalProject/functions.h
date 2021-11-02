extern FILE * AComplementFile, * BComplementFile, * PiComplementFile;
extern long double pStar, pStarComplement;
extern long double ** AComplement, ** BComplement;
extern long double * PiComplement;
extern int N, M, T;

bool compareAndUpdateModel()
{
	int i = 0, j = 0, t = 0, m = 0;
	if (pStarComplement > pStar)// && pStarComplement <= 1e-30)
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

		return true;
	}
	else return false;
}
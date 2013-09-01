/* BLAS: dgemm.f */

__declspec(dllexport) int dgemm(char TRANSA, char TRANSB, int M, int N, int K, double ALPHA, double * A, int LDA, double * B, int LDB, double BETA, double * C, int LDC)
{
	/* C = ALPHA*op(A)*op(B) + BETA*C
	 * 	where op(X) is one of
	 *		op(X) = X or op(X) = X**T.
	 *
	 *	op(A): M by K
	 *	op(B): K by N
	 *	    C: M by N
	 *
	 *	TRANSA = 'N' or 'n' --> op(A) = A,
	 *	TRANSA = 'T' or 't' --> op(A) = A**T.
	 *
	 *	TRANSB = 'N' or 'n' --> op(B) = B,
	 *	TRANSB = 'T' or 't' --> op(B) = B**T.
	 *
	 *	LDA is the leading dimension of A. It's the number of columns of A here.
	 *	LDB is the leading dimension of B. It's the number of columns of B herel
	 *	LDC is the leading dimension of C. It's the number of columns of C herel
	 */
	 
	int i, j, k, index;
	short NOTA, NOTB;	/* NO Trans A/B */
	int NROWA, NCOLA;
	int NROWB;
	double temp;
	 
	if( TRANSA == 'T' || TRANSA == 't' )
		NOTA = 0;
	else if( TRANSA == 'N' || TRANSA == 'n' )
		NOTA = 1;
	else
		return -1;
	
	if( TRANSB == 'T' || TRANSB == 't' )
		NOTB = 0;
	else if( TRANSB == 'N' || TRANSB == 'n' )
		NOTB = 1;
	else
		return -2;
	 
	if( NOTA ) { NROWA = M; NCOLA = K; }
		else { NROWA = K; NCOLA = M; }
		
	if( NOTB ) { NROWB = K; }
		else { NROWB = N; }
	
	/* Quick return without any calculation. */
	if( M == 0 || N == 0 || ( (ALPHA == 0.0 || K == 0) && (BETA == 1.0) ) )
		return 0;
	
	/* If ALPHA==0.0, nothing to do with A or B. */
	if( ALPHA == 0.0 )
	{
		if( BETA == 0.0 )	/* C(M)(N) is set to be 0.0 */
		{
			for( i = 0; i < M; i++ )
			{
				index = i * LDC;
				for( j = 0; j < N; j++ )
					C[index++] = 0.0;
			}
			return 0;
		}
		else				/* BETA == 1.0: already returned. */
		{
			for( i = 0; i < M; i++ )
			{
				index = i * LDC;
				for( j = 0; j < N; j++ )
					C[index++] *= BETA;
			}
			return 0;
		}
	}
	
	/* Start the operations. */
	if( NOTB )
	{
		/* B: K by N */
		if( NOTA )
		{
			/* A: M by K */
			for( i = 0; i < M; i++ )
			{
				if( BETA == 0.0 )
				{
					index = i * LDC;
					for( j = 0; j < N; j++ )
						C[index++] = 0.0;
				}
				else if( BETA != 1.0 )
				{
					index = i * LDC;
					for( j = 0; j < N; j++ )
						C[index++] *= BETA;
				}
				
				for( k = 0; k < K; k++ )
				{
					if( A[i*LDA+k] != 0 )
					{
						temp = ALPHA * A[i*LDA+k];
						index = i * LDC;
						for( j = 0; j < N; j++ )
							C[index++] += ( temp * B[k*LDB+j] );
					}
				}
			}
		}
		else
		{
			/* A: K by M */
			for( i = 0; i < M; i++ )
			{
				for( j = 0; j < N; j++ )
				{
					temp = 0.0;
					for( k = 0; k < K; k++ )
						temp += ( A[k*LDA+i] * B[k*LDB+j] );
						
					if( BETA == 0.0 )
						C[i*LDC+j] = ( ALPHA * temp );
					else
						C[i*LDC+j] = ALPHA * temp + BETA * C[i*LDC+j];
				}
			}
		}
	}
	else
	{
		/* B: N by K */
		if( NOTA )
		{
			/* A: M by K */
			for( i = 0; i < M; i++ )
			{
				if( BETA == 0.0 )
				{
					index = i * LDC;
					for( j = 0; j < N; j++ )
					{
						C[index++] = 0.0;
					}
				}
				else if( BETA != 1.0 )
				{
					index = i * LDC;
					for( j = 0; j < N; j++ )
					{
						C[index++] *= BETA;
					}
				}
				
				for( k = 0; k < K; k++ )
				{
					if( A[i*LDA+k] != 0.0 )
					{
						temp = ALPHA * A[i*LDA+k];
						index = i * LDC;
						for( j = 0; j < N; j++ )
						{
							C[index++] += ( temp * B[j*LDB+k] );
						}
					}
				}
			}
		}
		else
		{
			/* A: K by M */
			for( i = 0; i < M; i++ )
			{
				for( j = 0; j < N; j++ )
				{
					temp = 0.0;
					for( k = 0; k < K; k++ )
						temp += ( A[k*LDA+i] * B[j*LDB+k] );
						
					if( BETA == 0.0 )
						C[i*LDC+j] = ( ALPHA * temp );
					else
						C[i*LDC+j] = ALPHA * temp + BETA * C[i*LDC+j];
				}
			}
		}
	}
	
	return 0;
}

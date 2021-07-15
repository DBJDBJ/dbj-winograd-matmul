#ifndef strassen_winograd_naive_h
#define strassen_winograd_naive_h

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/*
* multiply two square matrices, by the book
*/
static inline void naive_mult(
	const unsigned N,
	double m1[static N][N], double m2[static N][N], double res[static N][N]
)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			double sum = 0.;
			for (int k = 0; k < N; k++)
			{
				sum += m1[i][k] * m2[k][j];
			}
			res[i][j] = sum;
		}
	}
}

#endif // strassen_winograd_naive_h
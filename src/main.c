//
//  Copyright (c) 2021 by dbj at dbj dot org. All rights reserved. https://dbj.org/license_dbj
//

#include "d:\MACHINE_WIDE\ubut\ubench.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef double matrix_data_type;

#include "naive.h"

#define winograd_api_implementation 1
#include "winograd.h"

// we take all matrices as square
// that simplifies things considerably
#define DBJ_MATRIX_SIDE 512

static struct app_data_struct
{
	unsigned int N;
	// each matrix is N x N
	matrix_data_type* m1;
	matrix_data_type* m2;
	matrix_data_type* naive_res;
	matrix_data_type* winograd_res;
	//matrix strassen_res;

	double* winograd_row;
	double* winograd_col;

	size_t winograd_row_size;
	size_t winograd_col_size;

} app_data = { .N = DBJ_MATRIX_SIDE };

static inline matrix_data_type* alloc_matrix(unsigned rows, unsigned cols) {
	matrix_data_type* rez = calloc(rows * cols, sizeof(matrix_data_type));
	assert(rez);
	return rez;
}

static inline void random_matrix(unsigned int count, matrix_data_type m[static count])
{
	for (unsigned j = 0; j < count; ++j)
	{
		m[j] = (matrix_data_type)(rand() / RAND_MAX);
	}
}

static void app_start(void)
{
	srand((unsigned int)time(NULL));

	app_data.m1 = alloc_matrix(app_data.N, app_data.N);
	app_data.m2 = alloc_matrix(app_data.N, app_data.N);
	app_data.naive_res = alloc_matrix(app_data.N, app_data.N);
	app_data.winograd_res = alloc_matrix(app_data.N, app_data.N);

	random_matrix(app_data.N, app_data.m1);
	random_matrix(app_data.N, app_data.m2);
	//zero_matrix(app_data.N, app_data.N, &app_data.naive_res);
	//zero_matrix(app_data.N, app_data.N, &app_data.winograd_res);

	app_data.winograd_row_size = app_data.N * sizeof(double);
	app_data.winograd_col_size = app_data.N * sizeof(double);

	app_data.winograd_row = calloc(1, app_data.N * sizeof(double));
	app_data.winograd_col = calloc(1, app_data.N * sizeof(double));

	// m1 and m2 are not changing between the calls to winograd_mult
	// so we do the preproc step only once and here
	winograd_preprocess(
		app_data.N, app_data.N, (void*)app_data.m1,
		app_data.N, app_data.N, (void*)app_data.m2,
		app_data.winograd_row_size, app_data.winograd_row,
		app_data.winograd_col_size, app_data.winograd_col);
}

static void app_end(void)
{
	free(app_data.m1);
	free(app_data.m2);
	free(app_data.naive_res);
	free(app_data.winograd_res);
	free(app_data.winograd_row);
	free(app_data.winograd_col);
}

UBENCH(dbj_winograd, mult_naive)
{
	naive_mult(app_data.N,
		(void*)app_data.m1,
		(void*)app_data.m2,
		(void*)app_data.naive_res);
}

UBENCH(dbj_winograd, winograd_mult)
{
	winograd_mult(
		app_data.N, app_data.N, (void*)app_data.m1,
		app_data.N, app_data.N, (void*)app_data.m2,
		app_data.N, app_data.N, (void*)app_data.winograd_res,
		app_data.winograd_row_size, app_data.winograd_row,
		app_data.winograd_col_size, app_data.winograd_col);
}

UBENCH_STATE();
int main(int argc, const char* const argv[])
{
	system(" ");
	app_start();
	return ubench_main(argc, argv);
	app_end();
}

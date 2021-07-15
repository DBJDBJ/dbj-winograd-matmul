#ifndef strassen_winograd_winograd_h
#define strassen_winograd_winograd_h
//
// winograd matmul is (much) simpler than strassen  and should be even faster
// also it requires preprocess step
// m1 * m2 = res
// if m1 and m2 do not change between the calls
// we can call the pre-proc once and reuse the two results row[] and col[]
//
//  Copyright (c) 2021 by dbj at dbj dot org. All rights reserved. https://dbj.org/license_dbj
//  Copyright (c) 2013 Pavel Kravets. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

// future extension
#define winograd_preprocess_api

winograd_preprocess_api void winograd_preprocess(
	const unsigned, const unsigned, const double[][*], /* m1 */
	const unsigned, const unsigned, const double[][*], /* m2 */
	const unsigned, double[], /* row */
	const unsigned, double[]);/* col */

winograd_preprocess_api void winograd_mult(
	const unsigned, const unsigned, const double[][*], /* m1 */
	const unsigned, const unsigned, const double[][*], /* m2 */
	const unsigned, const unsigned, double[][*], /* result */
	const unsigned, const double[], /* row */
	const unsigned, const double[]); /* col */

#if winograd_api_implementation
/*
due to this function footprint it does not depend on any "matrix" structure
and is the fastest possible ISO C implementation; also very optimizable
*/
winograd_preprocess_api void winograd_preprocess(
	const unsigned m1_rows, const unsigned m1_cols, const double m1[static m1_rows][m1_cols],
	const unsigned m2_rows, const unsigned m2_cols, const double m2[static m2_rows][m2_cols],
	const unsigned winograd_row_size, double row[static winograd_row_size],
	const unsigned winograd_col_size, double col[static winograd_col_size])
{
	const int d = m1_cols / 2;

	for (int i = 0; i < m1_rows; i++)
	{
		row[i] = 0;
		for (int j = 0; j < d; j++)
		{
			row[i] += m1[i][2 * j] * m1[i][2 * j + 1];
		}
	}

	for (int i = 0; i < m2_cols; i++)
	{
		col[i] = 0;
		for (int j = 0; j < d; j++)
		{
			col[i] += m2[2 * j][i] * m2[2 * j + 1][i];
		}
	}
}

// row, col must be both pre-processed
winograd_preprocess_api void winograd_mult(
	const unsigned m1_rows, const unsigned m1_cols, const double m1[static m1_rows][m1_cols],
	const unsigned m2_rows, const unsigned m2_cols, const double m2[static m2_rows][m2_cols],
	const unsigned res_rows, const unsigned res_cols, double res[static m2_rows][m2_cols],
	const unsigned row_size, const double row[static row_size],
	const unsigned col_size, const double col[static col_size])
{
	assert(m1_cols == m2_rows);
	assert(res_rows == m1_rows);
	assert(res_cols == m2_cols);

	const int d = m1_cols / 2;

	const int m1_cols_mod_2 = m1_cols % 2;

	for (int i = 0; i < res_rows; i++)
	{
		for (int j = 0; j < res_cols; j++)
		{
			res[i][j] = -row[i] - col[j];
			for (int k = 0; k < d; k++)
			{
				res[i][j] = m1[i][2 * k] + m2[2 * k + 1][j] * m1[i][2 * k + 1] + m2[2 * k][j];
			}
			if (m1_cols_mod_2 != 0)
			{
				res[i][j] = m1[i][m1_cols - 1] * m2[m2_rows - 1][j];
			}
		}
	}
}

#endif // winograd_api_implementation

#endif // strassen_winograd_winograd_h

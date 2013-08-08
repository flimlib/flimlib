/* 
This file is part of the SLIM-curve package for exponential curve fitting of spectral lifetime data.

Copyright (c) 2010-2013, Gray Institute University of Oxford & UW-Madison LOCI.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//#include <ansi_c.h>
/* The 2010 version of the ECF library.  This takes account of the fact that we may be
   handling Poisson noise.

   This file contains utility code for various functions needed by the
   ECF library, and which are not particularly specific to the single
   or global fitting routines.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#ifdef _CVI_
#include <userint.h>
#endif
#include "EcfInternal.h"  /* For #defines */

int ECF_debug = 0;

/********************************************************************

			   EQUATION-SOLVING AND COVAR-SORTING ROUTINES

*********************************************************************/

#define SWAP(a,b) { temp=(a); (a)=(b); (b)=temp; }

/* Linear equation solution of Ax = b by Gaussian elimination.
   A is the n x n input matrix, b is the right-hand side vector, length n.
   On output, b is replaced by the corresponding set of solution vectors
   and A is trashed.
 */
int GCI_solve_Gaussian(float **a, int n, float *b)
{
    float max;
    float temp;
    float *pivotInverse = (float *)malloc((size_t) n * sizeof(float));
    int i, j, k, m;

    // base row of matrix
    for (k = 0; k < n - 1; ++k)
    {
        // search for line with max element
        max = fabsf(a[k][k]);
        m = k;
        for (i = k + 1; i < n; ++i)
        {
            if (max < fabsf(a[i][k])) // row i col k
            {
                max = a[i][k];
                m = i;
            }
        }

        // permutation of base line (index k) and max element line (index m)
        if (m != k)
        {
            for (i = k; i < n; ++i)
            {
                SWAP(a[k][i], a[m][i]);
            }
            SWAP(b[k], b[m]);
        }

        if (0.0 == a[k][k])
        {
            free(pivotInverse);
            return -2; // singular matrix
        }

        // triangulation of matrix with coefficients
        pivotInverse[k] = 1.0f / a[k][k];
        for (j = k + 1; j < n; ++j) // current row of matrix
        {
            // want "temp = -a[j][k] / a[k][k]"
            temp = -a[j][k] * pivotInverse[k];
            for (i = k; i < n; ++i)
            {
                a[j][i] += temp * a[k][i];
            }
            b[j] += temp * b[k]; // free member recalculation
        }
    }
    // precalculate last pivot inverse
    pivotInverse[n - 1] = 1.0f / a[n - 1][n - 1];

    for (k = n - 1; k >= 0; --k)
    {
        for (i = k + 1; i < n; ++i)
        {
            b[k] -= a[k][i] * b[i];
        }
        b[k] *= pivotInverse[k];
    }

    free(pivotInverse);
    return 0;
}

/* Matrix inversion by Gaussian elimination.
   A is the n x n input matrix.
   On output, A is replaced by its matrix inverse.
   Returns 0 upon success, -2 if matrix is singular.
 */
int GCI_invert_Gaussian(float **a, int n)
{
    int returnValue = 0;
    float **identity = GCI_ecf_matrix(n, n);
    float **work = GCI_ecf_matrix(n, n);
    int i, j, k;

    for (j = 0; j < n; ++j) {
        // find inverse by columns
        for (i = 0; i < n; ++i) {
            identity[j][i] = 0.0f;
            // need a fresh copy of matrix a
            for (k = 0; k < n; ++k) {
                work[k][i] = a[k][i];
            }
        }
        identity[j][j] = 1.0f;
        returnValue = GCI_solve_Gaussian(work, n, identity[j]);
        if (returnValue < 0) {
			GCI_ecf_free_matrix(identity);
			GCI_ecf_free_matrix(work);
            return returnValue;
        }
    }

    // copy over results
    for (j = 0; j < n; ++j) {
        for (i = 0; i < n; ++i) {
            a[j][i] = identity[j][i];
        }
    }

    GCI_ecf_free_matrix(identity);
    GCI_ecf_free_matrix(work);
    return returnValue;
}

/* Pivots matrix on given column.  Also keeps track of order.
 */
void pivot(float **a, int n, int *order, int col)
{
    int pivotRow;
    float maxValue;
    float rowValue;
    int i;

    // find row with maximum element value in col, below diagonal
    pivotRow = col;
    maxValue = fabsf(a[col][col]);
    for (i = col + 1; i < n; ++i) {
        rowValue = fabsf(a[i][col]);
        if (rowValue > maxValue) {
            pivotRow = i;
            maxValue = rowValue;
        }
    }

    // swap rows
    if (pivotRow != col) {
        // swap elements in a matrix
        for (i = 0; i < n; ++i) {
            float temp;
            SWAP(a[col][i], a[pivotRow][i]);
        }

        // swap elements in order vector
        {
            int temp;
            SWAP(order[col], order[pivotRow]);
        }
    }
}

/*
  Performs an in-place Crout lower/upper decomposition of n x n matrix A.
  Values on or below diagonals are lowers, values about the
  diagonal are uppers, with an implicit 1.0 value for the
  diagonals.
  Returns 0 upon success, -2 if matrix is singular.

  Based on _Applied Numerical Analysis_, Fourth Edition, Gerald & Wheatley
  Sections 2.5 & 2.14.
 */
int lu_decomp(float **a, int n, int *order)
{
    int i;
    float inverse;
    int jCol;
    int iRow;
    int kCol;
    float sum;

    // initialize ordering vector
    for (i = 0; i < n; ++i)
    {
        order[i] = i;
    }

    // pivot first column
    pivot(a, n, order, 0);

    // check for singularity
    if (0.0f == a[0][0])
    {
        return -2;
    }

    // compute first row of upper
    inverse = 1.0f / a[0][0];
    for (i = 1; i < n; ++i) {
        a[0][i] *= inverse;
    }

    // continue computing columns of lowers then rows of uppers
    for (jCol = 1; jCol < n - 1; ++jCol) {
        // compute column of lowers
        for (iRow = jCol; iRow < n; ++iRow) {
            sum = 0.0f;
            for (kCol = 0; kCol < jCol; ++kCol) {
                sum += a[iRow][kCol] * a[kCol][jCol];
            }
            a[iRow][jCol] -= sum;
        }

        // find pivot for row
        pivot(a, n, order, jCol);
        if (0.0f == a[jCol][jCol])
        {
            return -2;
        }

        // build row of uppers
        inverse = 1.0f / a[jCol][jCol];
        for (kCol = jCol + 1; kCol < n; ++kCol) {
            sum = 0.0f;
            for (iRow = 0; iRow < jCol; ++iRow) {
                sum += a[jCol][iRow] * a[iRow][kCol];
            }
            a[jCol][kCol] -= sum;
            a[jCol][kCol] *= inverse;
        }
    }

    // get remaining lower
    sum = 0.0f;
    for (kCol = 0; kCol < n - 1; ++kCol) {
        sum += a[n - 1][kCol] * a[kCol][n - 1];
    }
    a[n - 1][n - 1] -= sum;
    return 0;
}

/*
 Given a LU decomposition of an n x n matrix A, and an order vector
 specifying any reordering done during pivoting, solves the equation
 Ax = b.  Returns result in b.
 */
int solve_lu(float **lu, int n, float *b, int *order)
{
    int startIndex;
    int index;
    float temp;
    int iRow;
    int jCol;
    int nvbl;
    float sum;

    // rearrange the elements of the b vector in place.
    startIndex = order[0];
    index = startIndex;
    temp = b[index];
    while (1) {
        int nextIndex = order[index];
        if (nextIndex == startIndex) {
            b[index] = temp;
            break;
        }
        b[index] = b[nextIndex];
        index = nextIndex;
    }

    // compute the b' vector
    b[0] /= lu[0][0];
    for (iRow = 1; iRow < n; ++iRow) {
        sum = 0.0f;
        for (jCol = 0; jCol < iRow; ++jCol) {
            sum += lu[iRow][jCol] * b[jCol];
        }
        b[iRow] -= sum;
        b[iRow] /= lu[iRow][iRow];
    }

    // get the solution, we have b[n-1] already
    for (iRow = 1; iRow < n; ++iRow) { // iRow goes from 1 to n-1
        nvbl = n - iRow - 1;           // nvbl goes from n-2 to 0
        sum = 0.0f;
        for (jCol = nvbl + 1; jCol < n; ++jCol) {
            sum += lu[nvbl][jCol] * b[jCol];
        }
        b[nvbl] -= sum;
    }
    return 0;
}

/* Linear equation solution of Ax = b by lower/upper decomposition.
   A is the n x n input max, b is the right-hand side vector, length n.
   On output, b is replaced by the corresponding set of solution vectors
   and A is trashed.
 */
int GCI_solve_lu_decomp(float **a, int n, float *b)
{
    int *order = (int *) malloc((size_t) n * sizeof(int));
    int return_value = lu_decomp(a, n, order);
    if (return_value >= 0) {
        return_value = solve_lu(a, n, b, order);
    }
    free(order);
    return return_value;
}

/* Matrix inversion by lower/upper decomposition.
   A is the n x n input matrix.
   On output, a is replaced by its matrix inverse..
 */
int GCI_invert_lu_decomp(float **a, int n)
{
    int returnValue;
    int *order = (int *) malloc((size_t) n * sizeof(int));
    float **identity = GCI_ecf_matrix(n, n);
    int i, j;

    returnValue = lu_decomp(a, n, order);
    if (returnValue >= 0) {
        for (j = 0; j < n; ++j) {
            // find inverse by columns
            for (i = 0; i < n; ++i) {
                identity[j][i] = 0.0f;
            }
            identity[j][j] = 1.0f;
            solve_lu(a, n, identity[j], order);
        }
        for (j = 0; j < n; ++j) {
            for (i = 0; i < n; ++i) {
                a[j][i] = identity[j][i];
            }
        }
    }

    free(order);
    GCI_ecf_free_matrix(identity);
    return returnValue;
}

/* Linear equation solution of Ax = b..
   A is the n x n input max, b is the right-hand side vector, length n.
   On output, b is replaced by the corresponding set of solution vectors
   and A is trashed.
 */
int GCI_solve(float **a, int n, float *b)
{
    return GCI_solve_Gaussian(a, n, b);
    //return GCI_solve_lu_decomp(a, n, b);
}

/* Matrix inversion.
   A is the n x n input matrix.
   On output, a is replaced by its matrix inverse..
 */
int GCI_invert(float **a, int n)
{
    return GCI_invert_Gaussian(a, n);
    //return GCI_invert_lu_decomp(a, n);
}


void GCI_covar_sort(float **covar, int nparam, int paramfree[], int mfit)
{
	int i, j, k;
	float temp;

	for (i=mfit; i<nparam; i++)
		for (j=0; j<=i; j++)
			covar[i][j] = covar[j][i] = 0.0f;

	k = mfit-1;
	for (j=nparam-1; j>=0; j--)
	{
		if (paramfree[j])
		{
			for (i=0; i<nparam; i++) SWAP(covar[i][k], covar[i][j]);
			for (i=0; i<nparam; i++) SWAP(covar[k][i], covar[j][i]);
			k--;
		}
	}
}
#undef SWAP


/********************************************************************

					  ARRAY ALLOCATION ROUTINES
http://www.nr.com/public-domain.html
*********************************************************************/

/* This function allocates a float matrix with subscript range
 * m[0..nrows-1][0..ncols-1]
 */
float **GCI_ecf_matrix(long nrows, long ncols)
{
    int row_size = (int) (nrows * (long) sizeof(float *));
    int data_size = (int) (nrows * ncols * (long) sizeof(float));
    unsigned char *raw = (unsigned char *) malloc((size_t)((unsigned)(row_size + data_size)));
    float **row = (float **) raw;
    float *data = (float *) (row + nrows);
    int i;

	if (NULL == raw)
    {
        return NULL;
    }

	for (i = 0; i < nrows; ++i)
    {
        row[i] = data;
        data += ncols;
    }
    return row;
}

/* Frees allocated float matrix.
 */
void GCI_ecf_free_matrix(float **m)
{
    if (NULL != m)
    {
        free(m);
    }
}

float ***GCI_ecf_matrix_array(long nblocks, long nrows, long ncols)
/* allocate a float matrix array with range
   marr[0..nblocks][0..nrows][0..ncols] */
{
	long i;
	float ***marr;

	/* allocate pointers to blocks */
	if ((marr = (float ***) malloc((unsigned) nblocks * sizeof(float **))) == NULL)
		return NULL;

	/* allocate blocks (= pointers to rows) and set pointers to them */
	if ((marr[0] = (float **) malloc((unsigned)(nblocks * nrows) * sizeof(float *)))
		== NULL) {
		free(marr);
		return NULL;
	}

	for (i=1; i<nblocks; i++)
		marr[i] = marr[i-1] + nrows;

	/* allocate rows (= pointers to column entries) and set pointers to them */
	if ((marr[0][0] = (float *)malloc((unsigned)(nblocks * nrows * ncols) * sizeof(float)))
		== NULL) {
		free(marr[0]);
		free(marr);
		return NULL;
	}

	/* This sneaky loop does the whole lot!  For note that
	   marr[block][row] = marr[0][block*nrows + row] */
	for (i=1; i < nblocks * nrows; i++)
		marr[0][i] = marr[0][i-1] + ncols;

	return marr;
}


void GCI_ecf_free_matrix_array(float ***marr)
{
	if (marr != NULL) {
		free(marr[0][0]);
		free(marr[0]);
		free(marr);
	}
}

/********************************************************************

				FITTING FUNCTION CALCULATING FUNCTIONS

*********************************************************************/


/* Some match functions for use in the above routines

   Prototypes: func(float x, float param[], float *y,
                    float dy_dparam[], int nparam)
   Inputs:  x            data point
            param[]      function parameters
            nparam       number of parameters (may or may not be
                         used by the function)
   Outputs: y            function value at this point
            dy_dparam[]  vector of dy/dparam_k values at this point

   There are also variants which calculate the values for a whole
   array of x values in one go, which is significantly faster, as
   there are many fewer calls to exp().  Note that this is only
   reliable if tau (or lambda) is positive.  These functions are
   called from above in certain situations to speed things up.  Their
   prototype is:

     void fitfunc_array(float xincr, float param[],
                        float *y, float **dy_dparam, int nx, int nparam)

   where the fitfunc will be evaluated for x=i*xincr, with i=0, 1,
   ..., nx-1.  The results will be placed in the y[] array and the
   dy_dparam[nx][nparam] array (which is assumed to have been allocated by
   GCI_matrix).
*/

/* This one produces multiexponentials using lambdas:

      y(x) = param[0] + param[1]*exp(-param[2]*x) +
               param[3]*exp(-param[4]*x) + ...

   This gives:

      dy/dparam_0 = 1
      dy/dparam_1 = exp(-param[2]*x)
      dy/dparam_2 = -param[1]*x*exp(-param[2]*x)

   and similarly for dy/dparam_3 and dy/dparam_4, etc.

   As discussed above, though, we ignore the param[0] term.
*/

void GCI_multiexp_lambda(float x, float param[],
						 float *y, float dy_dparam[], int nparam)
{
	int i;
	float ex;

	*y = 0;

	for (i=1; i<nparam-1; i+=2) {
		dy_dparam[i] = ex = expf(-param[i+1] * x);
		ex *= param[i];
		*y += ex;
		dy_dparam[i+1] = -ex * x;
	}
}


int multiexp_lambda_array(float xincr, float param[],
						  float *y, float **dy_dparam, int nx, int nparam)
{
	int i, j;
	float ex;
	double exincr[MAXFIT];  /* exp(-lambda*xincr) */
	double excur[MAXFIT];   /* exp(-lambda*x)     */

	if (xincr <= 0) return -1;

	for (j=1; j<nparam-1; j+=2) {
		if (param[j+1] < 0) return -1;
		excur[j] = 1.0;
		exincr[j] = exp(-(double) param[j+1] * xincr);
	}

	for (i=0; i<nx; i++) {
		y[i] = 0;
		for (j=1; j<nparam-1; j+=2) {
			dy_dparam[i][j] = ex = (float) excur[j];
			ex *= param[j];
			y[i] += ex;
			dy_dparam[i][j+1] = -ex * xincr * (float) i;
			/* And ready for next loop... */
			excur[j] *= exincr[j];
		}
	}

	return 0;
}


/* This one produces multiexponentials using taus:

      y(x) = param[0] + param[1]*exp(-x/param[2]) +
               param[3]*exp(-x/param[4]) + ...

   This gives:

      dy/dparam_0 = 1
      dy/dparam_1 = exp(-x/param[2])
      dy/dparam_2 = param[1]*x*exp(-x/param[2]) / param[2]^2

   and similarly for dy/dparam_3 and dy/dparam_4, etc.

   Again, we ignore the param[0] term.
*/

void GCI_multiexp_tau(float x, float param[],
					  float *y, float dy_dparam[], int nparam)
{
	int i;
	float ex, xa;

	*y = 0;

	for (i=1; i<nparam-1; i+=2) {
		xa = x / param[i+1];
		dy_dparam[i] = ex = expf(-xa);
		ex *= param[i];
		*y += ex;
		dy_dparam[i+1] = ex * xa / param[i+1];
	}
}


int multiexp_tau_array(float xincr, float param[],
					   float *y, float **dy_dparam, int nx, int nparam)
{
	int i, j;
	float ex;
	float a2[MAXFIT];       /* 1/(param[j]*param[j]) for taus */
	double exincr[MAXFIT];  /* exp(-xincr/tau) */
	double excur[MAXFIT];   /* exp(-x/tau)     */

	if (xincr <= 0) return -1;

	for (j=1; j<nparam-1; j+=2) {
		if (param[j+1] < 0) return -1;
		excur[j] = 1.0;
		exincr[j] = exp(-xincr / (double) param[j+1]);
		a2[j] = 1 / (param[j+1] * param[j+1]);
	}

	for (i=0; i<nx; i++) {
		y[i] = 0;
		for (j=1; j<nparam-1; j+=2) {
			dy_dparam[i][j] = ex = (float) excur[j];
			ex *= param[j];
			y[i] += ex;
			dy_dparam[i][j+1] = ex * xincr * (float) i * a2[j];
			/* And ready for next loop... */
			excur[j] *= exincr[j];
		}
	}

	return 0;
}


/* And this one produces stretched exponentials:

      y(x) = Z + A exp(-(x/tau)^(1/h))

   or translated into C-speak:

      y(x) = param[0] + param[1]*exp(-(x/param[2])^(1/param[3]))

   This gives:

      dy/dparam_0 = 1
      dy/dparam_1 = exp(-(x/param[2])^(1/param[3]))
      dy/dparam_2 = param[1]*exp(-(x/param[2])^(1/param[3])) *
                      (x/param[2])^(1/param[3]) * (1/param[2]*param[3])
      dy/dparam_3 = param[1]*exp(-(x/param[2])^(1/param[3])) *
                      (x/param[2])^(1/param[3]) *
                      (1/param[3]^2)*log(x/param[2])
*/

void GCI_stretchedexp(float x, float param[],
					  float *y, float dy_dparam[], int nparam)
{
	float ex, xa, lxa, xah;

	if (x > 0) {
		xa = x / param[2];         /* xa = x/param[2] */
		lxa = logf(xa);            /* lxa = log(x/param[2]) */
		xah = expf(lxa / param[3]); /* xah = exp(log(x/param[2])/param[3])
		                                  = (x/param[2])^(1/param[3]) */
		dy_dparam[1] = ex = expf(-xah);
		                           /* ex = exp(-(x/param[2])^(1/param[3])) */
		ex *= param[1];            /* ex = param[1] *
		                                     exp(-(x/param[2])^(1/param[3])) */
		*y = ex;                   /* y is now correct */
		ex *= xah / param[3];	   /* ex = param[1] * exp(...) *
		                                      (x/param[2])^(1/param[3]) *
		                                      1/param[3]   */
		dy_dparam[2] = ex / param[2];
		dy_dparam[3] = ex * lxa / param[3];
	} else if (x > -1e-10) {  // Almost zero
		*y = param[1];
		dy_dparam[1] = 1;
		dy_dparam[2] = dy_dparam[3] = 0;
	} else {
		/* Ouch: x<0 */
		fprintf(stderr, "Can't have x < 0 in stretched exponential!!\n");
	}
}

/* This is actually essentially the same as the other version; because
   of the power (x/tau)^(1/h), we cannot use the same tricks as for
   the multiexponential case, unfortunately. */

int stretchedexp_array(float xincr, float param[],
					   float *y, float **dy_dparam, int nx, int nparam)
{
	int i;
	float ex, lxa, xah, a2inv, a3inv;
	double xa, xaincr;

	if (xincr == 0) return -1;
	xa=0;
	xaincr = xincr / param[2];
	a2inv = 1/param[2];
	a3inv = 1/param[3];

	/* When x=0 */
	y[0] = param[1];
	dy_dparam[0][1] = 1;
	dy_dparam[0][2] = dy_dparam[0][3] = 0;

	for (i=1; i<nx; i++) {
		xa += xaincr;       /* xa = (xincr*i)/param[2] */
		lxa = logf((float) xa);     /* lxa = log(x/param[2]) */
		xah = expf(lxa * a3inv);  /* xah = exp(log(x/param[2])/param[3])
		                                = (x/param[2])^(1/param[3]) */
		dy_dparam[i][1] = ex = expf(-xah);
		                    /* ex = exp(-(x/param[2])^(1/param[3])) */
		ex *= param[1];     /* ex = param[1]*exp(-(x/param[2])^(1/param[3])) */
		y[i] = ex;          /* y is now correct */
		ex *= xah * a3inv;  /* ex = param[1] * exp(...) *
		                              (x/param[2])^(1/param[3]) * 1/param[3] */
		dy_dparam[i][2] = ex * a2inv;
		dy_dparam[i][3] = ex * lxa * a3inv;
	}

	return 0;
}


/********************************************************************

			   CHECKING AND RESTRICTED FITTING ROUTINES

*********************************************************************/


/* This is the default routine */
#define MIN_Z -1e5f
#define MIN_Z_FACTOR 0.4f
#define MAX_Z 1e10f  /* Silly */
#define MIN_A 0
#define MAX_A 1e10f  /* Silly */
#define MIN_TAU 0.001f  /* Works for lambda too */
#define MAX_TAU 1000
#define MIN_H 1
#define MAX_H 10

int check_ecf_params (float param[], int nparam,
					void (*fitfunc)(float, float [], float *, float [], int))
{
	if (fitfunc == GCI_multiexp_lambda || fitfunc == GCI_multiexp_tau) {
		switch (nparam) {
		case 3:
			if (param[0] < MIN_Z ||
				param[0] < -MIN_Z_FACTOR * fabs(param[1]) ||
				param[0] > MAX_Z)
				return -21;
			if (param[1] < MIN_A || param[1] > MAX_A)
				return -22;
			if (param[2] < MIN_TAU || param[2] > MAX_TAU)
				return -23;
			break;

		case 5:
			if (param[0] < MIN_Z ||
				param[0] < -MIN_Z_FACTOR * (fabs(param[1]) + fabs(param[3])) ||
				param[0] > MAX_Z)
				return -21;
			if (param[1] < MIN_A || param[1] > MAX_A)
				return -22;
			if (param[2] < MIN_TAU || param[2] > MAX_TAU)
				return -23;
			if (param[3] < MIN_A || param[3] > MAX_A)
				return -24;
			if (param[4] < MIN_TAU || param[4] > MAX_TAU)
				return -25;
			break;

		case 7:
			if (param[0] < MIN_Z ||
				param[0] < -MIN_Z_FACTOR * (fabs(param[1]) + fabs(param[3]) +
											fabs(param[5])) ||
				param[0] > MAX_Z)
				return -21;
			if (param[1] < MIN_A || param[1] > MAX_A)
				return -22;
			if (param[2] < MIN_TAU || param[2] > MAX_TAU)
				return -23;
			if (param[3] < MIN_A || param[3] > MAX_A)
				return -24;
			if (param[4] < MIN_TAU || param[4] > MAX_TAU)
				return -25;
			if (param[5] < MIN_A || param[5] > MAX_A)
				return -26;
			if (param[6] < MIN_TAU || param[6] > MAX_TAU)
				return -27;
			break;
		}
	} else if (fitfunc == GCI_stretchedexp) {
		if (param[0] < MIN_Z || param[0] < -MIN_Z_FACTOR * fabs(param[1]) ||
			param[0] > MAX_Z)
			return -21;
		if (param[1] < MIN_A || param[1] > MAX_A)
			return -22;
		if (param[2] < MIN_TAU || param[2] > MAX_TAU)
			return -23;
		if (param[3] < MIN_H || param[3] > MAX_H)
			return -24;
	}
	return 0;
}


/* For the user-specified version, we have some global variables to
   store the settings between invocations. */

static int restrain_nparam=0;  /* How many parameters have we set up? */
static int restrain_restraining[MAXFIT];  /* Do we check parameter i? */
static float restrain_minval[MAXFIT]; /* Minimum acceptable parameter value */
static float restrain_maxval[MAXFIT]; /* Maximum acceptable parameter value */

int GCI_set_restrain_limits(int nparam, int restrain[],
							float minval[], float maxval[])
{
	int i;

	if (nparam < 0 || nparam > MAXFIT)
		return -1;

	/* We're going to be doing something, so clear the memory */
	restrain_nparam = 0;

	for (i=0; i<nparam; i++) {
		if (restrain[i]) {
			restrain_restraining[i] = 1;

			if (minval[i] > maxval[i])
				return -2;
			restrain_minval[i] = minval[i];
			restrain_maxval[i] = maxval[i];
		} else
			restrain_restraining[i] = 0;
	}

	restrain_nparam = nparam;
	return 0;
}

/* original version
int check_ecf_user_params (float param[], int nparam,
					void (*fitfunc)(float, float [], float *, float [], int))
{
	int i;

	if (restrain_nparam != nparam) {
		dbgprintf(0, "Using user-defined parameter restraining with "
				  "wrong number of parameters:\n"
				  "actual nparam = %d, user restraining nparam = %d\n"
				  "Defaulting to standard tests\n", nparam, restrain_nparam);
		return check_ecf_params(param, nparam, fitfunc);
	}

	for (i=0; i<nparam; i++) {
		if (restrain_restraining[i]) {
			if ((param[i] < restrain_minval[i]) ||
				(param[i] > restrain_maxval[i]))
				return -21;
		}
	}

	return 0;
}
*/
// new version from J Gilbey 31.03.03
int check_ecf_user_params (float param[], int nparam,
                                        void (*fitfunc)(float, float [], float *, float [], int))
{
        int i;


        if (restrain_nparam != nparam) {
                dbgprintf(0, "Using user-defined parameter restraining with "
                                  "wrong number of parameters:\n"
                                  "actual nparam = %d, user restraining nparam = %d\n"
                                  "Defaulting to standard tests\n", nparam, restrain_nparam);
                return check_ecf_params(param, nparam, fitfunc);
        }


        for (i=0; i<nparam; i++) {
                if (restrain_restraining[i]) {
                        if (param[i] < restrain_minval[i])
                                param[i] = restrain_minval[i];
                        else if (param[i] > restrain_maxval[i])
                                param[i] = restrain_maxval[i];
                }
        }


        return 0;
}


/********************************************************************

			   ERROR BOUNDS AND CHI-SQUARED CALCULATION

*********************************************************************/

static float chisq50[MAXFIT+1] = { 0.0f, 0.45f, 1.39f, 2.37f, 3.36f,
                                   4.35f, 5.35f, 6.35f, 7.34f, 8.34f,
                                   9.34f, 10.34f, 11.34f, 12.34f, 13.34f,
                                   14.34f, 15.34f, 16.34f, 17.34f, 18.34f,
                                   19.34f };
static float chisq68[MAXFIT+1] = { 0.0f, 0.99f, 2.28f, 3.51f, 4.70f,
                                   5.86f, 7.01f, 8.14f, 9.27f, 10.39f,
                                   11.50f, 12.60f, 13.70f, 14.80f, 15.89f,
                                   16.98f, 18.07f, 19.15f, 20.23f, 21.31f,
                                   22.38f };
static float chisq90[MAXFIT+1] = { 0.0f, 2.71f, 4.61f, 6.25f, 7.78f,
                                   9.24f, 10.64f, 12.02f, 13.36f, 14.68f,
                                   15.99f, 17.27f, 18.55f, 19.81f, 21.06f,
                                   22.31f, 23.54f, 24.77f, 25.99f, 27.20f,
                                   28.41f };
static float chisq95[MAXFIT+1] = { 0.0f, 3.84f, 5.99f, 7.81f, 9.49f,
                                   11.07f, 12.59f, 14.07f, 15.51f, 16.92f,
                                   18.31f, 19.68f, 21.03f, 22.36f, 23.68f,
                                   25.00f, 26.30f, 27.59f, 28.87f, 30.14f,
                                   31.41f };

/**
 * Estimate errors for fitted parameters.
 * 
 * @param alpha    the inverse of the covariance matrix (with xero rows and
 *                 columns according to paramfree[])
 * @param nparam   the number of fitting parameters
 * @param mfit     the number of free parameters
 * @param d        output array for the eigenvalues
 * @param v        output matrix for the eigenvectors (produced with GCI_matrix)
 * @param interval chisquare percentage
 * @return         0 on success, < 0 on error
 *
 * This is based on "The Jacobi Method for Real Symmetric Matrices" by H.
 * Rutishauser, _Handbook for Automatic Computation, Volume II, Linear Algebra_,
 * Wilkinson & Reinsch, Springer-Verlag, 1971.
 */
int GCI_marquardt_estimate_errors(float **alpha, int nparam, int mfit,
								  float d[], float **v, float interval)
{
	float sm, c, s, t, h, g, tau, theta, thresh;
	float b[MAXFIT], z[MAXFIT], mult, chisq;
	int p, q, i, j, ret;
	
	switch ((int) (interval + 0.01f)) {
	case 50:
		chisq = chisq50[mfit];
		break;

	case 68:
		chisq = chisq68[mfit];
		break;

	case 90:
		chisq = chisq90[mfit];
		break;

	case 95:
		chisq = chisq95[mfit];
		break;

	default:
		ret = GCI_chisq(mfit, interval/100, &chisq);
		if (ret) return -1;
		break;
	}

	if (nparam > MAXFIT)
		return -2;

    // initialization	
	for (p = 0; p < nparam; ++p) {
		for (q = 0; q < nparam; ++q) {
			// identity matrix
			v[p][q] = (p == q) ? 1.0f : 0.0f;
		}
		b[p] = d[p] = alpha[p][p];
		z[p] = 0.0f;
	}
	
    // loop until convergence	
	for (i = 1; i <= 50; ++i) {
		
		// sum super-diagonal elements
		sm = 0.0f;
		for (p = 0; p < nparam; ++p) {
			for (q = p + 1; q < nparam; ++q) {
				sm += fabsf(alpha[p][q]);
			}
		}

		// check if zeroed
		if (0.0f == sm) {
			
			// Restore the alpha matrix before returning
			for (p = 0; p < nparam - 1; p++)
				for (q = p + 1; q < nparam; q++)
					alpha[p][q] = alpha[q][p];
			
			// Use the chisq values to find the semi-major axes
			for (p = 0; p < nparam; p++) {
				if (d[p] != 0.0f) {
					mult = sqrtf(chisq/d[p]);
					for (q = 0; q < nparam; q++)
						v[q][p] *= mult;
				}
			}

			return 0;
		}
		
        // first three sweeps use non-zero threshold		
		thresh = (i < 4) ? 0.2f * sm / (float)(nparam * nparam) : 0.0f;
		
		for (p = 0; p < nparam - 1; ++p) {
			for (q = p + 1; q < nparam; ++q) {
				g = 100.0f * fabsf(alpha[p][q]);
				if (i > 4
				    && fabsf(d[p]) + g == fabsf(d[p])
				    && fabsf(d[q]) + g == fabsf(d[q]))
				{
					alpha[p][q] = 0.0f;
				}
				else if (fabsf(alpha[p][q]) > thresh)
				{
					// rotate to zero diagonals
					h = d[q] - d[p];
					if (fabsf(h) + g == fabsf(h)) {
						t = alpha[p][q] / h;
					}
					else {
						// computing tangent of rotation angle
						theta = 0.5f * h / alpha[p][q];
						t = 1.0f / (fabsf(theta) + sqrtf(1.0f + theta*theta));
						if (theta < 0.0f) {
							t = -t;
						}
					}
					c = 1.0f / sqrtf(1.0f + t * t);
					s = t * c;
					tau = s / (1.0f + c);
					h = t * alpha[p][q];

					z[p] -= h;
					z[q] += h;
					d[p] -= h;
					d[q] += h;
					alpha[p][q] = 0.0f;
					
					for (j = 0; j < p; ++j) {
					    g = alpha[j][p];
					    h = alpha[j][q];
					    alpha[j][p] = g - s * (h + g * tau);
					    alpha[j][q] = h + s * (g - h * tau);
					}
					
					for (j = p + 1; j < q; ++j) {
						g = alpha[p][j];
						h = alpha[j][q];
						alpha[p][j] = g - s * (h + g * tau);
						alpha[j][q] = h + s * (g - h * tau);
					}

					for (j = q + 1; j < nparam; ++j) {
						g = alpha[p][j];
						h = alpha[q][j];
						alpha[p][j] = g - s * (h + g * tau);
						alpha[q][j] = h + s * (g - h * tau);
					}

					for (j = 0; j < nparam; ++j) {
						g = v[j][p];
						h = v[j][q];
						v[j][p] = g - s * (h + g * tau);
						v[j][q] = h + s * (g - h * tau);
					}
				}
			}
		}
		for (p = 0; p < nparam; ++p) {
			b[p] += z[p];
			d[p] = b[p];
			z[p] = 0.0f;
		}
    }
	// too many iterations
	return -1;
}

/**
 * Calculates lower incomplete gamma function by iterative power series expansion.
 */
float GCI_incomplete_gamma(float a, float x) {
	float returnValue = 0.0f;
	float multiplicand = powf(x, a) * expf(-x);
	int iter;
	float sum = 0.0f;
	float factor = 1.0f;
	float powerOfX = 1.0f;
	float value;
	
	for (iter = 0; iter < 1000; ++iter) {
		// this is "a * (a + 1) * (a + 2) * ... * (a + iter)"
		factor *= (a + (float) iter);
		
		// keep a sum
		sum += powerOfX / factor;
		
		// calculate value this iteration
		value = multiplicand * sum;
		if (value == returnValue) {
			// stopped changing, done
			return returnValue;
		}
		returnValue = value;
		
		// this will be "x ** iter" next time
		powerOfX *= x;
	}
	return returnValue;
}

/**
 * Log of Gamma function.
 * 
 * Based on LogGamma at http://www.johndcook.com/stand_alone_code.html:
 * "All code here is in the public domain. Do whatever you want with it, no
 * strings attached. Use at your own risk."
 * 
 * GCI_log_gamma and GCI_gamma are interdependent.
 * 
 * @param x must be positive
 * @return 
 */
float GCI_log_gamma(float x) {
	// Abramowitz and Stegun 6.1.41
    // Asymptotic series should be good to at least 11 or 12 figures
    // For error analysis, see Whittiker and Watson
    // A Course in Modern Analysis (1927), page 252
	
    static const double c[8] =
    {
		1.0/12.0,
		-1.0/360.0,
		1.0/1260.0,
		-1.0/1680.0,
		1.0/1188.0,
		-691.0/360360.0,
		1.0/156.0,
		-3617.0/122400.0
    };
    static const double halfLogTwoPi = 0.91893853320467274178032973640562;
	int i;
	double z, sum, series, logGamma;
	
    if (x < 12.0f)
    {
        return logf(fabsf(GCI_gamma(x)));
    }
	

    z = 1.0/(x*x);
    sum = c[7];
    for (i=6; i >= 0; i--)
    {
        sum *= z;
        sum += c[i];
    }
    series = sum/x;
	

    logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;    
	return (float) logGamma;
}

/**
 * Gamma function.
 * 
 * Based on Gamma at http://www.johndcook.com/stand_alone_code.html:
 * "All code here is in the public domain. Do whatever you want with it, no
 * strings attached. Use at your own risk."
 * 
 * GCI_gamma and GCI_log_gamma are interdependent.
 *
 * @param x must be positive
 * @return 
 */
float GCI_gamma(float x) {
	const double gamma = 0.577215664901532860606512090; // Euler's gamma constant
	
	// numerator coefficients for approximation over the interval (1,2)
	static const double p[] =
	{
		-1.71618513886549492533811E+0,
		2.47656508055759199108314E+1,
		-3.79804256470945635097577E+2,
		6.29331155312818442661052E+2,
		8.66966202790413211295064E+2,
		-3.14512729688483675254357E+4,
		-3.61444134186911729807069E+4,
		6.64561438202405440627855E+4
	};
		
	// denominator coefficients for approximation over the interval (1,2)
	static const double q[] =
	{
		-3.08402300119738975254353E+1,
		3.15350626979604161529144E+2,
		-1.01515636749021914166146E+3,
		-3.10777167157231109440444E+3,
		2.25381184209801510330112E+4,
		4.75584627752788110767815E+3,
		-1.34659959864969306392456E+5,
		-1.15132259675553483497211E+5
	};
	
	int n, i, arg_was_less_than_one;
	double y, z, num, den, result;
		
		
	// Split the function domain into three intervals:
	// (0, 0.001), [0.001, 12), and (12, infinity)
	
	///////////////////////////////////////////////////////////////////////////
	// First interval: (0, 0.001)
	//
	// For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
	// So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
	// The relative error over this interval is less than 6e-7.
	
	if (x < 0.001f)
		return (float)(1.0/(x*(1.0 + gamma*x)));
	
	///////////////////////////////////////////////////////////////////////////
	// Second interval: [0.001, 12)
    
	if (x < 12.0f)
	{
		// The algorithm directly approximates gamma over (1,2) and uses
		// reduction identities to reduce other arguments to this interval.
		
		y = x;
		n = 0;
		arg_was_less_than_one = (y < 1.0);
		
		// Add or subtract integers as necessary to bring y into (1,2)
		// Will correct for this below
		if (arg_was_less_than_one)
		{
			y += 1.0;
		}
		else
		{
			//n = static_cast<int> (floor(y)) - 1;  // will use n later
			n = (int) (floor(y)) - 1;
			y -= n;
		}
		
		num = 0.0;
		den = 1.0;
		
		z = y - 1;
		for (i = 0; i < 8; i++)
		{
			num = (num + p[i])*z;
			den = den*z + q[i];
		}
		result = num/den + 1.0;
		
		// Apply correction if argument was not initially in (1,2)
		if (arg_was_less_than_one)
		{
			// Use identity gamma(z) = gamma(z+1)/z
			// The variable "result" now holds gamma of the original y + 1
			// Thus we use y-1 to get back the orginal y.
			result /= (y-1.0);
		}
		else
		{
			// Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
			for (i = 0; i < n; i++)
				result *= y++;
		}
		
		return (float) result;
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Third interval: [12, infinity)
	
	if (x > 171.624f)
	{
		//printf("correct answer too large to display");
		return 0.0f; //TODO s/b +infinity
		// Correct answer too large to display. Force +infinity.
		//double temp = DBL_MAX;
		//return temp*2.0;
	}
	
	return expf(GCI_log_gamma(x));
}

float GCI_gammap(float a, float x) {
	return GCI_incomplete_gamma(a, x) / GCI_gamma(a);
}

/**
 * Find chi-square using bisection.
 * 
 * @param nu
 * @param chisq
 * @param root
 * @return 
 */
#define MAX_ITERS 40
#define ACC 0.0002f
int GCI_chisq(int nu, float chisq, float *root)
{
	int j;
	float x1, x2, val1, val2, mid, mid_val;
	
	x1 = 0.1f;
	x2 = 40.0f;
	if (chisq <= 0.0f || chisq >= 1.0f)
		return -1;

	if ((val1 = GCI_gammap(0.5f * (float) nu, 0.5f * x1)) < 0)
		return -3;

	j = 0;
	while ((val1 > chisq) && (++j <= MAX_ITERS)) {
		x1 /= 2;
		if ((val1 = GCI_gammap(0.5f * (float) nu, 0.5f * x1)) < 0)
			return -3;
	}
	if (j > MAX_ITERS)
		return -1;

	if ((val2 = GCI_gammap(0.5f * (float) nu, 0.5f * x2)) < 0)
		return -4;

	j = 0;
	while ((val2 < chisq) && (++j <= MAX_ITERS)) {
		x2 *= 2;
		if ((val2 = GCI_gammap(0.5f * (float) nu, 0.5f * x2)) < 0)
			return -4;
	}
	if (j > MAX_ITERS)
		return -2;

    // find value by bisection	
	for (j = 0; j < MAX_ITERS; ++j) {
		// compute bisection and value
		mid = (x1 + x2) / 2;
		mid_val = GCI_gammap(0.5f * (float) nu, 0.5f * mid);
		if (mid_val < 0.0f) {
			return -4;
		}
		
		// found the root?
		if (mid_val == chisq) {
			*root = mid;
			return 0;
		}
		
		// assign to endpoint of same sign
		if (mid_val * val1 > 0.0f) {
			x1 = mid;
			val1 = mid_val;
		}
		else {
			x2 = mid;
			val2 = mid_val;
		}
		
		// reached limit on interval size?
		if (fabsf(x2 - x1) < ACC && fabsf((x2 - x1) / mid) < ACC) {
			*root = mid;
			return 0;
		}
	}
	// too many iterations		
	return -5;
}
#undef MAX_ITERS
#undef ACC


/********************************************************************

							  MISCELLANY

*********************************************************************/


/* From filters */
int ECF_Find_Float_Max (float data[], int np, float *max_val)
{
	int i, maxi;
	float *ptr=data;

	if (np < 1)	   /* daft input */
		return -1;

	maxi = 0;
	*max_val = *ptr;

	for (i=1; i<np; i++)
		if (*++ptr>*max_val) {
			*max_val = *ptr;
			maxi = i;
		}

	return maxi;
}


/* Debugging version of printf */
#ifdef _CVI_
int dbgprintf(int dbg_level, const char *format, ...)
{
	int ret = 0;
	char msg[1000];
	va_list va;

	if (ECF_debug >= dbg_level) {
		va_start(va, format);
		ret = vsprintf(msg, format, va);  // Here's hoping, no vsnprintf ...
		MessagePopup("Debug info", msg);
		va_end(va);
	}

	return ret;
}
#else
int dbgprintf(int dbg_level, const char *format, ...)
{
	int ret = 0;
	va_list va;

	if (ECF_debug >= dbg_level) {
		va_start(va, format);
		ret = vfprintf(stderr, format, va);
		va_end(va);
	}

	return ret;
}
#endif

#ifdef TESTCHISQ
int main(int ac, char **av)
{
	int i, nu, ret;
	int debug = 0;
	float root, chisq;
	float percents[4] = { 0.50, 0.68, 0.90, 0.95 };

	for (i=0; i<4; i++) {
		chisq = percents[i];
		if (debug)
			fprintf(stderr, "Finding x for chisq=%.2f:\n", chisq);
		printf("static float chisq%d[MAXFIT+1] = { 0", (int)(100*chisq + 0.1));
		for (nu=1; nu<=20; nu++) {
			ret = GCI_chisq(nu, chisq, &root);
			if (ret != 0)
				fprintf(stderr, "error %d with nu=%d, chisq=%.2f\n",
						ret, nu, chisq);
			if (debug)
				fprintf(stderr, "nu=%2d  x=%6.2f  chisq(%d,%.3f)=%6.4f\n",
						nu, root, nu, root, GCI_gammp(0.5*nu, 0.5*root));
			printf(",%s%.2f", nu % 5 == 0 ? "\n                                   " : " ", root);
		}
		if (debug) fprintf(stderr, "\n");
		printf(" };\n");
	}
}
#endif

//***************************************** ExportParams ***********************************************/

void ECF_ExportParams_start (char path[])
{
	ecf_exportParams = 1;
	if (path) strcpy(ecf_exportParams_path, path);
}

void ECF_ExportParams_stop (void)
{
	ecf_exportParams = 0;
}

static FILE *ecf_exportFileStream;

void ecf_ExportParams_OpenFile (void)
{
	ecf_exportFileStream = fopen(ecf_exportParams_path, "a");
}

void ecf_ExportParams_CloseFile (void)
{
	fprintf(ecf_exportFileStream, "\n");
	fclose(ecf_exportFileStream);
}

void ecf_ExportParams (float param[], int nparam, float chisq)
{
	int i;

	for (i=0; i<nparam; i++) fprintf(ecf_exportFileStream, "%g, ", param[i]);
	fprintf(ecf_exportFileStream, "%g\n", chisq);
}

// Emacs settings:
// Local variables:
// mode: c
// c-basic-offset: 4
// tab-width: 4
// End:

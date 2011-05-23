//#include <ansi_c.h>
/* The 2010 version of the ECF library.  This has basically been
   completely rewritten to avoid license issues.  
   Also, this takes account of the fact that we may be
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

//TODO ARG
//#define SPEEDUP1 1
//#define SPEEDUP2 1
//#define SPEEDUP3 1

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
    float pivotInverse[n];
    int i, j, k, m;

    /*int q, w;
    printf("----------\n");
    for (q = 0; q < n; ++q) {
        for (w = 0; w < n; ++w) {
            printf("%f ", a[q][w]);
        }
        printf("  %f\n", b[q]);
    } */

    // base row of matrix
    for (k = 0; k < n - 1; ++k)
    {
        // search for line with max element
        max = fabs(a[k][k]);
        m = k;
        for (i = k + 1; i < n; ++i)
        {
            if (max < fabs(a[i][k])) // row i col k
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
            return -2; // singular matrix
        }

        // triangulation of matrix with coefficients
        pivotInverse[k] = 1.0 / a[k][k];
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
    pivotInverse[n - 1] = 1.0 / a[n - 1][n - 1];

    for (k = n - 1; k >= 0; --k)
    {
        for (i = k + 1; i < n; ++i)
        {
            b[k] -= a[k][i] * b[i];
        }
        b[k] *= pivotInverse[k];
    }

    /*printf("====>\n");
    for (q = 0; q < n; ++q) {
        for (w = 0; w < n; ++w) {
            printf("%f ", a[q][w]);
        }
        printf("  %f\n", b[q]);
    }*/

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
    float identity[n][n];
    float **work = GCI_ecf_matrix(n, n);
    int i, j, k;

    for (j = 0; j < n; ++j) {
        // find inverse by columns
        for (i = 0; i < n; ++i) {
            identity[j][i] = 0.0;
            // need a fresh copy of matrix a
            for (k = 0; k < n; ++k) {
                work[k][i] = a[k][i];
            }
        }
        identity[j][j] = 1.0;
        returnValue = GCI_solve_Gaussian(work, n, identity[j]);
        if (returnValue < 0) {
            return returnValue;
        }
    }
    GCI_ecf_free_matrix(work);

    // copy over results
    for (j = 0; j < n; ++j) {
        for (i = 0; i < n; ++i) {
            a[j][i] = identity[j][i];
        }
    }
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
    maxValue = fabs(a[col][col]);
    for (i = col + 1; i < n; ++i) {
        rowValue = fabs(a[i][col]);
        if (rowValue > maxValue) {
            pivotRow = i;
            maxValue = rowValue;
        }
    }

    #ifdef SPEEDUP1
    //TODO ARG this is actually slower!
    // swap rows
    float *ptr1;
    float *ptr2;
    if (pivotRow != col) {
        // swap elements in a matrix
        ptr1 = &a[col][0];
        ptr2 = &a[pivotRow][0];
        for (i = 0; i < n; ++i) {
            float temp;
            SWAP(*ptr1, *ptr2);
            ++ptr1;
            ++ptr2;
        }

        // swap elements in order vector
        {
            int temp;
            SWAP(order[col], order[pivotRow]);
        }
    }
    #else
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
    #endif
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
    #ifdef SPEEDUP2
    //TODO ARG this is *slightly* slower
    int *order_ptr = order;
    for (i = 0; i < n; ++i)
    {
        *order_ptr++ = i;
    }
    #else
    for (i = 0; i < n; ++i)
    {
        order[i] = i;
    }
    #endif

    // pivot first column
    pivot(a, n, order, 0);

    // check for singularity
    if (0.0 == a[0][0])
    {
        return -2;
    }

    // compute first row of upper
    inverse = 1.0 / a[0][0];
    #ifdef SPEEDUP3
    //TODO ARG this is *much* slower!!!
    //  Note compiler probably realizes a[0] is a constant anyway
    float *a_0_ptr = a[0];
    for (i = 1; i < n; ++i) {
        *a_0_ptr++ *= inverse;
    }
    #else
    for (i = 1; i < n; ++i) {
        a[0][i] *= inverse;
    }
    #endif

    // continue computing columns of lowers then rows of uppers
    for (jCol = 1; jCol < n - 1; ++jCol) {
        // compute column of lowers
        for (iRow = jCol; iRow < n; ++iRow) {
            sum = 0.0;
            for (kCol = 0; kCol < jCol; ++kCol) {
                sum += a[iRow][kCol] * a[kCol][jCol];
            }
            a[iRow][jCol] -= sum;
        }

        // find pivot for row
        pivot(a, n, order, jCol);
        if (0.0 == a[jCol][jCol])
        {
            return -2;
        }

        // build row of uppers
        inverse = 1.0 / a[jCol][jCol];
        for (kCol = jCol + 1; kCol < n; ++kCol) {
            sum = 0.0;
            for (iRow = 0; iRow < jCol; ++iRow) {
                sum += a[jCol][iRow] * a[iRow][kCol];
            }
            a[jCol][kCol] -= sum;
            a[jCol][kCol] *= inverse;
        }
    }

    // get remaining lower
    sum = 0.0;
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
//TODO check for lu[i][i] != 0?
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
        sum = 0.0;
        int jCol;
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
    int order[n];
    int return_value = lu_decomp(a, n, order);
    if (return_value >= 0) {
        return_value = solve_lu(a, n, b, order);
    }
    return return_value;
}

/* Matrix inversion by lower/upper decomposition.
   A is the n x n input matrix.
   On output, a is replaced by its matrix inverse..
 */
int GCI_invert_lu_decomp(float **a, int n)
{
    int returnValue;
    int order[n];
    float identity[n][n];
    int i, j;

    returnValue = lu_decomp(a, n, order);
    if (returnValue >= 0) {
        for (j = 0; j < n; ++j) {
            // find inverse by columns
            for (i = 0; i < n; ++i) {
                identity[j][i] = 0.0;
            }
            identity[j][j] = 1.0;
            solve_lu(a, n, identity[j], order);
        }
        for (j = 0; j < n; ++j) {
            for (i = 0; i < n; ++i) {
                a[j][i] = identity[j][i];
            }
        }
    }
    return returnValue;
}

/* Linear equation solution of Ax = b..
   A is the n x n input max, b is the right-hand side vector, length n.
   On output, b is replaced by the corresponding set of solution vectors
   and A is trashed.
 */
int GCI_solve(float **a, int n, float *b)
{
    //return GCI_solve_Gaussian(a, n, b);
    return GCI_solve_lu_decomp(a, n, b);
}

/* Matrix inversion.
   A is the n x n input matrix.
   On output, a is replaced by its matrix inverse..
 */
int GCI_invert(float **a, int n)
{
    //return GCI_invert_Gaussian(a, n);
    return GCI_invert_lu_decomp(a, n);
}


void GCI_covar_sort(float **covar, int nparam, int paramfree[], int mfit)
{
	int i, j, k;
	float temp;

	for (i=mfit; i<nparam; i++)
		for (j=0; j<=i; j++)
			covar[i][j] = covar[j][i] = 0.0;

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
    int row_size = nrows * sizeof(float *);
    int data_size = nrows * ncols * sizeof(float);
    unsigned char *raw = malloc(row_size + data_size);
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
	if ((marr = (float ***) malloc(nblocks * sizeof(float **))) == NULL)
		return NULL;

	/* allocate blocks (= pointers to rows) and set pointers to them */
	if ((marr[0] = (float **) malloc(nblocks * nrows * sizeof(float *)))
		== NULL) {
		free(marr);
		return NULL;
	}

	for (i=1; i<nblocks; i++)
		marr[i] = marr[i-1] + nrows;

	/* allocate rows (= pointers to column entries) and set pointers to them */
	if ((marr[0][0] = (float *)malloc(nblocks * nrows * ncols * sizeof(float)))
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
		dy_dparam[i] = ex = exp(-param[i+1] * x);
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
			dy_dparam[i][j] = ex = excur[j];
			ex *= param[j];
			y[i] += ex;
			dy_dparam[i][j+1] = -ex * xincr * i;
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
		dy_dparam[i] = ex = exp(-xa);
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
			dy_dparam[i][j] = ex = excur[j];
			ex *= param[j];
			y[i] += ex;
			dy_dparam[i][j+1] = ex * xincr * i * a2[j];
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
		lxa = log(xa);             /* lxa = log(x/param[2]) */
		xah = exp(lxa / param[3]); /* xah = exp(log(x/param[2])/param[3])
		                                  = (x/param[2])^(1/param[3]) */
		dy_dparam[1] = ex = exp(-xah);
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
		lxa = log(xa);      /* lxa = log(x/param[2]) */
		xah = exp(lxa * a3inv);  /* xah = exp(log(x/param[2])/param[3])
		                                = (x/param[2])^(1/param[3]) */
		dy_dparam[i][1] = ex = exp(-xah);
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
#define MIN_Z -1e5
#define MIN_Z_FACTOR 0.4
#define MAX_Z 1e10  /* Silly */
#define MIN_A 0
#define MAX_A 1e10  /* Silly */
#define MIN_TAU 0.001  /* Works for lambda too */
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

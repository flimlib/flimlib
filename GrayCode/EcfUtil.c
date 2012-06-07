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
#include <math.h>
#ifdef _CVI_
#include <userint.h>
#endif
#include "EcfInternal.h"  /* For #defines */

int ECF_debug = 0;


/********************************************************************

			   GAUSS-JORDAN AND COVAR-SORTING ROUTINES

*********************************************************************/


/* Linear equation solution by Gauss-Jordan elimination.
   a[0..n-1][0..n-1] is the input matrix,
   b[0..n] is input containing the single right-hand side vector.
   On output, a is replaced by its matrix inverse and b is replaced by
   the corresponding set of solution vectors. */

#define SWAP(a,b) { temp=(a); (a)=(b); (b)=temp; }

int GCI_gauss_jordan(float **a, int n, float *b)
{
	int indxc[MAXFIT], indxr[MAXFIT], ipiv[MAXFIT];
	int i, icol, irow, j, k, l, ll;
	float big, dum, pivinv, temp;

	for (j=0; j<n; j++)
		ipiv[j] = 0;
	for (i=0; i<n; i++) {
		big = 0.0;
		for (j=0; j<n; j++)
			if (ipiv[j] != 1)
				for (k=0; k<n; k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
				}

		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0; l<n; l++)
				SWAP(a[irow][l], a[icol][l]);
			if (b != NULL)
				SWAP(b[irow], b[icol]);
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0) {
                    printf("EcfUtil 71: Singular matrix\n");
			return -2; /* Singular Matrix */
        }

		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;
		for (l=0; l<n; l++)
			a[icol][l] *= pivinv;
		if (b != NULL)
			b[icol] *= pivinv;
		for (ll=0; ll<n; ll++)
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l=0; l<n; l++)
					a[ll][l] -= a[icol][l] * dum;
				if (b != NULL)
					b[ll] -= b[icol] * dum;
			}
	}

	for (l=n-1; l>=0; l--) {
		if (indxr[l] != indxc[l])
			for (k=0; k<n; k++)
				SWAP(a[k][indxr[l]], a[k][indxc[l]]);
	}
	return 0;
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

*********************************************************************/


/* This function, based on the public domain NR matrix function, allocates a float
   matrix with subscript range m[0..nrows-1][0..ncols-1] */

float **GCI_ecf_matrix(long nrows, long ncols)
{
	long i;
	float **m;

	/* allocate pointers to rows */
	if ((m = (float **) malloc(nrows * sizeof(float *))) == NULL)
		return NULL;

	/* allocate rows and set pointers to them */
	if ((m[0] = (float *) malloc(nrows * ncols * sizeof(float))) == NULL) {
		free(m);
		return NULL;
	}

	for (i=1; i<nrows; i++)
		m[i] = m[i-1] + ncols;

	/* return pointer to array of pointers to rows */
	return m;
}


void GCI_ecf_free_matrix(float **m)
{
	if (m != NULL) {
		free(m[0]);
		free(m);
	}
}


/* The next two functions are vaguely based on the NR f3tensor functions */
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

//printf("in multiexp_lambda_array xincr %g nx is %d nparam %d\n", xincr, nx, nparam);

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

        dbgprintf(0, "IN GCI SET RESTRAIN LIMITS");

	if (nparam < 0 || nparam > MAXFIT)
		return -1;

	/* We're going to be doing something, so clear the memory */
	restrain_nparam = 0;

	for (i=0; i<nparam; i++) {
		if (restrain[i]) {
			restrain_restraining[i] = 1;

			if (minval[i] > maxval[i]) {
                            dbgprintf(0, "minVal %g maxval %g at %d\n", minval[i], maxval[i], i);
				return -2;
                }
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


#define ROTATE(a, i, j, k, l) \
	g = a[i][j]; \
	h = a[k][l]; \
	a[i][j] = g - s * (h + g*tau); \
	a[k][l] = h + s * (g - h*tau);

/* This is jacobi.c and xjacobi.c suitably modified. */
/* Input:
     alpha     the inverse of the covariance matrix (with zero rows
               and columns according to paramfree[])
     nparam    the number of fitting parameters
     d[]       output array for the eigenvalues
     v[][]     output matrix for the eigenvectors (produced with
               GCI_matrix)
   Return:
      0  on success
     <0  on error
*/

static float chisq50[MAXFIT+1] = { 0, 0.45, 1.39, 2.37, 3.36,
                                   4.35, 5.35, 6.35, 7.34, 8.34,
                                   9.34, 10.34, 11.34, 12.34, 13.34,
                                   14.34, 15.34, 16.34, 17.34, 18.34,
                                   19.34 };
static float chisq68[MAXFIT+1] = { 0, 0.99, 2.28, 3.51, 4.70,
                                   5.86, 7.01, 8.14, 9.27, 10.39,
                                   11.50, 12.60, 13.70, 14.80, 15.89,
                                   16.98, 18.07, 19.15, 20.23, 21.31,
                                   22.38 };
static float chisq90[MAXFIT+1] = { 0, 2.71, 4.61, 6.25, 7.78,
                                   9.24, 10.64, 12.02, 13.36, 14.68,
                                   15.99, 17.27, 18.55, 19.81, 21.06,
                                   22.31, 23.54, 24.77, 25.99, 27.20,
                                   28.41 };
static float chisq95[MAXFIT+1] = { 0, 3.84, 5.99, 7.81, 9.49,
                                   11.07, 12.59, 14.07, 15.51, 16.92,
                                   18.31, 19.68, 21.03, 22.36, 23.68,
                                   25.00, 26.30, 27.59, 28.87, 30.14,
                                   31.41 };

int GCI_marquardt_estimate_errors(float **alpha, int nparam, int mfit,
								  float d[], float **v, float interval)
{
	int j, q, p, i, ret;
	float tresh, theta, tau, t, sm, s, h, g, c;
	float b[MAXFIT], z[MAXFIT], mult, chisq;

	switch ((int) (interval + 0.01)) {
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

	for (p=0; p<nparam; p++) {
		for (q=0; q<nparam; q++)
			v[p][q] = 0;
		v[p][p] = 1;
	}
	for (p=0; p<nparam; p++) {
		b[p] = d[p] = alpha[p][p];
		z[p] = 0;
	}
	for (i=1; i<=50; i++) {  /* iteration counter */
		sm = 0;
		for (p=0; p<nparam-1; p++) {
			for (q=p+1; q<nparam; q++)
				sm += fabs(alpha[p][q]);
		}
		if (sm == 0) {
			// Restore the alpha matrix before returning
			for (p=0; p<nparam-1; p++)
				for (q=p+1; q<nparam; q++)
					alpha[p][q] = alpha[q][p];
			// Use the chisq values to find the semimajor axes
			for (p=0; p<nparam; p++) {
				if (d[p] != 0) {
					mult = sqrt(chisq/d[p]);
					for (q=0; q<nparam; q++)
						v[q][p] *= mult;
				}
			}
			return 0;
		}
		if (i < 4)
			tresh = 0.2 * sm / (nparam * nparam);
		else
			tresh = 0.0;
		for (p=0; p<nparam-1; p++) {
			for (q=p+1; q<nparam; q++) {
				g = 100.0 * fabs(alpha[p][q]);
				if (i > 4 && (float) (fabs(d[p]) + g) == (float) fabs(d[p])
					&& (float) (fabs(d[q]) + g) == (float) fabs(d[q]))
					alpha[p][q] = 0;
				else if (fabs(alpha[p][q]) > tresh) {
					h = d[q] - d[p];
					if ((float) (fabs(h) + g) == (float) fabs(h))
						t = alpha[p][q] / h;
					else {
						theta = 0.5 * h / alpha[p][q];
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
						if (theta < 0) t = -t;
					}
					c = 1.0 / sqrt(1 + t*t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * alpha[p][q];
					z[p] -= h;
					z[q] += h;
					d[p] -= h;
					d[q] += h;
					alpha[p][q] = 0.0;
					for (j=0; j<p; j++) {
						ROTATE(alpha, j, p, j, q);
					}
					for (j=p+1; j<q; j++) {
						ROTATE(alpha, p, j, j, q);
					}
					for (j=q+1; j<nparam; j++) {
						ROTATE(alpha, p, j, q, j);
					}
					for (j=0; j<nparam; j++) {
						ROTATE(v, j, p, j, q);
					}
				}
			}
		}
		for (p=0; p<nparam; p++) {
			b[p] += z[p];
			d[p] = b[p];
			z[p] = 0;
		}
	}
	return -1;
}
#undef ROTATE

/* Find chisq inverse values.  This part is straight out of NR, with
   only minor modifications. */

float GCI_gammln(float xx)
{
	double x, y, tmp, ser;
	static double cof[6] = { 76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155,
		0.1208650973866179e-2, -0.5395239384953e-5 };
	int j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x+0.5) * log(tmp);
	ser = 1.000000000190015;
	for (j=0; j<=5; j++)
		ser += cof[j] / ++y;
	return (float) (-tmp + log(2.5066282746310005 * ser / x));
}

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

int GCI_gcf(float *gammcf, float a, float x, float *gln)
{
	int i;
	float an, b, c, d, del, h;

	*gln = GCI_gammln(a);
	b = x + 1.0 - a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h = d;
	for (i=1; i<=ITMAX; i++) {
		an = -i * (i-a);
		b += 2.0;
		d = an*d + b;
		if (fabs(d) < FPMIN) d = FPMIN;
		c = b + an/c;
		if (fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) return -1;
	/* was: nrerror("a too large, ITMAX too small in gcf"); */
	*gammcf = exp(-x + a*log(x) - (*gln)) * h;
	return 0;
}
#undef FPMIN

int GCI_gser(float *gamser, float a, float x, float *gln)
{
	int n;
	float sum, del, ap;

	*gln = GCI_gammln(a);
	if (x <= 0.0) {
		if (x < 0.0)
			return -1;
		/* was: nrerror("x less than 0 in routine gser"); */
		*gamser = 0.0;
		return 0;
	} else {
		ap = a;
		del = sum = 1.0 / a;
		for (n=0; n<ITMAX; n++) {
			++ap;
			del *= x / ap;
			sum += del;
			if (fabs(del) < fabs(sum) * EPS) {
				*gamser = sum * exp(-x + a*log(x) - (*gln));
				return 0;
			}
		}
		return -1;
		/* was: nrerror("a too large, ITMAX too small in routine gser"); */
	}
}
#undef ITMAX
#undef EPS

float GCI_gammp(float a, float x)
{
	float gamser, gammcf, gln;

	if (x < 0.0 || a <= 0.0) return -1;
	/* was: nrerror("Invalid arguments in routine gammp"); */
	if (x < (a+1.0)) {
		if (GCI_gser(&gamser, a, x, &gln) != 0)
			return -1;
		else
			return gamser;
	} else {
		if (GCI_gcf(&gammcf, a, x, &gln) != 0)
			return -1;
		else
			return 1.0 - gammcf;
	}
}

/* Modified from rtbis for our needs: we replace the general func by a
   chisq-evaluating funtion. */
#define JMAX 40
#define ACC 0.0002
int GCI_chisq(int nu, float chisq, float *root)
{
	int j;
	float x1, x2, dx, f, fmid, xmid, rtb;

	x1 = 0.1;
	x2 = 40;
	if (chisq <= 0 || chisq >= 1)
		return -1;

	if ((f = GCI_gammp(0.5*nu, 0.5*x1)) < 0)
		return -3;

	j = 0;
	while ((f > chisq) && (++j <= JMAX)) {
		x1 /= 2;
		if ((f = GCI_gammp(0.5*nu, 0.5*x1)) < 0)
			return -3;
	}
	if (j > JMAX)
		return -1;

	if ((fmid = GCI_gammp(0.5*nu, 0.5*x2)) < 0)
		return -4;

	j = 0;
	while ((fmid<chisq) && (++j <= JMAX)) {
		x2 *= 2;
		if ((fmid = GCI_gammp(0.5*nu, 0.5*x2)) < 0)
			return -4;
	}
	if (j > JMAX)
		return -2;

	dx = x2 - x1;
	rtb = x1;

	for (j=0; j<JMAX; j++) {
		xmid = rtb + (dx /= 2);
		if ((fmid = GCI_gammp(0.5*nu, 0.5*xmid)) < 0)
			return -4;
		if (fmid <= chisq) rtb = xmid;
		if ( ((dx < ACC) && (dx/rtb < ACC)) || fmid == chisq) {
			*root = rtb;
			return 0;
		}
	}
	return -5;
	/* was: nrerror("Too many bisections in rtbis"); */
}
#undef JMAX
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

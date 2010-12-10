#include <ansi_c.h>
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

/*
 * Based on LinearEquationSolving from linearEquations.c by Henry Guennadi Levkin.
 */
int GCI_gauss_jordan(float **a, int n, float *b)
{
    float max;
    float tmp;

    int i, j, k, m;

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
                tmp = a[k][i];
                a[k][i] = a[m][i];
                a[m][i] = tmp;
            }
            tmp = b[k];
            b[k] = b[m];
            b[m] = tmp;
        }

        if (0.0 == a[k][k])
        {
            return -2; // singular matrix
        }

        // triangulation of matrix with coefficients
        for (j = k + 1; j < n; ++j) // current row of matrix
        {
            tmp = -a[j][k] / a[k][k];
            for (i = k; i < n; ++i)
            {
                a[j][i] += tmp * a[k][i];
            }
            b[j] += tmp * b[k]; // free member recalculation
        }
    }

    for (k = n - 1; k >= 0; --k)
    {
        for (i = k + 1; i < n; ++i)
        {
            b[k] -= a[k][i] * b[i];
        }
        b[k] /= a[k][k];
    }
    return 0;
}

//==============================================================================
// return 1 if system not solving
// nDim - system dimension
// pfMatr - matrix with coefficients
// pfVect - vector with free members
// pfSolution - vector with system solution
// pfMatr becames trianglular after function call
// pfVect changes after function call
//
// Developer: Henry Guennadi Levkin
//
//==============================================================================
int LinearEquationsSolving(int nDim, double* pfMatr, double* pfVect, double* pfSolution)
{
  double fMaxElem;
  double fAcc;

  int i , j, k, m;


  for(k=0; k<(nDim-1); k++) // base row of matrix
  {
    // search of line with max element
    fMaxElem = fabs( pfMatr[k*nDim + k] );
    m = k;
    for(i=k+1; i<nDim; i++)
    {
      if(fMaxElem < fabs(pfMatr[i*nDim + k]) )
      {
        fMaxElem = pfMatr[i*nDim + k];
        m = i;
      }
    }

    // permutation of base line (index k) and max element line(index m)
    if(m != k)
    {
      for(i=k; i<nDim; i++)
      {
        fAcc               = pfMatr[k*nDim + i];
        pfMatr[k*nDim + i] = pfMatr[m*nDim + i];
        pfMatr[m*nDim + i] = fAcc;
      }
      fAcc = pfVect[k];
      pfVect[k] = pfVect[m];
      pfVect[m] = fAcc;
    }

    if( pfMatr[k*nDim + k] == 0.) return 1; // needs improvement !!!

    // triangulation of matrix with coefficients
    for(j=(k+1); j<nDim; j++) // current row of matrix
    {
      fAcc = - pfMatr[j*nDim + k] / pfMatr[k*nDim + k];
      for(i=k; i<nDim; i++)
      {
        pfMatr[j*nDim + i] = pfMatr[j*nDim + i] + fAcc*pfMatr[k*nDim + i];
      }
      pfVect[j] = pfVect[j] + fAcc*pfVect[k]; // free member recalculation
    }
  }

  for(k=(nDim-1); k>=0; k--)
  {
    pfSolution[k] = pfVect[k];
    for(i=(k+1); i<nDim; i++)
    {
      pfSolution[k] -= (pfMatr[k*nDim + i]*pfSolution[i]);
    }
    pfSolution[k] = pfSolution[k] / pfMatr[k*nDim + k];
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

/* This function allocates a float matrix with subscript range
 * m[0..nrows-1][0..ncols-1]
 */
float **GCI_ecf_matrix(long nrows, long ncols)
{
    int row_size = nrows * sizeof(float *);
    int data_size = nrows * ncols * sizeof(float);
    unsigned char *raw = malloc(row_size + data_size);
    if (NULL == raw)
    {
        return NULL;
    }
    float **row = (float **) raw;
    float *data = (float *) (row + nrows);
    int i;
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

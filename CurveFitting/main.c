// INCLUDES & DEFINES FOR GSL

/*#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
 
#define N 40
 */

// INCLUDES & DEFINES FOR NR
/*
#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define NPT 100
#define MA 3
#define SPREAD 0.001
*/

// COMBINED FOR GSL & NR
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "nr.h"
#include "nrutil.h"

// for timing
#include "time.h"

// for levmar
#include "levmar.h" // <levmar.h>

#ifndef LM_DBL_PREC
#error Example program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif
// end levmar

#define N 40

#define NPT 40
#define MA 3
#define SPREAD 0.001

// Paul Barber code
typedef enum { NOISE_CONST, NOISE_GIVEN, NOISE_POISSON_DATA,
	NOISE_POISSON_FIT, NOISE_GAUSSIAN_FIT, NOISE_MLE } noise_type;

double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*10)/CLOCKS_PER_SEC;
	return diffms;
}

void generateSingleExponentialDataFloat(float *xPtr, float *yPtr, float *sigmaPtr, int m, float xOffset, float xInc, float A, float lambda, float b) {
	printf("generate %d float exponential data, x %f %f, A %f lambda %f b %f\n", m, xOffset, xInc, A, lambda, b);
	float x = xOffset;
	for (int i = 0; i < m; ++i) {
		*(xPtr++) = x;
		*(yPtr++) = A * expf(-lambda*x) + b;
		*(sigmaPtr++) = 0.1;
		x += xInc;
	}
}

void generateSingleExponentialDataDouble(double *xPtr, double *yPtr, float *sigmaPtr, int m, double xOffset, float xInc, double A, double lambda, double b) {
	printf("generate %d double exponential data, x %f %f, A %f lambda %f b %f\n", m, xOffset, xInc, A, lambda, b);
	float x = xOffset;
	for (int i = 0; i < m; ++i) {
		*(xPtr++) = x;
		*(yPtr++) = A * exp(-lambda*x) + b;
		*(sigmaPtr++) = 0.1;
		x += xInc;
	}
}


// GSL CODE:

/*
 *  expfit.c
 *  CurveFitting
 *
 *  Created by Aivar Grislis on 3/17/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
/* expfit.c -- model functions for exponential + background */

struct data {
	size_t n;
	double * y;
	double * sigma;
};

int
expb_f (const gsl_vector * x, void *data, 
		gsl_vector * f)
{
	size_t n = ((struct data *)data)->n;
	double *y = ((struct data *)data)->y;
	double *sigma = ((struct data *) data)->sigma;
	
	double A = gsl_vector_get (x, 0);
	double lambda = gsl_vector_get (x, 1);
	double b = gsl_vector_get (x, 2);
	
	size_t i;
	
	for (i = 0; i < n; i++) // loops over data; doesn't handle multi-component
	{
		/* Model Yi = A * exp(-lambda * i) + b */
		double t = i + 1;  // GSL is indexed starting with 0 & NR at 1, compensate
		double Yi = A * exp (-lambda * t) + b;
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	
	return GSL_SUCCESS;
}

int
expb_df (const gsl_vector * x, void *data, 
		 gsl_matrix * J)
{
	size_t n = ((struct data *)data)->n;
	double *sigma = ((struct data *) data)->sigma;
	
	double A = gsl_vector_get (x, 0);
	double lambda = gsl_vector_get (x, 1);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double t = i + 1; // GSL is indexed starting with 0 & NR at 1, compensate
		double s = sigma[i];
		double e = exp(-lambda * t);
		gsl_matrix_set (J, i, 0, e/s); 
		gsl_matrix_set (J, i, 1, -t * A * e/s);
		gsl_matrix_set (J, i, 2, 1/s);
	}
	return GSL_SUCCESS;
}

int
expb_fdf (const gsl_vector * x, void *data,
		  gsl_vector * f, gsl_matrix * J)
{
	expb_f (x, data, f);
	expb_df (x, data, J);
	
	return GSL_SUCCESS;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s);

int
do_gsl_fit(void)
{
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int i, iter = 0;
	const size_t n = N;
	const size_t p = 3;
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	double y[N], sigma[N];
	struct data d = { n, y, sigma};
	gsl_multifit_function_fdf f;
	//double x_init[3] = { 1.0, 0.0, 0.0 };
	double x_init[3] = { 4.5, 2.2, 2.8 };
	gsl_vector_view x = gsl_vector_view_array (x_init, p);
	const gsl_rng_type * type;
	gsl_rng * r;
	
	printf("========START GSL FIT==========");
	
	gsl_rng_env_setup();
	
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);
	
	f.f = &expb_f;
	f.df = &expb_df;
	f.fdf = &expb_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;
	
	printf("=====GSL FIT=====\n");
	
	/* This is the data to be fitted */
	
	for (i = 0; i < n; i++)
	{
		double t = i + 1;  // GSL is indexed starting with 0 & NR at 1, compensate
		y[i] = 1.0 + 5.0 * exp (-0.1 * t); // NO NOISE
		//NO NOISE + gsl_ran_gaussian (r, 0.1);
		sigma[i] = 0.1;
		printf ("%f %f %f\n", i, y[i], sigma[i]);
	};
	printf("\n");
	
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	
	print_state (iter, s);
	
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		
		printf ("status = %s\n", gsl_strerror (status));
		
		print_state (iter, s);
		
		if (status)
			break;
		
		status = gsl_multifit_test_delta (s->dx, s->x,
										  1e-4, 1e-4);
	}
	while (status == GSL_CONTINUE && iter < 500);
	
	gsl_multifit_covar (s->J, 0.0, covar);
	
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	
	{ 
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
		
		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
		
		printf ("A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("b      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
	}
	
	printf ("status = %s\n", gsl_strerror (status));
	
	printf("=====END GSL FIT=====\n");
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	gsl_rng_free (r);
	return 0;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
	printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
			"|f(x)| = %g\n",
			iter,
			gsl_vector_get (s->x, 0), 
			gsl_vector_get (s->x, 1),
			gsl_vector_get (s->x, 2), 
			gsl_blas_dnrm2 (s->f));
}

// LEVMAR FIT:

/* model to be fitted to measurements: x_i = p[0]*exp(-p[1]*i) + p[2], i=0...n-1 */
void expfunc(double *p, double *x, int m, int n, void *data)
{
	register int i;
	
	for(i=0; i<n; ++i){
		x[i]=p[0]*exp(-p[1]*(1.0 + i)) + p[2];
	}
}

/* Jacobian of expfunc() */
void jacexpfunc(double *p, double *jac, int m, int n, void *data)
{   
	register int i, j;
	printf("LEVMAR GET JACOB FOR PARAMS %f %f %f\n", p[0], p[1], p[2]);
	
	/* fill Jacobian row by row */
	for(i=j=0; i<n; ++i){
		float e = exp(-p[i]*(1.0 + i));
		jac[j++]=e;
		jac[j++]=-p[0]*(1.0 + i) * e;
		jac[j++]=1.0;
	}
}
		
int
do_levmar_fit(void)
{
	const int n=40, m=3; // 40 measurements, 3 parameters
	double p[m], x[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	register int i;
	int ret;
	
	printf("========= BEGIN LEVMAR FIT ============");
	
	for(i=0; i<n; ++i) {
		x[i]=(5.0*exp(-0.1*(1.0+i)) + 1.0);
		printf("%f %f\n", 1.0 + i, x[i]);
	}
	printf("\n");
	
	/* initial parameters estimate: (1.0, 0.0, 0.0) */
	//p[0]=1.0; p[1]=0.0; p[2]=0.0;
	p[0]=4.5; p[1]=2.2; p[2]=2.8;
	
	/* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 
	
	/* invoke the optimization function */
	ret=dlevmar_der(expfunc, jacexpfunc, p, x, m, n, 1000, opts, info, NULL, NULL, NULL); // with analytic Jacobian
	//ret=dlevmar_dif(expfunc, p, x, m, n, 1000, opts, info, NULL, NULL, NULL); // without Jacobian
	printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
	printf("Best fit parameters: %.7g %.7g %.7g\n", p[0], p[1], p[2]);
	printf("========== END LEVMAR FIT ===========");
	return 0;
}

// NR CODE:

//void fmodel(x,a,y,dyda,na)   
//float x,a[],*y,dyda[];   
//int na;   
void fmodel(float x, float a[], float *y,float dyda[] , int na)   
{
int i;
	float A, lambda, b;
	
	*y = 0.0f;
	
	// loop for multi-component fits
	//JUST IN CASE THIS IS BUGGY for (i = 1; i < na - 1; i+=3) {
	{
		A      = a[i];
		lambda = a[i+1];
		b      = a[i+2];
		//printf("A is %f lambda is %f b is %f", A, lambda, b);
		
        float e = exp(-lambda * x);
		
		*y += A * e + b;
		
		//printf("x is %f y is %f\n", x, *y);

		float sigma = 0.1;
		
		dyda[i]   = e;         // note, divide all three by sigma, like GSL????
		dyda[i+1] = -x * A * e;
		dyda[i+2] = 1.0;
    }
}

do_nr_fit(void)
{
	long idum=(-911);
	int i,*ia,iter,itst,j,k,mfit=MA;
	float alamda,chisq,ochisq,*x,*y,*sig,**covar,**alpha;
	static float a[MA+1]=
	{0.0,5.0,2.0,3.0};
	static float gues[MA+1]=
	{0.0,4.5,2.2,2.8};
	
	ia=ivector(1,MA);
	x=vector(1,NPT);
	y=vector(1,NPT);
	sig=vector(1,NPT);
	covar=matrix(1,MA,1,MA);
	alpha=matrix(1,MA,1,MA);

// OLD DATA (SUM OF TWO GAUSSIANS)
	/* First try a sum of two Gaussians */
/*	for (i=1;i<=NPT;i++) {
		x[i]=0.1*i;
		y[i]=0.0;
		for (j=1;j<=MA;j+=3) {
			y[i] += a[j]*exp(-SQR((x[i]-a[j+1])/a[j+2]));
		}
		y[i] *= (1.0+SPREAD*gasdev(&idum));
		sig[i]=SPREAD*y[i];
	}
*/
// NEW DATA (EXPONENTIAL)
	const gsl_rng_type * type;
	gsl_rng * r;
	
	gsl_rng_env_setup();
	
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);
	
	printf("=====NR FIT=====\n");
	float expectedA = 1.0; //5.0;
	float expectedLambda = 2.0; //0.1;
	float expectedb = 1.0;
	generateSingleExponentialDataFloat(&x[1], &y[1], &sig[1], 40, 1.0, 1.0, expectedA, expectedLambda, expectedb);
	for (int i = 1; i <= 40; ++i) {
		printf("%d %f %f %f\n", i, x[i], y[i], sig[i]);
	}
	printf("\n");
	
	for (i=1;i<=mfit;i++) ia[i]=1;
	for (i=1;i<=MA;i++) a[i]=gues[i];
	//ARG for (iter=1;iter<=2;iter++) {
	for (iter=1; iter <=1; ++iter) {
		alamda = -1; // signal to mrqmin to initialize
		mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fmodel,&alamda);
		k=1;
		itst=0;
		for (;;) {
			printf("\n%s %2d %17s %10.4f %10s %9.2e\n","Iteration #",k,
				   "chi-squared:",chisq,"alamda:",alamda);
			printf("%8s %8s %8s\n",
				   "a[1]","a[2]","a[3]");
			for (i=1;i<=MA;i++) printf("%e",a[i]);
			printf("\n");
			k++;
			ochisq=chisq;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fmodel,&alamda);
			if (chisq > ochisq)
				itst=0;
			else if (fabs(ochisq-chisq) < 0.1) {
				printf("CHISQ IS CLOSE old %e new %e\n", ochisq, chisq);
				itst++;
		    }
			if (itst < 4) continue;
			alamda=0.0;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fmodel,&alamda);
			printf("\nUncertainties:\n");
			for (i=1;i<=MA;i++) printf("%9.4f",sqrt(covar[i][i]));
			printf("\n");
			printf("\nExpected results:\n");
			printf(" %8.2f %8.2f %8.2f\n",
				   expectedA, expectedLambda, expectedb);
			break;
		}
		if (iter == 1) {
			printf("ITER IS ONE!!!\n");
		}
		break; //ARG
		if (0 && iter == 1) {
			printf("press return to continue with constraint\n");
			(void) getchar();
			printf("holding a[2] /*and a[5]*/ constant\n");
			for (j=1;j<=MA;j++) a[j] += 0.1;
			a[2]=2.0;
			ia[2]=0;
			//a[5]=5.0; reduced MA from 6 to 3 to handle only one component
			//ia[5]=0;
		}
	}
	
	printf("=====END NR FIT=====\n");
	
	free_matrix(alpha,1,MA,1,MA);
	free_matrix(covar,1,MA,1,MA);
	free_vector(sig,1,NPT);
	free_vector(y,1,NPT);
	free_vector(x,1,NPT);
	free_ivector(ia,1,MA);
	return 0;

	// THE OLD VERSION
	// used MA 6
	/*
	long idum=(-911);
	int i,*ia,iter,itst,j,k,mfit=MA;
	float alamda,chisq,ochisq,*x,*y,*sig,**covar,**alpha;
	static float a[MA+1]=
	{0.0,5.0,2.0,3.0,2.0,5.0,3.0};
	static float gues[MA+1]=
	{0.0,4.5,2.2,2.8,2.5,4.9,2.8};
	
	ia=ivector(1,MA);
	x=vector(1,NPT);
	y=vector(1,NPT);
	sig=vector(1,NPT);
	covar=matrix(1,MA,1,MA);
	alpha=matrix(1,MA,1,MA);
	/ * First try a sum of two Gaussians * /
	for (i=1;i<=NPT;i++) {
		x[i]=0.1*i;
		y[i]=0.0;
		for (j=1;j<=MA;j+=3) {
			y[i] += a[j]*exp(-SQR((x[i]-a[j+1])/a[j+2]));
		}
		y[i] *= (1.0+SPREAD*gasdev(&idum));
		sig[i]=SPREAD*y[i];
	}
	for (i=1;i<=mfit;i++) ia[i]=1;
	for (i=1;i<=MA;i++) a[i]=gues[i];
	for (iter=1;iter<=2;iter++) {
		alamda = -1;
		mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fgauss,&alamda);
		k=1;
		itst=0;
		for (;;) {
			printf("\n%s %2d %17s %10.4f %10s %9.2e\n","Iteration #",k,
				   "chi-squared:",chisq,"alamda:",alamda);
			printf("%8s %8s %8s %8s %8s %8s\n",
				   "a[1]","a[2]","a[3]","a[4]","a[5]","a[6]");
			for (i=1;i<=6;i++) printf("%9.4f",a[i]);
			printf("\n");
			k++;
			ochisq=chisq;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fgauss,&alamda);
			if (chisq > ochisq)
				itst=0;
			else if (fabs(ochisq-chisq) < 0.1)
				itst++;
			if (itst < 4) continue;
			alamda=0.0;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fgauss,&alamda);
			printf("\nUncertainties:\n");
			for (i=1;i<=6;i++) printf("%9.4f",sqrt(covar[i][i]));
			printf("\n");
			printf("\nExpected results:\n");
			printf(" %7.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
				   5.0,2.0,3.0,2.0,5.0,3.0);
			break;
		}
		if (iter == 1) {
			printf("press return to continue with constraint\n");
			(void) getchar();
			printf("holding a[2] and a[5] constant\n");
			for (j=1;j<=MA;j++) a[j] += 0.1;
			a[2]=2.0;
			ia[2]=0;
			a[5]=5.0;
			ia[5]=0;
		}
	}
	free_matrix(alpha,1,MA,1,MA);
	free_matrix(covar,1,MA,1,MA);
	free_vector(sig,1,NPT);
	free_vector(y,1,NPT);
	free_vector(x,1,NPT);
	free_ivector(ia,1,MA);
	return 0;
*/
}

/*
 
// MY VERSION; GET ALL OF J AT ONCE
 
//
// Build an approximate Hessian by multiplying
// J' * J, then augment the diagonal.
//
float lambda = 0.1;
float dotProduct;

// for all columns
for (j = 0; j < nParams; ++j) {
	// will build (lower) triangle
	for (i = 0; i <= j; ++i) {
		// get dot product of row i of J'
 	    // and column j of J
		dotProduct = 0.0;
		for (k = 0; k < nData; ++k) {
			dotProduct += J[k][i] * J[k][j];
		}
		if (i == j) {
			// augment the diagonal
			A[i][j] = dotProduct * (1 + lambda);
		}
		else {
			// matrix is symmetric
			A[i][j] = A[j][i] = dotProduct;
		}
	}
}
 
// MY VERSION; GET J ROW AT A TIME
 
// zero out lower triangle
for (j = 0; j < nParams; ++j) {
    for (i = 0; i <= j; ++i) {
	    A[i][j] = 0.0;
    }
}
 
for (k = 0; k < nData; ++k) {
    // Jrow[] = ???
    for (j = 0; j < nParams; ++j) {
       for (i = 0; i <= j; ++i) {
           A[i][j] += Jrow[i]*Jrow[j];
           }
    }
 }
 // once you optimize this, add B vector,
 // add sig[i] weighting, compute chisq,
 // this will look about the same.
 

// ANY SENSE IN A COL AT A TIME?
 
// same as 1st except
// get Jcol[] (replaces J[][j]) for j loop
// get Jcol[] (replaces J[][i]) for i loop
// Getting the derivative for one fit param at a time is a bit weird
 
*/

/* FOR REFERENCE MRQCOF
 
 void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[],
			int ma, float **alpha, float beta[], float *chisq,
			void (*funcs)(float, float [], float *, float [], int))
{
	int i,j,k,l,m,mfit=0;
	float ymod,wt,sig2i,dy,*dyda;
	
	dyda=vector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				beta[j] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i;
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_vector(dyda,1,ma);
}
*/

// Odd bug:
//
// This was declared in the following style:
// void fmodel(x,a,y,dyda,na)   
// float x,a[],*y,dyda[];   
// int na;   
//
// x was always 0.0, other parameters worked fine!?

void nrSolutionModel(float x, float a[], float *y, float dyda[], int na) {
	float A, lambda, b;
	
	A      = a[1];
	lambda = a[2];
	b      = a[3];
	
	//printf("A %f lambda %f b %f ", A, lambda, b);
	
	float e = exp(-lambda * x);
	
	//*y += A * e + b; // WHY PLUS EQUALS?? (for components!!!)
	*y = A * e + b;
	
	dyda[1] = e;
	dyda[2] = -A * x * e;
	dyda[3] = 1.0;
	//printf("x is %f", x);
	//printf("Jacobean row is %f %f %f\n", dyda[1], dyda[2], dyda[3]);
}


void nrSolution(float x[], float y[], float sig[], int ndata, float a[], int ma, float **alpha, float beta[], float *chisq, void (*funcs)(float, float [], float *, float [], int)) {
	int i,j,k,l,m;
	float ymod,wt,sig2i,dy,*dyda;
	
	dyda=vector(1,ma);
	for (j=1;j<=ma;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy = y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			{
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++) {
					alpha[j][++k] += wt*dyda[m];
				}
				beta[j] += dy*wt;
			}
		}
		*chisq+=dy*dy*sig2i;
	}
	
	for (j=2;j<=ma;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_vector(dyda,1,ma);
}

//int onetime = 1;
void altSolutionModel(float x[], int nData, float a[], int nParams, float y[], float **J) {
	float A, lambda, b;
	
	A      = a[0];
	lambda = a[1];
	b      = a[2];
	
    for (int i = 0; i < nData; ++i) {
		
		float tmpX = x[i];
		
		float e = exp(-lambda * tmpX);
		
		y[i] = A * e + b;
		
		//if (onetime) {
		//	printf("A %f lambda %f b %f\n", A, lambda, b);
		//	printf("i %d x %f y %f\n", i, x[i], y[i]);
		//}
		
		J[i][0] = e;
		J[i][1] = -A * tmpX * e;
		J[i][2] = 1.0;
		//printf("Jacobean row %d is %f %f %f\n", i, J[i][1], J[i][2], J[i][3]);
	}
	//onetime = 0;
}

void altSolution(float x[], float y[], float sigma[], int nData, float a[], int nParams, float **alpha, float beta[], float *chisq, void (*modelFunction)(float[], int, float[], int, float[], float **)) {
	int i, j, k;
	float betaSum;
	float lambda = 0.1;
	float dotProduct;
	float weight;
	float yDiff;
	
	float *yCalc = vector(0, nData-1);
	float **J = matrix(0, nData-1, 0, nParams-1);
	
    (*modelFunction)(x, nData, a, nParams, yCalc, J);
	
	//
	// Build an approximate Hessian by multiplying
	// J^T * J, then augment the diagonal.
	//
	// Note H or Hessian w/b a better name than alpha, A usually means J.
	//
    *chisq = 0.0;
	
	// for all columns
	for (j = 0; j < nParams; ++j) {
		betaSum = 0.0;
		// row loop, only need to consider lower triangle
		for (i = 0; i <= j; ++i) {
			dotProduct = 0.0;
			// for all data
			for (k = 0; k < nData; ++k) {
				weight = 1.0 / (sigma[k] * sigma[k]);
				
				// get dot product of row i of J^T and col j of J
				dotProduct += (J[k][i] * J[k][j]) * weight;
				
				// once per column and data index
				if (0 == i) {
					// compute different between data and model
					yDiff = y[k] - yCalc[k];
					
					// accumulate approximate gradient
					betaSum += (yDiff * J[k][j]) * weight;
					
					// once per data index
					if (0 == j) {
						// calculate chi square
						*chisq += (yDiff * yDiff) * weight;
					}
					//printf("k %d j %d y diffs %f J %f betaSum %f\n", k, j, (y[k] - yCalc[k]), J[k][j], betaSum);
				}
			}
			alpha[i][j] = dotProduct;
			if (i != j) {
				// matrix is symmetric
				alpha[j][i] = dotProduct;
			}
		}
        beta[j] = betaSum;
		//printf("BETASUM %f\n", betaSum);
	}
	
	free_vector(yCalc, 0, nData-1);
	free_matrix(J, 0, nData-1, 0, nParams-1);
}
/*
void altSolution2(float x[], float y[], float sigma[], int nData, float a[], int fixed[], int nParams, float **alpha, float beta[], float *chisq, void (*modelFunction)(float[], int, float[], int, float[], float **)) {
	int i, j, k;
	int i_fit, j_fit;
	float betaSum;
	float lambda = 0.1;
	float dotProduct;
	float weight;
	float yDiff;
	
	float *yCalc = vector(0, nData-1);
	float **J = matrix(0, nData-1, 0, nParams-1);
	
    (*modelFunction)(x, nData, a, nParams, yCalc, J);
	
	//
	// Build an approximate Hessian by multiplying
	// J^T * J, then augment the diagonal.
	//
	// Note H or Hessian w/b a better name than alpha, A usually means J.
	//
    *chisq = 0.0;
	
	// keep track of indices for fitted parameters
	i_fit = 0;
	j_fit = 0;
	
	// initiliaze sums
	betaSum = 0.0;
	dotProduct = 0.0;
	
	// for all columns
	for (i = 0; i < nParams; ++i) {
		// row loop, only need to consider lower triangle
		for (j = 0; j < i; ++j) {
			// for all data
			for (k = 0; k < nData; ++k) {
				weight = 1.0 / (sigma[k] * sigma[k]); // c/b: if ((fit[i] && fit[j]) || (0 == i && 0 == j))
				if (fit[i] && fit[j]) {
					dotProduct += (J[k][j] * J[k][i]) * weight;
				}
				
				// once per column and data index
				if (0 == j) {
					// computer difference between data and model
					yDiff = y[k] - yCalc[k];
					
					// once per data index
					if (0 == i) {
						// accumulate chi square
						*chisq += (yDiff * yDiff) * weight; //TODO sure this is the NR way?
					}
					
					// accumulate approprimate gradient
					if (fit[i]) {
					    betaSum += (yDiff * J[k][i]) * weight);
					}
				}
			}

			if (fit[i] && fit[j]) {
				alpha[j][i] = dotProduct;
				if (i != j) {
					// matrix is symmetric
					alpha[i][j] = dotProduct;
				}
				dotProduct = 0.0;
				++j_fit;
			}
		}
		
		if (fit[i]) {			
			beta[i] = betaSum;
			betaSum = 0.0;
			++i_fit;
		}
	}
	
	for (k = 0; k < nData; ++k) {
		weight = 1.0 / (sigma[k] * sigma[k]);
        yDiff = y[k] - yCalc[k];       //TODO 'Error' is inappropriate term; residual; diff
		
		i_fit = 0;
		for (i = 0; i < nParams; ++i) {
			if (fit[i] ) {
				j_fit = 0;
				for (j = 0; j < i; ++j) {
					if (fit[j]) {
						alpha[j_fit][i_fit] += (J[k][j] * J[k][i]) * weight;						
						++j_fit;
					}
				}
				beta[i_fit] += (yDiff * J[k][i]) * weight;	
				++i_fit;
			}
		}				
	    *chisq += (yDiff * yDiff) * weight;
	}
	
	// for all columns //TODO WHY J I K??
	for (j = 0; j < nParams; ++j) {
		betaSum = 0.0;
		// row loop, only need to consider lower triangle
		for (i = 0; i <= j; ++i) {
			dotProduct = 0.0;
			// for all data
			for (k = 0; k < nData; ++k) {
				weight = 1.0 / (sigma[k] * sigma[k]);
				
				// get dot product of row i of J^T and col j of J
				dotProduct += (J[k][i] * J[k][j]) * weight; // only if fit i and fit j
				
				// once per column and data index
				if (0 == i) {
					// compute different between data and model
					yDiff = y[k] - yCalc[k];
					
					// once per data index
					if (0 == j) {
						// calculate chi square
						*chisq += (yDiff * yDiff) * weight;
					}
					
					// accumulate approximate gradient
					betaSum += (yDiff * J[k][j]) * weight; // only if !fixed[j]
					
					//printf("k %d j %d y diffs %f J %f betaSum %f\n", k, j, (y[k] - yCalc[k]), J[k][j], betaSum);
				}
			}
			alpha[i][j] = dotProduct;  // only if !fixed[i] && !fixed[j]
			if (i != j) {
				// matrix is symmetric
				alpha[j][i] = dotProduct; // only if !fixed[i] && !fixed[j]
			}
		}
        beta[j] = betaSum; // only if !fixed[j]
		//printf("BETASUM %f\n", betaSum);
	}
	
	free_vector(yCalc, 0, nData-1);
	free_matrix(J, 0, nData-1, 0, nParams-1);
}
*/

float setUpMatrices(float x[], float y[], float sigma[], int nData, float a[], int nParams, float **alpha, float beta[],
					void (*modelFunction)(float[], int, float[], int, float[], float **),
					int noise, int version) {
	int i, j, k;
	float betaSum;
	float dotProduct;
	float weight;
	float yDiff;
	float chiSquare = 0.0;	
	float *yCalc = vector(0, nData-1);
	float **J = matrix(0, nData-1, 0, nParams-1);

	// fake the fitted (should this param be fitted?) vector
	int *fitted = vector(0, nParams-1);
	for (int i = 0; i < nParams; ++i) {
		fitted[i] = 1;
	}
	
	// get model fit
    (*modelFunction)(x, nData, a, nParams, yCalc, J);

	if (0) 
	for (int k = 0; k < nData; ++k) {
		for (int i = 0; i < nParams; ++i) {
			printf(" %f ", J[k][i]);
		}
		printf("\n");
	}
	
    switch (version) {
		case 0:
		{
			//
			// Build an approximate Hessian by multiplying
			// J^T * J and an approximate Gradient by multiplying
			// J^T * f.
			//
			// Note H or Hessian w/b a better name than alpha, A is usually J.
			//
			
			for (j = 0; j < nParams; ++j) {                                       // for all columns
				betaSum = 0.0;
				for (i = 0; i <= j; ++i) {                                        // row loop, need only consider lower triangle
					dotProduct = 0.0;
					for (k = 0; k < nData; ++k) {                                 // for all data
						weight = 1.0 / (sigma[k] * sigma[k]);                       // NOTE inefficiency here!  c/b computed one time only!
						dotProduct += J[k][i] * J[k][j] * weight;                 // get dot product of row i of J^T and col j of J
						if (0 == i) {                                             // once per column for every data index
							yDiff = y[k] - yCalc[k];                              // compute difference between data and model
							betaSum += yDiff * J[k][j] * weight;                 // accumulate approximate gradient
							if (0 == j) {                                         // once per data index
								chiSquare += (yDiff * yDiff) * weight;          // calculate chi square
							}
						}
					}
					alpha[i][j] = dotProduct;
					if (i != j) {
						alpha[j][i] = dotProduct;                                 // matrix is symmetrical
					}
				}
				beta[j] = betaSum;
			}
		}
		break;
			
		case 1:
		{
			// Same as above approach but introduces the 'fitted[]' array.
			// keep track of indices for fitted parameters
			int i_fit = 0;
			int j_fit = 0;
			
			// initiliaze sums
			betaSum = 0.0;
			dotProduct = 0.0;
			
			// for all columns
			for (i = 0; i < nParams; ++i) {
				// row loop, only need to consider lower triangle
				for (j = 0; j <= i; ++j) {
					// for all data
					for (k = 0; k < nData; ++k) {
						weight = 1.0 / (sigma[k] * sigma[k]); // c/b: if ((fit[i] && fit[j]) || (0 == i && 0 == j))
						if (fitted[i] && fitted[j]) {
							dotProduct += J[k][j] * J[k][i] * weight;						}
						
						// once per column and data index
						if (0 == j) {
							// compute difference between data and model
							yDiff = y[k] - yCalc[k];
							
							// once per data index
							if (0 == i) {
								// accumulate chi square
								chiSquare += (yDiff * yDiff) * weight;
							}
							
							// accumulate approprimate gradient
							if (fitted[i]) {
								betaSum += (yDiff * J[k][i]) * weight;
							}
						}
					}
					
					if (fitted[i] && fitted[j]) {
						alpha[j][i] = dotProduct;
						if (i != j) {
							// matrix is symmetric
							alpha[i][j] = dotProduct;
						}
						dotProduct = 0.0;
						++j_fit;
					}
				}
				
				if (fitted[i]) {			
					beta[i] = betaSum;
					betaSum = 0.0;
					++i_fit;
				}
			}
		}
		break;
			
			
		case 2:
		{
			int i_fit, j_fit;
			
			// zero the matrices
			for (i = 0; i < nParams; ++i) {
				for (j = 0; j <= i; ++j) {
					alpha[j][i] = 0.0;
				}
				beta[i] = 0.0;
			}
			chiSquare = 0.0;
			
			// Loop through the data
            for (k = 0; k < nData; ++k) {
				weight = 1.0 / (sigma[k] * sigma[k]);
				yDiff = y[k] - yCalc[k];
				
				i_fit = 0;
				for (i = 0; i < nParams; ++i) {
					if (fitted[i] ) {
						j_fit = 0;
						for (j = 0; j <= i; ++j) {
							if (fitted[j]) {
								alpha[j_fit][i_fit] += (J[k][j] * J[k][i]) * weight;						
								++j_fit;
							}
						}
						beta[i_fit] += (yDiff * J[k][i]) * weight;	
						++i_fit;
					}
				}				
				chiSquare += (yDiff * yDiff) * weight;
			}
			
			// copy over triangle
			for (i = 0; i < nParams; ++i) {
				for (j = 0; j < i; ++j) {
					if (i != j) {
					    alpha[i][j] = alpha[j][i];
					}
				}
			}
		}
		break;
	}
	
	free_vector(yCalc, 0, nData-1);
	free_matrix(J, 0, nData-1, 0, nParams-1);
	
	// free fake fitted
	free_vector(fitted, 0, nParams-1);
	
	return chiSquare;
}
	

void solveMatrices(float lambda, float solution[], float **alpha, float beta[], int nParams) {
	float **oneda = matrix(1, nParams, 1, 1); // NR KLUDGE SO I CAN TEMPORARILY CALL NR'S GAUSS JORDAN; NR makes this a static
	float **alphaNR = matrix(1, nParams, 1, nParams);

	printf("%e lambda is\n", lambda);

	// augment the diagonal; this is the key to the Levenberg-Marquardt algorithm
	++lambda; // (1 + lambda)
	for (int i = 0; i < nParams; ++i) {
		alpha[i][i] *= lambda;
		oneda[i+1][1] = beta[i]; // KLUDGE
		for (int j = 0; j < nParams; ++j) {
			alphaNR[i+1][j+1] = alpha[i][j];
		}
	}
	for (int i = 0; i < nParams; ++i) {
		for (int j = 0; j < nParams; ++j) {
			printf(" %e ", alpha[i][j]);
		}
		printf(" X  %e\n", oneda[i+1][1]);
	}
	gaussj(alphaNR, nParams, oneda, 1); // FROM NR
	for (int i = 0; i < nParams; ++i) {
		solution[i] = oneda[i+1][1]; // KLUDGE
		for (int j = 0; j < nParams; ++j) {
			alpha[i][j] = alphaNR[i+1][j+1];
		}
	}
}

void altSolver(float x[], float y[], float sigma[], int nData, float a[], int nParams,
			   void (*modelFunction)(float[], int, float[], int, float[], float **),
			   int noise, int version) {
	int iterations = 0;
	float chiSquare;
	float solvedChiSquare;
	float lambda = 0.001; // NR 15.5
	int approximated = 0;
	int decreaseCounter = 0;
	float **alpha = matrix(0, nParams-1, 0, nParams-1);
	float *beta = vector(0, nParams-1);
	float *solution = vector(0, nParams-1);
	
	chiSquare = setUpMatrices(x, y, sigma, nData, a, nParams, alpha, beta, modelFunction, noise, version);
	while (!approximated) {
		++iterations;
		solveMatrices(lambda, solution, alpha, beta, nParams);
		for (int i = 0; i < nParams; ++i) {
			printf(" %e ", solution[i]);
			a[i] += solution[i];
			printf(" --> %e ", a[i]);
		}
		printf("\n");
		if (abs(a[0] + a[1] + a[2]) < 0.00001) { // KLUDGE quick, dirty and inaccurate! //NOT WORKING
			approximated = 1;
		}
		// are proposed changes trivial?
		//if (0) { //(norm(solution)  <= 2) {
		//	approximated = 1;
		//}
		//else {
		//}

        // solving |lambda * I + A^T * A| solution = -A ^ T * f
		solvedChiSquare = setUpMatrices(x, y, sigma, nData, a, nParams, alpha, beta, modelFunction, noise, version);
		printf("chi sq is %e\n", solvedChiSquare);
		
		if (solvedChiSquare <= chiSquare) {  // <= solved a problem when chiSquare became zero!  o'wise approximated never got set
			// successful step; favor Gauss-Newton
			lambda /= 10;
			printf("successful step; decrease lambda\n");
			
			// NR 15.5: "Stop 1st or 2nd time chi square decreases < 0.01"
			//printf("chi square decreases %f\n", chiSquare - solvedChiSquare);
			if (chiSquare - solvedChiSquare < 0.001) {
				printf("decreased < 0.01\n");
				++decreaseCounter;
				if (decreaseCounter > 1) {
					approximated = 1;
				}
			}
			else {
				decreaseCounter = 0;
			}
			
			chiSquare = solvedChiSquare;
		}
		else {
			// unsuccessful step; favor steepest descent
			lambda *= 10;
			decreaseCounter = 0;
			printf("unsuccessful step; increate lambda\n");
			
			// back out changes
			for (int i = 0; i < nParams; ++i) {
				a[i] -= solution[i];
			}
		}
	}
	printf("%d iterations\n", iterations);
	
	free_matrix(alpha, 0, nParams-1, 0, nParams-1);
	free_vector(beta, 0, nParams-1);
}
	
/* 
 * This code was used to compare NR with my version initially.
 *
 * It compares one cycle of NR with mine.
 */
void checkStrategies() {
	printf("CHECK STRATEGIES\n");

	int nParams = 3;
	int nData = 40;
	float chisq;
	float **alphaNR = matrix(1, nParams, 1, nParams);
	float *xNR = vector(1, nData);
	float *yNR = vector(1, nData);
	float *sigmaNR = vector(1, nData);
	float *aNR = vector(1, nParams);
	float *betaNR = vector(1, nParams);
	for (int i = 1; i <= nParams; ++i) {
		aNR[i] = i * 0.5;
	}
	for (int i = 1; i <= nData; ++i) {
		xNR[i] = i * 0.5;
		yNR[i] = i * 0.5;
		sigmaNR[i] = i * 0.1;
	}
	clock_t beginNR = clock();
	nrSolution(xNR, yNR, sigmaNR, nData, aNR, nParams, alphaNR, betaNR, &chisq, nrSolutionModel);
	clock_t endNR = clock();
	printf("NR took %f\n", diffclock(endNR, beginNR));
	for (int i = 1; i <= nParams; ++i) {
		for (int j = 1; j <= nParams; ++j) {
			printf(" %e ", alphaNR[i][j]);
		}
		printf("-     %e\n", betaNR[i]);
	}
	printf("chiSquare is %e\n", chisq);
	
	float **alphaALT = matrix(0, nParams-1, 0, nParams-1);
	float *xALT = vector(0, nData-1);
	float *yALT = vector(0, nData-1);
	float *sigmaALT = vector(0, nData-1);
	float *aALT = vector(0, nParams-1);
	float *betaALT = vector(0, nParams-1);
	for (int i = 0; i < nParams; ++i) {
		aALT[i] = (i + 1) * 0.5;
	}
	for (int i = 0; i < nData; ++i) {
		xALT[i] = (i + 1) * 0.5;
		yALT[i] = (i + 1) * 0.5;
		sigmaALT[i] = (i + 1) * 0.1;
	}
	clock_t beginALT = clock();
	altSolution(xALT, yALT, sigmaALT, nData, aALT, nParams, alphaALT, betaALT, &chisq, NOISE_CONST);
	clock_t endALT = clock();
	printf("ALT took %f\n", diffclock(endALT, beginALT));
	for (int i = 0; i < nParams; ++i) {
		for (int j = 0; j < nParams; ++j) {
			printf(" %e ", alphaALT[i][j]);
		}
		printf("-     %e\n", betaALT[i]);
	}
	printf("chiSquare is %e\n", chisq);
}

void do_alt_fit(int noise, int version) {
	printf("=====ALT FIT=====%d\n", version);

	int nParams = 3;
	int nData = 40;
	float *x = vector(0, nData-1);
	float *y = vector(0, nData-1);
	float *sigma = vector(0, nData-1);
	
	float expectedA = 1.0; //ARG 5.0;
	float expectedLambda = 2.0; //ARG 0.1;
	float expectedb = 1.0;

	generateSingleExponentialDataFloat(&x[0], &y[0], &sigma[0], 40, 1.0, 1.0, expectedA, expectedLambda, expectedb);
	for (int i = 0; i < 40; ++i) {
		printf("DATA: %d %f %f %f\n", i, x[i], y[i], sigma[i]);
	}
	
	// guess initial fit params
	float *a = vector(0, nParams-1);
	a[0] = 4.5;
	a[1] = 2.2;
	a[2] = 2.8; // same as NR uses
	
	altSolver(x, y, sigma, nData, a, nParams, altSolutionModel, noise, version);
	
	printf("solution is:\n");
	for (int i = 0; i < nParams; ++i) {
		printf(" %f\n", a[i]);
	}
	printf("=====END ALT FIT=====%d\n", version);
}
	
void gaussianModel(float x[], int nData, float a[], int nParams, float y[], float **J) {
	float A, B, C;
	
	A = a[0];
	B = a[1];
	C = a[2];
	
	for (int i = 0; i < nData; ++i) {
		float e = exp(-((x[i] - B)*(x[i] - B))/(2*C*C));
		
		y[i] = A * e;
		
		J[i][0] = e;
		J[i][1] = A * (x[i] - B)/(C*C);
		J[i][2] = A * (x[i] - B)*(x[i] - B)/(C*C*C);
	}
}

void fitGaussian() {
	printf("=====GAUSSIAN FIT=====\n");
	int nData = 4;
	float *x = vector(0, nData-1);
	x[0] = 1.0;
	x[1] = 2.0;
	x[2] = 3.0;
	x[3] = 4.0;
	float *y = vector(0, nData-1);
	y[0] = 56610.0;
	y[1] = 283550.0;
	y[2] = 2431860.0;
	y[3] = 3543000.0;
	float *sigma = vector(0, nData-1);
	sigma[0] = 0.1;
	sigma[1] = 0.1;
	sigma[2] = 0.1;
	sigma[3] = 0.1;
	
	// guess initial fit params
	int nParams = 3;
	float *a = vector(0, nParams-1);
	a[0] = 100.0;
	a[1] = 100.0;
	a[2] = 100.0;
	
	altSolver(x, y, sigma, nData, a, nParams, gaussianModel, NOISE_CONST, 0);
	
	printf("Gaussian solution is:\n");
	for (int i = 0; i < nParams; ++i) {
		printf(" %f\n", a[i]);
	}
	//TODO free vectors, etc.
	printf("===GAUSSIAN FIT=======\n");
}
	
	
	
	int
	main(void)
	{
		// fitGaussian was an attempt to fit a Gaussian to the leading edge of the decay curve.
		//
		// seemed really flakey, dependent on the initial guesses, chi square was huge
		// the equation I had didn't give any vertical adjustment of the Gaussian.
		//
		// SPCImage says "The SPCImage software provides an 'estimation' of this response function
		// by calculating the first derivative of the rising part of the fluoreoscence."
		//
		//
		//fitGaussian();
		//checkStrategies();
		//do_levmar_fit();
		clock_t beginALT = clock();
		do_alt_fit(NOISE_CONST, 0);
		printf("%f\n", diffclock(clock(), beginALT));
		clock_t beginALTnoise = clock();
		do_alt_fit(NOISE_CONST, 1);
		printf("%f\n", diffclock(clock(), beginALTnoise));
		clock_t beginALT2 = clock();
		do_alt_fit(NOISE_CONST, 2);
		printf("%f\n", diffclock(clock(), beginALT2));
		//do_gsl_fit();
		//do_levmar_fit();
		clock_t beginNR = clock();
        do_nr_fit();
		printf("%f\n", diffclock(clock(), beginNR));
}

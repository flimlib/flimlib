/*
 *  tests.c
 *  Untitled
 *
 *  Created by Aivar Grislis on 9/22/11.
 *  Copyright 2011 LOCI. All rights reserved.
 *
 */

#include "EcfInternal.h"
#include "json.h"
#include "output.h"
#include "parser.h"
#include "tests.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Predeclarations
// the next fn uses GCI_triple_integral_*() to fit repeatedly until chisq_target is met
int Old_GCI_triple_integral_fitting_engine(float xincr, float y[], int fit_start, int fit_end,
										   float instr[], int ninstr, noise_type noise, float sig[],
										   float *Z, float *A, float *tau, float *fitted, float *residuals,
										   float *chisq, float chisq_target);

// the next fn uses GCI_marquardt_instr() to fit repeatedly until chisq_target is met
int Old_GCI_marquardt_fitting_engine(float xincr, float *trans, int ndata, int fit_start, int fit_end, 
									 float prompt[], int nprompt,
									 noise_type noise, float sig[],
									 float param[], int paramfree[],
									 int nparam, restrain_type restrain,
									 void (*fitfunc)(float, float [], float *, float [], int),
									 float *fitted, float *residuals, float *chisq,
									 float **covar, float **alpha, float **erraxes,
									 float chisq_target, float chisq_delta, int chisq_percent);

int Old_GCI_marquardt_step(float x[], float y[], int ndata,
						   noise_type noise, float sig[],
						   float param[], int paramfree[], int nparam,
						   restrain_type restrain,
						   void (*fitfunc)(float, float [], float *, float [], int),
						   float yfit[], float dy[],
						   float **covar, float **alpha, float *chisq,
						   float *alambda, int *pmfit, float *pochisq, float *paramtry, float *beta, float *dparam);

int Old_GCI_marquardt_step_instr(float xincr, float y[],
								 int ndata, int fit_start, int fit_end,
								 float instr[], int ninstr,
								 noise_type noise, float sig[],
								 float param[], int paramfree[], int nparam,
								 restrain_type restrain,
								 void (*fitfunc)(float, float [], float *, float [], int),
								 float yfit[], float dy[],
								 float **covar, float **alpha, float *chisq,
								 float *alambda, int *pmfit, float *pochisq, float *paramtry, float *beta, float *dparam,
								 float **pfnvals, float ***pdy_dparam_pure, float ***pdy_dparam_conv,
								 int *pfnvals_len, int *pdy_dparam_nparam_size);

int GCI_marquardt_step(float x[], float y[], int ndata,
					   noise_type noise, float sig[],
					   float param[], int paramfree[], int nparam,
					   restrain_type restrain,
					   void (*fitfunc)(float, float [], float *, float [], int),
					   float yfit[], float dy[],
					   float **covar, float **alpha, float *chisq,
					   float *alambda, int *pmfit, float *pochisq, float *paramtry, float *beta, float *dparam);

int GCI_marquardt_step_instr(float xincr, float y[],
							 int ndata, int fit_start, int fit_end,
							 float instr[], int ninstr,
							 noise_type noise, float sig[],
							 float param[], int paramfree[], int nparam,
							 restrain_type restrain,
							 void (*fitfunc)(float, float [], float *, float [], int),
							 float yfit[], float dy[],
							 float **covar, float **alpha, float *chisq,
							 float *alambda, int *pmfit, float *pochisq, float *paramtry, float *beta, float *dparam,
							 float **pfnvals, float ***pdy_dparam_pure, float ***pdy_dparam_conv,
							 int *pfnvals_len, int *pdy_dparam_nparam_size);

void dump1D(float *data, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        printf("%f ", data[i]);
    }
    printf("\n");
}

void dump2D(float **data, int n) {
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            printf("%f ", data[i][j]);
        }
        printf("\n");
    }
}

void dump1Dint(int *data, int n) {
	int i;
	for (i = 0; i < n; ++i) {
		printf("%d ", data[i]);
	}
	printf("\n");
}

/*
 * Compares a test generated integer output with the output value specified in the JSON document.
 * Must be an exact match.
 * If no output value is specified, this test generated integer output is added to the JSON document
 */
void compareIntExact(json_t *outputs, char *name, int value, int *changed, int *success) {
	int unlikely_value = 12345678;
	int goldenValue = getInt(outputs, name, unlikely_value);
	
	if (goldenValue == unlikely_value) { //TODO camelcase vs underscore naming?
		json_insert_child(outputs, putInt(name, value));
		*changed = TRUE;
	}
	else {
		printf("%s got %d expected %d\n", name, value, goldenValue, value);
		if (value != goldenValue) {
			printf("FAIL %d\n", *changed);
			*changed = FALSE;
			printf("success was %d\n", *success);
			*success = FALSE;
		}
	}
}

/*
 * Compares a test generated integer output with the output value specified in the JSON document.
 * Must match within tolerance.
 * If no output value is specified, this test generated integer output is added to the JSON document
 */
//TODO ARG not used; not tested
void compareInt(json_t *outputs, char *name, int value, float tolerance, int *changed, int *success) {
	int unlikely_value = 12345678;
	int goldenValue = getInt(outputs, name, unlikely_value);
	
	if (goldenValue == unlikely_value) {
		json_insert_child(outputs, putInt(name, value));
		*changed = TRUE;
	}
	else {
		float percent = 100.0f * (value - goldenValue) / goldenValue;
		printf("%s percentage difference %f\n", name, percent);
		if (fabs(percent) > tolerance) {
			printf("FAIL, tolerance is %f\n", tolerance);
			*success = FALSE;
		}
	}
}

/*
 * Compares a test generated integer array output with the output value specified in the JSON document.
 * Must match within tolerance.
 * If no output value is specified, this test generated integer output is added to the JSON document
 */
//TODO ARG not used; not tested
void compareInt1D(json_t *outputs, char *name, int *values, int size, float tolerance, int *changed, int *success) {
	int success_local = TRUE;
	int goldenSize;
	int *goldenValues = getInt1D(outputs, name, &goldenSize);
	int i;
	
	if (NULL == goldenValues) {
		// if no outputs specified, the current values are written to JSON
		json_insert_child(outputs, putInt1D(name, values, size));
		*changed = TRUE;
	}
	else {
		// check whether within tolerance
		float *percents = malloc(size * sizeof(float));
		
		assert(size == goldenSize);
		
		for (i = 0; i < size; ++i) {
			percents[i] = 100.0f * (values[i] - goldenValues[i]) / goldenValues[i];
			if (fabs(percents[i]) > tolerance) {
				success_local = FALSE;
			}
		}
		printf("%s percentage difference:\n", name);
		dump1D(percents, size);
		if (!success_local) {
			printf("FAIL, tolerance is %f\n", tolerance);
			*success = FALSE;
		}
	}
}

/*
 * Compares a test generated float output with the output value specified in the JSON document.
 * Must match within tolerance.
 * If no output value is specified, this test generated integer output is added to the JSON document
 */
void compareFloat(json_t *outputs, char *name, float value, float tolerance, int *changed, int *success) {
	float goldenValue = getFloat(outputs, name, 1234.5678f); //TODO bad!
	
	if (1234.5678f == goldenValue) {
		// if no outputs specified, the current values are written to JSON
		json_insert_child(outputs, putFloat(name, value));
		*changed = TRUE;
	}
	else {
		// check whether within tolerance
		float percent = 100.0f * (value - goldenValue) / goldenValue;
		printf("%s percentage difference %f\n", name, percent);
		if (fabs(percent) > tolerance) {
			printf("FAIL, tolerance is %f\n", tolerance);
			*success = FALSE;
		}
	}
}

/*
 * Compares a test generated float array output with the output value specified in the JSON document.
 * Must match within tolerance.
 * If no output value is specified, this test generated float output is added to the JSON document
 */
void compareFloat1D(json_t *outputs, char *name, float *values, int size, float tolerance, int *changed, int *success) {
	int success_local = TRUE;
	int goldenSize;
	float *goldenValues = getFloat1D(outputs, name, &goldenSize);
	int i;
	
	if (NULL == goldenValues) {
		// if no outputs specified, the current values are written to JSON document
		json_insert_child(outputs, putFloat1D(name, values, size));
		*changed = TRUE;
	}
	else {
		// check whether within tolerance
		float *percents = malloc(size * sizeof(float));
		
		assert(size == goldenSize);
		
		for (i = 0; i < size; ++i) {
			percents[i] = 100.0f * (values[i] - goldenValues[i]) / goldenValues[i];
			if (fabs(percents[i]) > tolerance) {
				success_local = FALSE;
			}
		}
		printf("%s percentage difference:\n", name);
		dump1D(percents, size);
		if (!success_local) {
			printf("FAIL, tolerance is %f\n", tolerance);
			*success = FALSE;
		}
	}
}

/*
 * Compares a test generated float 2D array output with the output value specified in the JSON file.
 * Must match within tolerance.
 * If no output value is specified, this test generated integer output is added to the JSON.
 */
void compareFloat2D(json_t *outputs, char *name, float **values, int rows, int cols, float tolerance, int *changed, int *success) {
	int success_local = TRUE;
	int goldenRows, goldenCols;
	float **goldenValues = getFloat2D(outputs, name, &goldenRows, &goldenCols);
	int i,j;
	
	if (NULL == goldenValues) {
		// if no outputs specified, the current values are written to JSON
		json_insert_child(outputs, putFloat2D(name, values, rows, cols));
		*changed = TRUE;
	}
	else {
		// check whether within tolerance
		float **percents = GCI_ecf_matrix(rows, cols);
		
		assert(rows == goldenRows);
		assert(cols == goldenCols);
		
		for (i = 0; i < rows; ++i) {
			for (j = 0; j < cols; ++j) {
				percents[i][j] = 100.0f * (values[i][j] - goldenValues[i][j]) / goldenValues[i][j];
				if (fabs(percents[i][j]) > tolerance) {
					success_local = FALSE;
				}
			}
		}
		printf("%s percentage difference:\n", name);
		dump2D(percents, rows);
		if (!success_local) {
			printf("FAIL, tolerance is %f\n", tolerance);
			*success = FALSE;
		}
	}
}

/*
 * Custom SLIMCurve fitting tests.
 */
void do_fit(int exponents, json_t *inputs, json_t *outputs, float tolerance, int *changed, int *success)
{
	int old;
	int ndata, size;
	int *counts;
	float x_inc;
	float *y;
	int fit_start, fit_stop;
	float *instr;
	int n_instr;
	int noise;
	float *sig;
	float *fitted;
	float *residuals;
	float chi_square;
	float chi_square_target;
	float a1, a2, a3, t1, t2, t3, z;
	int i;
	int return_value;
	int start;
	int restrain;
	int nparam;
	float *param;
	int *paramfree;
	float **covar;
	float **alpha;
	float **erraxes;
	float chi_square_delta;
	int chi_square_percent;
	
	// if no outputs specified, run old versions to get reference value outputs
	old = 0;
	if (NULL == outputs->child) {
		printf("old\n");
		old = 1;
	}
	
	start = getInt(inputs, "start", 0);
	ndata = getInt(inputs, "ndata", 0);
	counts = getInt1D(inputs, "counts", &size); //TODO ARG s/b int or float here???
	assert(size == ndata);
	y = (float *)malloc(ndata * sizeof(float));
	sig = (float *)malloc(ndata * sizeof(float));
	fitted = (float *)malloc(ndata * sizeof(float));
	residuals = (float *)malloc(ndata * sizeof(float));
	for (i = 0; i < ndata; ++i) {
		y[i] = (float) counts[start + i];
		//printf("y[%d] is %f\n", i, y[i]);
		sig[i] = 1.0; //TODO grab sigma from JSON file
		//printf("sig[%d] is %f\n", i, sig[i]);
	}
	
	instr = getFloat1D(inputs, "instr", &n_instr);

	noise = getInt(inputs, "noise", 0);
	
	a1 = 1.0;
	z = 2.0;
	t1 = 3.0;
	a1 = getFloat(inputs, "a", a1);
	z = getFloat(inputs, "z", z);
	t1 = getFloat(inputs, "t", t1);
	
	x_inc = 0.1f;
	x_inc = getFloat(inputs, "xinc", x_inc);
	fit_start = ndata/4;
	fit_start = getInt(inputs, "fitstart", fit_start);
	fit_stop = 7*ndata/8;					 
	fit_stop = getInt(inputs, "fitstop", fit_stop);
		   
	chi_square_target = 1.0; // make something up for now
	chi_square_target = getFloat(inputs, "chisquaretarget", chi_square_target);
	
	// do a RLD fit		
	if (old) {
		return_value =  Old_GCI_triple_integral_fitting_engine(
															   x_inc,
															   y,
															   fit_start,
															   fit_stop,
															   instr,
															   n_instr,
															   noise,
															   sig,
															   &z,
															   &a1,
															   &t1,
															   fitted,
															   residuals,
															   &chi_square,
															   chi_square_target
															   );		}
	else {
		return_value =  GCI_triple_integral_fitting_engine(
														   x_inc,
														   y,
														   fit_start,
														   fit_stop,
														   instr,
														   n_instr,
														   noise,
														   sig,
														   &z,
														   &a1,
														   &t1,
														   fitted,
														   residuals,
														   &chi_square,
														   chi_square_target
														   );
	}
	
	printf("Triple Integral chisq is %f\n", chi_square);
	printf("Triple Integral return value %d a %f t %f z %f\n", return_value, a1, t1, z);
	
	restrain = 0;
	
	switch (exponents) {
		case 1:
			nparam = 3;
			param = malloc(nparam * sizeof(float));
			param[0] = z;
			param[1] = a1;
			param[2] = t1;
			break;
		case 2:
			nparam = 5;
			param = malloc(nparam * sizeof(float));
			param[0] = z;
			param[1] = a1;
			param[2] = t1;
			param[3] = a1/2;
			param[4] = t1/2;
			break;
		case 3:
			nparam = 7;
			param = malloc(nparam * sizeof(float));
			param[0] = z;
			param[1] = a1;
			param[2] = t1;
			param[3] = a1/2;
			param[4] = t1/2;
			param[5] = a1/4;
			param[6] = t1/4;
			break;
	}
	paramfree = malloc(nparam * sizeof(int));
	for (i = 0; i < nparam; ++i) {
		paramfree[i] = 1;
	}
	
	covar = GCI_ecf_matrix(nparam, nparam);
	alpha = GCI_ecf_matrix(nparam, nparam);
	erraxes = GCI_ecf_matrix(nparam, nparam);

	chi_square_delta = 0.1f;
	chi_square_delta = getFloat(inputs, "chisquaredelta", chi_square_delta);
	chi_square_percent = 500;	
	chi_square_percent = getInt(inputs, "chisquarepercent", chi_square_percent);
	
	// do a LMA fit
	if (old) {
		return_value = Old_GCI_marquardt_fitting_engine(
														x_inc,
														y,
														ndata,
														fit_start,
														fit_stop,
														instr,
														n_instr,
														noise,
														sig,
														param,
														paramfree,
														nparam,
														restrain,
														GCI_multiexp_tau,
														fitted,
														residuals,
														&chi_square,
														covar,
														alpha,
														erraxes,
														chi_square_target,
														chi_square_delta,
														chi_square_percent);
	}
	else {
		return_value = GCI_marquardt_fitting_engine(
													x_inc,
													y,
													ndata,
													fit_start,
													fit_stop,
													instr,
													n_instr,
													noise,
													sig,
													param,
													paramfree,
													nparam,
													restrain,
													GCI_multiexp_tau,
													fitted,
													residuals,
													&chi_square,
													covar,
													alpha,
													erraxes,
													chi_square_target,
													chi_square_delta,
													chi_square_percent);
	}

	printf("Marquardt chisq is %f\n", chi_square);
	
	switch (exponents) {
		case 1:
			z  = param[0];
			a1 = param[1];
			t1 = param[2];

            printf("Marquardt return value %d a %f t %f z %f\n", return_value, a1, t1, z);

			compareFloat(outputs, "a", a1, tolerance, changed, success);
			compareFloat(outputs, "t", t1, tolerance, changed, success);
			compareFloat(outputs, "z", z, tolerance, changed, success);
			break;
		case 2:
			z  = param[0];
			a1 = param[1];
			t1 = param[2];
			a2 = param[3];
			t2 = param[4];

            printf("Marquardt return value %d a1 %f t1 %f a2 %f t2 %f z %f\n", return_value, a1, t1, a2, t2, z);

			compareFloat(outputs, "a1", a1, tolerance, changed, success);
			compareFloat(outputs, "t1", t1, tolerance, changed, success);
			compareFloat(outputs, "a2", a2, tolerance, changed, success);
			compareFloat(outputs, "t2", t2, tolerance, changed, success);
			compareFloat(outputs, "z", z, tolerance, changed, success);
			break;
		case 3:
			z  = param[0];
			a1 = param[1];
			t1 = param[2];
			a2 = param[3];
			t2 = param[4];
			a3 = param[5];
			t3 = param[6];

            printf("Marquardt return value %d a1 %f t1 %f a2 %f t2 %f a3 %f t3 %f z %f\n", return_value, a1, t1, a2, t2, a3, t3, z);

			compareFloat(outputs, "a1", a1, tolerance, changed, success);
			compareFloat(outputs, "t1", t1, tolerance, changed, success);
			compareFloat(outputs, "a2", a2, tolerance, changed, success);
			compareFloat(outputs, "t2", t2, tolerance, changed, success);
			compareFloat(outputs, "a3", a3, tolerance, changed, success);
			compareFloat(outputs, "t3", t3, tolerance, changed, success);
			compareFloat(outputs, "z", z, tolerance, changed, success);
			break;
	}
	
	free(counts);
	free(y);
	free(sig);
	free(fitted);
	free(residuals);
	if (NULL != instr) {
		free(instr);
	}
	free(param);
	free(paramfree);
	GCI_ecf_free_matrix(covar);
	GCI_ecf_free_matrix(alpha);
	GCI_ecf_free_matrix(erraxes);	
}

/*
 * Processes a single test.  This method is customized for particular code being tested.
 */
void do_test(char *test, json_t *inputs, json_t *outputs, float tolerance, int *changed, int *success)
{
	if (strcmp("mono-exp", test) == 0) {
		do_fit(1, inputs, outputs, tolerance, changed, success);
	}
	else if (strcmp("bi-exp", test) == 0) {
		do_fit(2, inputs, outputs, tolerance, changed, success);
	}
	else if (strcmp("tri-exp", test) == 0) {
		do_fit(3, inputs, outputs, tolerance, changed, success);
	}
	else if (strcmp("gci_solve", test) == 0) {
		float **a;
		int n;
		float *b;
		int rows, cols, size;
			
		a = getFloat2D(inputs, "a", &rows, &cols);
		n = getInt(inputs, "n", rows);
		b = getFloat1D(inputs, "b", &size);
			
		printf("a\n");
		dump2D(a, rows);
		printf("b\n");
		dump1D(b, size);
			
		if (n != rows || n != cols || n != size) {
			printf("mismatch n %d rows %d cols %d size %d\n", n, rows, cols, size);
		}
			
		GCI_solve(a, n, b);
			
		compareFloat1D(outputs, "b", b, n, tolerance, changed, success);
			
		// free up heap space
		GCI_ecf_free_matrix(a);
		free(b);
	}
}

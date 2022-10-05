/*
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2022 University of Oxford and Board of Regents of the
 * University of Wisconsin-Madison.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

/* The 2003 version of the ECF library.  This takes account of the fact that we may be
   handling Poisson noise.

   This file contains functions for global analysis, referring to a
   little code from the single transient analysis functions.  Utility
   code is found in EcfUtil.c and single transient analysis code in
   EcfSingle.c.
*/


/**

						   GLOBAL ANALYSIS


   We only work with the case of multiple transients, all with the
   same instrument/prompt response, the same xincr, the same number of
   points, and so on.  The recommended fitting algorithm is a
   three-step process:

   (1) Sum the transients and use this to get initial estimates for
       the global parameters we are estimating.

   (2) Fixing these parameters, perform a Marquardt fit on each of the
       transients.  In our cases, this will be fairly efficient, as we
       will not need to repeatedly calculate the exponential decay.

   (3) Finally, perform the global fit.  There's lots of interesting
       maths here which we'll discuss when we get to it.  Note that
       again we will not need to repeatedly calculate the exponential
       decays for each transient, which will hopefully make the
       process significantly faster.

   We provide special code to perform these steps in the case of
   multiexponential and stretched exponential fits where we are aiming
   to globally fit all of the taus and the h parameter (in the
   stretched exponential case); these can be performed far more
   efficiently than the general case.  We also provide code to perform
   step (3) in the general case; steps (1) and (2) will have to be
   handled on a case-by-case basis by the calling code.  We also
   provide a version of step (3) to handle the case of arbitrary
   x data with no instrument response.

 \file EcfGlobal.c
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EcfInternal.h"
#include "EcfGlobal.h" 

typedef struct {
	float **P, **Q, ***S;
} global_matrix;

typedef struct {
	float *global, *local;
} global_vector;

/* Predeclarations */

int GCI_alloc_global_matrix(global_matrix *m,
							int global, int local, int ntrans);
void GCI_free_global_matrix(global_matrix *m);
void GCI_copy_global_matrix(global_matrix dest, global_matrix src,
							int global, int local, int ntrans);
int GCI_alloc_global_vector(global_vector *v,
							int global, int local, int ntrans);
void GCI_free_global_vector(global_vector *v);
void GCI_copy_global_vector(global_vector dest, global_vector src,
							int global, int local, int ntrans);
int GCI_marquardt_global_exps_est_globals_instr(
			float xincr, float **trans, int ndata, int ntrans,
			int fit_start, int fit_end, float instr[], int ninstr,
			noise_type noise, float sig[], int ftype,
			float **param, int paramfree[], int nparam, float gparam[],
			restrain_type restrain, float chisq_delta,
			float fitted[], float residuals[],
			float **covar, float **alpha, float *chisq_global);
int GCI_marquardt_global_exps_est_params_instr(
			float xincr, float **trans, int ndata, int ntrans,
			int fit_start, int fit_end, float instr[], int ninstr,
			noise_type noise, float sig[], int ftype,
			float **param, int paramfree[], int nparam, restrain_type restrain, float chisq_delta,
			float exp_pure[], float *exp_conv[],
			float **fitted, float **residuals,
			float **covar, float **alpha, float chisq_trans[],
			int drop_bad_transients);
int GCI_marquardt_global_exps_calculate_exps_instr(
			float xincr, int ndata, float instr[], int ninstr,
			int ftype, float param[], int nparam,
			float exp_pure[], float *exp_conv[]);
int GCI_marquardt_global_exps_do_fit_single(
		float xincr, float y[], int ndata, int fit_start, int fit_end,
		noise_type noise, float sig[], int ftype,
		float param[], int paramfree[], int nparam, restrain_type restrain, float chisq_delta,
		float *exp_conv[], float *fitted, float *residuals,
		float **covar, float **alpha, float *chisq_trans);
int GCI_marquardt_global_exps_single_step(
				float xincr, float y[],
				int ndata, int fit_start, int fit_end,
				noise_type noise, float sig[], int ftype,
				float param[], int paramfree[], int nparam,
				restrain_type restrain,
				float *exp_conv[],
				float yfit[], float dy[],
				float **covar, float **alpha, float *chisq,
				float *alambda);
int GCI_marquardt_global_compute_exps_fn(
			float xincr, float y[],
			int ndata, int fit_start, int fit_end,
			noise_type noise, float sig[], int ftype,
			float param[], int paramfree[], int nparam,
			float *exp_conv[],
			float yfit[], float dy[],
			float **alpha, float *beta, float *chisq, float old_chisq);
int GCI_marquardt_global_compute_exps_fn_final(
			float xincr, float y[],
			int ndata, int fit_start, int fit_end,
			noise_type noise, float sig[], int ftype,
			float param[], int paramfree[], int nparam,
			float *exp_conv[],
			float yfit[], float dy[], float *chisq);
int GCI_marquardt_global_exps_do_fit_instr(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[], int ftype,
		float **param, int paramfree[], int nparam, restrain_type restrain, float chisq_delta,
		float exp_pure[], float *exp_conv[],
		float **fitted, float **residuals,
		float **covar_scratch, float **alpha_scratch,
		float *chisq_trans, float *chisq_global,
		int drop_bad_transients);
int GCI_marquardt_global_exps_global_step(
				float xincr, float **trans,
				int ndata, int ntrans, int fit_start, int fit_end,
				float instr[], int ninstr,
				noise_type noise, float sig[], int ftype,
				float **param, int paramfree[], int nparam,
				restrain_type restrain,
				float exp_pure[], float *exp_conv[],
				float **yfit, float **dy,
				float *chisq_trans, float *chisq_global,
				float **alpha_scratch, float *alambda,
				int drop_bad_transients);
int GCI_marquardt_global_compute_global_exps_fn(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[], int ftype,
		float **param, int paramfree[], int nparam,
		int mfit_global, int mfit_local, int gindex[], int lindex[],
		float exp_pure[], float *exp_conv[],
		float **yfit, float **dy, global_matrix alpha, global_vector beta,
		float **alpha_scratch, float *chisq_trans, float *chisq_global,
		int drop_bad_transients);
int GCI_marquardt_global_compute_global_exps_fn_final(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[], int ftype,
		float **param, int paramfree[], int nparam,
		int mfit_global, int mfit_local, int gindex[], int lindex[],
		float exp_pure[], float *exp_conv[],
		float **yfit, float **dy,
		float *chisq_trans, float *chisq_global, int drop_bad_transients);

int GCI_marquardt_global_solve_eqn(global_matrix A, global_vector b,
			int mfit_global, int mfit_local, int ntrans);

int GCI_marquardt_global_generic_do_fit_instr(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[],
		float **param, int paramfree[], int nparam, int gparam[],
		restrain_type restrain, float chisq_delta,
		void (*fitfunc)(float, float [], float *, float [], int),
		float **fitted, float **residuals,
		float **covar_scratch, float **alpha_scratch,
		float *chisq_trans, float *chisq_global);
int GCI_marquardt_global_generic_global_step(
				float xincr, float **trans,
				int ndata, int ntrans, int fit_start, int fit_end,
				float instr[], int ninstr,
				noise_type noise, float sig[],
				float **param, int paramfree[], int nparam, int gparam[],
				restrain_type restrain, float chisq_delta,
				void (*fitfunc)(float, float [], float *, float [], int),
				float **yfit, float **dy,
				float *chisq_trans, float *chisq_global,
				float **alpha_scratch, float *alambda);
int GCI_marquardt_global_compute_global_generic_fn(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[],
		float **param, int paramfree[], int nparam, int gparam[],
		int mfit_global, int mfit_local, int gindex[], int lindex[],
		void (*fitfunc)(float, float [], float *, float [], int),
		float **yfit, float **dy, global_matrix alpha, global_vector beta,
		float **alpha_scratch, float *chisq_trans, float *chisq_global,
		float alambda,
		float **pfnvals, float ***pdy_dparam_pure, float ***pdy_dparam_conv,
		int *pfnvals_len, int *pdy_dparam_nparam_size);
int GCI_marquardt_global_compute_global_generic_fn_final(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[],
		float **param, int paramfree[], int nparam, int gparam[],
		int mfit_global, int mfit_local, int gindex[], int lindex[],
		void (*fitfunc)(float, float [], float *, float [], int),
		float **yfit, float **dy,
		float *chisq_trans, float *chisq_global,
		float **pfnvals, float ***pdy_dparam_pure, float ***pdy_dparam_conv,
		int *pfnvals_len, int *pdy_dparam_nparam_size);

/*         *****             UTILITY CODE             *****        */

/* First, we define functions to handle the data structures we will be
   using later on to store the alpha and covar matrices and the beta
   etc. vectors for global analysis */

int GCI_alloc_global_matrix(global_matrix *m,
							int global, int local, int ntrans)
{
	if (global <= 0 || local < 0 || ntrans <= 0)
		return -2;

	if ((m->P = GCI_ecf_matrix(global, global)) == NULL)
		return -1;
	if (local > 0) {
		if ((m->Q = GCI_ecf_matrix(global, ntrans*local)) == NULL) {
			GCI_ecf_free_matrix(m->P);
			return -1;
		}
		if ((m->S = GCI_ecf_matrix_array(ntrans, local, local)) == NULL) {
			GCI_ecf_free_matrix(m->P);
			GCI_ecf_free_matrix(m->Q);
			return -1;
		}
	} else {
		m->Q = NULL;
		m->S = NULL;
	}

	return 0;
}

void GCI_free_global_matrix(global_matrix *m)
{
	GCI_ecf_free_matrix(m->P);
	if (m->Q != NULL) {
		GCI_ecf_free_matrix(m->Q);
		GCI_ecf_free_matrix_array(m->S);
	}
}

void GCI_copy_global_matrix(global_matrix dest, global_matrix src,
							int global, int local, int ntrans)
{
	int i, j, k;

	for (i=0; i<global; i++)
		for (j=0; j<global; j++)
			dest.P[i][j] = src.P[i][j];

	if (local > 0) {
		for (i=0; i<global; i++)
			for (j=0; j<ntrans*local; j++)
				dest.Q[i][j] = src.Q[i][j];

		for (i=0; i<ntrans; i++)
			for (j=0; j<local; j++)
				for (k=0; k<local; k++)
					dest.S[i][j][k] = src.S[i][j][k];
	}
}


int GCI_alloc_global_vector(global_vector *v,
							int global, int local, int ntrans)
{
	if (global <= 0 || local < 0 || ntrans <= 0)
		return -2;

	if ((v->global = (float *) malloc((unsigned) global * sizeof(float))) == NULL)
		return -1;
	if (local > 0) {
		if ((v->local =
			 (float *) malloc((unsigned)(ntrans * local) * sizeof(float))) == NULL) {
			free(v->global);
			return -1;
		}
	} else {
		v->local = NULL;
	}

	return 0;
}

void GCI_free_global_vector(global_vector *v)
{
	free(v->global);
	if (v->local != NULL)
		free(v->local);
}

void GCI_copy_global_vector(global_vector dest, global_vector src,
							int global, int local, int ntrans)
{
	int i;

	for (i=0; i<global; i++)
		dest.global[i] = src.global[i];

	if (local > 0) {
		for (i=0; i<ntrans*local; i++)
			dest.local[i] = src.local[i];
	}
}


/** EXPONENTIALS CODE 

   Now the code for performing a global fit for multi-exponential taus
   and stretched exponentials.  This is the function which is called
   from external programs. */

/* I don't even want to _contemplate_ error axes for this!  It would
   be computationally very messy, as there would be very large numbers
   of large vectors involved. */

int GCI_marquardt_global_exps_instr(float xincr, float **trans,
					int ndata, int ntrans, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[], int ftype,
					float **param, int paramfree[], int nparam,
					restrain_type restrain, float chisq_delta,
					float **fitted, float **residuals,
					float chisq_trans[], float *chisq_global, int *df,
					int drop_bad_transients)
{
	float **covar, **alpha, *scaled_instr, instrsum;
	int i, j, ret;
	int mlocal, mglobal;
	float gparam[MAXFIT];
	float *exp_pure, *exp_conv[MAXFIT];
//	double time;

//	time=Timer();

	/* Some basic parameter checks */
	if (xincr <= 0) return -1;
	if (ntrans < 1) return -1;
	if (ndata < 1)  return -1;
	if (fit_start < 0 || fit_end > ndata) return -1;
//	if (ninstr < 1) return -1;
	if (nparam < 1) return -1;

	switch (ftype) {
	case FIT_GLOBAL_MULTIEXP:
		if (nparam % 2 != 1) {
			dbgprintf(1, "global fitting: "
					  "multiexp needs odd number of parameters\n");
			return -1;
		}
		break;

	case FIT_GLOBAL_STRETCHEDEXP:
		if (nparam != 4) {
			dbgprintf(1, "global fitting: "
					  "stretched exp needs precisely 4 parameters\n");
			return -1;
		}
		break;

	default:
		dbgprintf(1, "global fitting: unknown fitting type\n");
		return -1;
	}

	if ((covar = GCI_ecf_matrix(nparam, nparam)) == NULL)
		return -2;

	if ((alpha = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		GCI_ecf_free_matrix(covar);
		return -3;
	}

	if ((scaled_instr = (float *) malloc((unsigned) ninstr * sizeof(float))) == NULL && ninstr>1) {
		GCI_ecf_free_matrix(covar);
		GCI_ecf_free_matrix(alpha);
		return -4;
	}

	/* Also allocate space for the exp_pure and exp_conv arrays */
	if ((exp_conv[0] = (float *) malloc((unsigned)(nparam * ndata) * sizeof(float)))
		== NULL) {
		GCI_ecf_free_matrix(covar);
		GCI_ecf_free_matrix(alpha);
		free(scaled_instr);
		return -5;
	}

	exp_pure = exp_conv[0];
	for (i=1; i<nparam; i++)
		exp_conv[i] = exp_conv[0] + i * ndata;

	/* Scale the instrument response */
	for (i=0, instrsum=0; i<ninstr; i++)
		instrsum += instr[i];
	if (instrsum == 0) instrsum=1.0;  //return -6;
	for (i=0; i<ninstr; i++)
		scaled_instr[i] = instr[i] / instrsum;

//	printf("that took: %f secs.\n", Timer()-time); time=Timer();
	dbgprintf(2, "About to enter step (1)\n");

	/* Step (1): estimate the global taus */
	ret = GCI_marquardt_global_exps_est_globals_instr(
			xincr, trans, ndata, ntrans, fit_start, fit_end,
			instr != NULL ? scaled_instr : NULL, ninstr, noise, sig,
			ftype, param, paramfree, nparam, gparam, restrain, chisq_delta,
			fitted[0], residuals[0], covar, alpha, chisq_global);

	if (ret != 0) {
		dbgprintf(1, "Step (1) failed, ret = %d\n", ret);
		GCI_ecf_free_matrix(covar);
		GCI_ecf_free_matrix(alpha);
		free(scaled_instr);
		free(exp_conv[0]);
		return -10 + ret;
	}

//	printf("that took: %f secs.\n", Timer()-time); time=Timer();
	/* Copy the estimated global taus to the parameters array */

	switch (ftype)
	{
		case FIT_GLOBAL_MULTIEXP:
			for (i=2; i<nparam; i+=2)
				for (j=0; j<ntrans; j++)
					param[j][i] = gparam[i];
		break;

		case FIT_GLOBAL_STRETCHEDEXP:
			for (i=2; i<nparam; i++)  /* param 2 = tau, param 3 = h */
				for (j=0; j<ntrans; j++)
					param[j][i] = gparam[i];
		break;

		default:
			dbgprintf(1, "global_exps_instr: please update me!\n");
			GCI_ecf_free_matrix(covar);
			GCI_ecf_free_matrix(alpha);
			free(scaled_instr);
			free(exp_conv[0]);
			return -1;
	}

//	printf("that took: %f secs.\n", Timer()-time); time=Timer();
	dbgprintf(2, "About to enter step (2)\n");

	/* Step (2): use these taus to estimate initial values for all of
	   the individual transient parameters */
	ret = GCI_marquardt_global_exps_est_params_instr(
			xincr, trans, ndata, ntrans, fit_start, fit_end,
			instr != NULL ? scaled_instr : NULL, ninstr, noise, sig, ftype,
			param, paramfree, nparam, restrain, chisq_delta,
			exp_pure, exp_conv, fitted, residuals, covar, alpha,
			chisq_trans, drop_bad_transients);

	if (ret != 0) {
		dbgprintf(1, "Step (2) failed, ret = %d\n", ret);
		GCI_ecf_free_matrix(covar);
		GCI_ecf_free_matrix(alpha);
		free(scaled_instr);
		free(exp_conv[0]);
		return -20 + ret;
	}

//	printf("that took: %f secs.\n", Timer()-time); time=Timer();
	dbgprintf(2, "About to enter step (3)\n");

	/* Step (3): now that we have estimates for initial values for all
	   parameters, we can do the global Marquardt fitting.  Note that
	   covar and alpha are only provided as scratch space. */
	ret = GCI_marquardt_global_exps_do_fit_instr(
			xincr, trans, ndata, ntrans, fit_start, fit_end,
			instr != NULL ? scaled_instr : NULL, ninstr, noise, sig, ftype,
			param, paramfree, nparam, restrain, chisq_delta,
			exp_pure, exp_conv, fitted, residuals, covar, alpha,
			chisq_trans, chisq_global, drop_bad_transients);

	GCI_ecf_free_matrix(covar);
	GCI_ecf_free_matrix(alpha);
	free(scaled_instr);
	free(exp_conv[0]);

	if (ret < 0) {
		dbgprintf(1, "Step (3) failed, ret = %d\n", ret);
		return -30 + ret;
	}

//	printf("that took: %f secs.\n", Timer()-time); time=Timer();
	dbgprintf(2, "Step (3) succeeded, ret = %d\n", ret);

	/* Before we return, calculate the number of degrees of freedom */
	/* The number of degrees of freedom is given by:
	      d.f. = ntrans * ((fit_end - fit_start) - # free local parameters)
	                 - # free global parameters
	*/

	if (drop_bad_transients) {
		*df = 0;
		for (i=0; i<ntrans; i++) {
			if (chisq_trans[i] > 0)
				(*df)++;
		}
	} else
		*df = ntrans;

	mglobal = mlocal = 0;

	switch (ftype)
	{
	case FIT_GLOBAL_MULTIEXP:
		for (i=2; i<nparam; i+=2)
		if (paramfree[i]) mglobal++;
		for (i=1; i<nparam; i+=2)
			if (paramfree[i]) mlocal++;
		if (paramfree[0]) mlocal++;
		break;

	case FIT_GLOBAL_STRETCHEDEXP:
		if (paramfree[0]) mlocal++;
		if (paramfree[1]) mlocal++;
		if (paramfree[2]) mglobal++;
		if (paramfree[3]) mglobal++;
		break;

	default:
		dbgprintf(1, "global_exps_instr: please update me!\n");
		return -1;
	}

	*df *= ((fit_end - fit_start) - mlocal);
	*df -= mglobal;

//	printf("mlocal %d\nmglobal %d\ndf %d\n", mlocal, mglobal, *df);

//	printf("that took: %f secs.\n", Timer()-time); time=Timer();
	return ret;
}


int GCI_marquardt_global_exps_est_globals_instr(
			float xincr, float **trans, int ndata, int ntrans,
			int fit_start, int fit_end, float instr[], int ninstr,
			noise_type noise, float sig[], int ftype,
			float **param, int paramfree[], int nparam, float gparam[],
			restrain_type restrain, float chisq_delta,
			float fitted[], float residuals[],
			float **covar, float **alpha, float *chisq_global)
{
	int i, j, ret, nparamfree;
	float *summed, *tptr;
	int data_start;
	float Z, A, tau;
	void (*fitfunc)(float, float [], float *, float [], int);

	if ((summed = (float *) calloc((unsigned) ndata, sizeof(float))) == NULL)
		return -1;

	for (i=0; i<ntrans; i++) {
		tptr = trans[i];
		for (j=0; j<ndata; j++)
			summed[j] += tptr[j];
	}

	/* This code is now lifted from fitting.c, appropriately modified */
	data_start = fit_start + ECF_Find_Float_Max(&summed[fit_start],
												fit_end - fit_start, &A);
//	ret = GCI_triple_integral_instr(xincr, summed, data_start, fit_end,
//									instr, ninstr, noise, sig,
//									&Z, &A, &tau, NULL, NULL, NULL);

	ret = GCI_triple_integral_fitting_engine(xincr, summed, data_start, fit_end,
							 				instr, ninstr, noise, sig,
							  				&Z, &A, &tau, NULL, NULL, NULL, 1.5f*(float)(fit_end-fit_start-3));

	dbgprintf(3, "In est_globals_instr, triple integral ret = %d\n", ret);

	if (ret < 0) {
		Z = 0;
		ECF_Find_Float_Max(&summed[fit_start], fit_end - fit_start, &A);
		tau = 2.0;
	}


	/* We set gparam[] to be an array which holds initial estimates of
	   the parameters for the _sum_ of all the transients, for those
	   parameters which are fixed and therefore "known"; the rest will
	   be estimated later.  It doesn't matter if we set a few other
	   values, as these will be overwritten later.  We could also
	   merge this switch() with the next one, but then the code would
	   possibly be a little less easy to follow, so we won't. */

	switch (ftype) {
	case FIT_GLOBAL_MULTIEXP:
		/* Z */
		if (! paramfree[0])
		{
			gparam[0] = 0;
			for (j=0; j<ntrans; j++) gparam[0] += param[j][0];
		}

		/* A's */
		for (i=1; i<nparam; i++)
			if (! paramfree[i])
			{
				gparam[i] = 0;
				for (j=0; j<ntrans; j++) gparam[i] += param[j][i];
			}

		/* taus last (was first) */
		for (i=2; i<nparam; i+=2)
			gparam[i] = param[0][i];

		break;

	case FIT_GLOBAL_STRETCHEDEXP:
		/* Z */
		if (! paramfree[0]) {
			gparam[0] = 0;
			for (j=0; j<ntrans; j++) gparam[0] += param[j][0];
		}

		/* A */
		if (! paramfree[1]) {
			gparam[1] = 0;
			for (j=0; j<ntrans; j++) gparam[1] += param[j][1];
		}

		/* tau and h last (were first) */
		for (i=2; i<nparam; i++)
			gparam[i] = param[0][i];

		break;

	default:
		dbgprintf(1, "global_exps_est_globals_instr: please update me!\n");
		free(summed);
		return -1;
	}


	/* Now we can set any non-fixed parameters to meaningful initial
	   estimates */
	switch (ftype) {
	case FIT_GLOBAL_MULTIEXP:
		fitfunc = GCI_multiexp_tau;

		switch (nparam) {
		case 3:
			if (paramfree[0]) gparam[0] = Z;
			if (paramfree[1]) gparam[1] = A;
			if (paramfree[2]) gparam[2] = tau;
			break;

		case 5:
			if (paramfree[0]) gparam[0] = Z;
			if (paramfree[1]) gparam[1] = A*3/4;
			if (paramfree[2]) gparam[2] = tau;
			if (paramfree[3]) gparam[3] = A*1/4;
			if (paramfree[4]) gparam[4] = tau*2/3;
			break;

		default:
			if (nparam<7) {
				free(summed);
				return -2;
			}
			if (paramfree[0]) gparam[0] = Z;
			if (paramfree[1]) gparam[1] = A*3/4;
			if (paramfree[2]) gparam[2] = tau;
			if (paramfree[3]) gparam[3] = A*1/6;
			if (paramfree[4]) gparam[4] = tau*2/3;
			if (paramfree[5]) gparam[5] = A*1/6;
			if (paramfree[6]) gparam[6] = tau/3;
			for (i=7; i<nparam; i+=2) {
				if (paramfree[i]) gparam[i] = A/(float)i;
				if (paramfree[i+1]) gparam[i+1] = tau/(float)i;
			}
			break;
		}
		break;

	case FIT_GLOBAL_STRETCHEDEXP:
		fitfunc = GCI_stretchedexp;

		if (paramfree[0]) gparam[0] = Z;
		if (paramfree[1]) gparam[1] = A;
		if (paramfree[2]) gparam[2] = tau;
		if (paramfree[3]) gparam[3] = 1.5;  /* h */
		break;

	default:
		dbgprintf(1, "est_globals_instr: please update me!\n");
		free(summed);
		return -1;
	}

	for (i=0, nparamfree=0; i<nparam; i++) if (paramfree[i]) nparamfree++;

	/* Note that the only values in the gparam array which are of
	   interest are the taus and h for stretched exp, so we don't need
	   to rescale Z and the A's back again */
//	ret = GCI_marquardt_instr(xincr, summed, ndata, fit_start, fit_end,
//							  instr, ninstr, noise, sig,
//							  gparam, paramfree, nparam, restrain, fitfunc,
//							  fitted, residuals, covar, alpha, chisq_global,
//							  0, NULL);

	ret = GCI_marquardt_fitting_engine(xincr, summed, ndata, fit_start, fit_end,
							  instr, ninstr, noise, sig,
							  gparam, paramfree, nparam, restrain, fitfunc,
							  fitted, residuals, chisq_global, covar, alpha,
							  NULL, 1.5f*(float)(fit_end-fit_start-nparamfree), chisq_delta, 0);

	dbgprintf(3, "In est_globals_instr, marquardt ret = %d\n", ret);

	free(summed);

	if (ret < 0)
		return -3;
	else
		return 0;
}


int GCI_marquardt_global_exps_est_params_instr(
			float xincr, float **trans, int ndata, int ntrans,
			int fit_start, int fit_end, float instr[], int ninstr,
			noise_type noise, float sig[], int ftype,
			float **param, int paramfree[], int nparam, restrain_type restrain, float chisq_delta,
			float exp_pure[], float *exp_conv[],
			float **fitted, float **residuals,
			float **covar, float **alpha, float chisq_trans[],
			int drop_bad_transients)
{
	int i, j, sortkey[MAXFIT], paramfree_local[MAXFIT], tempi, ret;
	int data_start;
	float Z, A, tau;

	/* We begin by estimating the non-tau parameters for each
	   transient.  We do this by performing a triple-integral fit and
	   using the Z and A resulting as a basis for our initial
	   parameters.  In the case of multiexponential fitting, we also
	   sort the taus into decreasing order, and assume that the
	   largest tau is the most significant component. */

	// PRB 03/07 Although **fitted and **residuals are provided only one "transient" is required and used, fitted[0] and residuals[0]

	if (ftype ==  FIT_GLOBAL_MULTIEXP) {
		/* Initialise */
		for (i=2; i<nparam; i+=2)
			sortkey[i] = i;

		/* Bubblesort :-) */
		for (i=2; i<nparam; i+=2)
			for (j=2*(nparam/2); j>i; j-=2)  /* nparam is odd */
				if (param[0][sortkey[j]] > param[0][sortkey[j-2]]) {
					tempi = sortkey[j];
					sortkey[j] = sortkey[j-2];
					sortkey[j-2] = tempi;
				}
	}

	dbgprintf(3, "In est_params_instr, parameters initialised to:\n");

	for (i=0; i<ntrans; i++) {
		/* This code is now lifted from fitting.c, appropriately modified */
		data_start = fit_start + ECF_Find_Float_Max(&trans[i][fit_start],
													fit_end - fit_start, &A);
//		ret = GCI_triple_integral_instr(xincr, trans[i], data_start, fit_end,
//										instr, ninstr, noise, sig,
//										&Z, &A, &tau, NULL, NULL, NULL);

		ret = GCI_triple_integral_fitting_engine(xincr, trans[i], data_start, fit_end,
							 				instr, ninstr, noise, sig,
							  				&Z, &A, &tau, NULL, NULL, NULL, 1.5f*(float)(fit_end-fit_start-3));
		if (ret < 0) {
			Z = 0;
			ECF_Find_Float_Max(&trans[i][fit_start], fit_end - fit_start, &A);
		}

		switch (ftype) {
		case FIT_GLOBAL_MULTIEXP:
			switch (nparam) {
			case 3:
				if (paramfree[0]) param[i][0] = Z;
				if (paramfree[1]) param[i][1] = A;
				break;

			case 5:
				if (paramfree[0]) param[i][0] = Z;
				if (paramfree[sortkey[2]-1]) param[i][sortkey[2]-1] = A*3/4;
				if (paramfree[sortkey[4]-1]) param[i][sortkey[4]-1] = A*1/4;
				break;

			default:
				if (nparam<7) { /* only actually need to do this once */
					return -1;
				}
				if (paramfree[0]) param[i][0] = Z;
				if (paramfree[sortkey[2]-1]) param[i][sortkey[2]-1] = A*3/4;
				if (paramfree[sortkey[4]-1]) param[i][sortkey[4]-1] = A*1/6;
				if (paramfree[sortkey[6]-1]) param[i][sortkey[6]-1] = A*1/6;
				/* this is all pretty meaningless from here on, anyway */
				for (j=8; j<nparam; j+=2) {
					if (paramfree[sortkey[j]-1]) param[i][sortkey[j]-1] = A/(float)j;
				}
				break;
			}
			break;

		case FIT_GLOBAL_STRETCHEDEXP:
			if (paramfree[0]) param[i][0] = Z;
			if (paramfree[1]) param[i][1] = A;
			break;

		default:
			dbgprintf(1, "est_params_instr: please update me!\n");
			return -1;
		}

		if (ECF_debug >= 3) {
			for (j=0; j<nparam; j++)
				dbgprintf(3, "param[%d][%d] = %.4g\n", i, j, param[i][j]);
		}
	}


	/* OK, the initial parameters are set up, we now do a Marquardt
	   fit on all of the transients simultaneously to get more decent
	   initial A and Z values.  But we do this manually without
	   recalculating the exponentials repeatedly.  Furthermore, since
	   the instrument response is convolved linearly with the
	   exponentials, we can get by doing the convolution once only as
	   well, making further major time savings. */

	if (GCI_marquardt_global_exps_calculate_exps_instr(
			xincr, ndata, instr, ninstr, ftype, param[0], nparam,
			exp_pure, exp_conv) != 0)
		return -2;

	/* Now we can do a Marquardt fit on each of the transients */

	/* Create a paramfree[] array which fixes all of the taus */
	switch (ftype) {
	case FIT_GLOBAL_MULTIEXP:
		paramfree_local[0] = paramfree[0];  /* Z */
		for (i=1; i<nparam; i+=2)
			paramfree_local[i] = paramfree[i];  /* the A's */
		for (i=2; i<nparam; i+=2)
			paramfree_local[i] = 0;   /* the taus */
		break;

	case FIT_GLOBAL_STRETCHEDEXP:
		paramfree_local[0] = paramfree[0];  /* Z */
		paramfree_local[1] = paramfree[1];  /* A */
		paramfree_local[2] = 0;  /* tau */
		paramfree_local[3] = 0;  /* h */
		break;

	default:
		dbgprintf(1, "est_params_instr: please update me!\n");
		return -1;
	}

	dbgprintf(3, "In est_params_instr, after do_fit_single, "
			  "parameters initialised to:\n");
	for (i=0; i<ntrans; i++) {
		ret = GCI_marquardt_global_exps_do_fit_single(
					xincr, trans[i], ndata, fit_start, fit_end,
					noise, sig, ftype,
					param[i], paramfree_local, nparam, restrain, chisq_delta, exp_conv,
					fitted[0], residuals[0], covar, alpha, &chisq_trans[i]);

		if (ret < 0) {
			if (drop_bad_transients) {
				dbgprintf(2, "In est_params_instr, transient %d gave "
						  "do_fit_single return val %d; dropping it\n",
						  i, ret);
				chisq_trans[i] = -1;
				continue;
			} else {
				dbgprintf(1, "In est_params_instr, transient %d gave "
						  "do_fit_single return val %d\n", i, ret);
				return -10 + ret;
			}
		}

		/* We try a second time with these parameters if we got
		   nonsense */
		if (chisq_trans[i] > 20 * (fit_end - fit_start)) {
			ret = GCI_marquardt_global_exps_do_fit_single(
					xincr, trans[i], ndata, fit_start, fit_end,
					noise, sig, ftype,
					param[i], paramfree_local, nparam, restrain, chisq_delta, exp_conv,
					fitted[0], residuals[0], covar, alpha, &chisq_trans[i]);

			/* Improved? */
			if (ret < 0 || chisq_trans[i] > 20 * (fit_end - fit_start)) {
				if (drop_bad_transients) {
					dbgprintf(2, "In est_params_instr, transient %d gave "
							  "do_fit_single return val %d, chisq value %.3f; "
							  "dropping it\n", i, ret, chisq_trans[i]);
					chisq_trans[i] = -1;
					continue;
				} else {
					dbgprintf(1, "In est_params_instr, transient %d gave "
							  "do_fit_single return val %d, "
							  "chisq value %.3f\n", i, ret, chisq_trans[i]);
					return -10 + ret;
				}
			}

			if (ECF_debug >= 3) {
				for (j=0; j<nparam; j++)
					dbgprintf(3, "param[%d][%d] = %.4g\n", i, j, param[i][j]);
			}
		}
	}

	return 0;
}


/** This finds values of exp(-x/tau) and x*exp(-x/tau)/tau^2, which are
   needed later for the multiexponential case, and finds the
   equivalents in the stretched exponential case. */
// should also now handle no instrument response, i.e. instr=NULL

int GCI_marquardt_global_exps_calculate_exps_instr(
			float xincr, int ndata, float instr[], int ninstr,
			int ftype, float param[], int nparam,
			float exp_pure[], float *exp_conv[])
{
	int i, j, k;
	int convpts;
	double excur;  /* exp(-x/tau) */
	double exincr; /* exp(-xincr/tau) */
	float *expi;
	float ex, lxa, xah, a2inv, a3inv;  /* for stetched exp */
	double xa, xaincr;

	switch (ftype) {
	case FIT_GLOBAL_MULTIEXP:
		/* Not quite the most efficient way to do this, but not bad */
		/* First we calculate the exp(-x/tau) array */
		for (i=2; i<nparam; i+=2) {
			expi = exp_conv[i];

			excur = 1.0;
			exincr = exp(-xincr/(double)param[i]);
			for (j=0; j<ndata; j++) {
				exp_pure[j] = (float) excur;
				excur *= exincr;

				/* And convolve the exponentials with the instrument response if possible */
				expi[j] = 0;
				convpts = (ninstr <= j) ? ninstr-1 : j;
				if (convpts<=0 || instr==NULL) expi[j] = exp_pure[j];
				else for (k=0; k<=convpts; k++) expi[j] += exp_pure[j-k]*instr[k];
			}
		}

		/* Now we repeat the exercise for x*exp(-x/tau) / tau^2 */
		for (i=2; i<nparam; i+=2) {
			expi = exp_conv[i-1];

			excur = 1.0 / (param[i]*param[i]);  /* 1/tau^2 */
			exincr = exp(-xincr/(double)param[i]);
			for (j=0; j<ndata; j++) {
				exp_pure[j] = (xincr*(float)i) * (float)excur; /* x*exp(-x/tau) / tau^2 */
				excur *= exincr;

				/* And convolve the exponentials with the instrument response if possible */
				expi[j] = 0;
				convpts = (ninstr <= j) ? ninstr-1 : j;
				if (convpts<=0 || instr==NULL) expi[j] = exp_pure[j];
				else for (k=0; k<=convpts; k++) expi[j] += exp_pure[j-k]*instr[k];
			}
		}
		break;

	case FIT_GLOBAL_STRETCHEDEXP:
		/* We start by essentially repeating the stretchedexp_array
		   function */
		xa=0;
		xaincr = xincr / param[2];
		a2inv = 1/param[2];
		a3inv = 1/param[3];

		/* When x=0 */
		exp_conv[1][0] = 1.0;
		exp_conv[2][0] = exp_conv[3][0] = 0.0;

		for (i=1; i<ndata; i++) {
			xa += xaincr;       /* xa = (xincr*i)/param[2] */
			lxa = (float)log(xa);      /* lxa = log(x/param[2]) */
			xah = expf(lxa * a3inv);  /* xah = exp(log(x/param[2])/param[3])
		                                    = (x/param[2])^(1/param[3]) */
			exp_conv[1][i] = ex = expf(-xah);
			                    /* ex = exp(-(x/param[2])^(1/param[3])) */
			ex *= xah * a3inv;  /* ex = exp(...) * (x/param[2])^(1/param[3]) *
			                              1/param[3] */
			exp_conv[2][i] = ex * a2inv;
			exp_conv[3][i] = ex * lxa * a3inv;
		}

		if (ninstr>0 && instr!=NULL)   // else exp_conv already contains the pure data
		{
			/* Now convolve these with the instrument response */
			for (i=1; i<4; i++)
			{
				expi = exp_conv[i];

				for (j=0; j<ndata; j++)
					exp_pure[j] = expi[j];  /* save the array in temp storage */

				for (j=0; j<ndata; j++)
				{
					expi[j] = 0;
					convpts = (ninstr <= j) ? ninstr-1 : j;
					for (k=0; k<=convpts; k++)
						expi[j] += exp_pure[j-k]*instr[k];
				}
			}
		}
		break;

	default:
		dbgprintf(1, "calculate_exps_instr: please update me!\n");
		return -1;
	}

	return 0;
}


/** This is just like the normal GCI_marquardt_instr, except that it
   is designed for the multiexp case where we provide the exponential
   decays in advance, and where we don't care about error axes */
int GCI_marquardt_global_exps_do_fit_single(
		float xincr, float y[], int ndata, int fit_start, int fit_end,
		noise_type noise, float sig[], int ftype,
		float param[], int paramfree[], int nparam, restrain_type restrain,
		float chisq_delta, float *exp_conv[], float *fitted, float *residuals,
		float **covar, float **alpha, float *chisq)
{
	float alambda, ochisq;
	int mfit;
	int i, k, itst, itst_max;

	itst_max = (restrain == ECF_RESTRAIN_DEFAULT) ? 4 : 6;

	mfit = 0;
	for (i=0; i<nparam; i++) {
		if (paramfree[i])
			mfit++;
	}

	alambda = -1;
	if (GCI_marquardt_global_exps_single_step(
				xincr, y, ndata, fit_start, fit_end,
				noise, sig, ftype, param, paramfree, nparam, restrain,
				exp_conv, fitted, residuals, covar, alpha,
				chisq, &alambda) != 0) {
		return -1;
	}

	k = 1;  /* Iteration counter */
	itst = 0;
	for (;;) {
		k++;
		if (k > MAXITERS) {
			return -2;
		}

		ochisq = *chisq;
		if (GCI_marquardt_global_exps_single_step(
					xincr, y, ndata, fit_start, fit_end,
					noise, sig, ftype, param, paramfree, nparam, restrain,
					exp_conv, fitted, residuals, covar, alpha,
					chisq, &alambda) != 0) {
			return -3;
		}

		if (*chisq > ochisq)
			itst = 0;
		else if (ochisq - *chisq < chisq_delta)
			itst++;

		if (itst < itst_max) continue;

		/* Endgame; this also handles correcting the chi-squared values */
		alambda=0.0;
		if (GCI_marquardt_global_exps_single_step(
					xincr, y, ndata, fit_start, fit_end,
					noise, sig, ftype, param, paramfree, nparam, restrain,
					exp_conv, fitted, residuals, covar, alpha,
					chisq, &alambda) != 0) {
			return -4;
		}

		return k;  /* We're done now */
	}
}


/** And this one is basically a specialised GCI_marquardt_instr_step */
int GCI_marquardt_global_exps_single_step(
				float xincr, float y[],
				int ndata, int fit_start, int fit_end,
				noise_type noise, float sig[], int ftype,
				float param[], int paramfree[], int nparam,
				restrain_type restrain,
				float *exp_conv[],
				float yfit[], float dy[],
				float **covar, float **alpha, float *chisq,
				float *alambda)
{
	int j, k, l, ret;
	static int mfit;
	static float ochisq, paramtry[MAXFIT], beta[MAXFIT], dparam[MAXFIT];
	static void (*fitfunc)(float, float [], float *, float [], int);

	if (nparam > MAXFIT)
		return -10;
	if (xincr <= 0)
		return -11;
	if (fit_start < 0 || fit_start > fit_end || fit_end > ndata)
		return -12;

	/* Initialisation */
	/* We assume we're given sensible starting values for param[] */
	if (*alambda < 0.0) {
		switch (ftype) {
		case FIT_GLOBAL_MULTIEXP:
			fitfunc = GCI_multiexp_tau;
			break;

		case FIT_GLOBAL_STRETCHEDEXP:
			fitfunc = GCI_stretchedexp;
			break;

		default:
			dbgprintf(1, "exps_single_step: please update me!\n");
			return -1;
		}

		for (mfit=0, j=0; j<nparam; j++)
			if (paramfree[j])
				mfit++;

		if (GCI_marquardt_global_compute_exps_fn(
					xincr, y, ndata, fit_start, fit_end, noise, sig,
					ftype, param, paramfree, nparam, exp_conv,
					yfit, dy, alpha, beta, chisq, 0.0) != 0)
			return -2;

		*alambda = 0.001f;
		ochisq = *chisq;
		for (j=0; j<nparam; j++)
			paramtry[j] = param[j];
	}

	/* Alter linearised fitting matrix by augmenting diagonal elements */
	for (j=0; j<mfit; j++) {
		for (k=0; k<mfit; k++)
			covar[j][k] = alpha[j][k];
		covar[j][j] = alpha[j][j] * (1.0f + (*alambda));
		dparam[j] = beta[j];
	}

	/* Matrix solution; GCI_solve solves Ax=b rather than AX=B */
	if (GCI_solve(covar, mfit, dparam) != 0)
		return -1;

	/* Once converged, calculate corrected chi-squared values */
	if (*alambda == 0) {
		if (GCI_marquardt_global_compute_exps_fn_final(
					xincr, y, ndata, fit_start, fit_end, noise, sig,
					ftype, param, paramfree, nparam, exp_conv,
					yfit, dy, chisq) != 0)
			return -4;
		/* Irrelevant */
		// if (mfit < nparam) {  /* no need to do this otherwise */
		// 	GCI_covar_sort(covar, nparam, paramfree, mfit);
		// 	GCI_covar_sort(alpha, nparam, paramfree, mfit);
		// }
		return 0;
	}

	/* Did the trial succeed? */
	for (j=0, l=0; l<nparam; l++)
		if (paramfree[l])
			paramtry[l] = param[l] + dparam[j++];

	if (restrain == ECF_RESTRAIN_DEFAULT)
		ret = check_ecf_params (paramtry, nparam, fitfunc);
	else
		ret = check_ecf_user_params (paramtry, nparam, fitfunc);

	if (ret != 0) {
		/* Bad parameters, increase alambda and return */
		*alambda *= 10.0f;
		return 0;
	}

	if (GCI_marquardt_global_compute_exps_fn(
				xincr, y, ndata, fit_start, fit_end, noise, sig,
				ftype, paramtry, paramfree, nparam, exp_conv,
				yfit, dy, covar, dparam, chisq, ochisq) != 0)
		return -2;

	/* Success, accept the new solution */
	if (*chisq < ochisq) {
		*alambda *= 0.1f;
		ochisq = *chisq;
		for (j=0; j<mfit; j++) {
			for (k=0; k<mfit; k++)
				alpha[j][k] = covar[j][k];
			beta[j] = dparam[j];
		}
		for (l=0; l<nparam; l++)
			param[l] = paramtry[l];
	} else { /* Failure, increase alambda and return */
		*alambda *= 10.0f;
		*chisq = ochisq;
	}

	return 0;
}


/** This is a streamlined GCI_marquardt_compute_fn_instr */
int GCI_marquardt_global_compute_exps_fn(
			float xincr, float y[],
			int ndata, int fit_start, int fit_end,
			noise_type noise, float sig[], int ftype,
			float param[], int paramfree[], int nparam,
			float *exp_conv[],
			float yfit[], float dy[],
			float **alpha, float *beta,
			float *chisq, float old_chisq)
{
	int i, j, k, mfit;
	float dy_dparam[MAXBINS][MAXFIT];
	float alpha_weight[MAXBINS];
	float beta_weight[MAXBINS];
	float weight;
	int i_free;
	int j_free;
	float dot_product;
	float beta_sum;
	float dy_dparam_k_i;

	for (j=0, mfit=0; j<nparam; j++)
		if (paramfree[j])
			mfit++;

	*chisq = 0.0f;

	switch (ftype) {
		case FIT_GLOBAL_MULTIEXP:
			switch (noise) {
				case NOISE_CONST:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// multi-exponential fit
						yfit[i] = param[0];  /* Z */
						dy_dparam[i][0] = 1.0f;

						for (j=1; j<nparam; j+=2) {
							yfit[i] += param[j] * exp_conv[j+1][i];
							/* A_j . exp(-x/tau_j) */
							dy_dparam[i][j] = exp_conv[j+1][i];
							/* exp(-x/tau_j) */
							dy_dparam[i][j+1] = param[j] * exp_conv[j][i];
							/* A_j * x * exp(-x/tau_j) / tau_j^2 */
						}
						dy[i] = y[i] - yfit[i];

						// constant noise
						weight = 1.0f / sig[0];
						alpha_weight[i] = weight; // 1 / (sig[0] * sig[0]);
						weight *= dy[i];
						beta_weight[i] = weight; // dy[i] / (sig[0] * sig[0]);
						weight *= dy[i];
						*chisq += weight; // (dy[i] * dy[i]) / (sig[0] * sig[0]);
					}
					break;
				case NOISE_GIVEN:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// multi-exponential fit
						yfit[i] = param[0];  /* Z */
						dy_dparam[i][0] = 1.0f;

						for (j=1; j<nparam; j+=2) {
							yfit[i] += param[j] * exp_conv[j+1][i];
							/* A_j . exp(-x/tau_j) */
							dy_dparam[i][j] = exp_conv[j+1][i];
							/* exp(-x/tau_j) */
							dy_dparam[i][j+1] = param[j] * exp_conv[j][i];
							/* A_j * x * exp(-x/tau_j) / tau_j^2 */
						}
						dy[i] = y[i] - yfit[i];

						// given noise
						weight = 1.0f / (sig[i] * sig[i]);
						alpha_weight[i] = weight; // 1 / (sig[i] * sig[i])
						weight *= dy[i];
						beta_weight[i] = weight; // dy[i] / (sig[i] * sig[i])
						weight *= dy[i];
						*chisq += weight; // (dy[i] * dy[i]) / (sig[i] * sig[i])
					}
					break;
				case NOISE_POISSON_DATA:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// multi-exponential fit
						yfit[i] = param[0];  /* Z */
						dy_dparam[i][0] = 1.0f;

						for (j=1; j<nparam; j+=2) {
							yfit[i] += param[j] * exp_conv[j+1][i];
							/* A_j . exp(-x/tau_j) */
							dy_dparam[i][j] = exp_conv[j+1][i];
							/* exp(-x/tau_j) */
							dy_dparam[i][j+1] = param[j] * exp_conv[j][i];
							/* A_j * x * exp(-x/tau_j) / tau_j^2 */
						}
						dy[i] = y[i] - yfit[i];

						// poisson noise based on data
						weight = (y[i] > 15 ? 1.0f / y[i] : 1.0f / 15);
						alpha_weight[i] = weight; // 1 / sig(i)
						weight *= dy[i];
						beta_weight[i] = weight; // dy[i] / sig(i)
						weight *= dy[i];
						*chisq += weight; // (dy[i] * dy[i]) / sig(i)
					}
					break;
				case NOISE_POISSON_FIT:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// multi-exponential fit
						yfit[i] = param[0];  /* Z */
						dy_dparam[i][0] = 1.0f;

						for (j=1; j<nparam; j+=2) {
							yfit[i] += param[j] * exp_conv[j+1][i];
							/* A_j . exp(-x/tau_j) */
							dy_dparam[i][j] = exp_conv[j+1][i];
							/* exp(-x/tau_j) */
							dy_dparam[i][j+1] = param[j] * exp_conv[j][i];
							/* A_j * x * exp(-x/tau_j) / tau_j^2 */
						}
						dy[i] = y[i] - yfit[i];

						// poisson noise based on fit
						weight = (yfit[i] > 15 ? 1.0f / yfit[i] : 1.0f / 15);
						alpha_weight[i] = weight; // 1 / sig(i)
						weight *= dy[i];
						beta_weight[i] = weight; // dy(i) / sig(i)
						weight *= dy[i];
						*chisq += weight; // (dy(i) * dy(i)) / sig(i)
					}
					break;
				case NOISE_GAUSSIAN_FIT:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// multi-exponential fit
						yfit[i] = param[0];  /* Z */
						dy_dparam[i][0] = 1.0f;

						for (j=1; j<nparam; j+=2) {
							yfit[i] += param[j] * exp_conv[j+1][i];
							/* A_j . exp(-x/tau_j) */
							dy_dparam[i][j] = exp_conv[j+1][i];
							/* exp(-x/tau_j) */
							dy_dparam[i][j+1] = param[j] * exp_conv[j][i];
							/* A_j * x * exp(-x/tau_j) / tau_j^2 */
						}
						dy[i] = y[i] - yfit[i];

						// gaussian noise based on fit
						weight = (yfit[i] > 1.0f ? 1.0f / yfit[i] : 1.0f);
						alpha_weight[i] = weight; // 1 / sig(i)
						weight *= dy[i];
						beta_weight[i] = weight; // dy[i] / sig(i)
						weight *= dy[i];
						*chisq += weight; // dy[i] / sig(i)
					}
					break;
				case NOISE_MLE:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// multi-exponential fit
						yfit[i] = param[0];  /* Z */
						dy_dparam[i][0] = 1.0f;

						for (j=1; j<nparam; j+=2) {
							yfit[i] += param[j] * exp_conv[j+1][i];
							/* A_j . exp(-x/tau_j) */
							dy_dparam[i][j] = exp_conv[j+1][i];
							/* exp(-x/tau_j) */
							dy_dparam[i][j+1] = param[j] * exp_conv[j][i];
							/* A_j * x * exp(-x/tau_j) / tau_j^2 */
						}
						dy[i] = y[i] - yfit[i];


						// maximum likelihood estimation noise
						weight = (yfit[i] > 1 ? 1.0f / yfit[i] : 1.0f);
						alpha_weight[i] = weight * y[i] / yfit[i];
						beta_weight[i] = dy[i] * weight;
						if (yfit[i] > 0.0f) {
							*chisq += (0.0f == y[i])
									? 2.0f * yfit[i]
									: 2.0f * (yfit[i] - y[i]) - 2.0f * y[i] * logf(yfit[i] / y[i]);
						}
					}
					if (*chisq <= 0.0f) {
						*chisq = 1.0e38f; // don't let chisq=0 through yfit being all -ve
					}
					break;
				default:
					return -3;
			}
			break;
		case FIT_GLOBAL_STRETCHEDEXP:
			switch (noise) {
				case NOISE_CONST:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// stretched exponential fit
						yfit[i] = param[0] + param[1] * exp_conv[1][i];
						dy[i] = y[i] - yfit[i];

						dy_dparam[i][0] = 1.0f;
						dy_dparam[i][1] = exp_conv[1][i];
						dy_dparam[i][2] = param[1] * exp_conv[2][i];
						dy_dparam[i][3] = param[1] * exp_conv[3][i];

						// constant noise
						weight = 1.0f / sig[0];
						alpha_weight[i] = weight; // 1 / (sig[0] * sig[0]);
						weight *= dy[i];
						beta_weight[i] = weight; // dy[i] / (sig[0] * sig[0]);
						weight *= dy[i];
						*chisq += weight; // (dy[i] * dy[i]) / (sig[0] * sig[0]);
					}
					break;
				case NOISE_GIVEN:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// stretched exponential fit
						yfit[i] = param[0] + param[1] * exp_conv[1][i];
						dy[i] = y[i] - yfit[i];

						dy_dparam[i][0] = 1.0f;
						dy_dparam[i][1] = exp_conv[1][i];
						dy_dparam[i][2] = param[1] * exp_conv[2][i];
						dy_dparam[i][3] = param[1] * exp_conv[3][i];

						// given noise
						weight = 1.0f / (sig[i] * sig[i]);
						alpha_weight[i] = weight; // 1 / (sig[i] * sig[i])
						weight *= dy[i];
						beta_weight[i] = weight; // dy[i] / (sig[i] * sig[i])
						weight *= dy[i];
						*chisq += weight; // (dy[i] * dy[i]) / (sig[i] * sig[i])
					}
					break;
				case NOISE_POISSON_DATA:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// stretched exponential fit
						yfit[i] = param[0] + param[1] * exp_conv[1][i];
						dy[i] = y[i] - yfit[i];

						dy_dparam[i][0] = 1.0f;
						dy_dparam[i][1] = exp_conv[1][i];
						dy_dparam[i][2] = param[1] * exp_conv[2][i];
						dy_dparam[i][3] = param[1] * exp_conv[3][i];

						// poisson noise based on data
						weight = (y[i] > 15 ? 1.0f / y[i] : 1.0f / 15);
						alpha_weight[i] = weight; // 1 / sig(i)
						weight *= dy[i];
						beta_weight[i] = weight; // dy[i] / sig(i)
						weight *= dy[i];
						*chisq += weight; // (dy[i] * dy[i]) / sig(i)
					}
					break;
				case NOISE_POISSON_FIT:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// stretched exponential fit
						yfit[i] = param[0] + param[1] * exp_conv[1][i];
						dy[i] = y[i] - yfit[i];

						dy_dparam[i][0] = 1.0f;
						dy_dparam[i][1] = exp_conv[1][i];
						dy_dparam[i][2] = param[1] * exp_conv[2][i];
						dy_dparam[i][3] = param[1] * exp_conv[3][i];

						// poisson noise based on fit
						weight = (yfit[i] > 15 ? 1.0f / yfit[i] : 1.0f / 15);
						alpha_weight[i] = weight; // 1 / sig(i)
						weight *= dy[i];
						beta_weight[i] = weight; // dy(i) / sig(i)
						weight *= dy[i];
						*chisq += weight; // (dy(i) * dy(i)) / sig(i)
					}
					break;
				case NOISE_GAUSSIAN_FIT:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// stretched exponential fit
						yfit[i] = param[0] + param[1] * exp_conv[1][i];
						dy[i] = y[i] - yfit[i];

						dy_dparam[i][0] = 1.0f;
						dy_dparam[i][1] = exp_conv[1][i];
						dy_dparam[i][2] = param[1] * exp_conv[2][i];
						dy_dparam[i][3] = param[1] * exp_conv[3][i];

						// gaussian noise based on fit
						weight = (yfit[i] > 1.0f ? 1.0f / yfit[i] : 1.0f);
						alpha_weight[i] = weight; // 1 / sig(i)
						weight *= dy[i];
						beta_weight[i] = weight; // dy[i] / sig(i)
						weight *= dy[i];
						*chisq += weight; // dy[i] / sig(i)
					}
					break;
				case NOISE_MLE:
					// loop over all data
					for (i = fit_start; i < fit_end; ++i) {
						// stretched exponential fit
						yfit[i] = param[0] + param[1] * exp_conv[1][i];
						dy[i] = y[i] - yfit[i];

						dy_dparam[i][0] = 1.0f;
						dy_dparam[i][1] = exp_conv[1][i];
						dy_dparam[i][2] = param[1] * exp_conv[2][i];
						dy_dparam[i][3] = param[1] * exp_conv[3][i];

						// maximum likelihood estimation noise
						weight = (yfit[i] > 1 ? 1.0f / yfit[i] : 1.0f);
						alpha_weight[i] = weight * y[i] / yfit[i];
						beta_weight[i] = dy[i] * weight;
						if (yfit[i] > 0.0f) {
							*chisq += (0.0f == y[i])
									? 2.0f * yfit[i]
									: 2.0f * (yfit[i] - y[i]) - 2.0f * y[i] * logf(yfit[i] / y[i]);
						}
					}
					if (*chisq <= 0.0f) {
						*chisq = 1.0e38f; // don't let chisq=0 through yfit being all -ve
					}
					break;
				default:
					return -3;
			}
			break;
		default:
			dbgprintf(1, "compute_exps_fn: please update me!\n");
			return -1;
	}

	// Check if chi square worsened:
	if (0.0f != old_chisq && *chisq >= old_chisq) {
		// don't bother to set up the matrices for solution
		return 0;
	}

	i_free = 0;
	// for all columns
	for (i = 0; i < nparam; ++i) {
		if (paramfree[i]) {
			j_free = 0;
			beta_sum = 0.0f;
			// row loop, only need to consider lower triangle
			for (j = 0; j <= i; ++j) {
				if (paramfree[j]) {
					dot_product = 0.0f;
					if (0 == j_free) { // true only once for each outer loop i
						// for all data
						for (k = fit_start; k < fit_end; ++k) {
							dy_dparam_k_i = dy_dparam[k][i];
							dot_product += dy_dparam_k_i * dy_dparam[k][j] * alpha_weight[k]; //TODO ARG make it [i][k] and just *dy_dparam++ it.
							beta_sum += dy_dparam_k_i * beta_weight[k];
						}
					}
					else {
						// for all data
						for (k = fit_start; k < fit_end; ++k) {
							dot_product += dy_dparam[k][i] * dy_dparam[k][j] * alpha_weight[k];
						}
					} // k loop

					alpha[j_free][i_free] = alpha[i_free][j_free] = dot_product;
					// if (i_free != j_free) {
					//     // matrix is symmetric
					//     alpha[i_free][j_free] = dot_product; //TODO dotProduct s/b including fixed parameters????!!!
					// }
					++j_free;
				}
			} // j loop
			beta[i_free] = beta_sum;
			++i_free;
		}
	} // i loop

	return 0;
}


/** And this is a final variant which computes the true chi-squared
   values and the full fit, as in EcfSingle.c */
int GCI_marquardt_global_compute_exps_fn_final(
			float xincr, float y[],
			int ndata, int fit_start, int fit_end,
			noise_type noise, float sig[], int ftype,
			float param[], int paramfree[], int nparam,
			float *exp_conv[],
			float yfit[], float dy[], float *chisq)
{
	int i, j, mfit;
	float sig2i;

	for (j=0, mfit=0; j<nparam; j++)
		if (paramfree[j])
			mfit++;

	/* Calculation of the fitting data will depend upon the type of
	   noise.  Since there's no convolution involved here, this is
	   very easy. */

	switch (ftype) {
	case FIT_GLOBAL_MULTIEXP:
		switch (noise) {
		case NOISE_CONST:
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0];  /* Z */
				for (j=1; j<nparam; j+=2) {
					yfit[i] += param[j] * exp_conv[j+1][i];
				        /* A_j . exp(-x/tau_j) */
				}
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end)
					*chisq += dy[i] * dy[i];
			}

			/* Now divide by sigma^2 */
			sig2i = 1.0f / (sig[0] * sig[0]);
			*chisq *= sig2i;
			break;

		case NOISE_GIVEN:  /* This is essentially the NR version */
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0];  /* Z */
				for (j=1; j<nparam; j+=2) {
					yfit[i] += param[j] * exp_conv[j+1][i];
				        /* A_j . exp(-x/tau_j) */
				}
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
					sig2i = 1.0f / (sig[i] * sig[i]);
					*chisq += dy[i] * dy[i] * sig2i;
				}
			}
			break;

		case NOISE_POISSON_DATA:
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0];  /* Z */
				for (j=1; j<nparam; j+=2) {
					yfit[i] += param[j] * exp_conv[j+1][i];
					    /* A_j . exp(-x/tau_j) */
				}
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
					/* we still don't let the variance drop below 1 */
					sig2i = (y[i] > 1 ? 1.0f/y[i] : 1.0f);
					*chisq += dy[i] * dy[i] * sig2i;
				}
			}
			break;

		case NOISE_POISSON_FIT:
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0];  /* Z */
				for (j=1; j<nparam; j+=2) {
					yfit[i] += param[j] * exp_conv[j+1][i];
					    /* A_j . exp(-x/tau_j) */
				}
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
					/* we still don't let the variance drop below 1 */
					sig2i = (yfit[i] > 1 ? 1.0f/yfit[i] : 1.0f);
					*chisq += dy[i] * dy[i] * sig2i;
				}
			}
			break;

		case NOISE_MLE:
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0];  /* Z */
				for (j=1; j<nparam; j+=2) {
					yfit[i] += param[j] * exp_conv[j+1][i];
					    /* A_j . exp(-x/tau_j) */
				}
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
//					sig2i = (yfit[i] > 1 ? 1.0/yfit[i] : 1.0);
//					*chisq += dy[i] * dy[i] * sig2i;
					if (yfit[i]<=0.0)
						; // do nothing
					else if (y[i]==0.0)
						*chisq += 2.0f*yfit[i];   // to avoid NaN from log
					else
						*chisq += 2.0f*(yfit[i]-y[i]) - 2.0f*y[i]*logf(yfit[i]/y[i]); // was dy[i] * dy[i] * sig2i;
				}
			}
			if (*chisq <= 0.0f) *chisq = 1.0e38f; // don't let chisq=0 through yfit being all -ve
			break;


		case NOISE_GAUSSIAN_FIT:
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0];  /* Z */
				for (j=1; j<nparam; j+=2) {
					yfit[i] += param[j] * exp_conv[j+1][i];
					    /* A_j . exp(-x/tau_j) */
				}
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
					sig2i = (yfit[i] > 1 ? 1.0f/yfit[i] : 1.0f);
					*chisq += dy[i] * dy[i] * sig2i;
				}
			}
			break;

		default:
			return -3;
			/* break; */   // (unreachable code)
		}
		break;

	case FIT_GLOBAL_STRETCHEDEXP:
		switch (noise) {
		case NOISE_CONST:
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0] + param[1] * exp_conv[1][i];
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end)
					*chisq += dy[i] * dy[i];
			}

			/* Now divide by sigma^2 */
			sig2i = 1.0f / (sig[0] * sig[0]);
			*chisq *= sig2i;
			break;

		case NOISE_GIVEN:  /* This is essentially the NR version */
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0] + param[1] * exp_conv[1][i];
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
					sig2i = 1.0f / (sig[i] * sig[i]);
					*chisq += dy[i] * dy[i] * sig2i;
				}
			}
			break;

		case NOISE_POISSON_DATA:
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0] + param[1] * exp_conv[1][i];
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
					/* we still don't let the variance drop below 1 */
					sig2i = (y[i] > 1 ? 1.0f/y[i] : 1.0f);
					*chisq += dy[i] * dy[i] * sig2i;
				}
			}
			break;

		case NOISE_POISSON_FIT:
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0] + param[1] * exp_conv[1][i];
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
					/* we still don't let the variance drop below 1 */
					sig2i = (yfit[i] > 1 ? 1.0f/yfit[i] : 1.0f);
					*chisq += dy[i] * dy[i] * sig2i;
				}
			}
			break;

		case NOISE_MLE:
			*chisq = 0.0f;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0] + param[1] * exp_conv[1][i];
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
//					sig2i = (yfit[i] > 1 ? 1.0/yfit[i] : 1.0);
//					*chisq += dy[i] * dy[i] * sig2i;
					if (yfit[i]<=0.0f)
						; // do nothing
					else if (y[i]==0.0)
						*chisq += 2.0f*yfit[i];   // to avoid NaN from log
					else
						*chisq += 2.0f*(yfit[i]-y[i]) - 2.0f*y[i]*logf(yfit[i]/y[i]); // was dy[i] * dy[i] * sig2i;
				}
			}
			if (*chisq <= 0.0) *chisq = 1.0e38f; // don't let chisq=0 through yfit being all -ve
			break;

		case NOISE_GAUSSIAN_FIT:
			*chisq = 0.0;
			/* Summation loop over all data */
			for (i=0; i<ndata; i++) {
				yfit[i] = param[0] + param[1] * exp_conv[1][i];
				dy[i] = y[i] - yfit[i];

				/* And find chi^2 */
				if (i >= fit_start && i < fit_end) {
					/* we still don't let the variance drop below 1 */
					sig2i = (yfit[i] > 1 ? 1.0f/yfit[i] : 1.0f);
					*chisq += dy[i] * dy[i] * sig2i;
				}
			}
			break;

		default:
			return -3;
			/* break; */   // (unreachable code)
		}
		break;

	default:
		dbgprintf(1, "compute_exps_fn: please update me!\n");
		return -1;
	}

	return 0;
}


/** This one does the actual global fitting for multiexponential taus.
   It is basically similar to the above do_fit_single function, except
   that now we handle the extra intricacies involved in global
   fitting, in particular, the much larger alpha matrix is handled in
   a special way. */

int GCI_marquardt_global_exps_do_fit_instr(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[], int ftype,
		float **param, int paramfree[], int nparam, restrain_type restrain,
		float chisq_delta, float exp_pure[], float *exp_conv[],
		float **fitted, float **residuals,
		float **covar_scratch, float **alpha_scratch,
		float *chisq_trans, float *chisq_global,
		int drop_bad_transients)
{
	float alambda, ochisq_global, *ochisq_trans;
	int i, k, itst, itst_max;
	int ret;

	itst_max = (restrain == ECF_RESTRAIN_DEFAULT) ? 4 : 6;

	/* If there are no global parameters being fitted, we simply fit
	   each local set. */
	switch (ftype) {
	case FIT_GLOBAL_MULTIEXP:
		for (i=2; i<nparam; i+=2) {
			if (paramfree[i]) {
				i = 0; /* sentinel value */
				break;
			}
		}
		break;

	case FIT_GLOBAL_STRETCHEDEXP:
		for (i=2; i<nparam; i++) {
			if (paramfree[i]) {
				i = 0; /* sentinel value */
				break;
			}
		}
		break;

	default:
		dbgprintf(1, "exps_do_fit_instr: please update me!\n");
		return -1;
	}

	if (i > 0) { /* no globals to fit */
		if (GCI_marquardt_global_exps_calculate_exps_instr(
				xincr, ndata, instr, ninstr, ftype, param[0], nparam,
				exp_pure, exp_conv) != 0)
			return -2;

		*chisq_global = 0;

		for (i=0; i<ntrans; i++) {
			if (drop_bad_transients && chisq_trans[i] < 0)
				continue;

			ret = GCI_marquardt_global_exps_do_fit_single(
					xincr, trans[i], ndata, fit_start, fit_end,
					noise, sig, ftype, param[i], paramfree, nparam, restrain,
////					exp_conv, fitted[i], residuals[i],
					chisq_delta, exp_conv, fitted[0], residuals[0],
					covar_scratch, alpha_scratch, &chisq_trans[i]);
			if (ret < 0) {
				if (drop_bad_transients) {
					dbgprintf(1, "In do_fit_instr, do_fit_single returned %d "
							  "for transient %d, dropping it\n", ret, i);
					chisq_trans[i] = -1;
				} else {
					dbgprintf(1, "In do_fit_instr, do_fit_single returned %d "
							  "for transient %d\n", ret, i);
					return -10 + ret;
				}
			} else {
				*chisq_global += chisq_trans[i];
			}
		}
		return 0;
	}

	/* If there are no free local variables to fit, we still do the
	   global fitting, but we have to be a little careful in some of
	   the later routines */

	/* Now allocate all of the arrays we will need. */

	if ((ochisq_trans = (float *) malloc((unsigned) ntrans * sizeof(float))) == NULL)
		return -1;

	/* We now begin our standard Marquardt loop, with several
	   modifications */
	alambda = -1;
	ret = GCI_marquardt_global_exps_global_step(
				xincr, trans, ndata, ntrans, fit_start, fit_end,
				instr, ninstr, noise, sig, ftype,
				param, paramfree, nparam, restrain, exp_pure, exp_conv,
				fitted, residuals, chisq_trans, chisq_global,
				alpha_scratch, &alambda, drop_bad_transients);
	if (ret != 0) {
		dbgprintf(1, "In do_fit_instr, first global_step returned %d\n", ret);
		if (ret != -1) {
			/* Wasn't a memory error, so unallocate arrays */
			alambda = 0.0;
			GCI_marquardt_global_exps_global_step(
				xincr, trans, ndata, ntrans, fit_start, fit_end,
				instr, ninstr, noise, sig, ftype,
				param, paramfree, nparam, restrain, exp_pure, exp_conv,
				fitted, residuals, chisq_trans, chisq_global,
				alpha_scratch, &alambda, drop_bad_transients);
		}
		free(ochisq_trans);
		return ret;
	}

	k = 1;  /* Iteration counter */
	itst = 0;
	for (;;) {
		dbgprintf(3, "In do_fit_instr, beginning iteration %d:\n", k);
		dbgprintf(3, " itst = %d, chisq_global = %.4f\n", itst, *chisq_global);

		k++;
		if (k > MAXITERS) {
			free(ochisq_trans);
			return -2;
		}

		ochisq_global = *chisq_global;
		for (i=0; i<ntrans; i++)
			ochisq_trans[i] = chisq_trans[i];

		ret = GCI_marquardt_global_exps_global_step(
					xincr, trans, ndata, ntrans, fit_start, fit_end,
					instr, ninstr, noise, sig, ftype,
					param, paramfree, nparam, restrain, exp_pure, exp_conv,
					fitted, residuals, chisq_trans, chisq_global,
					alpha_scratch, &alambda, drop_bad_transients);
		if (ret != 0) {
			dbgprintf(1, "In do_fit_instr, second global_step returned %d\n",
					  ret);
			/* Unallocate arrays */
			alambda = 0.0f;
			GCI_marquardt_global_exps_global_step(
				xincr, trans, ndata, ntrans, fit_start, fit_end,
				instr, ninstr, noise, sig, ftype,
				param, paramfree, nparam, restrain, exp_pure, exp_conv,
				fitted, residuals, chisq_trans, chisq_global,
				alpha_scratch, &alambda, drop_bad_transients);
			free(ochisq_trans);
			return ret;
		}

		if (*chisq_global > ochisq_global)
			itst = 0;
		else {
			/* Let's try this approach; I really don't know what will
			   be best */
			float maxdiff;

			maxdiff = 0.0f;
			for (i=0; i<ntrans; i++)
				if (ochisq_trans[i] - chisq_trans[i] > maxdiff)
					maxdiff = ochisq_trans[i] - chisq_trans[i];

			if (maxdiff < chisq_delta)
				itst++;
			dbgprintf(3, "In do_fit_instr, maxdiff = %.3f:\n", maxdiff);
		}

		if (itst < itst_max) continue;

		/* Endgame */
		alambda = 0.0f;
		ret = GCI_marquardt_global_exps_global_step(
					xincr, trans, ndata, ntrans, fit_start, fit_end,
					instr, ninstr, noise, sig, ftype,
					param, paramfree, nparam, restrain, exp_pure, exp_conv,
					fitted, residuals, chisq_trans, chisq_global,
					alpha_scratch, &alambda, drop_bad_transients);
		if (ret != 0) {
			dbgprintf(1, "In do_fit_instr, final global_step returned %d\n",
					  ret);
			free(ochisq_trans);
			return ret;
		}

		free(ochisq_trans);
		return k;  /* We're done now */
	}
}


/** And this one is basically a specialised GCI_marquardt_instr_step
   for the global fitting setup. */
int GCI_marquardt_global_exps_global_step(
				float xincr, float **trans,
				int ndata, int ntrans, int fit_start, int fit_end,
				float instr[], int ninstr,
				noise_type noise, float sig[], int ftype,
				float **param, int paramfree[], int nparam,
				restrain_type restrain,
				float exp_pure[], float *exp_conv[],
				float **yfit, float **dy,
				float *chisq_trans, float *chisq_global,
				float **alpha_scratch, float *alambda,
				int drop_bad_transients)
{
	int i, j, ret;
	static global_matrix alpha, covar;
	static global_vector beta, dparam;
	static float **paramtry;
	static int mfit_local, mfit_global;
	static int gindex[MAXFIT], lindex[MAXFIT];
	static float ochisq_global, *ochisq_trans;
	static void (*fitfunc)(float, float [], float *, float [], int);
	static int initialised=0;

	if (nparam > MAXFIT)
		return -10;
	if (xincr <= 0)
		return -11;
	if (fit_start < 0 || fit_start > fit_end || fit_end > ndata)
		return -12;

	/* Initialisation */
	/* We assume we're given sensible starting values for param[] */
	if (*alambda < 0.0f) {
		/* Start by allocating lots of variables we will need */
		mfit_local = mfit_global = 0;

		switch (ftype) {
		case FIT_GLOBAL_MULTIEXP:
			fitfunc = GCI_multiexp_tau;

			/* We know that all of param[2i], the taus, are the global
			   variables, and that the param[2i+1], the A's, are the
			   local variables, along with param[0] = Z.  We stored
			   the indices of the local and global free variables in
			   lindex and gindex respectively. */
			if (paramfree[0]) {
				lindex[mfit_local++] = 0;
			}
			for (i=1; i<nparam; i+=2) {
				if (paramfree[i])
					lindex[mfit_local++] = i;
			}
			for (i=2; i<nparam; i+=2) {
				if (paramfree[i])
					gindex[mfit_global++] = i;
			}
			break;

		case FIT_GLOBAL_STRETCHEDEXP:
			fitfunc = GCI_stretchedexp;

			if (paramfree[0])
				lindex[mfit_local++] = 0;  /* Z */
			if (paramfree[1])
				lindex[mfit_local++] = 1;  /* A */
			if (paramfree[2])
				gindex[mfit_global++] = 2;  /* tau */
			if (paramfree[3])
				gindex[mfit_global++] = 3;  /* h */
			break;

		default:
			dbgprintf(1, "exps_global_step: please update me!\n");
			return -1;
		}

		if (initialised) {
			GCI_free_global_matrix(&alpha);   GCI_free_global_matrix(&covar);
			GCI_ecf_free_matrix(paramtry);        GCI_free_global_vector(&beta);
			GCI_free_global_vector(&dparam);  free(ochisq_trans);
			initialised = 0;
		}

		if (GCI_alloc_global_matrix(&alpha, mfit_global, mfit_local, ntrans)
			!= 0)
			return -1;

		if (GCI_alloc_global_matrix(&covar, mfit_global, mfit_local, ntrans)
			!= 0) {
			GCI_free_global_matrix(&alpha);
			return -1;
		}

		if ((paramtry = GCI_ecf_matrix(ntrans, nparam)) == NULL) {
			GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
			return -1;
		}

		if (GCI_alloc_global_vector(&beta, mfit_global, mfit_local, ntrans)
			!= 0) {
			GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
			GCI_ecf_free_matrix(paramtry);
			return -1;
		}

		if (GCI_alloc_global_vector(&dparam, mfit_global, mfit_local, ntrans)
			!= 0) {
			GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
			GCI_ecf_free_matrix(paramtry);       GCI_free_global_vector(&beta);
			return -1;
		}

		if ((ochisq_trans = (float *) malloc((unsigned) ntrans * sizeof(float)))
			== NULL) {
			GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
			GCI_ecf_free_matrix(paramtry);       GCI_free_global_vector(&beta);
			GCI_free_global_vector(&dparam);
			return -1;
		}

		initialised = 1;

		if (GCI_marquardt_global_compute_global_exps_fn(
					xincr, trans, ndata, ntrans,
					fit_start, fit_end, instr, ninstr, noise, sig, ftype,
					param, paramfree, nparam, mfit_global, mfit_local,
					gindex, lindex, exp_pure, exp_conv,
					yfit, dy, alpha, beta, alpha_scratch,
					chisq_trans, chisq_global, drop_bad_transients) != 0)
			return -2;

		*alambda = 0.001f;
		ochisq_global = *chisq_global;
		for (i=0; i<ntrans; i++)
			ochisq_trans[i] = chisq_trans[i];

		/* Initialise paramtry to param */
		for (i=0; i<ntrans; i++) {
			for (j=0; j<nparam; j++)
				paramtry[i][j] = param[i][j];
		}
	}

	/* Once converged, evaluate covariance matrix */
	if (*alambda == 0) {
		if (GCI_marquardt_global_compute_global_exps_fn_final(
					xincr, trans, ndata, ntrans,
					fit_start, fit_end, instr, ninstr, noise, sig, ftype,
					param, paramfree, nparam, mfit_global, mfit_local,
					gindex, lindex, exp_pure, exp_conv, yfit, dy,
					chisq_trans, chisq_global, drop_bad_transients) != 0)
			return -3;
		/* Don't need to do this here; if we wished to, we'd have to
		   move this code (the "if (*alambda == 0)" block) to after
		   the Gauss-Jordan call.  We'd also need to rewrite it for
		   our situation.... */
		//	if (mfit < nparam) {  /* no need to do this otherwise */
		//		GCI_covar_sort(covar, nparam, paramfree, mfit);
		//		GCI_covar_sort(alpha, nparam, paramfree, mfit);
		//	}
		GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
		GCI_ecf_free_matrix(paramtry);       GCI_free_global_vector(&beta);
		GCI_free_global_vector(&dparam); free(ochisq_trans);
		initialised = 0;
		return 0;
	}

	/* Alter linearised fitting matrix by augmenting diagonal
	   elements. */
	GCI_copy_global_matrix(covar, alpha, mfit_global, mfit_local, ntrans);
	GCI_copy_global_vector(dparam, beta, mfit_global, mfit_local, ntrans);
	for (j=0; j<mfit_global; j++)
		covar.P[j][j] *= 1.0f + (*alambda);
	for (i=0; i<ntrans; i++)
		for (j=0; j<mfit_local; j++)
			covar.S[i][j][j] *= 1.0f + (*alambda);

	/* Matrix solution; GCI_solve solves Ax=b rather than AX=B */
	if (GCI_marquardt_global_solve_eqn(covar, dparam,
									   mfit_global, mfit_local, ntrans) != 0)
		return -3;

	/* Did the trial succeed?  Modify param by dparam... */
	for (i=0; i<ntrans; i++) {
		for (j=0; j<mfit_global; j++)
			paramtry[i][gindex[j]] = param[i][gindex[j]] + dparam.global[j];
		for (j=0; j<mfit_local; j++)
			paramtry[i][lindex[j]] =
				param[i][lindex[j]] + dparam.local[i*mfit_local + j];
	}

	for (i=0; i<ntrans; i++) {
		if (drop_bad_transients && chisq_trans[i] < 0)
			continue;

		if (restrain == ECF_RESTRAIN_DEFAULT)
			ret = check_ecf_params (paramtry[i], nparam, fitfunc);
		else
			ret = check_ecf_user_params (paramtry[i], nparam, fitfunc);

		if (ret != 0) {
			/* Bad parameters, increase alambda and return */
			*alambda *= 10.0f;
			return 0;
		}
	}

	if (GCI_marquardt_global_compute_global_exps_fn(
				xincr, trans, ndata, ntrans,
				fit_start, fit_end, instr, ninstr, noise, sig, ftype,
				paramtry, paramfree, nparam, mfit_global, mfit_local,
				gindex, lindex, exp_pure, exp_conv,
				yfit, dy, covar, dparam, alpha_scratch,
				chisq_trans, chisq_global, drop_bad_transients) != 0)
		return -2;

	/* Success, accept the new solution */
	if (*chisq_global < ochisq_global) {
		*alambda *= 0.1f;
		ochisq_global = *chisq_global;
		for (i=0; i<ntrans; i++)
			ochisq_trans[i] = chisq_trans[i];
		GCI_copy_global_matrix(alpha, covar, mfit_global, mfit_local, ntrans);
		GCI_copy_global_vector(beta, dparam, mfit_global, mfit_local, ntrans);
		for (i=0; i<ntrans; i++) {
			for (j=0; j<nparam; j++)
				param[i][j] = paramtry[i][j];
		}
	} else { /* Failure, increase alambda and return */
		*alambda *= 10.0f;
		*chisq_global = ochisq_global;
		for (i=0; i<ntrans; i++)
			chisq_trans[i] = ochisq_trans[i];
	}

	return 0;
}


/* Here we use alpha only for scratch space */
int GCI_marquardt_global_compute_global_exps_fn(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[], int ftype,
		float **param, int paramfree[], int nparam,
		int mfit_global, int mfit_local, int gindex[], int lindex[],
		float exp_pure[], float *exp_conv[],
		float **yfit, float **dy, global_matrix alpha, global_vector beta,
		float **alpha_scratch, float *chisq_trans, float *chisq_global,
		int drop_bad_transients)
{
	int i, j, k, ret;
	float beta_scratch[MAXFIT];  /* scratch space */

	/* Calculate the exponential array once only */
	if (GCI_marquardt_global_exps_calculate_exps_instr(
			xincr, ndata, instr, ninstr, ftype, param[0], nparam,
			exp_pure, exp_conv) != 0)
		return -1;

	/* We initialise P and beta_global to zero; the others don't
	   matter, as they will be totally overwritten */
	for (i=0; i<mfit_global; i++) {
		for (j=0; j<mfit_global; j++)
			alpha.P[i][j] = 0;
		beta.global[i] = 0;
	}
	*chisq_global = 0.0f;

	for (i=0; i<ntrans; i++) {
		if (drop_bad_transients && chisq_trans[i] < 0) {
			for (j=0; j<mfit_global; j++)
				for (k=0; k<mfit_local; k++)
					alpha.Q[j][i*mfit_local + k] = 0.0f;
			for (j=0; j<mfit_local; j++) {
				/* Make this component of S an identity matrix and of
				   beta zero */
				for (k=0; k<mfit_local; k++)
					alpha.S[i][j][k] = (j == k) ? 1.0f : 0.0f;
				beta.local[i*mfit_local + j] = 0;
			}
			continue;
		}

		/* This transient is fine! */
		ret = GCI_marquardt_global_compute_exps_fn(
					xincr, trans[i], ndata, fit_start, fit_end, noise, sig,
					ftype, param[i], paramfree, nparam,
////					exp_conv, yfit[i], dy[i],
					exp_conv, yfit[0], dy[0],
					alpha_scratch, beta_scratch, &chisq_trans[i], 0.0f);

		if (ret != 0) {
			if (drop_bad_transients) {
				dbgprintf(3, "In compute_global_exps_fn, "
						  "compute_exps_fn returned %d for transient %d\n",
						  ret, i);
				chisq_trans[i] = -1;
				continue;
			} else {
				dbgprintf(1, "In compute_global_exps_fn, "
						  "compute_exps_fn returned %d for transient %d\n",
						  ret, i);
				return -2;
			}
		}

		/* So now have to populate alpha and beta with the contents of
		   alpha_scratch and beta_scratch. */

		for (j=0; j<mfit_global; j++) {
			for (k=0; k<mfit_global; k++)
				alpha.P[j][k] += alpha_scratch[gindex[j]][gindex[k]];
			for (k=0; k<mfit_local; k++)
				alpha.Q[j][i*mfit_local + k] =
					alpha_scratch[gindex[j]][lindex[k]];
			beta.global[j] += beta_scratch[gindex[j]];
		}
		for (j=0; j<mfit_local; j++) {
			for (k=0; k<mfit_local; k++)
				alpha.S[i][j][k] = alpha_scratch[lindex[j]][lindex[k]];
			beta.local[i*mfit_local + j] = beta_scratch[lindex[j]];
		}

		*chisq_global += chisq_trans[i];
	}

	return 0;
}


/** The final variant */
int GCI_marquardt_global_compute_global_exps_fn_final(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[], int ftype,
		float **param, int paramfree[], int nparam,
		int mfit_global, int mfit_local, int gindex[], int lindex[],
		float exp_pure[], float *exp_conv[],
		float **yfit, float **dy,
		float *chisq_trans, float *chisq_global, int drop_bad_transients)
{
	int i, ret;

	/* Calculate the exponential array once only */
	if (GCI_marquardt_global_exps_calculate_exps_instr(
			xincr, ndata, instr, ninstr, ftype, param[0], nparam,
			exp_pure, exp_conv) != 0)
		return -1;

	*chisq_global = 0.0f;

	for (i=0; i<ntrans; i++) {
		if (drop_bad_transients && chisq_trans[i] < 0)
			continue;

		/* This transient is fine! */
		ret = GCI_marquardt_global_compute_exps_fn_final(
					xincr, trans[i], ndata, fit_start, fit_end, noise, sig,
					ftype, param[i], paramfree, nparam,
////					exp_conv, yfit[i], dy[i], &chisq_trans[i]);
					exp_conv, yfit[0], dy[0], &chisq_trans[i]);

		if (ret != 0) {
			if (drop_bad_transients) {
				dbgprintf(3, "In compute_global_exps_fn_final, "
						  "compute_exps_fn_final returned %d "
						  "for transient %d\n",
						  ret, i);
				chisq_trans[i] = -1;
				continue;
			} else {
				dbgprintf(1, "In compute_global_exps_fn_final, "
						  "compute_exps_fn_final returned %d "
						  "for transient %d\n",
						  ret, i);
				return -2;
			}
		}

		*chisq_global += chisq_trans[i];
	}

	return 0;
}


/** GCI_marquardt_global_solve_eqn.

   This function solves the equation Ax=b, where A is the alpha
   matrix, which has the form:

     A = (P Q)
         (R S)

   Here P is an mfit_global x mfit_global square matrix, S is a block
   diagonal matrix with ntrans blocks, each of size mfit_local x
   mfit_local, and Q and R are the right sizes to make the whole
   matrix square.  We solve it by inverting the matrix A using the
   formulae given in Numerical Recipes section 2.7, then multiplying
   the inverse by b to get x.  We are not too concerned about
   accuracy, as this does not affect the solution found, only the
   route taken to find it.

   Numerical Recipes, section 2.7, notes that A^{-1} is given by:

      (P' Q')
      (R' S')

   with:

      P' = (P - Q.S^{-1}.R)^{-1}
      Q' = - (P - Q.S^{-1}.R)^{-1} . (Q.S^{-1})
      R' = - (S^{-1}.R) . (P - Q.S^{-1}.R)^{-1}
      S' = S^{-1} + (S^{-1}.R) . (P - Q.S^{-1}.R)^{-1} . (Q.S^{-1})

   We also make use of the fact that A is symmetric, so in particular,
   (S^{-1}.R) = (Q.S^{-1})^T and R' = Q'^T.

   We are given A as a global_matrix and b as a global_vector.  This
   function destroys the original matrix and returns the solution in
   place of b.
*/

int GCI_marquardt_global_solve_eqn(global_matrix A, global_vector b,
			int mfit_global, int mfit_local, int ntrans)
{
	int row, col, block, i, j;
	float x_temp[MAXFIT], x_temp2[MAXFIT];
	static float **QS;
	static global_vector x;
	static int saved_global=0, saved_local=0, saved_ntrans=0;

	/* If no local parameters, just do a straight matrix solution */
	if (mfit_local == 0) {
		if (GCI_solve(A.P, mfit_global, b.global) != 0)
			return -2;
		return 0;
	}

	/* Allocate arrays if necessary */
	if ((saved_global != mfit_global) || (saved_local != mfit_local) ||
		(saved_ntrans != ntrans)) {
		if (saved_global > 0) {
			GCI_ecf_free_matrix(QS);
			GCI_free_global_vector(&x);
			saved_global = 0;
		}
		if ((QS = GCI_ecf_matrix(mfit_global, mfit_local*ntrans)) == NULL)
			return -1;
		if (GCI_alloc_global_vector(&x, mfit_global, mfit_local, ntrans)
			!= 0) {
			GCI_ecf_free_matrix(QS);
			return -1;
		}
		saved_global = mfit_global;
		saved_local = mfit_local;
		saved_ntrans = ntrans;
	}

	/* Start by inverting S */
	for (block=0; block<ntrans; block++)
		if (GCI_invert(A.S[block], mfit_local) != 0)
			return -2;

	/* Calculate Q.S^{-1} */
	for (row=0; row<mfit_global; row++)
		for (block=0; block<ntrans; block++)
			for (i=0; i<mfit_local; i++) {
				QS[row][block*mfit_local + i] = 0;
				for (j=0; j<mfit_local; j++)
					QS[row][block*mfit_local + i] +=
						A.Q[row][block*mfit_local + j] * A.S[block][j][i];
			}

	/* Now find P - Q.S^{-1}.R */
	for (row=0; row<mfit_global; row++)
		for (col=0; col<mfit_global; col++)
			for (i=0; i<ntrans*mfit_local; i++)
				A.P[row][col] -= QS[row][i] * A.Q[col][i];  /* Q = R^T */

	/* And invert it to get P' */
	if (GCI_invert(A.P, mfit_global) != 0)
		return -3;

	/* Now overwrite Q with Q' */
	for (row=0; row<mfit_global; row++)
		for (col=0; col<ntrans*mfit_local; col++) {
			A.Q[row][col] = 0;
			for (i=0; i<mfit_global; i++)
				A.Q[row][col] -= A.P[row][i] * QS[i][col];  /* P contains P' */
		}

	/* Finally, we can solve to find x */
	/* We have x.global = P'.(b.global) + Q'.(b.local)
	   and     x.local  = R'.(b.global) + S'.(b.local)
	   We do x.global first. */
	for (row=0; row<mfit_global; row++) {
		x.global[row] = 0;
		for (i=0; i<mfit_global; i++)
			x.global[row] += A.P[row][i] * b.global[i];
		/* Recall that Q now contains Q' */
		for (i=0; i < ntrans * mfit_local; i++)
			x.global[row] += A.Q[row][i] * b.local[i];
	}

	/* Now x_local; the R'.b_global component first, recalling that R'
	   is Q' transposed and that Q' is stored in Q. */
	for (row=0; row < ntrans * mfit_local; row++) {
		x.local[row] = 0;
		for (i=0; i<mfit_global; i++)
			x.local[row] += A.Q[i][row] * b.global[i];
	}

	/* Now S' = S^{-1} + (S^{-1}.R).(P-Q.S^{-1}.R)^{-1}.(Q.S^{-1}).
	   We first handle the S^{-1} term, then the remaining term.  */
	for (block=0; block<ntrans; block++)
		for (row=0; row<mfit_local; row++) {
			for (j=0; j<mfit_local; j++)
				x.local[block*mfit_local + row] +=
					A.S[block][row][j] * b.local[block*mfit_local + j];
		}

	/* For the remaining term, we have an x.local[row] contribution of

	      sum_j sum_k sum_l
	          (S^{-1}.R)_{row,j} . (P-Q.S^{-1}.R)^{-1}_{j,k} .
	          (Q.S^{-1})_{k,l} . b.local_{l}

	   In order to save computations, we calculate the matrices once
	   only.  We start with (Q.S^{-1}) . b_local, which is an
	   mfit_global x 1 array, and store this in x_temp; premultiplying
	   this by the square matrix P' again gives an mfit_global x 1
	   array, which goes in x_temp2, then premultiplying this by
	   (S^{-1}.R) gives an (ntrans * mfit_local) x 1 array which is
	   added directly onto x_local.  Recall also that S^{-1}.R is the
	   transpose of Q.S^{-1}, which is currently stored in QS, and
	   that the middle term is currently stored in P. */

	for (row=0; row<mfit_global; row++) {
		x_temp[row] = 0;
		for (i=0; i < ntrans*mfit_local; i++)
			x_temp[row] += QS[row][i] * b.local[i];
	}
	for (row=0; row<mfit_global; row++) {
		x_temp2[row] = 0;
		for (i=0; i<mfit_global; i++)
			x_temp2[row] += A.P[row][i] * x_temp[i];
	}
	/* Again, S^{-1}.R is the transpose of Q.S^{-1} */
	for (row=0; row < ntrans * mfit_local; row++) {
		for (i=0; i<mfit_global; i++)
			x.local[row] += QS[i][row] * x_temp2[i];
	}

	/* And we're done, once we've copied x into b */
	GCI_copy_global_vector(b, x, mfit_global, mfit_local, ntrans);

	return 0;
}


/** GENERIC FUNCTION CODE

   These functions are essentially the same as the above functions
   GCI_marquardt_global_exps_instr and
   GCI_marquardt_global_exps_do_fit_instr and the latter's dependents,
   except that this version takes an arbitrary function and a list of
   global parameters.  This function is designed to be called from
   external code.  It is nowhere near as efficient as the streamlined
   code above for globally fitting taus for multi-exponential models.
   Also, this function must be provided with meaningful starting
   estimates for all parameters.
   */

int GCI_marquardt_global_generic_instr(float xincr, float **trans,
					int ndata, int ntrans, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[],
					float **param, int paramfree[], int nparam, int gparam[],
					restrain_type restrain, float chisq_delta,
					void (*fitfunc)(float, float [], float *, float [], int),
					float **fitted, float **residuals,
					float chisq_trans[], float *chisq_global, int *df)
{
	float **covar, **alpha, *scaled_instr, instrsum;
	int i, ret;
	int mlocal, mglobal;

	/* Some basic parameter checks */
	if (xincr <= 0) return -1;
	if (ntrans < 1) return -1;
	if (ndata < 1)  return -1;
	if (fit_start < 0 || fit_end > ndata) return -1;
	if (ninstr < 1) return -1;
	if (nparam < 1) return -1;

	if ((covar = GCI_ecf_matrix(nparam, nparam)) == NULL)
		return -2;

	if ((alpha = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		GCI_ecf_free_matrix(covar);
		return -3;
	}

	if ((scaled_instr = (float *) malloc((unsigned) ninstr * sizeof(float))) == NULL) {
		GCI_ecf_free_matrix(covar);
		GCI_ecf_free_matrix(alpha);
		return -4;
	}

	/* Scale the instrument response */
	for (i=0, instrsum=0; i<ninstr; i++)
		instrsum += instr[i];
	if (instrsum == 0) {
		GCI_ecf_free_matrix(covar);
		GCI_ecf_free_matrix(alpha);
		free(scaled_instr);
		return -6;
	}
	
	for (i=0; i<ninstr; i++)
		scaled_instr[i] = instr[i] / instrsum;

	/* Now call the global fitting function. */
	ret = GCI_marquardt_global_generic_do_fit_instr(
			xincr, trans, ndata, ntrans, fit_start, fit_end,
			scaled_instr, ninstr, noise, sig,
			param, paramfree, nparam, gparam, restrain, chisq_delta,
			fitfunc, fitted, residuals, covar, alpha,
			chisq_trans, chisq_global);

	GCI_ecf_free_matrix(covar);
	GCI_ecf_free_matrix(alpha);
	free(scaled_instr);
	GCI_marquardt_cleanup();

	if (ret < 0) {
		dbgprintf(1, "Fit failed, ret = %d\n", ret);
		GCI_ecf_free_matrix(covar);
		GCI_ecf_free_matrix(alpha);
		free(scaled_instr);
		return -10 + ret;
	}

	dbgprintf(2, "Fit succeeded, ret = %d\n", ret);

	/* Before we return, calculate the number of degrees of freedom */
	/* The number of degrees of freedom is given by:
	      d.f. = ntrans * ((fit_end - fit_start) - # free local parameters)
	                 - # free global parameters
	*/

	mglobal = mlocal = 0;
	for (i=0; i<nparam; i++)
		if (paramfree[i]) {
			if (gparam[i]) mglobal++;
			else mlocal++;
		}

	*df = ntrans * ((fit_end - fit_start) - mlocal) - mglobal;

	GCI_ecf_free_matrix(covar);
	GCI_ecf_free_matrix(alpha);
	free(scaled_instr);

	return ret;
}


int GCI_marquardt_global_generic_do_fit_instr(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[],
		float **param, int paramfree[], int nparam, int gparam[],
		restrain_type restrain, float chisq_delta,
		void (*fitfunc)(float, float [], float *, float [], int),
		float **fitted, float **residuals,
		float **covar_scratch, float **alpha_scratch,
		float *chisq_trans, float *chisq_global)
{
	// PRB 03/07 Although **fitted and **residuals are provided only one "transient" is required and used, fitted[0] and residuals[0]

	float alambda, ochisq_global, *ochisq_trans;
	int i, k, itst, itst_max;
	int ret;

	itst_max = (restrain == ECF_RESTRAIN_DEFAULT) ? 4 : 6;

	/* If there are no global parameters being fitted, we simply fit
	   each local set. */
	for (i=0; i<nparam; i++) {
		if (gparam[i] && paramfree[i]) {
			i = -1; /* sentinel value */
			break;
		}
	}

	if (i >= 0) { /* no globals to fit */
		*chisq_global = 0;

		for (i=0; i<ntrans; i++) {
			ret = GCI_marquardt_instr(xincr, trans[i],
					ndata, fit_start, fit_end, instr, ninstr, noise, sig,
					param[i], paramfree, nparam, restrain, fitfunc,
					fitted[0], residuals[0], covar_scratch, alpha_scratch,
					&chisq_trans[i], chisq_delta, 0, NULL);
			if (ret < 0) {
				dbgprintf(1, "In do_fit_instr, marquardt_instr returned %d "
						  "for transient %d\n", ret, i);
				return -10 + ret;
			} else {
				*chisq_global += chisq_trans[i];
			}
		}
		return 0;
	}

	/* If there are no free local variables to fit, we still do the
	   global fitting, but we have to be a little careful in some of
	   the later routines */

	/* Now allocate all of the arrays we will need. */

	if ((ochisq_trans = (float *) malloc((unsigned) ntrans * sizeof(float))) == NULL)
		return -1;

	/* We now begin our standard Marquardt loop, with several
	   modifications */
	alambda = -1;
	ret = GCI_marquardt_global_generic_global_step(
				xincr, trans, ndata, ntrans, fit_start, fit_end,
				instr, ninstr, noise, sig,
				param, paramfree, nparam, gparam, restrain, chisq_delta, fitfunc,
				fitted, residuals, chisq_trans, chisq_global,
				alpha_scratch, &alambda);
	if (ret != 0) {
		dbgprintf(1, "In do_fit_instr, first global_step returned %d\n", ret);
		if (ret != -1) {
			/* Wasn't a memory error, so unallocate arrays */
			alambda = 0.0f;
			GCI_marquardt_global_generic_global_step(
				xincr, trans, ndata, ntrans, fit_start, fit_end,
				instr, ninstr, noise, sig,
				param, paramfree, nparam, gparam, restrain, chisq_delta, fitfunc,
				fitted, residuals, chisq_trans, chisq_global,
				alpha_scratch, &alambda);
		}
		free(ochisq_trans);
		return ret;
	}

	k = 1;  /* Iteration counter */
	itst = 0;
	for (;;) {
		dbgprintf(3, "In do_fit_instr, beginning iteration %d:\n", k);
		dbgprintf(3, " itst = %d, chisq_global = %.4f\n", itst, *chisq_global);

		k++;
		if (k > MAXITERS) {
			free(ochisq_trans);
			return -2;
		}

		ochisq_global = *chisq_global;
		for (i=0; i<ntrans; i++)
			ochisq_trans[i] = chisq_trans[i];

		ret = GCI_marquardt_global_generic_global_step(
					xincr, trans, ndata, ntrans, fit_start, fit_end,
					instr, ninstr, noise, sig,
					param, paramfree, nparam, gparam, restrain, chisq_delta, fitfunc,
					fitted, residuals, chisq_trans, chisq_global,
					alpha_scratch, &alambda);
		if (ret != 0) {
			dbgprintf(1, "In do_fit_instr, second global_step returned %d\n",
					  ret);
			/* Unallocate arrays */
			alambda = 0.0f;
			GCI_marquardt_global_generic_global_step(
				xincr, trans, ndata, ntrans, fit_start, fit_end,
				instr, ninstr, noise, sig,
				param, paramfree, nparam, gparam, restrain, chisq_delta, fitfunc,
				fitted, residuals, chisq_trans, chisq_global,
				alpha_scratch, &alambda);
			free(ochisq_trans);
			return ret;
		}

		if (*chisq_global > ochisq_global)
			itst = 0;
		else {
			/* Let's try this approach; I really don't know what will
			   be best */
			float maxdiff;

			maxdiff = 0.0f;
			for (i=0; i<ntrans; i++)
				if (ochisq_trans[i] - chisq_trans[i] > maxdiff)
					maxdiff = ochisq_trans[i] - chisq_trans[i];

			if (maxdiff < chisq_delta)
				itst++;
			dbgprintf(3, "In do_fit_instr, maxdiff = %.3f:\n", maxdiff);
		}

		if (itst < itst_max) continue;

		/* Endgame */
		alambda = 0.0f;
		ret = GCI_marquardt_global_generic_global_step(
					xincr, trans, ndata, ntrans, fit_start, fit_end,
					instr, ninstr, noise, sig,
					param, paramfree, nparam, gparam, restrain, chisq_delta, fitfunc,
					fitted, residuals, chisq_trans, chisq_global,
					alpha_scratch, &alambda);
		if (ret != 0) {
			dbgprintf(1, "In do_fit_instr, final global_step returned %d\n",
					  ret);
			free(ochisq_trans);
			return ret;
		}

		free(ochisq_trans);
		return k;  /* We're done now */
	}
}


/* And this one is basically a specialised GCI_marquardt_instr_step
   for the global fitting setup. */

#define do_frees \
	if (fnvals) free(fnvals);\
	if (dy_dparam_pure) GCI_ecf_free_matrix(dy_dparam_pure);\
	if (dy_dparam_conv) GCI_ecf_free_matrix(dy_dparam_conv);

int GCI_marquardt_global_generic_global_step(
				float xincr, float **trans,
				int ndata, int ntrans, int fit_start, int fit_end,
				float instr[], int ninstr,
				noise_type noise, float sig[],
				float **param, int paramfree[], int nparam, int gparam[],
				restrain_type restrain, float chisq_delta,
				void (*fitfunc)(float, float [], float *, float [], int),
				float **yfit, float **dy,
				float *chisq_trans, float *chisq_global,
				float **alpha_scratch, float *alambda)
{
	int i, j, ret;
	static global_matrix alpha, covar;
	static global_vector beta, dparam;
	static float **paramtry;
	static int mfit_local, mfit_global;
	static int gindex[MAXFIT], lindex[MAXFIT];
	static float ochisq_global, *ochisq_trans;
	static int initialised=0;

	// The following are declared here to retain some optimisation by not repeatedly mallocing
	// (only once per transient), but to remain thread safe.
	// They are malloced by lower fns but at the end, freed by this fn.
	// These vars were global or static before thread safety was introduced.
	float *fnvals=NULL, **dy_dparam_pure=NULL, **dy_dparam_conv=NULL;
	int fnvals_len=0, dy_dparam_nparam_size=0;

	if (nparam > MAXFIT)
		return -10;
	if (xincr <= 0)
		return -11;
	if (fit_start < 0 || fit_start > fit_end || fit_end > ndata)
		return -12;

	/* Initialisation */
	/* We assume we're given sensible starting values for param[] */
	if (*alambda < 0.0f) {
		/* Start by allocating lots of variables we will need */
		mfit_local = mfit_global = 0;

		for (i=0; i<nparam; i++) {
			if (paramfree[i]) {
				if (gparam[i])
					gindex[mfit_global++] = i;
				else
					lindex[mfit_local++] = i;
			}
		}

		if (initialised) {
			GCI_free_global_matrix(&alpha);   GCI_free_global_matrix(&covar);
			GCI_ecf_free_matrix(paramtry);        GCI_free_global_vector(&beta);
			GCI_free_global_vector(&dparam);  free(ochisq_trans);
			initialised = 0;
		}

		if (GCI_alloc_global_matrix(&alpha, mfit_global, mfit_local, ntrans)
			!= 0)
			return -1;

		if (GCI_alloc_global_matrix(&covar, mfit_global, mfit_local, ntrans)
			!= 0) {
			GCI_free_global_matrix(&alpha);
			return -1;
		}

		if ((paramtry = GCI_ecf_matrix(ntrans, nparam)) == NULL) {
			GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
			return -1;
		}

		if (GCI_alloc_global_vector(&beta, mfit_global, mfit_local, ntrans)
			!= 0) {
			GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
			GCI_ecf_free_matrix(paramtry);
			return -1;
		}

		if (GCI_alloc_global_vector(&dparam, mfit_global, mfit_local, ntrans)
			!= 0) {
			GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
			GCI_ecf_free_matrix(paramtry);       GCI_free_global_vector(&beta);
			return -1;
		}

		if ((ochisq_trans = (float *) malloc((unsigned) ntrans * sizeof(float)))
			== NULL) {
			GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
			GCI_ecf_free_matrix(paramtry);       GCI_free_global_vector(&beta);
			GCI_free_global_vector(&dparam);
			return -1;
		}

		initialised = 1;

		if (GCI_marquardt_global_compute_global_generic_fn(
					xincr, trans, ndata, ntrans,
					fit_start, fit_end, instr, ninstr, noise, sig,
					param, paramfree, nparam, gparam,
					mfit_global, mfit_local, gindex, lindex, fitfunc,
					yfit, dy, alpha, beta, alpha_scratch,
					chisq_trans, chisq_global, *alambda,
					&fnvals, &dy_dparam_pure, &dy_dparam_conv,
					&fnvals_len, &dy_dparam_nparam_size) != 0)
			return -2;

		*alambda = 0.001f;
		ochisq_global = *chisq_global;
		for (i=0; i<ntrans; i++)
			ochisq_trans[i] = chisq_trans[i];

		/* Initialise paramtry to param */
		for (i=0; i<ntrans; i++) {
			for (j=0; j<nparam; j++)
				paramtry[i][j] = param[i][j];
		}
	}

	/* Once converged, evaluate covariance matrix */
	if (*alambda == 0) {
		if (GCI_marquardt_global_compute_global_generic_fn_final(
					xincr, trans, ndata, ntrans,
					fit_start, fit_end, instr, ninstr, noise, sig,
					param, paramfree, nparam, gparam,
					mfit_global, mfit_local, gindex, lindex, fitfunc,
					yfit, dy, chisq_trans, chisq_global,
					&fnvals, &dy_dparam_pure, &dy_dparam_conv,
					&fnvals_len, &dy_dparam_nparam_size) != 0)
			return -3;
		/* Don't need to do this here; if we wished to, we'd have to
		   move this code (the "if (*alambda == 0)" block) to after
		   the Gauss-Jordan call.  We'd also need to rewrite it for
		   our situation.... */
		//	if (mfit < nparam) {  /* no need to do this otherwise */
		//		GCI_covar_sort(covar, nparam, paramfree, mfit);
		//		GCI_covar_sort(alpha, nparam, paramfree, mfit);
		//	}
		GCI_free_global_matrix(&alpha);  GCI_free_global_matrix(&covar);
		GCI_ecf_free_matrix(paramtry);       GCI_free_global_vector(&beta);
		GCI_free_global_vector(&dparam); free(ochisq_trans);
		initialised = 0;
		return 0;
	}

	/* Alter linearised fitting matrix by augmenting diagonal
	   elements. */
	GCI_copy_global_matrix(covar, alpha, mfit_global, mfit_local, ntrans);
	GCI_copy_global_vector(dparam, beta, mfit_global, mfit_local, ntrans);
	for (j=0; j<mfit_global; j++)
		covar.P[j][j] *= 1.0f + (*alambda);
	for (i=0; i<ntrans; i++)
		for (j=0; j<mfit_local; j++)
			covar.S[i][j][j] *= 1.0f + (*alambda);

	/* Matrix solution; GCI_solve solves Ax=b rather than AX=B */
	if (GCI_marquardt_global_solve_eqn(covar, dparam,
									   mfit_global, mfit_local, ntrans) != 0)
		return -3;

	/* Did the trial succeed?  Modify param by dparam... */
	for (i=0; i<ntrans; i++) {
		for (j=0; j<mfit_global; j++)
			paramtry[i][gindex[j]] = param[i][gindex[j]] + dparam.global[j];
		for (j=0; j<mfit_local; j++)
			paramtry[i][lindex[j]] =
				param[i][lindex[j]] + dparam.local[i*mfit_local + j];
	}

	for (i=0; i<ntrans; i++) {
		if (restrain == ECF_RESTRAIN_DEFAULT)
			ret = check_ecf_params (paramtry[i], nparam, fitfunc);
		else
			ret = check_ecf_user_params (paramtry[i], nparam, fitfunc);

		if (ret != 0) {
			/* Bad parameters, increase alambda and return */
			*alambda *= 10.0f;
			return 0;
		}
	}

	if (GCI_marquardt_global_compute_global_generic_fn(
				xincr, trans, ndata, ntrans,
				fit_start, fit_end, instr, ninstr, noise, sig,
				paramtry, paramfree, nparam, gparam,
				mfit_global, mfit_local, gindex, lindex, fitfunc,
				yfit, dy, covar, dparam, alpha_scratch,
				chisq_trans, chisq_global, *alambda,
				&fnvals, &dy_dparam_pure, &dy_dparam_conv,
				&fnvals_len, &dy_dparam_nparam_size) != 0)
		return -2;

	/* Success, accept the new solution */
	if (*chisq_global < ochisq_global) {
		*alambda *= 0.1f;
		ochisq_global = *chisq_global;
		for (i=0; i<ntrans; i++)
			ochisq_trans[i] = chisq_trans[i];
		GCI_copy_global_matrix(alpha, covar, mfit_global, mfit_local, ntrans);
		GCI_copy_global_vector(beta, dparam, mfit_global, mfit_local, ntrans);
		for (i=0; i<ntrans; i++) {
			for (j=0; j<nparam; j++)
				param[i][j] = paramtry[i][j];
		}
	} else { /* Failure, increase alambda and return */
		*alambda *= 10.0f;
		*chisq_global = ochisq_global;
		for (i=0; i<ntrans; i++)
			chisq_trans[i] = ochisq_trans[i];
	}

	return 0;
}


/* Here we use alpha only for scratch space */
int GCI_marquardt_global_compute_global_generic_fn(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[],
		float **param, int paramfree[], int nparam, int gparam[],
		int mfit_global, int mfit_local, int gindex[], int lindex[],
		void (*fitfunc)(float, float [], float *, float [], int),
		float **yfit, float **dy, global_matrix alpha, global_vector beta,
		float **alpha_scratch, float *chisq_trans, float *chisq_global,
		float alambda,
		float **pfnvals, float ***pdy_dparam_pure, float ***pdy_dparam_conv,
		int *pfnvals_len, int *pdy_dparam_nparam_size)
{
	int i, j, k, ret;
	float beta_scratch[MAXFIT];  /* scratch space */

	/* We initialise P and beta_global to zero; the others don't
	   matter, as they will be totally overwritten */
	for (i=0; i<mfit_global; i++) {
		for (j=0; j<mfit_global; j++)
			alpha.P[i][j] = 0.0f;
		beta.global[i] = 0.0f;
	}
	*chisq_global = 0.0f;

	for (i=0; i<ntrans; i++) {
		/* Only pass the true alambda, used for initialisation, for
		   the first transient */
		ret = GCI_marquardt_compute_fn_instr(
					xincr, trans[i], ndata, fit_start, fit_end,
					instr, ninstr, noise, sig,
					param[i], paramfree, nparam, fitfunc,
////					yfit[i], dy[i], alpha_scratch, beta_scratch,
					yfit[0], dy[0], alpha_scratch, beta_scratch,
					&chisq_trans[i], 0.0f, (i == 0) ? alambda : 0.0f, //TODO ARG added 0.0f here for new old_chisq parameter
					pfnvals, pdy_dparam_pure, pdy_dparam_conv,
					pfnvals_len, pdy_dparam_nparam_size);

		if (ret != 0) {
			dbgprintf(1, "In compute_global_generic_fn, "
					  "compute_fn_instr returned %d for transient %d\n",
					  ret, i);
			return -2;
		}

		/* So now have to populate alpha and beta with the contents of
		   alpha_scratch and beta_scratch. */

		for (j=0; j<mfit_global; j++) {
			for (k=0; k<mfit_global; k++)
				alpha.P[j][k] += alpha_scratch[gindex[j]][gindex[k]];
			for (k=0; k<mfit_local; k++)
				alpha.Q[j][i*mfit_local + k] =
					alpha_scratch[gindex[j]][lindex[k]];
			beta.global[j] += beta_scratch[gindex[j]];
		}
		for (j=0; j<mfit_local; j++) {
			for (k=0; k<mfit_local; k++)
				alpha.S[i][j][k] = alpha_scratch[lindex[j]][lindex[k]];
			beta.local[i*mfit_local + j] = beta_scratch[lindex[j]];
		}

		*chisq_global += chisq_trans[i];
	}

	return 0;
}


/* And the final variant */
int GCI_marquardt_global_compute_global_generic_fn_final(
		float xincr, float **trans, int ndata, int ntrans,
		int fit_start, int fit_end, float instr[], int ninstr,
		noise_type noise, float sig[],
		float **param, int paramfree[], int nparam, int gparam[],
		int mfit_global, int mfit_local, int gindex[], int lindex[],
		void (*fitfunc)(float, float [], float *, float [], int),
		float **yfit, float **dy,
		float *chisq_trans, float *chisq_global,
		float **pfnvals, float ***pdy_dparam_pure, float ***pdy_dparam_conv,
		int *pfnvals_len, int *pdy_dparam_nparam_size)
{
	int i, ret;

	*chisq_global = 0.0f;

	for (i=0; i<ntrans; i++) {
		/* Only pass the true alambda, used for initialisation, for
		   the first transient */
		ret = GCI_marquardt_compute_fn_final_instr(
					xincr, trans[i], ndata, fit_start, fit_end,
					instr, ninstr, noise, sig,
					param[i], paramfree, nparam, fitfunc,
//					yfit[i], dy[i], &chisq_trans[i]);
					yfit[0], dy[0], &chisq_trans[i],
					pfnvals, pdy_dparam_pure, pdy_dparam_conv,
					pfnvals_len, pdy_dparam_nparam_size);

		if (ret != 0) {
			dbgprintf(1, "In compute_global_generic_fn_final, "
					  "compute_fn_final_instr returned %d for transient %d\n",
					  ret, i);
			return -2;
		}

		*chisq_global += chisq_trans[i];
	}

	return 0;
}


// Emacs settings:
// Local variables:
// mode: c
// c-basic-offset: 4
// tab-width: 4
// End:

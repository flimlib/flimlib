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

   This file contains functions for support plane analysis.
*/

#include <stdio.h>
#include <stdlib.h>
#include "EcfInternal.h"


/** The support plane analysis code. 

   Calls the appropriate Marquardt
   function to perform a fit for each parameter value set we are
   interested in.
   
   * \file EcfSPA.h
 */

/**
   GCI_SPA_1D_*  performs a 1-dimensional support plane analysis.  It
                 takes the same fitting paramters as the corresponding
                 standard Marquardt function, without a few of the
                 transient-specific parameters, and with the following
                 extra parameters:

                 int spa_param      which parameter are we analysing?
				 int spa_nvalues    how many different parameter
                                    values will we calculate
                                    chi-squared values for?
				 float spa_low      the lowest parameter value
                 float spa_high     the highest parameter value
				 float chisq[]      the resulting chisq values

				 In the case of the global analysis for exponentials
				 function, there is an extra parameter, int df[],
				 which gives the number of degrees of freedom for each
				 value; this is relevant if the drop_bad_transients
				 parameter is set.  If an individual fit fails, the
				 corresponding chisq[] value is set to -1.  Also, the
				 chisq[] array runs from chisq[0], corresponding to
				 spa_low to chisq[spa_nvalues-1], corresponding to
				 spa_high.
*/
int GCI_SPA_1D_marquardt(
				float x[], float y[], int ndata,
				noise_type noise, float sig[],
				float param[], int paramfree[], int nparam,
				restrain_type restrain, float chisq_delta,
				void (*fitfunc)(float, float [], float *, float [], int),
				int spa_param, int spa_nvalues,
				float spa_low, float spa_high,
				float chisq[], void (*progressfunc)(float))
{
	int i, j, ret;
	float *fitted, *residuals, **covar, **alpha;
	float param_copy[MAXFIT];
	int paramfree_copy[MAXFIT];
	
	if (spa_param < 0 || spa_param >= nparam)  /* and so nparam > 0, too */
		return -1;
	if (spa_nvalues < 2)
		return -2;
	
	if ((fitted = (float *) malloc((size_t)ndata * sizeof(float))) == NULL)
		return -3;
	if ((residuals = (float *) malloc((size_t)ndata * sizeof(float))) == NULL) {
		free(fitted);
		return -3;
	}
	if ((covar = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		free(fitted);
		free(residuals);
		return -3;
	}
	if ((alpha = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		free(fitted);
		free(residuals);
		GCI_ecf_free_matrix(covar);
		return -3;
	}

	for (j=0; j<nparam; j++)
		paramfree_copy[j] = paramfree[j];
	paramfree_copy[spa_param] = 0;  /* we fix the parameter we are
									   analysing */

	for (i=0; i<spa_nvalues; i++) {
		/* Initialise parameter array each time */
		for (j=0; j<nparam; j++)
			param_copy[j] = param[j];
		/* Set the parameter we are analysing */
		param_copy[spa_param] =
			spa_low + (spa_high - spa_low) * (float)i / (float)(spa_nvalues - 1);

		ret = GCI_marquardt(x, y, ndata, noise, sig,
							param_copy, paramfree_copy, nparam, restrain,
							fitfunc, fitted, residuals, covar, alpha,
							&chisq[i], chisq_delta, 0, NULL);
			
		if (progressfunc)
            progressfunc ((float)i/(float)(spa_nvalues-1));

		if (ret < 0)
			chisq[i] = -1;
	}

	free(fitted);
	free(residuals);
	GCI_ecf_free_matrix(covar);
	GCI_ecf_free_matrix(alpha);

	return 0;
}

/**
   GCI_SPA_2D_*  performs a 2-dimensional support plane analysis.  It
                 is the same as the 1-dimensional variant, except that
                 now the extra parameters are:

                 int spa_param1, int spa_nvalues1,
				 float spa_low1, float spa_high1,
                 int spa_param1, int spa_nvalues1,
				 float spa_low1, float spa_high1,
				 float chisq[][], (int df[][])

				 Here chisq[i][j] will correspond to the i-th value of
				 param1 and the j-th value of param2; chisq[][] itself
				 will have been allocated by the GCI_ecf_matrix function.

   In all cases where the original function must be provided with
   reasonable initial estimates for the parameters, so must these
   functions.
*/

int GCI_SPA_2D_marquardt(
				float x[], float y[], int ndata,
				noise_type noise, float sig[],
				float param[], int paramfree[], int nparam,
				restrain_type restrain, float chisq_delta,
				void (*fitfunc)(float, float [], float *, float [], int),
				int spa_param1, int spa_nvalues1,
				float spa_low1, float spa_high1,
				int spa_param2, int spa_nvalues2,
				float spa_low2, float spa_high2,
				float **chisq, void (*progressfunc)(float))
{
	int i1, i2, j, ret, progress, total;
	float *fitted, *residuals, **covar, **alpha;
	float param_copy[MAXFIT];
	int paramfree_copy[MAXFIT];
	
	if (spa_param1 < 0 || spa_param1 >= nparam)  /* and so nparam > 0, too */
		return -1;
	if (spa_param2 < 0 || spa_param2 >= nparam)
		return -1;
	if (spa_param1 == spa_param2)
		return -1;
	if (spa_nvalues1 < 2 || spa_nvalues2 < 2)
		return -2;
	
	if ((fitted = (float *) malloc((size_t)ndata * sizeof(float))) == NULL)
		return -3;
	if ((residuals = (float *) malloc((size_t)ndata * sizeof(float))) == NULL) {
		free(fitted);
		return -3;
	}
	if ((covar = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		free(fitted);
		free(residuals);
		return -3;
	}
	if ((alpha = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		free(fitted);
		free(residuals);
		GCI_ecf_free_matrix(covar);
		return -3;
	}

	for (j=0; j<nparam; j++)
		paramfree_copy[j] = paramfree[j];
	paramfree_copy[spa_param1] = 0;
	paramfree_copy[spa_param2] = 0;  /* we fix the parameters we are
										analysing */

	progress = 0;
	total    = spa_nvalues1*spa_nvalues2; 
	for (i1=0; i1<spa_nvalues1; i1++) {
		for (i2=0; i2<spa_nvalues2; i2++) {
			/* Initialise parameter array each time */
			for (j=0; j<nparam; j++)
				param_copy[j] = param[j];
			/* Set the parameters we are analysing */
			param_copy[spa_param1] =
				spa_low1 + (spa_high1 - spa_low1) * (float)i1 / (float)(spa_nvalues1 - 1);
			param_copy[spa_param2] =
				spa_low2 + (spa_high2 - spa_low2) * (float)i2 / (float)(spa_nvalues2 - 1);

			ret = GCI_marquardt(x, y, ndata, noise, sig,
								param_copy, paramfree_copy, nparam, restrain,
								fitfunc, fitted, residuals, covar, alpha,
								&chisq[i1][i2], chisq_delta, 0, NULL);
			
			progress++;
            if (progressfunc)
                progressfunc ((float)progress/((float)(total-1)));
		
			if (ret < 0)
				chisq[i1][i2] = -1;
		}
	}

	free(fitted);
	free(residuals);
	GCI_ecf_free_matrix(covar);
	GCI_ecf_free_matrix(alpha);

	return 0;
}

/**
    1D Support plane analysis with instrument response. See GCI_SPA_1D_marquardt.

 */
int GCI_SPA_1D_marquardt_instr(
				float xincr, float y[],
				int ndata, int fit_start, int fit_end,
				float instr[], int ninstr,
				noise_type noise, float sig[],
				float param[], int paramfree[], int nparam,
				restrain_type restrain, float chisq_delta,
				void (*fitfunc)(float, float [], float *, float [], int),
				int spa_param, int spa_nvalues,
				float spa_low, float spa_high,
				float chisq[], float chisq_target, void (*progressfunc)(float))
{
	int i, j, ret;
	float *fitted, *residuals, **covar, **alpha;
	float param_copy[MAXFIT];
	int paramfree_copy[MAXFIT];
	
	if (spa_param < 0 || spa_param >= nparam)  /* and so nparam > 0, too */
		return -1;
	if (spa_nvalues < 2)
		return -2;
	
	if ((fitted = (float *) malloc((size_t)ndata * sizeof(float))) == NULL)
		return -3;
	if ((residuals = (float *) malloc((size_t)ndata * sizeof(float))) == NULL) {
		free(fitted);
		return -3;
	}
	if ((covar = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		free(fitted);
		free(residuals);
		return -3;
	}
	if ((alpha = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		free(fitted);
		free(residuals);
		GCI_ecf_free_matrix(covar);
		return -3;
	}

	for (j=0; j<nparam; j++)
		paramfree_copy[j] = paramfree[j];
	paramfree_copy[spa_param] = 0;  /* we fix the parameter we are
									   analysing */

	for (i=0; i<spa_nvalues; i++) {
		/* Initialise parameter array each time */
		for (j=0; j<nparam; j++)
			param_copy[j] = param[j];
		/* Set the parameter we are analysing */
		param_copy[spa_param] =
			spa_low + (spa_high - spa_low) * (float)i / (float)(spa_nvalues - 1);

		ret = GCI_marquardt_fitting_engine(xincr, y, ndata, fit_start, fit_end, instr, ninstr,
											noise, NULL, 
											param_copy, paramfree_copy, nparam, restrain,
							   				fitfunc, fitted, residuals, &chisq[i],
							   				covar, alpha, NULL, chisq_target, chisq_delta, 0);
			
		if (progressfunc)
            progressfunc ((float)i/(float)(spa_nvalues-1));
		
		if (ret < 0)
			chisq[i] = -1;
	}

	free(fitted);
	free(residuals);
	GCI_ecf_free_matrix(covar);
	GCI_ecf_free_matrix(alpha);
	GCI_marquardt_cleanup();

	return 0;
}


/**
    2D Support plane analysis with instrument response. See GCI_SPA_2D_marquardt.

 */
int GCI_SPA_2D_marquardt_instr(
				float xincr, float y[],
				int ndata, int fit_start, int fit_end,
				float instr[], int ninstr,
				noise_type noise, float sig[],
				float param[], int paramfree[], int nparam,
				restrain_type restrain, float chisq_delta,
				void (*fitfunc)(float, float [], float *, float [], int),
				int spa_param1, int spa_nvalues1,
				float spa_low1, float spa_high1,
				int spa_param2, int spa_nvalues2,
				float spa_low2, float spa_high2,
				float **chisq, float chisq_target, void (*progressfunc)(float))
{
	int i1, i2, j, ret, progress, total;
	float *fitted, *residuals, **covar, **alpha;
	float param_copy[MAXFIT];
	int paramfree_copy[MAXFIT];
	
	if (spa_param1 < 0 || spa_param1 >= nparam)  /* and so nparam > 0, too */
		return -1;
	if (spa_param2 < 0 || spa_param2 >= nparam)
		return -1;
	if (spa_param1 == spa_param2)
		return -1;
	if (spa_nvalues1 < 2 || spa_nvalues2 < 2)
		return -2;
	
	if ((fitted = (float *) malloc((size_t)ndata * sizeof(float))) == NULL)
		return -3;
	if ((residuals = (float *) malloc((size_t)ndata * sizeof(float))) == NULL) {
		free(fitted);
		return -3;
	}
	if ((covar = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		free(fitted);
		free(residuals);
		return -3;
	}
	if ((alpha = GCI_ecf_matrix(nparam, nparam)) == NULL) {
		free(fitted);
		free(residuals);
		GCI_ecf_free_matrix(covar);
		return -3;
	}

	for (j=0; j<nparam; j++)
		paramfree_copy[j] = paramfree[j];
	paramfree_copy[spa_param1] = 0;
	paramfree_copy[spa_param2] = 0;  /* we fix the parameters we are
										analysing */
	progress = 0;
	total    = spa_nvalues1*spa_nvalues2; 
	for (i1=0; i1<spa_nvalues1; i1++) {
		for (i2=0; i2<spa_nvalues2; i2++) {
			/* Initialise parameter array each time */
			for (j=0; j<nparam; j++)
				param_copy[j] = param[j];
			/* Set the parameters we are analysing */
			param_copy[spa_param1] =
				spa_low1 + (spa_high1 - spa_low1) * (float)i1 / (float)(spa_nvalues1 - 1);
			param_copy[spa_param2] =
				spa_low2 + (spa_high2 - spa_low2) * (float)i2 / (float)(spa_nvalues2 - 1);

			ret = GCI_marquardt_fitting_engine(xincr, y, ndata, fit_start, fit_end, instr, ninstr,
									noise, NULL,
									param_copy, paramfree_copy, nparam, restrain,
					   				fitfunc, fitted, residuals, &chisq[i1][i2],
					   				covar, alpha, NULL, chisq_target, chisq_delta, 0);
			
			progress++;
            if (progressfunc)
                progressfunc ((float)progress/((float)(total-1)));
		
			if (ret < 0)
				chisq[i1][i2] = -1;
	
		
		}
	}

	free(fitted);
	free(residuals);
	GCI_ecf_free_matrix(covar);
	GCI_ecf_free_matrix(alpha);
	GCI_marquardt_cleanup();

	return 0;
}


/**
    1D Global Support plane analysis with instrument response. See GCI_SPA_1D_marquardt.

 */
int GCI_SPA_1D_marquardt_global_exps_instr(
					float xincr, float **trans,
					int ndata, int ntrans, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[], int ftype,
					float **param, int paramfree[], int nparam,
					restrain_type restrain, float chisq_delta, int drop_bad_transients,
					int spa_param, int spa_nvalues,
					float spa_low, float spa_high,
					float chisq_global[], int df[], void (*progressfunc)(float))
{
	int i, j, k, ret;
	float **param_copy;
	float **fitted, **residuals, *chisq_trans;
	int paramfree_copy[MAXFIT];
	
	if (spa_param < 0 || spa_param >= nparam)  /* and so nparam > 0, too */
		return -1;
	if (spa_nvalues < 2)
		return -2;
	if (ntrans < 1)
		return -3;

	if ((param_copy = GCI_ecf_matrix(ntrans, nparam)) == NULL)
		return -4;
	if ((fitted = GCI_ecf_matrix(ntrans, ndata)) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		return -4;
	}
	if ((residuals = GCI_ecf_matrix(ntrans, ndata)) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		GCI_ecf_free_matrix(fitted);
		return -4;
	}
	if ((chisq_trans = (float *) malloc((size_t)ntrans * sizeof(float))) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		GCI_ecf_free_matrix(fitted);
		GCI_ecf_free_matrix(residuals);
		return -4;
	}

	/* We set up the paramfree array, and also count the number of
	   free parameters for the degrees of freedom calculation */
	for (j=0; j<nparam; j++)
		paramfree_copy[j] = paramfree[j];
	paramfree_copy[spa_param] = 0;  /* we fix the parameter we are
									   analysing */

	/* We only need to initialise this parameter array once at the
	   beginning, in case there are any global variables. */
	for (j=0; j<ntrans; j++) {
		for (k=0; k<nparam; k++)
			param_copy[j][k] = param[j][k];
	}

	for (i=0; i<spa_nvalues; i++) {
		/* Set the parameter we are analysing */
		for (j=0; j<ntrans; j++) {
			param_copy[j][spa_param] =
				spa_low + (spa_high - spa_low) * (float)i / (float)(spa_nvalues - 1);
		}

		ret = GCI_marquardt_global_exps_instr(
							xincr, trans, ndata, ntrans, fit_start, fit_end,
							instr, ninstr, noise, sig, ftype,
							param_copy, paramfree_copy, nparam, restrain, chisq_delta,
							fitted, residuals, chisq_trans,
							&chisq_global[i], &df[i], drop_bad_transients);
			
		if (progressfunc)
            progressfunc ((float)i/(float)(spa_nvalues-1));

		if (ret < 0)
			chisq_global[i] = -1;
	}

	GCI_ecf_free_matrix(param_copy);
	GCI_ecf_free_matrix(fitted);
	GCI_ecf_free_matrix(residuals);
	free(chisq_trans);
	GCI_marquardt_cleanup();

	return 0;
}


/**
    2D Global Support plane analysis with instrument response. See GCI_SPA_2D_marquardt.

 */
int GCI_SPA_2D_marquardt_global_exps_instr(
					float xincr, float **trans,
					int ndata, int ntrans, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[], int ftype,
					float **param, int paramfree[], int nparam,
					restrain_type restrain, float chisq_delta, int drop_bad_transients,
					int spa_param1, int spa_nvalues1,
					float spa_low1, float spa_high1,
					int spa_param2, int spa_nvalues2,
					float spa_low2, float spa_high2,
					float **chisq_global, int **df, void (*progressfunc)(float))
{
	int i1, i2, j, k, ret, progress, total;
	float **param_copy;
	float **fitted, **residuals, *chisq_trans;
	int paramfree_copy[MAXFIT];
	
	if (spa_param1 < 0 || spa_param1 >= nparam)  /* and so nparam > 0, too */
		return -1;
	if (spa_param2 < 0 || spa_param2 >= nparam)
		return -1;
	if (spa_param1 == spa_param2)
		return -1;
	if (spa_nvalues1 < 2 || spa_nvalues2 < 2)
		return -2;
	if (ntrans < 1)
		return -3;

	if ((param_copy = GCI_ecf_matrix(ntrans, nparam)) == NULL)
		return -4;
	if ((fitted = GCI_ecf_matrix(ntrans, ndata)) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		return -4;
	}
	if ((residuals = GCI_ecf_matrix(ntrans, ndata)) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		GCI_ecf_free_matrix(fitted);
		return -4;
	}
	if ((chisq_trans = (float *) malloc((size_t)ntrans * sizeof(float))) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		GCI_ecf_free_matrix(fitted);
		GCI_ecf_free_matrix(residuals);
		return -4;
	}

	/* We set up the paramfree array, and also count the number of
	   free parameters for the degrees of freedom calculation */
	for (j=0; j<nparam; j++)
		paramfree_copy[j] = paramfree[j];
	paramfree_copy[spa_param1] = 0;
	paramfree_copy[spa_param2] = 0;  /* we fix the parameters we are
										analysing */

	/* We only need to initialise this parameter array once at the
	   beginning, in case there are any global variables; all of the
	   other estimates are done by the ECF code itself */
	for (j=0; j<ntrans; j++) {
		for (k=0; k<nparam; k++)
			param_copy[j][k] = param[j][k];
	}

	progress = 0;
	total    = spa_nvalues1*spa_nvalues2; 
	for (i1=0; i1<spa_nvalues1; i1++) {
		for (i2=0; i2<spa_nvalues2; i2++) {
			/* Set the parameters we are analysing */
			for (j=0; j<ntrans; j++) {
				param_copy[j][spa_param1] = spa_low1 +
					(spa_high1 - spa_low1) * (float)i1 / (float)(spa_nvalues1 - 1);
				param_copy[j][spa_param2] = spa_low2 +
					(spa_high2 - spa_low2) * (float)i2 / (float)(spa_nvalues2 - 1);
			}

			ret = GCI_marquardt_global_exps_instr(
							xincr, trans, ndata, ntrans, fit_start, fit_end,
							instr, ninstr, noise, sig, ftype,
							param_copy, paramfree_copy, nparam, restrain, chisq_delta,
							fitted, residuals, chisq_trans,
							&chisq_global[i1][i2], &df[i1][i2],
							drop_bad_transients);
			
			progress++;
            if (progressfunc)
                progressfunc ((float)progress/((float)(total-1)));
		
			if (ret < 0)
				chisq_global[i1][i2] = -1;
		}
	}

	GCI_ecf_free_matrix(param_copy);
	GCI_ecf_free_matrix(fitted);
	GCI_ecf_free_matrix(residuals);
	free(chisq_trans);
	GCI_marquardt_cleanup();

	return 0;
}


/**
    1D Generic Global Support plane analysis with instrument response. See GCI_SPA_1D_marquardt.

 */
int GCI_SPA_1D_marquardt_global_generic_instr(
					float xincr, float **trans,
					int ndata, int ntrans, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[],
					float **param, int paramfree[], int nparam, int gparam[],
					restrain_type restrain, float chisq_delta,
					void (*fitfunc)(float, float [], float *, float [], int),
					int spa_param, int spa_nvalues,
					float spa_low, float spa_high,
					float chisq_global[], int df[], void (*progressfunc)(float))
{
	int i, j, k, ret;
	float **param_copy;
	float **fitted, **residuals, *chisq_trans;
	int paramfree_copy[MAXFIT];
	
	if (spa_param < 0 || spa_param >= nparam)  /* and so nparam > 0, too */
		return -1;
	if (spa_nvalues < 2)
		return -2;
	if (ntrans < 1)
		return -3;

	if ((param_copy = GCI_ecf_matrix(ntrans, nparam)) == NULL)
		return -4;
	if ((fitted = GCI_ecf_matrix(ntrans, ndata)) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		return -4;
	}
	if ((residuals = GCI_ecf_matrix(ntrans, ndata)) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		GCI_ecf_free_matrix(fitted);
		return -4;
	}
	if ((chisq_trans = (float *) malloc((size_t)ntrans * sizeof(float))) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		GCI_ecf_free_matrix(fitted);
		GCI_ecf_free_matrix(residuals);
		return -4;
	}

	/* We set up the paramfree array, and also count the number of
	   free parameters for the degrees of freedom calculation */
	for (j=0; j<nparam; j++)
		paramfree_copy[j] = paramfree[j];
	paramfree_copy[spa_param] = 0;  /* we fix the parameter we are
									   analysing */

	for (i=0; i<spa_nvalues; i++) {
		/* Initialise parameter array each time */
		for (j=0; j<ntrans; j++) {
			for (k=0; k<nparam; k++)
				param_copy[j][k] = param[j][k];
			/* Set the parameter we are analysing */
			param_copy[j][spa_param] =
				spa_low + (spa_high - spa_low) * (float)i / (float)(spa_nvalues - 1);
		}

		ret = GCI_marquardt_global_generic_instr(
						xincr, trans, ndata, ntrans, fit_start, fit_end,
						instr, ninstr, noise, sig,
						param_copy, paramfree_copy, nparam, gparam, restrain, chisq_delta,
						fitfunc, fitted, residuals, chisq_trans,
						&chisq_global[i], &df[i]);
			
		if (progressfunc)
            progressfunc ((float)i/(float)(spa_nvalues-1));

		if (ret < 0)
			chisq_global[i] = -1;
	}

	GCI_ecf_free_matrix(param_copy);
	GCI_ecf_free_matrix(fitted);
	GCI_ecf_free_matrix(residuals);
	free(chisq_trans);
	GCI_marquardt_cleanup();

	return 0;
}


/**
    2D Support plane analysis with instrument response. See GCI_SPA_2D_marquardt.

 */
int GCI_SPA_2D_marquardt_global_generic_instr(
					float xincr, float **trans,
					int ndata, int ntrans, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[],
					float **param, int paramfree[], int nparam, int gparam[],
					restrain_type restrain, float chisq_delta,
					void (*fitfunc)(float, float [], float *, float [], int),
					int spa_param1, int spa_nvalues1,
					float spa_low1, float spa_high1,
					int spa_param2, int spa_nvalues2,
					float spa_low2, float spa_high2,
					float **chisq_global, int **df, void (*progressfunc)(float))
{
	int i1, i2, j, k, ret, progress, total;
	float **param_copy;
	float **fitted, **residuals, *chisq_trans;
	int paramfree_copy[MAXFIT];
	
	if (spa_param1 < 0 || spa_param1 >= nparam)  /* and so nparam > 0, too */
		return -1;
	if (spa_param2 < 0 || spa_param2 >= nparam)
		return -1;
	if (spa_param1 == spa_param2)
		return -1;
	if (spa_nvalues1 < 2 || spa_nvalues2 < 2)
		return -2;
	if (ntrans < 1)
		return -3;

	if ((param_copy = GCI_ecf_matrix(ntrans, nparam)) == NULL)
		return -4;
	if ((fitted = GCI_ecf_matrix(ntrans, ndata)) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		return -4;
	}
	if ((residuals = GCI_ecf_matrix(ntrans, ndata)) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		GCI_ecf_free_matrix(fitted);
		return -4;
	}
	if ((chisq_trans = (float *) malloc((size_t)ntrans * sizeof(float))) == NULL) {
		GCI_ecf_free_matrix(param_copy);
		GCI_ecf_free_matrix(fitted);
		GCI_ecf_free_matrix(residuals);
		return -4;
	}

	/* We set up the paramfree array, and also count the number of
	   free parameters for the degrees of freedom calculation */
	for (j=0; j<nparam; j++)
		paramfree_copy[j] = paramfree[j];
	paramfree_copy[spa_param1] = 0;
	paramfree_copy[spa_param2] = 0;  /* we fix the parameters we are
										analysing */

	progress = 0;
	total    = spa_nvalues1*spa_nvalues2; 
	for (i1=0; i1<spa_nvalues1; i1++) {
		for (i2=0; i2<spa_nvalues2; i2++) {
			/* Initialise parameter array each time */
			for (j=0; j<ntrans; j++) {
				for (k=0; k<nparam; k++)
					param_copy[j][k] = param[j][k];
				/* Set the parameters we are analysing */
				param_copy[j][spa_param1] = spa_low1 +
					(spa_high1 - spa_low1) * (float)i1 / (float)(spa_nvalues1 - 1);
				param_copy[j][spa_param2] = spa_low2 +
					(spa_high2 - spa_low2) * (float)i2 / (float)(spa_nvalues2 - 1);
			}

			ret = GCI_marquardt_global_generic_instr(
						xincr, trans, ndata, ntrans, fit_start, fit_end,
						instr, ninstr, noise, sig,
						param_copy, paramfree_copy, nparam, gparam, restrain, chisq_delta,
						fitfunc, fitted, residuals, chisq_trans,
						&chisq_global[i1][i2], &df[i1][i2]);
			
			progress++;
            if (progressfunc)
                progressfunc ((float)progress/((float)(total-1)));
		
			if (ret < 0)
				chisq_global[i1][i2] = -1;
		}
	}

	GCI_ecf_free_matrix(param_copy);
	GCI_ecf_free_matrix(fitted);
	GCI_ecf_free_matrix(residuals);
	free(chisq_trans);
	GCI_marquardt_cleanup();

	return 0;
}


// Emacs settings:
// Local variables:
// mode: c
// c-basic-offset: 4
// tab-width: 4
// End:

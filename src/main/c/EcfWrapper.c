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

#include "EcfWrapper.h"
#include "Ecf.h"
#include "EcfInternal.h"
#include <stdio.h>
#include <stdlib.h>


noise_type getNoiseModels (int val)
{  // This is from TRI2 and must match the noise models, between TRI2 UI and FLIMLib - This is not a good thing to have! TRI2 should save the FLIMLib type

    switch (val)
    {
        case 1:
            return NOISE_GAUSSIAN_FIT;
            break;
        case 2:
            return NOISE_POISSON_DATA;
            break;
        case 3:
            return NOISE_POISSON_FIT;
            break;
        case 4:
            return NOISE_MLE;
            break;
        default:
            return NOISE_POISSON_FIT;
            break;
    }
}

/*
 * Rapid Lifetime Determination fit method.
 * 
 * Quickly estimates single exponential curve parameters fitted to decay curve.
 * 
 * See further documentation in EcfSingle.c.
 * 
 * x_inc time increment of decay histogram bins
 * y decay histogram values
 * fit_start start of data to be fitted on decay
 * fit_end end of data to be fitted on decay
 * instr instrument response decay curve (prompt, excitation, IRF, lamp function)
 * n_instr number of values in instrument response decay curve
 * noise noise model, see enum in Ecf.h
 * sig sigma array for decay values
 * z background fitted parameter
 * a amplitude fitted parameter
 * tau lifetime fitted parameter
 * fitted fitted decay values
 * chi_square fitted reduced chi square
 * chi_square_target targeted chi square
 */
int RLD_fit(
		double x_inc,
		double y[],
		int fit_start,
		int fit_end,
		double instr[],
		int n_instr,
		int noise,
		double sig[],
		double *z,
		double *a,
		double *tau,
		double fitted[],
		double *chi_square,
		double chi_square_target
		) {
	noise_type nt_noise = getNoiseModels(noise);
	int n_data = fit_end + 1;
	float x_inc_float = (float) x_inc;
	float *y_float = (float *)malloc((size_t) n_data * sizeof(float));
	float *sig_float = (float *)malloc((size_t) n_data * sizeof(float));
	float *instr_float = 0;
	float *fitted_float = 0;
	float *residuals_float = (float *)malloc((size_t) n_data * sizeof(float));
	float z_float          = (float) *z;
	float a_float          = (float) *a;
	float tau_float        = (float) *tau;
	float chi_square_float = (float) *chi_square;
	float chi_square_target_float = (float) chi_square_target;
	int return_value;

	int i;
	for (i = 0; i < n_data; i++) {
		y_float[i] = (float) y[i];
		sig_float[i] = ( sig ? (float) sig[i] : 1.0f );
	}

	if (instr) {
		instr_float = (float *)malloc((size_t) n_instr * sizeof(float));
		for (i = 0; i < n_instr; ++i) {
			instr_float[i] = (float) instr[i];
		}
	}

	if (fitted) {
		fitted_float = (float*)malloc((size_t) n_data * sizeof(float));
	}
	
	return_value =  GCI_triple_integral_fitting_engine(
			x_inc_float,
			y_float,
			fit_start,
			fit_end,
			instr_float,
			n_instr,
			nt_noise,
			sig_float,
			&z_float,
			&a_float,
			&tau_float,
			fitted_float,
			residuals_float,
			&chi_square_float,
			chi_square_target_float
			);
	*z          = (double) z_float;
	*a          = (double) a_float;
	*tau        = (double) tau_float;
	*chi_square = (double) chi_square_float;
   	if (y_float) 
		free(y_float);
	if (sig_float)
		free(sig_float);
	if (instr_float) {
		free(instr_float);
	}
	if (fitted_float) {
		free(fitted_float);
	}
	if (residuals_float)
		free(residuals_float);
	return return_value;
}

/*
 * Levenberg-Marquardt algorithm fit method.
 * 
 * Fits to decay curve.  Iteratively refines from an initial estimate.
 * 
 * See further documentation in EcfSingle.c.
 * 
 * x_inc time increment of decay histogram bins
 * y decay histogram values
 * fit_start start of data to be fitted on decay
 * fit_end end of data to be fitted on decay
 * instr instrument response decay curve (prompt, excitation, IRF, lamp function)
 * n_instr number of values in instrument response decay curve
 * noise noise model, see enum in Ecf.h
 * sig sigma array for decay values
 * param source/destination for fitted parameters
 * param_free indicates which parameters should be fixed vs. fitted
 * chi_square fitted reduced chi square
 * chi_square_target targeted chi square
 * chi_square_delta defines a trivial improvement in chi square
 */
int LMA_fit(
		double x_inc,
		double y[],
		int fit_start,
		int fit_end,
		double instr[],
		int n_instr,
		int noise,
		double sig[],
		double param[],
		int param_free[],
		int n_param,
		double fitted[],
		double *chi_square,
		double chi_square_target,
		double chi_square_delta
		) {

	noise_type nt_noise = getNoiseModels(noise);
	restrain_type restrain = ECF_RESTRAIN_DEFAULT;
	int chi_square_percent = 95;
	int n_data = fit_end + 1; 
	float *y_float = (float *)malloc((size_t) n_data * sizeof(float));
	float *sig_float = (float *)malloc((size_t) n_data * sizeof(float));
	float *instr_float = 0;
	float *fitted_float = 0;
	float *param_float;
	float *residuals_float = (float *)malloc((size_t) n_data * sizeof(float));
	float **covar_float    = GCI_ecf_matrix(n_data, n_data);
	float **alpha_float    = GCI_ecf_matrix(n_data, n_data);
	float **err_axes_float = GCI_ecf_matrix(n_data, n_data);
	float chi_square_float = (float) *chi_square;
	void (*fitfunc)(float, float [], float *, float[], int) = NULL;
	int return_value;
	int i;
	for (i = 0; i < n_data; ++i) {
		y_float[i] = (float) y[i];
		sig_float[i] = ( sig ? (float) sig[i] : 1.0f ); 
	}
	if (instr) {
		int i;
		instr_float = (float *)malloc((size_t) n_instr * sizeof(float));
		for (i = 0; i < n_instr; ++i) {
			instr_float[i] = (float) instr[i];
		}
	}

	if (fitted) {
		fitted_float = (float*)malloc((size_t) n_data * sizeof(float)); 
	}
	param_float = (float *)malloc((size_t) n_param * sizeof(float));
	switch (n_param) {
		case 3:
			// single exponential fit
			param_float[0] = (float) param[1]; // z
			param_float[1] = (float) param[2]; // a
			param_float[2] = (float) param[3]; // tau
			break;
		case 4: 
			// stretched exponential fit
			param_float[0] = (float) param[1]; // z
			param_float[1] = (float) param[2]; // a
			param_float[2] = (float) param[3]; // tau
			param_float[3] = (float) param[4]; // h
			break;
		case 5:
			// double exponential fit
			param_float[0] = (float) param[1]; // z
			param_float[1] = (float) param[2]; // a1
			param_float[2] = (float) param[3]; // tau1
			param_float[3] = (float) param[4]; // a2
			param_float[4] = (float) param[5]; // tau2
			break;
		case 7:
			// triple exponential fit
			param_float[0] = (float) param[1]; // z
			param_float[1] = (float) param[2]; // a1
			param_float[2] = (float) param[3]; // tau1
			param_float[3] = (float) param[4]; // a2
			param_float[4] = (float) param[5]; // tau2
			param_float[5] = (float) param[6]; // a3
			param_float[6] = (float) param[7]; // tau3
			break;
		default:
			break;
	}

	// choose appropriate fitting function
	if (4 == n_param) {
		fitfunc = GCI_stretchedexp;
	}
	else {
		fitfunc = GCI_multiexp_tau;
	}

	return_value = GCI_marquardt_fitting_engine(
			(float) x_inc,
			y_float,
			n_data,
			fit_start,
			fit_end,
			instr_float,
			n_instr,
			nt_noise,
			sig_float,
			param_float,
			param_free,
			n_param,
			restrain,
			fitfunc,
			fitted_float,
			residuals_float,
			&chi_square_float,
			covar_float,
			alpha_float,
			err_axes_float,
			(float) chi_square_target,
			(float) chi_square_delta,
			chi_square_percent);
	*chi_square = (double) chi_square_float;

	if (y_float)
		free(y_float);
	if (sig_float)
		free(sig_float);
	if (instr_float) {
		free(instr_float);
	}
	if (fitted_float) {
		free(fitted_float);
	}
	
	switch (n_param) {
		case 3:
			// single exponential fit
			param[0] = (double) chi_square_float;
			param[1] = (double) param_float[0]; // z
			param[2] = (double) param_float[1]; // a
			param[3] = (double) param_float[2]; // tau
			break;
		case 4:
			// stretched exponential fit
			
			// convert tau into <tau>
			param_float[2]
				= param_float[2] * param_float[3] * GCI_gamma(param_float[3]);
			
			param[0] = (double) chi_square_float;
			param[1] = (double) param_float[0]; // z
			param[2] = (double) param_float[1]; // a
			param[3] = (double) param_float[2]; // tau
			param[4] = (double) param_float[3]; // h
			break;
		case 5:
			// double exponential fit
			param[0] = (double) chi_square_float;
			param[1] = (double) param_float[0]; // z
			param[2] = (double) param_float[1]; // a1
			param[3] = (double) param_float[2]; // tau1
			param[4] = (double) param_float[3]; // a2
			param[5] = (double) param_float[4]; // tau2
			break;
		case 7:
			// triple exponential fit
			param[0] = (double) chi_square_float;
			param[1] = (double) param_float[0]; // z
			param[2] = (double) param_float[1]; // a1
			param[3] = (double) param_float[2]; // tau1
			param[4] = (double) param_float[3]; // a2
			param[5] = (double) param_float[4]; // tau2
			param[6] = (double) param_float[5]; // a3
			param[7] = (double) param_float[6]; // tau3
			break;
	}
	if (param_float)
		free(param_float);
	if (residuals_float)
		free(residuals_float);
	GCI_ecf_free_matrix(covar_float);
	GCI_ecf_free_matrix(alpha_float);
	GCI_ecf_free_matrix(err_axes_float);

	return return_value;
}


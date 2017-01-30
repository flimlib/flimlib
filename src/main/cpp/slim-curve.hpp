/*
 * #%L
 * SLIM Curve package for exponential curve fitting of spectral lifetime data.
 * %%
 * Copyright (C) 2010 - 2017 University of Oxford and Board of Regents of the
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

 /* 
  C++ interface for SLIM Curve.
*/

#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* calloc, exit, free */

#include "Ecf.h"
#include "GCI_Phasor.h"

#define SLIM_CURVE_SUCCESS 0
#define SLIM_CURVE_MEMORY_ERROR -1
#define SLIM_CURVE_PARAMETER_ERROR -2
#define SLIM_CURVE_FITTING_ERROR -3

// Set default values
float chisq_target = 1.0f;

class SLIMCurve {
public:
	float xincr;		///< The time increment inbetween the values in the y array.
	float *transient;   ///< The transient (time resolved) signal to be analysed, the 'data'.
	int ndata;          ///< The number of data points.
	int fit_start;      ///< The index into the y array marking the start to the data to be used in the fit.
	int fit_end;        ///< The index into the y array marking the end of the data to be used in the fit.
	float *instr=NULL;  ///< The instrument reponse (IRF) or prompt signal to be used (optional, can be NULL).
	int ninstr=0;       ///< The number of data points in the prompt (ignored if prompt = NULL).
	noise_type noise = NOISE_POISSON_FIT;   ///< The #noise_type to be used.
	float *sig;         ///< The standard deviation at each data point in y if #noise_type NOISE_GIVEN is used (optional, can pass NULL).
	float *param;       ///<
	int *paramfree;     ///< An array indicating which parameters are free (1), fixed (0)
	int nparam;         ///< The number of parameters.
	restrain_type restrain; ///< Parameter #restrain_type.Normally use ECF_RESTRAIN_DEFAULT.Use ECF_RESTRAIN_USER if restraining parameters has been setup via GCI_set_restrain_limits.
	float *fitted;      ///< An array containing values fitted to the data, the 'fit'. Fit points are coincident in time with the data points.
	float *residuals;   ///<
	float *chisq;       ///< The resulting raw chi squared value of the fit. To get the reduced chisq, divide by the degrees of freedom (fit_start - fit_end - nparam)
	float chisq_target=1.0f; ///< A raw chi squared value to aim for. If this value is reached fitting will stop. If you want to aim for a reduced chisq (say 1.1 or 1.0) you must multiply by the degree of freedom. (TRI2: "Try refits")
	float chisq_delta=0.000001f;  ///< An individual fit will continue if the chisq value changes by more then this amount.Try 1E-5. (TRI2: "Stopping Criterion")
	int chisq_percent=95;  ///< Defines the confidence interval when calculating the error axes, e.g. 95 % .
	float **covar;      ///< covar The covariance matrix. Allocate with a square matrix with #GCI_ecf_matrix(nparam, nparam).
	float **alpha;      ///< The alpha matrix.Allocate with a square matrix with #GCI_ecf_matrix(nparam, nparam).
	float **err_axes;   ///< The dimensions of the confidence ellipsoid of the chisq.See chisq_percent below.Allocate with a square matrix with #GCI_ecf_matrix(nparam, nparam).
	int iterations;     ///< Store of the number of iterations performed.

	void(*fitfunc)(float, float[], float *, float[], int);

	SLIMCurve() {}

	~SLIMCurve(){}

	int fitRLD() {

		float *Z = &(param[0]);
		float *A = &(param[1]);
		float *tau = &(param[2]);

		iterations = GCI_triple_integral_fitting_engine(xincr, transient, fit_start, fit_end,
			instr, ninstr, noise, sig, Z, A, tau, fitted, residuals, chisq, chisq_target);

		if (iterations >= 0) {
			return SLIM_CURVE_SUCCESS;
		}
		else {
			return SLIM_CURVE_FITTING_ERROR;
		}
	}



	int fitLMA() {

		iterations = GCI_marquardt_fitting_engine(xincr, transient, ndata, fit_start, fit_end,
			instr, ninstr, noise, sig, param, paramfree, nparam, restrain, fitfunc,
			fitted, residuals, chisq, covar, alpha, err_axes, chisq_target, chisq_delta, chisq_percent);

		if (iterations >= 0) {
			return SLIM_CURVE_SUCCESS;
		}
		else {
			return SLIM_CURVE_FITTING_ERROR;
		}
	}

  
};


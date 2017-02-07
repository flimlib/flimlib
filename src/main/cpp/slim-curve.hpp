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

  Caller _has_ to provide:
	transient array
	ndata
	param array

*/

#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* calloc, exit, free */

#include "Ecf.h"
#include "GCI_Phasor.h"

#define SLIM_CURVE_SUCCESS 0
#define SLIM_CURVE_MEMORY_ERROR -1
#define SLIM_CURVE_SETTINGS_ERROR -2
#define SLIM_CURVE_FITTING_ERROR -3

#define SLIM_CURVE_MONO 0
#define SLIM_CURVE_BI 1
#define SLIM_CURVE_TRI 2
#define SLIM_CURVE_STRETCHED 3


// Set default values
float chisq_target = 1.0f;

class SLIMCurve {
	// Local stores in case user does not provide them
	int * _paramfree = NULL;
	int _nparamfree = 0;
	float *_fitted = NULL;
	float *_residuals = NULL;
	float _chisq = 0.0f;
	float **_covar = NULL;
	float **_alpha = NULL;
	float **_err_axes = NULL;

public:
	float xincr;		///< The time increment inbetween the values in the y array.
	float *transient;   ///< The transient (time resolved) signal to be analysed, the 'data'.
	int ndata;          ///< The number of data points.
	int data_start;     ///< The index into the y array marking the start to the transient.
	int fit_start;      ///< The index into the y array marking the start to the data to be used in the fit. Data between data_start and fit_start is used for convolution with the IRF.
	int fit_end;        ///< The index into the y array marking the end of the data to be used in the fit.
	float *instr;       ///< The instrument reponse (IRF) or prompt signal to be used (optional, can be NULL).
	int ninstr;         ///< The number of data points in the prompt (ignored if prompt = NULL).
	noise_type noise;   ///< The #noise_type to be used.
	float *sig;         ///< The standard deviation at each data point in y if #noise_type NOISE_GIVEN is used (optional, can pass NULL).
	float *param;       ///<
	int *paramfree;     ///< An array indicating which parameters are free (1), fixed (0)
	int nparam;         ///< The number of parameters.
	restrain_type restrain; ///< Parameter #restrain_type.Normally use ECF_RESTRAIN_DEFAULT.Use ECF_RESTRAIN_USER if restraining parameters has been setup via GCI_set_restrain_limits.
	float *fitted;      ///< An array containing values fitted to the data, the 'fit'. Fit points are coincident in time with the data points.
	float *residuals;   ///< An array containing the difference between the fit and the data. Points are coincident in time with the data points.
	float *chisq;       ///< The resulting raw chi squared value of the fit. To get the reduced chisq, divide by the degrees of freedom (fit_end - fit_start - nparam), see getReducedChiSq()
	float chisq_target; ///< A raw chi squared value to aim for. If this value is reached fitting will stop. If you want to aim for a reduced chisq (say 1.1 or 1.0) you must multiply by the degree of freedom. (TRI2: "Try refits")
	float chisq_delta;  ///< An individual fit will continue if the chisq value changes by more then this amount.Try 1E-5. (TRI2: "Stopping Criterion")
	int chisq_percent;  ///< Defines the confidence interval when calculating the error axes, e.g. 95 % .
	float **covar;      ///< covar The covariance matrix. Allocate with a square matrix with #GCI_ecf_matrix(nparam, nparam).
	float **alpha;      ///< The alpha matrix.Allocate with a square matrix with #GCI_ecf_matrix(nparam, nparam).
	float **err_axes;   ///< The dimensions of the confidence ellipsoid of the chisq.See chisq_percent below.Allocate with a square matrix with #GCI_ecf_matrix(nparam, nparam).
	int iterations;     ///< Store of the number of iterations performed.

	void(*fitfunc)(float x, float param[], float *y, float dy_dparam[], int nparam);  ///< A function that executes the fitting function at pt x, given the parameters

	SLIMCurve() {
		// Default values
		xincr = 1.0;		
		transient = NULL;
		ndata = 0;          
		fit_start = 0;      
		fit_end = 0;        
		instr = NULL;
		ninstr = 0;         
		noise = NOISE_POISSON_FIT;
		sig = NULL;
		param = NULL;
		paramfree = NULL;
		nparam = 0;         
		restrain = ECF_RESTRAIN_DEFAULT;
		fitted = NULL;  
		residuals = NULL;
		chisq = NULL;
		chisq_target = 1.0;
		chisq_delta = 0.0000001f;
		chisq_percent = 95;  
		covar = NULL;
		alpha = NULL;
		err_axes = NULL;   
		iterations = 0; 
		fitfunc = GCI_multiexp_tau;
	}

	~SLIMCurve(){
		freePrivateVars();  // Should be done by fit functions, but just make sure.
	}

	int checkValues() {
		if (transient == NULL) return SLIM_CURVE_SETTINGS_ERROR;
		if (ndata == 0) return SLIM_CURVE_SETTINGS_ERROR;

		if (fit_start < 0) fit_start = 0;
		if (fit_end <= fit_start) fit_end = ndata-1;
		if (fit_end > ndata-1) return SLIM_CURVE_SETTINGS_ERROR;

		if (noise == NOISE_GIVEN && sig == NULL) return SLIM_CURVE_SETTINGS_ERROR;

		if (paramfree == NULL) {
			_paramfree = (int *)malloc(nparam * sizeof(int));
			if (_paramfree == NULL) return SLIM_CURVE_MEMORY_ERROR;
			for (int i = 0; i < nparam; i++) {
				_paramfree[i] = 1;
			}
			paramfree = _paramfree;
			_nparamfree = nparam;
		}
		else {
			for (int i = 0, _nparamfree = 0; i < nparam; i++) {
				if (paramfree[i]) _nparamfree++;
			}
		}

		if (fitted == NULL) {
			_fitted = (float *)malloc(ndata * sizeof(int));
			if (_fitted == NULL) return SLIM_CURVE_MEMORY_ERROR;
			fitted = _fitted;
		}

		if (residuals == NULL) {
			_residuals = (float *)malloc(ndata * sizeof(int));
			if (_residuals == NULL) return SLIM_CURVE_MEMORY_ERROR;
			residuals = _residuals;
		}

		if (chisq == NULL) chisq = &_chisq;

		if (covar == NULL) {
			_covar = GCI_ecf_matrix(ndata, ndata);
			if (_covar == NULL) return SLIM_CURVE_MEMORY_ERROR;
			covar = _covar;
		}

		if (alpha == NULL) {
			_alpha = GCI_ecf_matrix(ndata, ndata);
			if (_alpha == NULL) return SLIM_CURVE_MEMORY_ERROR;
			alpha = _alpha;
		}

		if (err_axes == NULL) {
			_err_axes = GCI_ecf_matrix(ndata, ndata);
			if (_err_axes == NULL) return SLIM_CURVE_MEMORY_ERROR;
			err_axes = _err_axes;
		}

		return SLIM_CURVE_SUCCESS;
	}

	void freePrivateVars() {
		if (_paramfree) {
			free(_paramfree);
			_paramfree = paramfree = NULL;
		}

		if (_fitted) {
			free(_fitted);
			_fitted = fitted = NULL;
		}

		if (_residuals) {
			free(_residuals);
			_residuals = residuals = NULL;
		}

		if (_covar) {
			GCI_ecf_free_matrix(_covar);
			_covar = covar = NULL;
		}

		if (_alpha) {
			GCI_ecf_free_matrix(_alpha);
			_alpha = alpha = NULL;
		}

		if (_err_axes) {
			GCI_ecf_free_matrix(_err_axes);
			_err_axes = err_axes = NULL;
		}
	}

	int fitRLD() {

		float *Z = &(param[0]);
		float *A = &(param[1]);
		float *tau = &(param[2]);

		int err = checkValues();
		if (err < 0) return err;

		iterations = GCI_triple_integral_fitting_engine(xincr, &transient[data_start], fit_start - data_start, fit_end - data_start,
			instr, ninstr, noise, sig, Z, A, tau, fitted, residuals, chisq, chisq_target*(fit_end - fit_start + ninstr - 3));

		freePrivateVars();

		return iterations;
	}

	int fitLMA(int type) {

		switch (type)
		{
		case SLIM_CURVE_MONO:
			nparam = 3;
			fitfunc = GCI_multiexp_tau;
			break;
		case SLIM_CURVE_BI:
			nparam = 5;
			fitfunc = GCI_multiexp_tau;
			break;
		case SLIM_CURVE_TRI:
			nparam = 7;
			fitfunc = GCI_multiexp_tau;
			break;
		case SLIM_CURVE_STRETCHED:
			fitfunc = GCI_stretchedexp;
			nparam = 4;
			break;
		default:
			break;
		}

		int err = checkValues();
		if (err < 0) return err;

		iterations = GCI_marquardt_fitting_engine(xincr, &transient[data_start], ndata - data_start, fit_start-data_start, fit_end-data_start,
			instr, ninstr, noise, sig, param, paramfree, nparam, restrain, fitfunc,
			fitted, residuals, chisq, covar, alpha, err_axes, chisq_target*(fit_end - fit_start + ninstr - _nparamfree), chisq_delta, chisq_percent);

		freePrivateVars();

		return iterations;
	}

	float getReducedChiSq(void) {
		return (*chisq / (fit_end - fit_start - nparam));
	}

};


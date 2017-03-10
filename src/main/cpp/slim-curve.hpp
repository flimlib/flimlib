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

 /** 
  C++ interface for SLIM Curve.

  Caller _has_ to at least use:
	SlimCurve.setupData(transient, ndata);
  but better is:
	SlimCurve.setupData(transient, ndata, xincr, transient_rise, transient_peak, transient_end);

  And should read results using e.g.:
	SlimCurve.param[SLIM_CURVE_BI_PARAM_Z]
	SlimCurve.param[SLIM_CURVE_BI_PARAM_A1]
	SlimCurve.param[SLIM_CURVE_BI_PARAM_TAU1]
	etc.

  Optionally can use: 
	SlimCurve.setupIRF(instr, ninstr);
	SlimCurve.restrainParameter(SLIM_CURVE_X_PARAM_Y, min, max);
	SlimCurve.fixParameter(SLIM_CURVE_X_PARAM_Y, value);

  Can also optionally get the fit and residuals; before fitting do this:
	float fitted[ndata];
	float residuals[ndata];
	SlimCurve.fitted = fitted;
	SlimCurve.residuals = residuals;
  and use fitted[] and residuals[] after the fit.
  Similarly for:
	SlimCurve.covar
	SlimCurve.alpha
	SlimCurve.err_axes

  A good sequence would be to use:
    SlimCurve.setupData(transient, ndata, xincr, transient_rise, transient_peak, transient_end);
	SlimCurve.setupIRF(instr, ninstr);
	SlimCurve.fitRLD();    // To get initial parameter estimates
	SlimCurve.noise_model = NOISE_MLE;  // To use maximum likelihood
	SlimCurve.fitLMA(SLIM_CURVE_MONO);  // To do a mono fit

	\file slim-curve.hpp
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

// parameter order as defined by RLD()
#define SLIM_CURVE_RLD_PARAM_Z 0
#define SLIM_CURVE_RLD_PARAM_A 1
#define SLIM_CURVE_RLD_PARAM_TAU 2

// parameter order as defined by GCI_multiexp_tau()
#define SLIM_CURVE_MONO_PARAM_Z 0
#define SLIM_CURVE_MONO_PARAM_A 1
#define SLIM_CURVE_MONO_PARAM_TAU 2
#define SLIM_CURVE_BI_PARAM_Z 0
#define SLIM_CURVE_BI_PARAM_A1 1
#define SLIM_CURVE_BI_PARAM_TAU1 2
#define SLIM_CURVE_BI_PARAM_A2 3
#define SLIM_CURVE_BI_PARAM_TAU2 4
#define SLIM_CURVE_TRI_PARAM_Z 0
#define SLIM_CURVE_TRI_PARAM_A1 1
#define SLIM_CURVE_TRI_PARAM_TAU1 2
#define SLIM_CURVE_TRI_PARAM_A2 3
#define SLIM_CURVE_TRI_PARAM_TAU2 4
#define SLIM_CURVE_TRI_PARAM_A3 5
#define SLIM_CURVE_TRI_PARAM_TAU3 6

// parameter order as defined by GCI_stretchedexp()
#define SLIM_CURVE_STRETCHED_PARAM_Z 0
#define SLIM_CURVE_STRETCHED_PARAM_A 1
#define SLIM_CURVE_STRETCHED_PARAM_TAU 2
#define SLIM_CURVE_STRETCHED_PARAM_H 3

// parameter order as defined by GCI_Phasor()
#define SLIM_CURVE_PHASOR_PARAM_Z 0
#define SLIM_CURVE_PHASOR_PARAM_TAU 2
#define SLIM_CURVE_PHASOR_PARAM_U 7
#define SLIM_CURVE_PHASOR_PARAM_V 8
#define SLIM_CURVE_PHASOR_PARAM_TAUP 9
#define SLIM_CURVE_PHASOR_PARAM_TAUM 10

// Set default values
float chisq_target = 1.0f;

/**
	Global class for a SLIM Curve fitter
*/
class SLIMCurve {
	// Fixing and restraining
	int _nparamfree;
	int _restrained[MAXFIT];
	float _restrained_min[MAXFIT];
	float _restrained_max[MAXFIT];

	// Local stores in case user does not provide them
	float *_fitted;
	float *_residuals;
	float **_covar;
	float **_alpha;
	float **_err_axes;

public:
	float time_incr;		///< The time increment inbetween the values in the y array.
	float *transient;   ///< The transient (time resolved) signal to be analysed, the 'data'.
	int ndata;          ///< The number of data points.
	int data_start;     ///< The index into the y array marking the start to the transient.
	int fit_start;      ///< The index into the y array marking the start to the data to be used in the fit. Data between data_start and fit_start is used for convolution with the IRF.
	int fit_end;        ///< The index into the y array marking the end of the data to be used in the fit.
	float *instr;       ///< The instrument reponse (IRF) or prompt signal to be used (optional, can be NULL).
	int ninstr;         ///< The number of data points in the prompt (ignored if prompt = NULL).
	noise_type noise_model;   ///< The #noise_type to be used.
	float *noise_sd;    ///< The standard deviation at each data point in y if #noise_type NOISE_GIVEN is used (optional, can pass NULL).
	float param[MAXFIT];  ///< The fitted parameters, also allows starting values to be passed.
	int paramfree[MAXFIT]; ///< An array indicating which parameters are free (1), fixed (0)
	int nparam;         ///< The number of parameters.
	restrain_type restrain; ///< Parameter #restrain_type.Normally use ECF_RESTRAIN_DEFAULT.Use ECF_RESTRAIN_USER if restraining parameters has been setup via GCI_set_restrain_limits.
	float *fitted;      ///< An array containing values fitted to the data, the 'fit'. Fit points are coincident in time with the data points.
	float *residuals;   ///< An array containing the difference between the fit and the data. Points are coincident in time with the data points.
	float chisq;       ///< The resulting raw chi squared value of the fit. To get the reduced chisq, divide by the degrees of freedom (fit_end - fit_start - nparam), see getReducedChiSq()
	float chisq_target; ///< A raw chi squared value to aim for. If this value is reached fitting will stop. If you want to aim for a reduced chisq (say 1.1 or 1.0) you must multiply by the degree of freedom. (TRI2: "Try refits")
	float chisq_delta;  ///< An individual fit will continue if the chisq value changes by more then this amount.Try 1E-5. (TRI2: "Stopping Criterion")
	int chisq_percent;  ///< Defines the confidence interval when calculating the error axes, e.g. 95 % .
	float **covar;      ///< covar The covariance matrix. Allocate with a square matrix with #GCI_ecf_matrix(nparam, nparam).
	float **alpha;      ///< The alpha matrix.Allocate with a square matrix with #GCI_ecf_matrix(nparam, nparam).
	float **err_axes;   ///< The dimensions of the confidence ellipsoid of the chisq.See chisq_percent below.Allocate with a square matrix with #GCI_ecf_matrix(nparam, nparam).
	int iterations;     ///< Store of the number of iterations performed.

	void(*fitfunc)(float x, float param[], float *y, float dy_dparam[], int nparam);  ///< A function that executes the fitting function at pt x, given the parameters

	/**
	 * Constructor: make sure all vars are assigned.
	 */
	SLIMCurve() {
		// Default values
		_nparamfree = 0;
		_fitted = NULL;
		_residuals = NULL;
		_covar = NULL;
		_alpha = NULL;
		_err_axes = NULL;
		time_incr = 1.0;		
		transient = NULL;
		ndata = 0;          
		data_start = 0;
		fit_start = 0;
		fit_end = 0;
		instr = NULL;
		ninstr = 0;         
		noise_model = NOISE_POISSON_FIT;
		noise_sd = NULL;
		memset(param, 0, MAXFIT * sizeof(int));
		memset(paramfree, 1, MAXFIT * sizeof(int));
		nparam = 0;
		restrain = ECF_RESTRAIN_DEFAULT;
		fitted = NULL;  
		residuals = NULL;
		chisq = 0;
		chisq_target = 1.0;
		chisq_delta = 0.0000001f;
		chisq_percent = 95;  
		covar = NULL;
		alpha = NULL;
		err_axes = NULL;   
		iterations = 0; 
		fitfunc = GCI_multiexp_tau;

		memset(_restrained, 0, MAXFIT * sizeof(int));
		memset(_restrained_min, 0, MAXFIT * sizeof(float));
		memset(_restrained_max, 0, MAXFIT * sizeof(float));
	}

	/**
	* Destructor: make sure all arrays are freed.
	*/
	~SLIMCurve(){
		freePrivateVars();  // Should be done by fit functions, but just make sure.
	}

	/**
		Check the input values and assign pointers to arrays if required.
		May use the private variables.
	*/
	int checkValues() {
		if (transient == NULL) return SLIM_CURVE_SETTINGS_ERROR;
		if (ndata == 0) return SLIM_CURVE_SETTINGS_ERROR;

		if (fit_start < 0) fit_start = 0;
		if (fit_end <= fit_start) fit_end = ndata-1;
		if (fit_end > ndata-1) return SLIM_CURVE_SETTINGS_ERROR;

		if (noise_model == NOISE_GIVEN && noise_sd == NULL) return SLIM_CURVE_SETTINGS_ERROR;

		for (int i = 0, _nparamfree = 0; i < nparam; i++) {
			if (paramfree[i]) _nparamfree++;
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

	/**
		If the private variables have been used to assign arrays, free them here.
	*/
	void freePrivateVars() {

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

	/**	
	* Convenience function to setup the data required for a fit
	* \param[in] ltransient The transient (time resolved) signal to be analysed, the 'data'.
	* \param[in] lndata The number of data points.
	* \param[in] lxincr The time increment in between the values in the y array.
	* \param[in] data_start The index into the y array marking the start to the data.
	* \param[in] fit_start The index into the y array marking the start to the data to be used in the fit.
	* \param[in] fit_end The index into the y array marking the end of the data to be used in the fit.
	*/
	void setupData(float* ltransient, int lndata, float lxincr=1.0, int ldata_start=0, int lfit_start=0, int lfit_end=0) {
		transient = ltransient;
		ndata = lndata;
		time_incr = lxincr;
		data_start = ldata_start;
		fit_start = lfit_start;
		fit_end = lfit_end;
	}

	/**
	* Convenience function to setup the IRF
	* \param[in] instr The instrument response (IRF) or prompt signal to be used.
	* \param[in] ninstr The number of data points in the prompt.
	*/
	void setupIRF(float* linstr, int lninstr) {
		instr = linstr;
		ninstr = lninstr;
	}

	/**
	* Wrapper for GCI_triple_integral_fitting_engine()
	* Uses the parameters already setup to perform the fit
	*/
	int fitRLD() {

		float *Z = &(param[SLIM_CURVE_RLD_PARAM_Z]);
		float *A = &(param[SLIM_CURVE_RLD_PARAM_A]);
		float *tau = &(param[SLIM_CURVE_RLD_PARAM_TAU]);

		int err = checkValues();
		if (err < 0) return err;

		iterations = GCI_triple_integral_fitting_engine(time_incr, &transient[data_start], fit_start - data_start, fit_end - data_start,
			instr, ninstr, noise_model, noise_sd, Z, A, tau, fitted, residuals, &chisq, chisq_target*(fit_end - fit_start + ninstr - 3));

		freePrivateVars();

		return iterations;
	}

	/**
	* Wrapper for GCI_marquardt_fitting_engine()
	* Uses the parameters already setup to perform the fit
	*/
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

		if (restrain==ECF_RESTRAIN_USER) GCI_set_restrain_limits(nparam, _restrained, _restrained_min, _restrained_max);

		iterations = GCI_marquardt_fitting_engine(time_incr, &transient[data_start], ndata - data_start, fit_start-data_start, fit_end-data_start,
			instr, ninstr, noise_model, noise_sd, param, paramfree, nparam, restrain, fitfunc,
			fitted, residuals, &chisq, covar, alpha, err_axes, chisq_target*(fit_end - fit_start + ninstr - _nparamfree), chisq_delta, chisq_percent);

		freePrivateVars();

		return iterations;
	}

	/**
	* Wrapper for GCI_Phasor()
	* Uses the parameters already setup to perform the fit
	*/
	int fitPhasor(float lZ=0.0f) {

		float *Z = &(param[SLIM_CURVE_PHASOR_PARAM_Z]);
		float *u = &(param[SLIM_CURVE_PHASOR_PARAM_U]);
		float *v = &(param[SLIM_CURVE_PHASOR_PARAM_V]);
		float *taup = &(param[SLIM_CURVE_PHASOR_PARAM_TAUP]);
		float *taum = &(param[SLIM_CURVE_PHASOR_PARAM_TAUM]);
		float *tau = &(param[SLIM_CURVE_PHASOR_PARAM_TAU]);

		*Z = lZ;  // This function requires a Z in order to subtract it

		int err = checkValues();
		if (err < 0) return err;

		GCI_Phasor(time_incr, &transient[data_start], fit_start - data_start, fit_end - data_start,
			Z, u, v, taup, taum, tau, fitted, residuals, &chisq);

		freePrivateVars();

		return SLIM_CURVE_SUCCESS;
	}

	/**
	* Convenience function for fixing a parameter in a LMA fit
	* /param[in] parameter Index into paramter array of the parameter you wish to fix, use e.g. SLIM_CURVE_MONO_TAU etc.
	* /param[in] value The value to fix the paramter to in future fits.
	*/
	void fixParameter(int parameter, float value) {
		paramfree[parameter] = 0;
		param[parameter] = value;
	}

	/**
	* Clear all fixing of paramters; make all parameters free
	*/
	void clearFixed(void) {
		memset(paramfree, 1, MAXFIT * sizeof(int));
	}

	/**
	* Convenience function for restraining a parameter in a LMA fit
	* /param[in] parameter Index into paramter array of the parameter you wish to fix, use e.g. SLIM_CURVE_MONO_TAU etc.
	* /param[in] min The min value to restrain the paramter to in future fits.
	* /param[in] max The max value to restrain the paramter to in future fits.
	*/
	void restrainParameter(int parameter, float min, float max) {
		_restrained[parameter] = 1;
		_restrained_min[parameter] = min;
		_restrained_max[parameter] = max;
		restrain = ECF_RESTRAIN_USER;
	}

	/**
	* Clear all restraining of paramters; make all parameters free
	*/
	void clearRestrained(void) {
		memset(_restrained, 0, MAXFIT * sizeof(int));
		memset(_restrained_min, 0, MAXFIT * sizeof(float));
		memset(_restrained_max, 0, MAXFIT * sizeof(float));
		restrain = ECF_RESTRAIN_DEFAULT;
	}

	/**
	* Get the reduced chi squared value for the last fit performed based in the chisq and the degrees of freedom
	* /return The reduce chisq value
	*/
	float getReducedChiSq(void) {
		return (chisq / (fit_end - fit_start - nparam));
	}

};


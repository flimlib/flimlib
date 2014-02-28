
/*
 * SLIM-curve package for exponential curve fitting of spectral lifetime data.
 * 
 * Copyright (C) 2010 - 2014 Gray Institute University of Oxford and Board of
 * Regents of the University of Wisconsin-Madison.
 * 
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
 * 
 */

/** SLIM Curve - Exponential Curve Fitting
 * \file Ecf.h
 */

/* This is Ecf.h, the public header file for the 2003 version of the
   ECF library. */

#ifndef _GCI_ECF
#define _GCI_ECF

/* #defines which are publically needed */

/** Noise Type.
'NOISE_CONST' - every data point is assumed to have the same supplied variance.<br>
'NOISE_GIVEN' - every data point can have an individual variance, given via a data array.<br>
'NOISE_GAUSSIAN_FIT' - Variance for Gaussian distribution is used at all data points
Variants based on the data or the fit.<br>
'NOISE_POISSON_DATA and NOISE_POISSON_FIT' - Variance for Gaussian distribution is used with a lower limit of 15, this being the point where the Gaussian approximation begins to break down with Poissonian data. 
Variants based on the data or the fit.<br>
'NOISE_MLE' - Maximum Likelihood Estimation through the use of the Poisson equation (Laurence and Chromy 2010).<br>
*/
typedef enum { NOISE_CONST, NOISE_GIVEN, NOISE_POISSON_DATA,
	       NOISE_POISSON_FIT, NOISE_GAUSSIAN_FIT, NOISE_MLE } noise_type;

/** Fit Type.
 * Chooses between a multi-exponential fit and the special 'stretched' exponential.
 */
typedef enum { FIT_GLOBAL_MULTIEXP, FIT_GLOBAL_STRETCHEDEXP } fit_type;

/** Restrain Type.
 * Chooses the restrain type, either the default with very wide sensible limits or 'user' defined limits. See GCI_set_restrain_limits.
 */
typedef enum { ECF_RESTRAIN_DEFAULT, ECF_RESTRAIN_USER } restrain_type;

/* Single transient analysis functions */

/**
 * The main entry point for Triple Integral or Rapid Lifetime Determination (RLD).
 * Uses GCI_triple_integral_*() to fit repeatedly until chisq_target is met or a maximum number of iterations are used.
 *
 * \param[in] xincr The time increment inbetween the values in the y array.
 * \param[in] y The transient (time resolved) signal to be analysed, the 'data'.
 * \param[in] fit_start The index into the y array marking the start to the data to be used in the fit.
 * \param[in] fit_end The index into the y array marking the end of the data to be used in the fit.
 * \param[in] instr The instrument reponse (IRF) or prompt signal to be used (optional, can pass NULL).
 * \param[in] noise The noise_type to be used.
 * \param[in] sig The standard deviation at each data point in y if noise_type NOISE_GIVEN is used (optional, can pass NULL).
 * \param[out] Z The returned background value from the fit.
 * \param[out] A The returned amplitude value from the fit.
 * \param[out] tau The returned lifetime value from the fit.
 * \param[out] fitted An array containing values fitted to the data, the 'fit'. Fit points are coincident in time with the data points.
 * \param[out] residuals An array containing the difference between the data and the fit.
 * \param[out] chisq The resulting reduced chi squared value of the fit.
 * \param[in] chisq_target A target chi squared value to aim for. If this value is reached fitting will stop.
 * \return An error code, 0 = success.
 */
int GCI_triple_integral_fitting_engine(float xincr, float y[], int fit_start, int fit_end,
							  float instr[], int ninstr, noise_type noise, float sig[],
							  float *Z, float *A, float *tau, float *fitted, float *residuals,
							  float *chisq, float chisq_target);

/**
 * The main entry point for LM and MLE fitting.
 * Uses GCI_marquardt_instr() to fit repeatedly until chisq_target is met or a maximum number of iterations are used.
 * This can be used to fit any function through the 'fitfunc', but is primarily used for single, multi and stretched exponential fits.
 * Predefined fitting models GCI_multiexp_tau and GCI_stretchedexp are provided.
 *
 * \param[in] xincr The time increment inbetween the values in the y array.
 * \param[in] trans The transient (time resolved) signal to be analysed, the 'data'.
 * \param[in] ndata The number of data points.
 * \param[in] fit_start The index into the y array marking the start to the data to be used in the fit. Some data before this start index is required for convolution with the prompt.
 * \param[in] fit_end The index into the y array marking the end of the data to be used in the fit.
 * \param[in] prompt The instrument reponse (IRF) or prompt signal to be used (optional, can pass NULL).
 * \param[in] nprompt The number of data points in the prompt (ignored if prompt = NULL).
 * \param[in] noise The noise_type to be used.
 * \param[in] sig The standard deviation at each data point in y if noise_type NOISE_GIVEN is used (optional, can pass NULL).
 * \param[out] param An array of parameters, the order of which must match fitfunc. 
 * \param[in] paramfree An array indicating which parameters are free (1), fixed (0)
 * \param[in] nparam The number of parameters.
 * \param[in] restrain Parameter restraining type. Normally use ECF_RESTRAIN_DEFAULT. Use ECF_RESTRAIN_USER if restraining parameters has been setup via GCI_set_restrain_limits.
 * \param[in] fitfunc Encodes the function to fit to the data, e.g. GCI_multiexp_tau.
 * \param[out] fitted An array containing values fitted to the data, the 'fit'. Fit points are coincident in time with the data points.
 * \param[out] residuals An array containing the difference between the data and the fit.
 * \param[out] chisq The resulting reduced chi squared value of the fit.
 * \param[in] chisq_target A target chi squared value to aim for. If this value is reached fitting will stop.
 * \return An error code, 0 = success.
 */
int GCI_marquardt_fitting_engine(float xincr, float *trans, int ndata, int fit_start, int fit_end, 
						float prompt[], int nprompt,
						noise_type noise, float sig[],
						float param[], int paramfree[],
					   int nparam, restrain_type restrain,
					   void (*fitfunc)(float, float [], float *, float [], int),
					   float *fitted, float *residuals, float *chisq,
					   float **covar, float **alpha, float **erraxes,
					   float chisq_target, float chisq_delta, int chisq_percent);

int GCI_triple_integral(float xincr, float y[],
						int fit_start, int fit_end,
						noise_type noise, float sig[],
						float *Z, float *A, float *tau,
						float *fitted, float *residuals,
						float *chisq, int division);
int GCI_triple_integral_instr(float xincr, float y[],
							  int fit_start, int fit_end,
							  float instr[], int ninstr,
							  noise_type noise, float sig[],
							  float *Z, float *A, float *tau,
							  float *fitted, float *residuals,
							  float *chisq, int division);

int GCI_marquardt(float x[], float y[], int ndata,
				  noise_type noise, float sig[],
				  float param[], int paramfree[], int nparam,
				  restrain_type restrain,
				  void (*fitfunc)(float, float [], float *, float [], int),
				  float *fitted, float *residuals,
				  float **covar, float **alpha, float *chisq,
				  float chisq_delta, float chisq_percent, float **erraxes);
int GCI_marquardt_instr(float xincr, float y[],
					int ndata, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[],
					float param[], int paramfree[], int nparam,
					restrain_type restrain,
					void (*fitfunc)(float, float [], float *, float [], int),
					float *fitted, float *residuals,
					float **covar, float **alpha, float *chisq,
					float chisq_delta, float chisq_percent, float **erraxes);
void GCI_marquardt_cleanup(void);

/* Global analysis analysis functions */

int GCI_marquardt_global_exps_instr(float xincr, float **trans,
					int ndata, int ntrans, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[], int ftype,
					float **param, int paramfree[], int nparam,
					restrain_type restrain, float chisq_delta, 
					float **fitted, float **residuals,
					float chisq_trans[], float *chisq_global, int *df,
					int drop_bad_transients);
int GCI_marquardt_global_generic_instr(float xincr, float **trans,
					int ndata, int ntrans, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[],
					float **param, int paramfree[], int nparam, int gparam[],
					restrain_type restrain, float chisq_delta, 
					void (*fitfunc)(float, float [], float *, float [], int),
					float **fitted, float **residuals,
					float chisq_trans[], float *chisq_global, int *df);

/* Support plane analysis functions */

int GCI_SPA_1D_marquardt(
				float x[], float y[], int ndata,
				noise_type noise, float sig[],
				float param[], int paramfree[], int nparam,
				restrain_type restrain, float chisq_delta,
				void (*fitfunc)(float, float [], float *, float [], int),
				int spa_param, int spa_nvalues,
				float spa_low, float spa_high,
				float chisq[], void (*progressfunc)(float));
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
				float **chisq, void (*progressfunc)(float));
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
				float chisq[], float chisq_target, void (*progressfunc)(float));
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
				float **chisq, float chisq_target, void (*progressfunc)(float));
int GCI_SPA_1D_marquardt_global_exps_instr(
					float xincr, float **trans,
					int ndata, int ntrans, int fit_start, int fit_end,
					float instr[], int ninstr,
					noise_type noise, float sig[], int ftype,
					float **param, int paramfree[], int nparam,
					restrain_type restrain, float chisq_delta, int drop_bad_transients,
					int spa_param, int spa_nvalues,
					float spa_low, float spa_high,
					float chisq_global[], int df[], void (*progressfunc)(float));
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
					float **chisq_global, int **df, void (*progressfunc)(float));
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
					float chisq_global[], int df[], void (*progressfunc)(float));
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
					float **chisq_global, int **df, void (*progressfunc)(float));

/* Setting restraints */

int GCI_set_restrain_limits(int nparam, int restrain[],
							float minval[], float maxval[]);

/* Predefined fitting models */

void GCI_multiexp_lambda(float x, float param[],
			 float *y, float dy_dparam[], int nparam);
void GCI_multiexp_tau(float x, float param[],
		      float *y, float dy_dparam[], int nparam);
void GCI_stretchedexp(float x, float param[],
		      float *y, float dy_dparam[], int nparam);

/* Utility functions */
float **GCI_ecf_matrix(long nrows, long ncols);
void GCI_ecf_free_matrix(float **m);

void ECF_ExportParams_start (char path[]);
void ECF_ExportParams_stop (void);

#endif /* _GCI_ECF */

// Emacs settings:
// Local variables:
// mode: c
// c-basic-offset: 4
// tab-width: 4
// End:

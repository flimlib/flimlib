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

 /** FLIMLib - Multi-Pixel Exponential Curve Fitting Header
  * \file EcfMultiple.h
  */

#ifndef _ECF_MULTIPLE_H
#define _ECF_MULTIPLE_H

#include "Ecf.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** A structure to represent strided 1D arrays */
struct array1d {
	float* data; /**<The pointer to the first element*/
	size_t sizes[1]; /**<The size of the array*/
	ptrdiff_t strides[1]; /**<the byte stride between elements*/
};

/** A structure to represent strided 2D arrays */
struct array2d {
	float* data; /**<The pointer to the first element*/
	size_t sizes[2]; /**<The number of rows, columns*/
	ptrdiff_t strides[2]; /**<The bytes between rows, columns*/
};

/** A structure to represent strided 3D arrays */
struct array3d {
	float* data; /**<The pointer to the first element*/
	size_t sizes[3]; /**<The number of layers, rows, columns*/
	ptrdiff_t strides[3]; /**<The bytes between layers, rows, columns*/
};

/** A structure to represent strided 1D arrays of booleans*/
struct array1d_int8 {
	int8_t* data; /**<The pointer to the first element*/
	size_t sizes[1]; /**<The size of the array*/
	ptrdiff_t strides[1]; /**<the byte stride between elements*/
};

/** A structure containing common flim parameters */
struct common_params {
	/**
	 * [in] The time between samples
	 */
	float xincr; 

	/**
	 * [in] The data to fit. The first axis is spatial and the second is 
	 * temporal
	 */
	struct array2d* trans; 

	/**
	 * [in] The index of the start of the fit. Some data before this start 
	 * index is required if convolving with the prompt. 
	 */
	int fit_start; 

	/**
	 * [in] The index of the end of the fit
	 */
	int fit_end; 

	/**
	 * [out] Fitted values coincident in time with data points 
	 * (optional, can pass NULL).
	 */
	struct array2d* fitted; 

	/**
	 * [out] The difference between the data and the fit 
	 * (optional, can pass NULL).
	 */
	struct array2d* residuals; 

	/**
	 * [out] The resulting raw (RLD and LMA) or reduced (Phasor) chi squared 
	 * value of the fit (optional, can pass NULL).
	 */
	struct array1d* chisq; 

	/**
	 * [in] A mask containing 1s or 0s to select which pixels to fit 
	 * (optional, can pass NULL to fit all pixels).
	 */
	struct array1d_int8* fit_mask; 
};

/** A structure containing parameters used in LMA */
struct marquardt_params {
	/**
	 * [in] instr The instrument reponse(IRF) or prompt signal to be used
	 * (optional, can pass NULL).
	 */
	struct array1d* instr;

	/**
	 * [in] The noise_type to be used
	 */
	noise_type noise;

	/**
	 * [in] The standard deviation at each data point in y if #noise_type 
	 * NOISE_GIVEN is used (optional, can pass NULL).
	 */
	struct array1d* sig; 

	/**
	 * [in,out] Provide parameter estimates, these are overridden with the
	 * fitted values. The first axis is spatial and the second must match the 
	 * parameters of `fitfunc`
	 */
	struct array2d* param; 

	/**
	 * [in] An array indicating which parameters are free (1), fixed (0) 
	 * (optional, can pass NULL to free all params).
	 */
	struct array1d_int8* paramfree; 

	/**
	 * [in] Restrain type to use. Normally use ECF_RESTRAIN_DEFAULT. 
	 * Use ECF_RESTRAIN_USER if restraining parameters has been setup via 
	 * GCI_set_restrain_limits
	 */
	restrain_type restrain; 

	/**
	 * [in] Encodes the function to fit to the data, e.g. #GCI_multiexp_tau
	 */
	void (*fitfunc)(float, float[], float*, float[], int); 

	/**
	 * [out] The covariance matrix (optional, can pass NULL).
	 */
	struct array3d* covar; 

	/**
	 * [out] The alpha matrix (optional, can pass NULL).
	 */
	struct array3d* alpha; 

	/**
	 * [out] The dimensions of the confidence ellipsoid of the chisq 
	 * (optional, can pass NULL).
	 */
	struct array3d* erraxes; 

	/**
	 * [in] A raw chi squared value to aim for. If this value is reached 
	 * fitting will stop. If you want to aim for a reduced chisq 
	 * (say 1.1 or 1.0) you must multiply by the degree of freedom. 
	 */
	float chisq_target; 

	/**
	 * [in] An individual fit will continue if the chisq value changes by 
	 * more then this amount. Try 1E-5.
	 */
	float chisq_delta; 

	/**
	 * [in] Defines the confidence interval when calculating the error axes, 
	 * e.g. 95 %.
	 */
	int chisq_percent; 
};

/** A structure containing parameters used in Phasor */
struct phasor_params {
	/**
	 * [in] Z must have been estimated previously so that it can be 
	 * subtracted from the data here.
	 */
	struct array1d* Z; 

	/**
	 * [out] u The 'horizontal' phasor coordinate.
	 */
	struct array1d* u; 

	/**
	 * [out] v The 'vertical' phasor coordinate.
	 */
	struct array1d* v; 

	/**
	 * [out] taup The lifetime calculated from the phase change.
	 */
	struct array1d* taup; 

	/**
	 * [out] taum The lifetime calculated from the amplitude change 
	 * (the demodulation).
	 */
	struct array1d* taum; 

	/**
	 * [out] tau The average of the other taus.
	 */
	struct array1d* tau; 
};

/** A structure containing parameters used in RLD */
struct triple_integral_params {
	/**
	 * [in] instr The instrument reponse(IRF) or prompt signal to be used 
	 * (optional, can pass NULL).
	 */
	struct array1d* instr; 

	/**
	 * [in] The noise_type to be used
	 */
	noise_type noise; 

	/**
	 * [in] The standard deviation at each data point in y if #noise_type 
	 * NOISE_GIVEN is used (optional, can pass NULL).
	 */
	struct array1d* sig; 

	/**
	 * [out] The returned background value from the fit.
	 */
	struct array1d* Z; 

	/**
	 * [out] A The returned amplitude value from the fit.
	 */
	struct array1d* A; 

	/**
	 * [out] tau The returned lifetime value from the fit.
	 */
	struct array1d* tau; 

	/**
	 * [in] A raw chi squared value to aim for. If this value is reached 
	 * fitting will stop. If you want to aim for a reduced chisq 
	 * (say 1.1 or 1.0) you must multiply by the degree of freedom. 
	 * A negative value will lead to only a single iteration
	 */
	float chisq_target; 
};

/** A structure containing parameters used in flim */
struct flim_params {
	/**
	 * [in,out] A structure to containing common flim parameters 
	 */
	struct common_params* common; 
	union {
		/**
		 * [in,out] A structure to containing parameters used in LMA
		 */
		struct marquardt_params *marquardt; 

		/**
		 * [in,out] A structure to containing parameters used in RLD
		 */
		struct triple_integral_params *triple_integral; 

		/**
		 * [in,out] A structure to containing parameters used in Phasor
		 */
		struct phasor_params *phasor; 
	};
};

/** 
 * Print a 1D array with a maximum number of elements in each dimension
 */
void print_array1d(struct array1d* arr, int max_print);

/** 
 * Print a 2D array with a maximum number of elements in each dimension
 */
void print_array2d(struct array2d* arr, int max_print);

/** 
 * Print a 3D array with a maximum number of elements in each dimension
 */
void print_array3d(struct array3d* arr, int max_print);

/** 
 * Print a 1D byte array with a maximum number of elements in each dimension
 */
void print_array1d_int8(struct array1d_int8* arr, int max_print);

/** 
 * Prints common flim parameters 
 */
void print_common_params(struct common_params *fit);

/** 
 * multi-pixel LMA fitting using strided arrays. Accesses parameters in-place 
 * when possible. Otherwise temporary copies will be made
 * \param[in,out] flim A flim_params structure containing parameters used in 
 * flim (must contain a marquardt_params structure)
 */
int GCI_marquardt_fitting_engine_many(struct flim_params* flim);

/**
 * multi-pixel RLD fitting using strided arrays. Accesses parameters in-place 
 * when possible. Otherwise temporary copies will be made. if chisq_target is 
 * negative, only one iteration of the algorithm will occur per pixel
 * \param[in,out] flim A flim_params structure containing parameters used in 
 * flim (must contain a triple_integral_params structure)
 */
int GCI_triple_integral_fitting_engine_many(struct flim_params* flim);

/**
 * multi-pixel Phasor analysis using strided arrays. Accesses parameters 
 * in-place when possible. Otherwise temporary copies will be made
 * \param[in,out] flim A flim_params structure containing parameters used in 
 * flim (must contain a phasor_params structure)
 */
int GCI_Phasor_many(struct flim_params* flim);


#ifdef __cplusplus
} // extern "C"
#endif

#endif /* _ECF_MULTIPLE_H */

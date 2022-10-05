/*-
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

/** FLIMLib - Multi-Pixel Exponential Curve Fitting
 * \file EcfMultiple.c
 */

#include "EcfMultiple.h"

#include "Ecf.h"
#include "GCI_Phasor.h"
#include "GCI_PhasorInternal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static inline float *array1d_float_ptr(struct array1d *arr, size_t i) {
	return (float *)((char *)arr->data + i * arr->strides[0]);
}

static inline int8_t *array1d_int8_ptr(struct array1d_int8 *arr, size_t i) {
	return (int8_t *)((char *)arr->data + i * arr->strides[0]);
}

static inline float *array2d_float_ptr(struct array2d *arr, size_t i, size_t j) {
	return (float *)((char *)arr->data + i * arr->strides[0] + j * arr->strides[1]);
}

static inline float* array3d_float_ptr(struct array3d* arr, size_t i, size_t j, size_t k) {
	return (float *)((char*)arr->data + i * arr->strides[0] + j * arr->strides[1] + k * arr->strides[2]);
}

// Allocate unstrided buffer of appropriate size only if source is
// strided; otherwise return NULL.
// Caller must later call free() on returned pointer.
static float* allocate_temp_1d(struct array1d* source) {
	if (source == NULL || source->strides[0] == sizeof(float))
		return NULL;
	return (float*)malloc(source->sizes[0] * sizeof(float));
}

// Allocate unstrided buffer of appropriate size only if source is
// strided; otherwise return NULL.
// Caller must later call free() on returned pointer.
static float* allocate_temp_2d_row(struct array2d* source) {
	if (source == NULL || source->strides[1] == sizeof(float))
		return NULL;
	return (float*)malloc(source->sizes[1] * sizeof(float));
}

// Allocate unstrided buffer of given size if source is strided or NULL.
// (Used for functions/parameters that require non-null buffers.)
// Caller must later call free() on returned pointer.
static float* allocate_required_temp_2d_row(struct array2d* source, size_t size) {
	if (source != NULL && source->strides[1] == sizeof(float))
		return NULL;
	return (float*)malloc(size * sizeof(float));
}

// Return pointer to unstrided data by copying into buffer if necessary.
static float* get_unstrided_1d_input(struct array1d* source, float* buffer) {
	if (source == NULL)
		return NULL;
	if (source->strides[0] == sizeof(float)) // source is unstrided
		return source->data;
	for (size_t col = 0; col < source->sizes[0]; col++)
		*(buffer + col) = *array1d_float_ptr(source, col);
	return buffer;
}

// Return pointer to unstrided data by copying into buffer if necessary.
static float* get_unstrided_2d_row_input(struct array2d* source, float* buffer, int row) {
	if (source == NULL)
		return NULL;
	if (source->strides[1] == sizeof(float)) // source is unstrided
		return array2d_float_ptr(source, row, 0);
	for (size_t col = 0; col < source->sizes[1]; col++) 
		*(buffer + col) = *array2d_float_ptr(source, row, col);
	return buffer;
}

// Return pointer to unstrided output buffer: either within dest if possible,
// or the provided temporary buffer.
// buffer must be return value of allocate_temp_2d_row.
static float* get_unstrided_2d_row_output(struct array2d* dest, float* buffer, int row) {
	if (dest == NULL)
		return NULL;
	if (dest->strides[1] == sizeof(float)) // dest is unstrided
		return array2d_float_ptr(dest, row, 0);
	return buffer;
}

// Used for functions/parameters that cannot be NULL.
static float* get_required_unstrided_2d_row_output(struct array2d* dest, float* buffer, int row) {
	if (dest != NULL && dest->strides[1] == sizeof(float)) // dest is unstrided
		return array2d_float_ptr(dest, row, 0);
	return buffer;
}

// Copy data in buffer to its final destination if necessary.
// buffer must be return value of get_unstrided_2d_row_output().
static void write_2d_row_output(float* buffer, struct array2d* dest, int row) {
	if (dest == NULL || dest->strides[1] == sizeof(float)) // buffer _is_ the unstrided destination
		return;
	for (size_t col = 0; col < dest->sizes[1]; col++)
		*array2d_float_ptr(dest, row, col) = *(buffer + col);
}

// Write a matrix allocated by GCI_ecf_matrix into an array3d
static void write_matrix(float** source, struct array3d* dest, int layer) {
	if (dest == NULL)
		return;
	for (size_t row = 0; row < dest->sizes[1]; row++)
		for (size_t col = 0; col < dest->sizes[2]; col++)
			*array3d_float_ptr(dest, layer, row, col) = *(source[0] + row * dest->sizes[2] + col);
}

// fill a specified row in an array2d with NaNs
static void write_2d_row_nan(struct array2d* dest, int row) {
	if (dest == NULL)
		return;
	for (size_t col = 0; col < dest->sizes[1]; col++)
		*array2d_float_ptr(dest, row, col) = NAN;
}


// fill a specified matrix in an array3d with NaNs
static void write_matrix_nan(struct array3d* dest, int layer) {
	if (dest == NULL)
		return;
	for (size_t row = 0; row < dest->sizes[1]; row++)
		for (size_t col = 0; col < dest->sizes[2]; col++)
			*array3d_float_ptr(dest, layer, row, col) = NAN;
}

void print_array1d(struct array1d* arr, int max_print) {
	if (arr == NULL) return;
	for (int col = 0; col < arr->sizes[0] && col <= max_print; col++) {
		printf("%f\t", *array1d_float_ptr(arr, col));
	}
	printf("\n");
}

void print_array2d(struct array2d* arr, int max_print) {
	if (arr == NULL) return;
	for (int row = 0; row < arr->sizes[0] && row <= max_print; row++) {
		for (int col = 0; col < arr->sizes[1] && col <= max_print; col++) {
			printf("%f\t", *array2d_float_ptr(arr, row, col));
		}
		printf("\n");
	}
}

void print_array3d(struct array3d* arr, int max_print) {
	if (arr == NULL) return;
	for (int layer = 0; layer < arr->sizes[0] && layer <= max_print; layer++) {
		printf("Layer %d\n", layer);
		for (int row = 0; row < arr->sizes[1] && row <= max_print; row++) {
			for (int col = 0; col < arr->sizes[2] && col <= max_print; col++) {
				printf("%f\t", *array3d_float_ptr(arr, layer, row, col));
			}
			printf("\n");
		}
	}
}

void print_array1d_int8(struct array1d_int8* arr, int max_print) {
	if (arr == NULL) return;
	for (int col = 0; col < arr->sizes[0] && col <= max_print; col++) {
		printf("%d\t", *array1d_int8_ptr(arr, col));
	}
	printf("\n");
}

void print_common_params(struct common_params *common) {
	printf("xincr: %f\n", common->xincr);
	printf("fit start: %d\n", common->fit_start);
	printf("fit end: %d\n", common->fit_end);
	printf("trans:\n");
	print_array2d(common->trans, 5);
	printf("fitted:\n");
	print_array2d(common->fitted, 5);
	printf("residuals:\n");
	print_array2d(common->residuals, 5);
	printf("chisq:\n");
	print_array1d(common->chisq, 5);
	printf("fit mask:\n");
	print_array1d_int8(common->fit_mask, 5);
}

// TODO return negative ints to report any errors
// TODO set all output values to be nan in the event of a failed fit
int GCI_marquardt_fitting_engine_many(struct flim_params* flim) {

	// ndata is the length of data to fit. will be sliced further using fit_start and fit_end
	// fit_start and fit_end are not redundant with the stride characteristics of trans since
	// on occasion data outside of the fit range is used for convolution with the "prompt" whatever that means
	int ndata = (int)(flim->common->trans->sizes[1]);
	int nparam = (int)(flim->marquardt->param->sizes[1]);
	int ninstr = (int)(flim->marquardt->instr == NULL ? 0 : flim->marquardt->instr->sizes[0]);

	// allocate inputs and outputs. they may be NULL (no memory used) if not needed
	float* temp_trans = allocate_temp_2d_row(flim->common->trans);
	float* temp_instr = allocate_temp_1d(flim->marquardt->instr);
	float* temp_sig = allocate_temp_1d(flim->marquardt->sig);
	float* temp_param = allocate_temp_2d_row(flim->marquardt->param);
	float* temp_fitted = allocate_required_temp_2d_row(flim->common->fitted, ndata);
	float* temp_residuals = allocate_required_temp_2d_row(flim->common->residuals, ndata);

	// these two don't get modified by the algorithm and are constant for all pixels
	float* unstrided_instr = get_unstrided_1d_input(flim->marquardt->instr, temp_instr);
	float* unstrided_sig = get_unstrided_1d_input(flim->marquardt->sig, temp_sig);

	// paramfree passed to marquardt must not be NULL! if NULL we want to leave all parameters free
	int* temp_paramfree = malloc(nparam * sizeof(int));
	if (flim->marquardt->paramfree == NULL && temp_paramfree != NULL) {
		for (int i = 0; i < nparam; i++)
			temp_paramfree[i] = 1;
	}
	else {
		for (int i = 0; i < nparam; i++)
			temp_paramfree[i] = *array1d_int8_ptr(flim->marquardt->paramfree, i);
	}

	// allocate the matrices. we will need these no matter what
	float** temp_covar = GCI_ecf_matrix(nparam, nparam);
	float** temp_alpha = GCI_ecf_matrix(nparam, nparam);
	// erraxes can be NULL. This avoids the estimate errors step
	float** temp_erraxes = flim->marquardt->erraxes == NULL ? NULL : GCI_ecf_matrix(nparam, nparam);

	for (int i = 0; i < flim->common->trans->sizes[0]; i++) { // iterate over rows of trans
		if (flim->common->fit_mask == NULL || *array1d_int8_ptr(flim->common->fit_mask, i)) { 
			// get input parameters trans and param in an unstrided format
			float* unstrided_trans = get_unstrided_2d_row_input(flim->common->trans, temp_trans, i);
			float* unstrided_param = get_unstrided_2d_row_input(flim->marquardt->param, temp_param, i);

			// get the outputs in unstrided format
			float* unstrided_fitted = get_required_unstrided_2d_row_output(flim->common->fitted, temp_fitted, i);
			float* unstrided_residuals = get_required_unstrided_2d_row_output(flim->common->residuals, temp_residuals, i);
			float* unstrided_chisq = flim->common->chisq == NULL ? NULL : array1d_float_ptr(flim->common->chisq, i);

			int error_code = 0;
			if(flim->common->fit_start >= flim->common->fit_end || flim->common->xincr <= 0.0){
				error_code = -1;
			}
			else{
				error_code = GCI_marquardt_fitting_engine(flim->common->xincr, unstrided_trans, ndata, flim->common->fit_start, 
					flim->common->fit_end, unstrided_instr, ninstr, flim->marquardt->noise, unstrided_sig, unstrided_param, 
					temp_paramfree, nparam, flim->marquardt->restrain, flim->marquardt->fitfunc, unstrided_fitted, 
					unstrided_residuals, unstrided_chisq, temp_covar, temp_alpha, temp_erraxes, flim->marquardt->chisq_target, 
					flim->marquardt->chisq_delta, flim->marquardt->chisq_percent);
			}
			
			if (error_code < 0){ // fit has failed
				// fill all outputs with NaN
				write_2d_row_nan(flim->marquardt->param, i);
				write_2d_row_nan(flim->common->fitted, i);
				write_2d_row_nan(flim->common->residuals, i);
				if (unstrided_chisq != NULL){
					*unstrided_chisq = NAN;
				}
				write_matrix_nan(flim->marquardt->covar, i);
				write_matrix_nan(flim->marquardt->alpha, i);
				write_matrix_nan(flim->marquardt->erraxes, i);
			}
			else{
				// time to handle outputs!
				write_2d_row_output(unstrided_param, flim->marquardt->param, i);
				write_2d_row_output(unstrided_fitted, flim->common->fitted, i);
				write_2d_row_output(unstrided_residuals, flim->common->residuals, i);

				// able to save memory if these array3d are NULL. This prevents copying over the results
				// it was not an option to modify these in-place since the function requires the ecf matrix format
				write_matrix(temp_covar, flim->marquardt->covar, i);
				write_matrix(temp_alpha, flim->marquardt->alpha, i);
				write_matrix(temp_erraxes, flim->marquardt->erraxes, i);
			}
		}
	}

	// free anything with temp_ in its name
	free(temp_trans);
	free(temp_instr);
	free(temp_sig);
	free(temp_param);
	free(temp_paramfree);
	free(temp_fitted);
	free(temp_residuals);
	GCI_ecf_free_matrix(temp_covar);
	GCI_ecf_free_matrix(temp_alpha);
	GCI_ecf_free_matrix(temp_erraxes);

	return 0;
}

int GCI_triple_integral_fitting_engine_many(struct flim_params* flim) {
	int ninstr = (int)(flim->marquardt->instr == NULL ? 0 : flim->marquardt->instr->sizes[0]);
	// set the target to be INFINITY since this will cause the algorithm to iterate only once
	float chisq_target_in = flim->triple_integral->chisq_target < 0 ? INFINITY : flim->triple_integral->chisq_target;
	
	// allocate inputs and outputs. they may be NULL (no memory used) if not needed
	float* temp_trans = allocate_temp_2d_row(flim->common->trans);
	float* temp_instr = allocate_temp_1d(flim->triple_integral->instr);
	float* temp_sig = allocate_temp_1d(flim->triple_integral->sig);
	float* temp_fitted = allocate_temp_2d_row(flim->common->fitted);
	float* temp_residuals = allocate_temp_2d_row(flim->common->residuals);
	
	// these two don't get modified by the algorithm and are constant for all pixels
	float* unstrided_instr = get_unstrided_1d_input(flim->triple_integral->instr, temp_instr);
	float* unstrided_sig = get_unstrided_1d_input(flim->triple_integral->sig, temp_sig);

	for (int i = 0; i < flim->common->trans->sizes[0]; i++) {
		if (flim->common->fit_mask == NULL || *array1d_int8_ptr(flim->common->fit_mask, i)) {

			float* unstrided_trans = get_unstrided_2d_row_input(flim->common->trans, temp_trans, i);
			float* unstrided_fitted = get_unstrided_2d_row_output(flim->common->fitted, temp_fitted, i);
			float* unstrided_residuals = get_unstrided_2d_row_output(flim->common->residuals, temp_residuals, i);
			float* unstrided_chisq = flim->common->chisq == NULL ? NULL : array1d_float_ptr(flim->common->chisq, i);
			int error_code = 0;
			if(flim->common->fit_start >= flim->common->fit_end || flim->common->xincr <= 0.0){
				error_code = -1;
			}
			else{
				error_code = GCI_triple_integral_fitting_engine(flim->common->xincr, unstrided_trans, flim->common->fit_start, flim->common->fit_end,
					unstrided_instr, ninstr, flim->triple_integral->noise, unstrided_sig,
					array1d_float_ptr(flim->triple_integral->Z, i), array1d_float_ptr(flim->triple_integral->A, i),
					array1d_float_ptr(flim->triple_integral->tau, i), unstrided_fitted, unstrided_residuals,
					unstrided_chisq, chisq_target_in);
			}

			if(error_code < 0){
				*array1d_float_ptr(flim->triple_integral->Z, i) = NAN;
				*array1d_float_ptr(flim->triple_integral->A, i) = NAN;
				*array1d_float_ptr(flim->triple_integral->tau, i) = NAN;
				write_2d_row_nan(flim->common->fitted, i);
				write_2d_row_nan(flim->common->residuals, i);
				if (unstrided_chisq != NULL){
					*unstrided_chisq = NAN;
				}
			}
			else{
				write_2d_row_output(unstrided_fitted, flim->common->fitted, i);
				write_2d_row_output(unstrided_residuals, flim->common->residuals, i);
			}
		}
	}

	free(temp_trans);
	free(temp_instr);
	free(temp_sig);
	free(temp_fitted);
	free(temp_residuals);

	return 0;
}

int GCI_Phasor_many(struct flim_params* flim) {

	float *cosine, *sine;
	int nBins = (flim->common->fit_end - flim->common->fit_start);
	if (nBins > 0) {
		cosine = malloc((size_t)nBins * sizeof(float));
		sine = malloc((size_t)nBins * sizeof(float));
		createSinusoids(nBins, cosine, sine);
	}
	else {
		cosine = sine = NULL;
	}

	float* temp_trans = allocate_temp_2d_row(flim->common->trans);
	float* temp_fitted = allocate_temp_2d_row(flim->common->fitted);
	float* temp_residuals = allocate_temp_2d_row(flim->common->residuals);

	for (int i = 0; i < flim->common->trans->sizes[0]; i++) {
		if (flim->common->fit_mask == NULL || *array1d_int8_ptr(flim->common->fit_mask, i)) {

			float* unstrided_trans = get_unstrided_2d_row_input(flim->common->trans, temp_trans, i);
			float* unstrided_fitted = get_unstrided_2d_row_output(flim->common->fitted, temp_fitted, i);
			float* unstrided_residuals = get_unstrided_2d_row_output(flim->common->residuals, temp_residuals, i);
			float* unstrided_chisq = flim->common->chisq == NULL ? NULL : array1d_float_ptr(flim->common->chisq, i);

			int error_code = 0;
			if(nBins <= 0 || flim->common->xincr <= 0.0){
				error_code = -1;
			}
			else{
				error_code = GCI_Phasor_compute(flim->common->xincr, unstrided_trans, flim->common->fit_start, flim->common->fit_end,
					array1d_float_ptr(flim->phasor->Z, i), cosine, sine, array1d_float_ptr(flim->phasor->u, i), 
					array1d_float_ptr(flim->phasor->v, i), array1d_float_ptr(flim->phasor->taup, i), 
					array1d_float_ptr(flim->phasor->taum, i), array1d_float_ptr(flim->phasor->tau, i),
					unstrided_fitted, unstrided_residuals, unstrided_chisq);
			}

			if (error_code < 0){
				*array1d_float_ptr(flim->phasor->u, i) = NAN;
				*array1d_float_ptr(flim->phasor->v, i) = NAN;
				*array1d_float_ptr(flim->phasor->taup, i) = NAN;
				*array1d_float_ptr(flim->phasor->taum, i) = NAN;
				*array1d_float_ptr(flim->phasor->tau, i) = NAN;
				write_2d_row_nan(flim->common->fitted, i);
				write_2d_row_nan(flim->common->residuals, i);
				if (unstrided_chisq != NULL){
					*unstrided_chisq = NAN;
				}
			}
			else{
				write_2d_row_output(unstrided_fitted, flim->common->fitted, i);
				write_2d_row_output(unstrided_residuals, flim->common->residuals, i);
			}
		}
	}
	free(temp_trans);
	free(temp_fitted);
	free(temp_residuals);

	free(cosine);
	free(sine);

	return 0;
}

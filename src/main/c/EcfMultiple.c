/*
* # % L
* FLIMLib package for exponential curve fitting of fluorescence lifetime data.
*%%
* Copyright(C) 2010 - 2015 University of Oxford and Board of Regents of the
* University of Wisconsin - Madison.
* %%
*This program is free software : you can redistribute it and /or modify
* it under the terms of the GNU General Public License as
* published by the Free Software Foundation, either version 3 of the
* License, or (at your option) any later version.
*
*This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public
* License along with this program.If not, see
* < http://www.gnu.org/licenses/gpl-3.0.html>.
*#L%
* /

/** FLIMLib - Multi-Pixel Exponential Curve Fitting

 * \file EcfMultiple.c
 */

#include <stdio.h>
#include "EcfMultiple.h" // header file
#include "Ecf.h" 
#include "GCI_Phasor.h"
#include <stdlib.h> 
#include <math.h> 


#define ARRAY1D_ELEM_PTR(arr, i) ((char *)(arr)->data + (i) * (arr)->strides[0])
#define ARRAY2D_ELEM_PTR(arr, i, j) ((char *)(arr)->data + (i) * (arr)->strides[0] + (j) * (arr)->strides[1])
#define ARRAY3D_ELEM_PTR(arr, i, j, k) ((char *)(arr)->data + (i) * (arr)->strides[0] + (j) * (arr)->strides[1] + (k) * (arr)->strides[2])

float* allocate_temp(struct array1d* source) {
	if (source == NULL || source->strides[0] == sizeof(float)) // we don't need this, don't allocate anything.
		return NULL;
	return (float*)malloc(source->sizes[0] * sizeof(float));
}

float* allocate_temp_row(struct array2d* source) {
	if (source == NULL || source->strides[1] == sizeof(float)) // we don't need this, don't allocate anything.
		return NULL;
	return (float*)malloc(source->sizes[1] * sizeof(float));
}

// marquardt fitting engine requires that fitted and residuals are not NULL
float* allocate_required_temp_row(struct array2d* source, size_t size) {
	if (source != NULL && source->strides[1] == sizeof(float)) // we don't need this, don't allocate anything.
		return NULL;
	return (float*)malloc(size * sizeof(float));
}

float* read_strided(struct array1d* source, float* destination) {
	if (source == NULL)
		return NULL;
	if (source->strides[0] == sizeof(float))
		return source->data; // strided source is actually unstrided! it will be modified in-place
	for (size_t col = 0; col < source->sizes[0]; col++) // iterate over columns
		*(destination + col) = *(float*)ARRAY1D_ELEM_PTR(source, col);
	return destination;
}

float* read_strided_row(struct array2d* source, float* destination, int row) {
	if (source == NULL)
		return NULL;
	if (source->strides[1] == sizeof(float))
		return ARRAY2D_ELEM_PTR(source, row, 0); // strided source is actually unstrided! it will be modified in-place
	for (size_t col = 0; col < source->sizes[1]; col++) // iterate over columns
		*(destination + col) = *(float*)ARRAY2D_ELEM_PTR(source, row, col);
	return destination;
}

float* prep_unstrided_output(struct array2d* strided_output, float* unstrided_output, int row) {
	if (strided_output == NULL) // the caller doesn't want this output so we pass the NULL onward
		return NULL;
	if (strided_output->strides[1] == sizeof(float)) // strided output is actually unstrided! it will be modified in-place
		return ARRAY2D_ELEM_PTR(strided_output, row, 0);
	return unstrided_output;
}

// marquardt fitting engine requires that fitted and residuals are not NULL
float* prep_required_unstrided_output(struct array2d* strided_output, float* unstrided_output, int row) {
	if (strided_output != NULL && strided_output->strides[1] == sizeof(float)) // strided output is actually unstrided! it will be modified in-place
		return ARRAY2D_ELEM_PTR(strided_output, row, 0);
	return unstrided_output;
}

void write_strided_row(float* source, struct array2d* destination, int row) {
	if (destination == NULL || destination->strides[1] == sizeof(float)) // if the output was NULL the user didn't want. If it was unstrided, it was already modified in-place
		return;
	for (size_t col = 0; col < destination->sizes[1]; col++) // iterate over columns
		*(float*)ARRAY2D_ELEM_PTR(destination, row, col) = *(source + col);
}

/* write a matrix allocated by GCI_ecf_matrix into an array3d */
void write_strided_matrix(float** source, struct array3d* destination, int layer) {
	if (destination == NULL) // even if the destination isn't strided, the data still needs to be copied since GCI_ecf_matrix has a pointer header
		return;
	for (size_t row = 0; row < destination->sizes[1]; row++) // iterate over columns
		for (size_t col = 0; col < destination->sizes[2]; col++) // iterate over columns
			*(float*)ARRAY3D_ELEM_PTR(destination, layer, row, col) = *(source[0] + row * destination->sizes[2] + col);
}

void print_array1d(struct array1d* arr, int max_print) {
	if (arr == NULL) return;
	for (int col = 0; col < arr->sizes[0] && col <= max_print; col++) {
		printf("%f\t", *(float*)ARRAY1D_ELEM_PTR(arr, col));
	}
	printf("\n");
}

void print_array2d(struct array2d* arr, int max_print) {
	if (arr == NULL) return;
	for (int row = 0; row < arr->sizes[0] && row <= max_print; row++) {
		for (int col = 0; col < arr->sizes[1] && col <= max_print; col++) {
			printf("%f\t", *(float*)ARRAY2D_ELEM_PTR(arr, row, col));
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
				printf("%f\t", *(float*)ARRAY3D_ELEM_PTR(arr, layer, row, col));
			}
			printf("\n");
		}
	}
}

void print_array1d_int8(struct array1d_int8* arr, int max_print) {
	if (arr == NULL) return;
	for (int col = 0; col < arr->sizes[0] && col <= max_print; col++) {
		printf("%d\t", *(int8_t*)ARRAY1D_ELEM_PTR(arr, col));
	}
	printf("\n");
}

void print_common(struct common_params *common) {
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
	size_t ndata = flim->common->trans->sizes[1];
	size_t nparam = flim->marquardt->param->sizes[1];
	int ninstr = flim->marquardt->instr == NULL ? 0 : flim->marquardt->instr->sizes[0];

	// allocate inputs and outputs. they may be NULL (no memory used) if not needed
	float* temp_trans = allocate_temp_row(flim->common->trans);
	float* temp_instr = allocate_temp(flim->marquardt->instr);
	float* temp_sig = allocate_temp(flim->marquardt->sig);
	float* temp_param = allocate_temp_row(flim->marquardt->param);
	float* temp_fitted = allocate_required_temp_row(flim->common->fitted, ndata);
	float* temp_residuals = allocate_required_temp_row(flim->common->residuals, ndata);
	float local_chisq;

	// these two don't get modified by the algorithm and are constant for all pixels
	float* unstrided_instr = read_strided(flim->marquardt->instr, temp_instr);
	float* unstrided_sig = read_strided(flim->marquardt->sig, temp_sig);

	//paramfree must not be NULL! if NULL we want to leave all parameters free
	int* temp_paramfree = malloc(nparam * sizeof(int));
	if (flim->marquardt->paramfree == NULL) {
		for (int i = 0; i < nparam; i++)
			temp_paramfree[i] = 1;
	}
	else
		for (int i = 0; i < nparam; i++) {
			int8_t pf = *ARRAY1D_ELEM_PTR(flim->marquardt->paramfree, i);
			temp_paramfree[i] = pf;
		}

	// allocate the matrices. we will need these no matter what
	float** temp_covar = GCI_ecf_matrix(nparam, nparam);
	float** temp_alpha = GCI_ecf_matrix(nparam, nparam);
	float** temp_erraxes = GCI_ecf_matrix(nparam, nparam);

	for (int i = 0; i < flim->common->trans->sizes[0]; i++) { // iterate over rows of trans
		if (flim->common->fit_mask == NULL || *ARRAY1D_ELEM_PTR(flim->common->fit_mask, i)) { 
			// get input parameters trans and param in an unstrided format
			float* unstrided_trans = read_strided_row(flim->common->trans, temp_trans, i);
			float* unstrided_param = read_strided_row(flim->marquardt->param, temp_param, i);

			// get the outputs in unstrided format
			float* unstrided_fitted = prep_required_unstrided_output(flim->common->fitted, temp_fitted, i);
			float* unstrided_residuals = prep_required_unstrided_output(flim->common->residuals, temp_residuals, i);
			float* unstrided_chisq = flim->common->chisq == NULL ? &local_chisq : ARRAY1D_ELEM_PTR(flim->common->chisq, i);

			GCI_marquardt_fitting_engine(flim->common->xincr, unstrided_trans, ndata, flim->common->fit_start, flim->common->fit_end,
				unstrided_instr, ninstr, flim->marquardt->noise, unstrided_sig, unstrided_param, temp_paramfree, nparam,
				flim->marquardt->restrain, flim->marquardt->fitfunc, unstrided_fitted, unstrided_residuals, unstrided_chisq,
				temp_covar, temp_alpha, temp_erraxes, flim->marquardt->chisq_target, flim->marquardt->chisq_delta, flim->marquardt->chisq_percent);
			
			// time to handle outputs!
			write_strided_row(unstrided_param, flim->marquardt->param, i);
			write_strided_row(unstrided_fitted, flim->common->fitted, i);
			write_strided_row(unstrided_residuals, flim->common->residuals, i);

			// able to save memory if these array3d are NULL. This prevents copying over the results
			// it was not an option to modify these in-place since the function requires the ecf matrix format
			write_strided_matrix(temp_covar, flim->marquardt->covar, i);
			write_strided_matrix(temp_alpha, flim->marquardt->alpha, i);
			write_strided_matrix(temp_erraxes, flim->marquardt->erraxes, i);
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
	int ninstr = flim->marquardt->instr == NULL ? 0 : flim->marquardt->instr->sizes[0];
	// set the target to be NAN since this will cause the algorithm to iterate only once
	float chisq_target_in = flim->triple_integral->chisq_target < 0 ? NAN : flim->triple_integral->chisq_target;
	
	// allocate inputs and outputs. they may be NULL (no memory used) if not needed
	float* temp_trans = allocate_temp_row(flim->common->trans);
	float* temp_instr = allocate_temp_row(flim->triple_integral->instr);
	float* temp_sig = allocate_temp_row(flim->triple_integral->sig);
	float* temp_fitted = allocate_temp_row(flim->common->fitted);
	float* temp_residuals = allocate_temp_row(flim->common->residuals);
	
	// these two don't get modified by the algorithm and are constant for all pixels
	float* unstrided_instr = read_strided(flim->triple_integral->instr, temp_instr);
	float* unstrided_sig = read_strided(flim->triple_integral->sig, temp_sig);

	for (int i = 0; i < flim->common->trans->sizes[0]; i++) {
		if (flim->common->fit_mask == NULL || *ARRAY1D_ELEM_PTR(flim->common->fit_mask, i)) {

			float* unstrided_trans = read_strided_row(flim->common->trans, temp_trans, i);
			float* unstrided_fitted = prep_unstrided_output(flim->common->fitted, temp_fitted, i);
			float* unstrided_residuals = prep_unstrided_output(flim->common->residuals, temp_residuals, i);
			float* unstrided_chisq = flim->common->chisq == NULL ? NULL : ARRAY1D_ELEM_PTR(flim->common->chisq, i);

			GCI_triple_integral_fitting_engine(flim->common->xincr, unstrided_trans, flim->common->fit_start, flim->common->fit_end,
				unstrided_instr, ninstr, flim->triple_integral->noise, unstrided_sig,
				ARRAY1D_ELEM_PTR(flim->triple_integral->Z, i), ARRAY1D_ELEM_PTR(flim->triple_integral->A, i),
				ARRAY1D_ELEM_PTR(flim->triple_integral->tau, i), unstrided_fitted, unstrided_residuals,
				unstrided_chisq, chisq_target_in);

			write_strided_row(unstrided_fitted, flim->common->fitted, i);
			write_strided_row(unstrided_residuals, flim->common->residuals, i);
		}
	}

	free(temp_trans);
	free(temp_instr);
	free(temp_sig);
	free(temp_fitted);
	free(temp_residuals);
}

int GCI_Phasor_many(struct flim_params* flim) {

	float* temp_trans = allocate_temp_row(flim->common->trans);
	float* temp_fitted = allocate_temp_row(flim->common->fitted);
	float* temp_residuals = allocate_temp_row(flim->common->residuals);

	for (int i = 0; i < flim->common->trans->sizes[0]; i++) {
		if (flim->common->fit_mask == NULL || *ARRAY1D_ELEM_PTR(flim->common->fit_mask, i)) {

			float* unstrided_trans = read_strided_row(flim->common->trans, temp_trans, i);
			float* unstrided_fitted = prep_unstrided_output(flim->common->fitted, temp_fitted, i);
			float* unstrided_residuals = prep_unstrided_output(flim->common->residuals, temp_residuals, i);
			float* unstrided_chisq = flim->common->chisq == NULL ? NULL : ARRAY1D_ELEM_PTR(flim->common->chisq, i);

			GCI_Phasor(flim->common->xincr, unstrided_trans, flim->common->fit_start, flim->common->fit_end,
				ARRAY1D_ELEM_PTR(flim->phasor->Z, i), ARRAY1D_ELEM_PTR(flim->phasor->u, i), ARRAY1D_ELEM_PTR(flim->phasor->v, i),
				ARRAY1D_ELEM_PTR(flim->phasor->taup, i), ARRAY1D_ELEM_PTR(flim->phasor->taum, i), ARRAY1D_ELEM_PTR(flim->phasor->tau, i),
				unstrided_fitted, unstrided_residuals, unstrided_chisq);

			write_strided_row(unstrided_fitted, flim->common->fitted, i);
			write_strided_row(unstrided_residuals, flim->common->residuals, i);
		}
	}
	free(temp_trans);
	free(temp_fitted);
	free(temp_residuals);
}
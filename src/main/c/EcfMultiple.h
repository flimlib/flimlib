/*
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2015 University of Oxford and Board of Regents of the
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

  /* This is EcfMultiple.h, the public header file for the 2021 version of EcfMultiple.c */

#include <stddef.h> // need this for ptrdiff_t
#include "Ecf.h" // need for noisetypes

#ifdef __cplusplus
extern "C" {
#endif

struct array1d {
	float* data;
	size_t size;
	ptrdiff_t stride;
};

struct array2d {
	float* data;
	size_t sizes[2]; // [num_rows, num_cols]
	ptrdiff_t strides[2]; // [bytes_between_rows, bytes_between_cols]
};

struct array3d {
	float* data;
	size_t sizes[3]; // [num_layers, num_rows, num_cols]
	ptrdiff_t strides[3]; // [bytes_between_layers, bytes_between_rows, bytes_between_cols]
};

struct common_params {
	float xincr;
	struct array2d* trans;
	int fit_start;
	int fit_end;
	struct array2d* fitted;
	struct array2d* residuals;
	struct array1d* chisq;
	struct array1d* fit_mask;
};

struct marquardt_params {
	struct array1d* instr;
	noise_type noise;
	struct array1d* sig;
	struct array2d* param;
	struct array1d* paramfree;
	restrain_type restrain;
	void (*fitfunc)(float, float[], float*, float[], int);
	struct array3d* covar;
	struct array3d* alpha;
	struct array3d* erraxes;
	float chisq_target;
	float chisq_delta;
	int chisq_percent;
};

struct phasor_params {
	struct array1d* Z;
	struct array1d* u;
	struct array1d* v;
	struct array1d* taup;
	struct array1d* taum;
	struct array1d* tau;
};

struct triple_integral_params {
	struct array1d* instr;
	noise_type noise;
	struct array1d* sig;
	struct array1d* Z;
	struct array1d* A;
	struct array1d* tau;
	float chisq_target;

};

struct flim_params {
	struct common_params* common;
	union {
		struct marquardt_params *marquardt;
		struct triple_integral_params *triple_integral;
		struct phasor_params *phasor;
	};
};


void print_array1d(struct array1d* arr, int max_print);

void print_array2d(struct array2d* arr, int max_print);

void print_array3d(struct array3d* arr, int max_print);

void print_common(struct common_params *fit);

/** multidimentional LMA fitting
* TODO write doc
*/
int GCI_marquardt_fitting_engine_many(struct flim_params* flim);

int GCI_triple_integral_fitting_engine_many(struct flim_params* flim);

int GCI_Phasor_many(struct flim_params* flim);
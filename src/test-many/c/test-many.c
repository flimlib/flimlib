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

/** 
 * This is a program to excersise the functions in EcfMultiple.h
 * Meant to be replaced with proper unit tests later
 * \file test-many.c
 */

#include "EcfMultiple.h"

#include <math.h>
#include <stdio.h>


#define NDATA 256
#define NPARAM 3
#define NROWS 2

void createExponential(float* time, float* photonCount, float a, float tau, float period, size_t count);

int main() {
	float a = 10.0;
	float tauIn = 1;
	float period = 0.04;
	float time[NDATA];
	float photonCount[NDATA];
	float scaledPhotonCount[NDATA];

	float Z_out;
	float A_out;
	float tau_out;
	createExponential(time, photonCount, a, tauIn, period, NDATA);
	GCI_triple_integral_fitting_engine(period, photonCount, 0, NDATA, NULL, 0, 3, NULL, &Z_out, &A_out, &tau_out, NULL, NULL, NULL, NAN);

	printf("%f\n", A_out);

	float covar_data[NPARAM * NPARAM * NROWS] = { 0 };
	float alpha_data[NPARAM * NPARAM * NROWS] = { 0 };
	float erraxes_data[NPARAM * NPARAM * NROWS] = { 0 };
	struct array3d covar = { covar_data, {NROWS, NPARAM, NPARAM}, {NPARAM * NPARAM * sizeof(float), NPARAM * sizeof(float), sizeof(float)} };
	struct array3d alpha = { NULL }; // passing null should be fine. skips copying step
	struct array3d erraxes = { NULL };

	float exponential2d[NDATA * NROWS] = { 0 };
	// fill with two experimental exponential curves
	createExponential(time, exponential2d, a, tauIn, period, NDATA);
	createExponential(time, exponential2d + NDATA, a * 2, tauIn * 2, period, NDATA); // double a and tau to see 2 different "pixels"
	struct array2d photonCount2d = { exponential2d, {NROWS, NDATA}, {NDATA * sizeof(float) , sizeof(float)} };

	float param2d_data[NPARAM * NROWS] = { 0.1, a + 0.1, tauIn + 0.1, 0.1, a + 0.1, tauIn + 0.1 }; // add 0.1 to make sure fitting works
	struct array2d param2d = { param2d_data, {NROWS, NPARAM}, {NPARAM * sizeof(float), sizeof(float)} };

	float fitted2d_data[NDATA * NROWS];
	struct array2d fitted2d = { fitted2d_data, {NROWS, NDATA}, {NDATA * sizeof(float) , sizeof(float)} };

	float residuals2d_data[NDATA * NROWS];
	struct array2d residuals2d = { residuals2d_data, {NROWS, NDATA}, {NDATA * sizeof(float) , sizeof(float)} };

	float *chisq_data[NROWS];
	struct array1d chisq = {chisq_data, NROWS, sizeof(float)};

	struct common_params common_in = {.xincr=period, .trans=&photonCount2d, .fit_start=0, .fit_end=NDATA, .fitted=&fitted2d, .residuals=&residuals2d, .chisq=&chisq, .fit_mask=NULL};

	struct marquardt_params marquardt_in = {.instr=NULL, .noise=NOISE_POISSON_FIT, .sig=NULL, .param=&param2d, .paramfree=NULL, 
		.restrain=ECF_RESTRAIN_DEFAULT, .fitfunc=GCI_multiexp_tau, .covar=&covar, .alpha=&alpha, .erraxes=&erraxes, .chisq_target=1.1, .chisq_delta=1E-5, .chisq_percent=95 };

	struct flim_params flim_in = { &common_in, &marquardt_in };

	GCI_marquardt_fitting_engine_many(&flim_in);

	printf("params\n");
	for (int i = 0; i < NPARAM * NROWS; i++)
		printf("%f\n", param2d.data[i]);

	printf("\n3d covariance matrix\n");
	//for (int i = 0; i < NPARAM * NPARAM * NROWS; i++)
	//	printf("%f\n", covar.data[i]);
	print_array3d(&covar, 5);

	common_in.chisq = NULL;
	common_in.residuals = NULL;

	// what happens if we only look at every other element in photonCount?
	photonCount2d.sizes[1] = NDATA/2;
	photonCount2d.strides[1] = 2 * sizeof(float);
	common_in.xincr = period * 2;

	// also lets reverse the order of param
	float param2d_data_rev[NPARAM * NROWS] = { tauIn + 0.1, a + 0.1, 0.1, tauIn + 0.1, a + 0.1, 0.1 };
	param2d.strides[1] = -sizeof(float);
	param2d.data = param2d_data_rev + 2;

	common_in.fit_end /= 2;

	GCI_marquardt_fitting_engine_many(&flim_in);

	printf("params\n");
	for (int i = 0; i < NPARAM * NROWS; i++)
		printf("%f\n", param2d_data_rev[i]);

	print_common_params(&common_in);

	float Z_data[NROWS] = { 0 };
	struct array1d Z = { .data = Z_data, .sizes = {NROWS}, .strides = {sizeof(float) } };
	float u_data[NROWS];
	struct array1d u = { .data = u_data, .sizes = {NROWS}, .strides = {sizeof(float) } };
	float v_data[NROWS];
	struct array1d v = { .data = v_data, .sizes = {NROWS}, .strides = {sizeof(float) } };
	float taup_data[NROWS];
	struct array1d taup = { .data = taup_data, .sizes = {NROWS}, .strides = {sizeof(float) } };
	float taum_data[NROWS];
	struct array1d taum = { .data = taum_data, .sizes = {NROWS}, .strides = {sizeof(float) } };
	float tau_data[NROWS];
	struct array1d tau = { .data = tau_data, .sizes = {NROWS}, .strides = {sizeof(float) } };

	
	struct phasor_params phasor_in = { .Z=&Z, .u=&u, .v=&v, .taup=&taup, .taum=&taum, .tau=&tau };

	flim_in.phasor = &phasor_in;

	print_common_params(&common_in);

	GCI_Phasor_many(&flim_in);

	print_common_params(&common_in);

	print_array1d(&tau, 5);


}

void createExponential(float* time, float* photonCount, float a, float tau, float period, size_t count) {
	for (int i = 0; i < count; i++) {
		float t = period * i;
		time[i] = t;
		photonCount[i] = a * exp(-t / tau);
	}
}




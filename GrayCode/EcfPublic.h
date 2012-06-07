//TODO ARG modified Paul Barber's ecf.h file, reduced to minimum set of curve fitting functions,
// replaced fitfunc with fit_type

/*
 Gray Institute's Exponential Curve Fitting (ECF) Library
 Copyright (C) 2003-@year@ Gray Institute for Radiation Oncology and Biology 
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* This is Ecf.h, the public header file for the  ECF library. */ //TODO ARG renamed EcfPublic.h

#ifndef _PUBLIC_GCI_ECF
#define _PUBLIC_GCI_ECF

/* #defines which are publically needed */

typedef enum { NOISE_CONST, NOISE_GIVEN, NOISE_POISSON_DATA,
	NOISE_POISSON_FIT, NOISE_GAUSSIAN_FIT, NOISE_MLE } noise_type;

typedef enum { FIT_GLOBAL_MULTIEXP, FIT_GLOBAL_STRETCHEDEXP } fit_type;

typedef enum { ECF_RESTRAIN_DEFAULT, ECF_RESTRAIN_USER } restrain_type;

/* Single transient analysis functions */

// the next fn uses GCI_triple_integral_*() to fit repeatedly until chisq_target is met //TODO ARG renamed nr_, these will initially be Numeric Recipes versions
int nr_GCI_triple_integral_fitting_engine(float xincr, float y[], int fit_start, int fit_end,
									   float instr[], int ninstr, noise_type noise, float sig[],
									   float *Z, float *A, float *tau, float *fitted, float *residuals,
									   float *chisq, float chisq_target);
			  
// the next fn uses GCI_marquardt_instr() to fit repeatedly until chisq_target is met //TODO ARG renamed nr_, these will initially be Numeric Recipes versions
int nr_GCI_marquardt_fitting_engine(float xincr, float *trans, int ndata, int fit_start, int fit_end, 
								 float prompt[], int nprompt, //TODO ARG is this actually instr[] & ninstr?
								 noise_type noise, float sig[],
								 float param[], int paramfree[],
								 int nparam, restrain_type restrain,
								 fit_type fit, //TODO ARG void (*fitfunc)(float, float [], float *, float [], int),
								 float *fitted, float *residuals, float *chisq); //,
								 //float **covar, float **alpha, float **erraxes,
								 //float chisq_target, int chisq_percent);

#endif /* _GCI_ECF */

// Emacs settings:
// Local variables:
// mode: c
// c-basic-offset: 4
// tab-width: 4
// End:

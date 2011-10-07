/* 
This file is part of the SLIM-curve package for exponential curve fitting of spectral lifetime data.

Copyright (c) 2010, 2011, Gray Institute University of Oxford & UW-Madison LOCI.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* This is EcfInternal.h, the header file for internal functions in
   the 2003 version of the ECF library. */

#ifndef _GCI_ECF_INTERNAL
#define _GCI_ECF_INTERNAL

#include "Ecf.h"  /* in case there's anything we need from there */

#define MAXFIT 20  /* The maximum number of parameters we'll ever try
		              to fit; saves dynamic allocation of small arrays.
				      If this is increased, then the arrays chisq50 etc.
				      in ecf.c will need to be extended.  The values can
				      be calculated by using the test function at the end
				      of the file. */

#define MAXITERS 80
#define MAXREFITS 10
#define MAXBINS 1024 /* Maximum number of lifetime bins; saves dynamic allocation of small arrays */

/* Functions from EcfSingle.c */

int GCI_marquardt_compute_fn(float x[], float y[], int ndata,
					 noise_type noise, float sig[],
					 float param[], int paramfree[], int nparam,
					 void (*fitfunc)(float, float [], float *, float [], int),
					 float yfit[], float dy[],
					 float **alpha, float beta[], float *chisq, float old_chisq,
					 float alambda);
int GCI_marquardt_compute_fn_instr(float xincr, float y[], int ndata,
				   int fit_start, int fit_end,
				   float instr[], int ninstr,
				   noise_type noise, float sig[],
				   float param[], int paramfree[], int nparam,
				   void (*fitfunc)(float, float [], float *, float [], int),
				   float yfit[], float dy[],
				   float **alpha, float beta[], float *chisq, float old_chisq,
				   float alambda,	
					float **pfnvals, float ***pdy_dparam_pure, float ***pdy_dparam_conv,
					int *pfnvals_len, int *pdy_dparam_nparam_size);
int GCI_marquardt_compute_fn_final(float x[], float y[], int ndata,
					 noise_type noise, float sig[],
					 float param[], int paramfree[], int nparam,
					 void (*fitfunc)(float, float [], float *, float [], int),
					 float yfit[], float dy[], float *chisq);
int GCI_marquardt_compute_fn_final_instr(float xincr, float y[], int ndata,
				   int fit_start, int fit_end,
				   float instr[], int ninstr,
				   noise_type noise, float sig[],
				   float param[], int paramfree[], int nparam,
				   void (*fitfunc)(float, float [], float *, float [], int),
				   float yfit[], float dy[], float *chisq,	
					float **pfnvals, float ***pdy_dparam_pure, float ***pdy_dparam_conv,
					int *pfnvals_len, int *pdy_dparam_nparam_size);

/* Functions from EcfGlobal.c */


/* Functions from EcfUtil.c */

int GCI_solve(float **a, int, float *b);
int GCI_invert(float **a, int);
void GCI_covar_sort(float **covar, int nparam, int paramfree[], int mfit);
float ***GCI_ecf_matrix_array(long nblocks, long nrows, long ncols);
void GCI_ecf_free_matrix_array(float ***marr);
int check_ecf_params (float param[], int nparam,
					void (*fitfunc)(float, float [], float *, float [], int));
int check_ecf_user_params (float param[], int nparam,
					void (*fitfunc)(float, float [], float *, float [], int));
int GCI_chisq(int nu, float chisq, float *root);
int multiexp_lambda_array(float xincr, float param[],
						  float *y, float **dy_dparam, int nx, int nparam);
int multiexp_tau_array(float xincr, float param[],
					   float *y, float **dy_dparam, int nx, int nparam);
int stretchedexp_array(float xincr, float param[],
					   float *y, float **dy_dparam, int nx, int nparam);
int ECF_Find_Float_Max (float data[], int np, float *max_val);

/* For debugging printing */
extern int ECF_debug;  /* defined in EcfUtil.c */
int dbgprintf(int dbg_level, const char *format, ...);

// Vars for the export of params at each iteration
int ecf_exportParams;
char ecf_exportParams_path[256];

void ecf_ExportParams_OpenFile (void);
void ecf_ExportParams_CloseFile (void);
void ecf_ExportParams (float param[], int nparam, float chisq);

#endif /* _GCI_ECF_INTERNAL */


// Emacs settings:
// Local variables:
// mode: c
// c-basic-offset: 4
// tab-width: 4
// End:

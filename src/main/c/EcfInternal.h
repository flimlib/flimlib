/* 
This file is part of the SLIM-curve package for exponential curve fitting of spectral lifetime data.

Copyright (c) 2010-2013, Gray Institute University of Oxford & UW-Madison LOCI.

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
#define MAXBINS 2048 /* Maximum number of lifetime bins; saves dynamic allocation of small arrays */

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
int GCI_solve_Gaussian(float **a, int n, float *b);
int GCI_invert_Gaussian(float **a, int n);
void pivot(float **a, int n, int *order, int col);
int lu_decomp(float **a, int n, int *order);
int solve_lu(float **lu, int n, float *b, int *order);
int GCI_solve_lu_decomp(float **a, int n, float *b);
int GCI_invert_lu_decomp(float **a, int n);
int GCI_solve(float **a, int n, float *b);
int GCI_invert(float **a, int n);
void GCI_covar_sort(float **covar, int nparam, int paramfree[], int mfit);
float **GCI_ecf_matrix(long nrows, long ncols);
void GCI_ecf_free_matrix(float **m);
float ***GCI_ecf_matrix_array(long nblocks, long nrows, long ncols);
void GCI_ecf_free_matrix_array(float ***marr);
void GCI_multiexp_lambda(float x, float param[],
						 float *y, float dy_dparam[], int nparam);
int multiexp_lambda_array(float xincr, float param[],
						  float *y, float **dy_dparam, int nx, int nparam);
void GCI_multiexp_tau(float x, float param[],
					  float *y, float dy_dparam[], int nparam);
int multiexp_tau_array(float xincr, float param[],
					   float *y, float **dy_dparam, int nx, int nparam);
void GCI_stretchedexp(float x, float param[],
					  float *y, float dy_dparam[], int nparam);
int stretchedexp_array(float xincr, float param[],
					   float *y, float **dy_dparam, int nx, int nparam);
int check_ecf_params (float param[], int nparam,
                      void (*fitfunc)(float, float [], float *, float [], int));
int GCI_set_restrain_limits(int nparam, int restrain[],
							float minval[], float maxval[]);
int check_ecf_user_params (float param[], int nparam,
                           void (*fitfunc)(float, float [], float *, float [], int));
int GCI_marquardt_estimate_errors(float **alpha, int nparam, int mfit,
								  float d[], float **v, float interval);
float GCI_incomplete_gamma(float a, float x);
float GCI_log_gamma(float x);
float GCI_gamma(float x);
float GCI_gammap(float a, float x);
int GCI_chisq(int nu, float chisq, float *root);
int ECF_Find_Float_Max (float data[], int np, float *max_val);

/* For debugging printing */
extern int ECF_debug;  /* defined in EcfUtil.c */
int dbgprintf(int dbg_level, const char *format, ...);
void ECF_ExportParams_start (char path[]);
void ECF_ExportParams_stop (void);
void ecf_ExportParams_OpenFile (void);
void ecf_ExportParams_CloseFile (void);
void ecf_ExportParams (float param[], int nparam, float chisq);

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

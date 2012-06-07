/*
 *  EcfPublic.c
 *  ECF
 *
 *  Created by Aivar Grislis on 7/2/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "EcfPublic.h"
#include "ecf.h"
#include "EcfInternal.h" //for MAXFIT
#include <stdio.h>

// the next fn uses GCI_triple_integral_*() to fit repeatedly until chisq_target is met
int nr_GCI_triple_integral_fitting_engine(float xincr, float y[], int fit_start, int fit_end,
									   float instr[], int ninstr, noise_type noise, float sig[],
									   float *Z, float *A, float *tau, float *fitted, float *residuals,
									   float *chisq, float chisq_target) {
//printf("HELLO\nWORLD!!!\n");
	return GCI_triple_integral_fitting_engine(xincr, y, fit_start, fit_end,
									   NULL /*instr*/, ninstr, noise, sig,
									   Z, A, tau, fitted, residuals,
									   chisq, chisq_target);
}
	

// the next fn uses GCI_marquardt_instr() to fit repeatedly until chisq_target is met
int nr_GCI_marquardt_fitting_engine(float xincr, float *trans, int ndata, int fit_start, int fit_end,
								 float prompt[], int nprompt, //TODO ARG is this actually instr[] & ninstr?
								 noise_type noise, float sig[],
								 float param[], int paramfree[],
								 int nparam, restrain_type restrain,
								 fit_type fit, //TODO ARG void (*fitfunc)(float, float [], float *, float [], int),
								 float *fitted, float *residuals, float *chisq) { //,
								 //float **covar, float **alpha, float **erraxes,
									//float chisq_target, int chisq_percent) {
	void (*fitfunc)(float, float [], float *, float[], int) = NULL;

//printf("Incoming params:\n");
//int i;
//for (i = 0; i < nparam; ++i) {
//printf(" %g ", param[i]);
//}
//printf("\n");
	
	switch (fit) {
		case FIT_GLOBAL_MULTIEXP:
 //printf("GCI_multiexp_lambda\n");
			fitfunc = GCI_multiexp_lambda; //TODO lambda VS tau
			break;
		case FIT_GLOBAL_STRETCHEDEXP:
//printf("GCI_stretchedexp\n");
			fitfunc = GCI_stretchedexp;
			break;
	}

//printf("before LMA\n");
        float **covarX = GCI_ecf_matrix(ndata,ndata);
        float **alphaX = GCI_ecf_matrix(ndata,ndata);
        float **erraxesX = GCI_ecf_matrix(ndata,ndata);
        float chisq_target = 100.0;
        int chisq_percent = 50;
//printf("about to LMA\n");

       // if (0)
        int     returnValue = GCI_marquardt_fitting_engine(xincr, trans, ndata, fit_start, fit_end,
								 prompt, nprompt,
								 noise, sig,
								 param, paramfree,
								 nparam, restrain,
								 fitfunc,
								 fitted, residuals, chisq,
								 covarX, alphaX, erraxesX,
								 chisq_target, chisq_percent);
//printf("back from LMA %d\n", returnValue);

//printf("Fitted params:\n");
//for (i = 0; i < nparam; ++i) {
//printf(" %g ", param[i]);
//}
//printf("\n");

        int n;
        float x, y;
        x = 0.0;
        float dy_dparam[MAXFIT];
        for (n = 0; n < ndata; ++n) {
            if (n < fit_start || n > fit_end) {
                fitted[n] = residuals[n] = 0.0;
            }
            else {
                printf("fitted %d\n", n);
                (*fitfunc)(x, param, &y, dy_dparam, nparam);
                fitted[n] = y;
                residuals[n] = trans[n] - y;
           }
            x += xincr;
        }


        return returnValue;
}


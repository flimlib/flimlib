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
/*=================================================================================*/
/* File:       bayes_MonoExpAnalysisBinned.h                                       */
/*---------------------------------------------------------------------------------*/
/* Purpose:    Header file for bayes_MonoExpAnalysisBinned.c.                      */
/*---------------------------------------------------------------------------------*/
/* References: [1] - "Bayesian data analysis and parameter inference and           */
/*                   applications to imaging based on photon counting",            */
/*                   ACC Coolen, KCL, December 8th 2006.                           */
/*---------------------------------------------------------------------------------*/
/* Notes:      None.                                                               */
/*=================================================================================*/
/* Revision history:                                                               */
/*---------------------------------------------------------------------------------*/
/* Date   | Modification                                                           */
/*---------------------------------------------------------------------------------*/
/* 090424 | Creation, mrowley.                                                     */
/*---------------------------------------------------------------------------------*/
/*        |                                                                        */
/*=================================================================================*/

#ifndef BAYES_MONO_EXP_ANALYSIS_BINNED_H
#define BAYES_MONO_EXP_ANALYSIS_BINNED_H


/***********************************************************************************/
/*                                                                                 */
/*                             DEFINES / SWITCHES                                  */
/*                                                                                 */
/***********************************************************************************/


#include "bayes_Types.h"
#include "bayes_RapidBayesDecayAnalysis.h"


struct BayesContainerMonoExpMinusLogProbW0W1
{
    /* Forward parameters... */
    int                           *data;
    int                            nbins;
    int                            fitstart;
    int                            nphotons;
    double                         interval;
    double                         modulationperiod;
    BayesInstrRsp_t               *instr;
    double                         hyperparam;
    /* Pre-calculated fixed instrument values... */
    double                        *upsilon1;
    BayesRapidMonoExpDiscreteValues_t *likelihoods;
    double                        *fluorescencedecayphotonlikelihoods;
    double                        *binwalls;
    /* Reverse parameters */
    double                         w0;
    /* Internal routine failure indicator... */
    int                            error;
    /* Hyperparameter optimisation... */
    double                         alphamin;
    BayesRapidMonoExpValueStore_t *grid;
    double                       **datalikelihoods;
};
typedef struct BayesContainerMonoExpMinusLogProbW0W1 MonoExpMinusLogProbW0W1_t;



struct BayesAveErrDistribution
{
    int      paramxid;
    int      paramyid;
    int      m;
    int      n;
    int      manualrange;
    float    xmin;
    float    xmax;
    float    ymin;
    float    ymax;
    float   *x;
    float   *y;
    float  **z;
};    
typedef struct BayesAveErrDistribution BayesAveErrDistn_t;



/***********************************************************************************/
/*                                                                                 */
/*                             FUNCTION PROTOTYPES                                 */
/*                                                                                 */
/***********************************************************************************/

#if 0
int direct_most_probable_wrapped_binned (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, int distr, float *val, float *params, float *fitted, float *residuals);

int indirect_most_probable_wrapped_binned (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, int distr, float *val, float *params, float *fitted, float *residuals);

int most_probable_freeinstrument_wrapped_binned (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, int distr, float *val, float *params, float *fitted, float *residuals);

int averages_and_errorbars_wrapped_binned (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
                                           int p, float interval, float alpha, float precision, int quick, int distr, float *val, float *params, float *fitted, float *residuals);

int gci_triple_integral_rld(int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, int distr, float *val, float *params, float *fitted, float *residuals);
#endif


int bayes_AveragesAndErrorBarsBinLikelihood(int                *data,
                                            int                 nbins,
                                            int                 fitstart,
                                            double             *binwalls,
                                            int                 nphotons,
		                                    double             *w0_mp,
			                                double             *w1_mp,
		                                    double             *w0_ave,
			                                double             *w1_ave,
                                            double             *dw0,
			                                double             *dw1,
		                                    BayesInstrRsp_t    *instr,
                                            float               interval,
                                            float               modulationperiod,
                                            float               alpha,
			                                float               precision,
			                                BayesAveErrDistn_t *distr,
                                            float              *minuslogprob);






int direct_most_probable_wrapped_bin_likelihood (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, BayesAveErrDistn_t *distr, float *val, float *params, float *fitted, float *residuals);

int indirect_most_probable_wrapped_bin_likelihood (float *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, BayesAveErrDistn_t *distr, float *val/*, float *params, float *fitted, float *residuals*/);

int most_probable_freeinstrument_wrapped_bin_likelihood (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, BayesAveErrDistn_t *distr, float *val, float *params, float *fitted, float *residuals);

int averages_and_errorbars_wrapped_bin_likelihood (float *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, BayesAveErrDistn_t *distr, float *val/*, float *params, float *fitted, float *residuals*/);

#if 0
int bayes_MonoExpDetermineMostProbAveErrParamValues(float                    *data,
                                                    int                       nbins,
                                                    int                       fitstart,
                                                    int                       fitend,
                                                    int                       nphotons,
                                                    ArrLikelihoodConstants_t *constants,
                                                    int                       nparams,
                                                    int                      *param_free,
                                                    float                    *param_mp,//w1,tau1,w2,tau2
                                                    float                    *param_ave,
                                                    float                    *param_err,
                                                    float                    *delay,
                                                    float                    *width,
                                                    float                     interval,
                                                    float                     modulationperiod,
                                                    float                     alpha,
                                                    float                     precision,
                                                    int                       quick,
                                                    BayesAveErrDistn_t       *distr,
                                                    float                    *val);
#else
int bayes_MonoExpDetermineMostProbAveErrParamValues(float                    *data,
                                                    int                       nbins,
                                                    int                       fitstart,
                                                    int                       fitend,
                                                    int                       nphotons,
                                                    ArrLikelihoodConstants_t *constants,
                                                    int                       nparams,
                                                    int                      *param_free,
                                                    float                    *param_mp,//w1,tau1,w2,tau2
                                                    float                    *param_ave,
                                                    float                    *param_err,
                                                    //float                    *delay,
                                                    //float                    *width,
                                                    BayesInstrRsp_t           *instr,
                                                    float                     interval,
                                                    float                     modulationperiod,
                                                    float                     alpha,
                                                    float                     precision,
                                                    int                       quick,
                                                    BayesAveErrDistn_t       *distr,
                                                    float                    *val);
#endif


#endif

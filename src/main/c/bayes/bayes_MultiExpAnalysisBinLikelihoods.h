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
#ifndef BAYES_MULTI_EXP_ANALYSIS_BINNED_H
#define BAYES_MULTI_EXP_ANALYSIS_BINNED_H


#include "bayes_Types.h"


struct BayesContainerMultiExpMinusLogProbParams
{
    /* Forward parameters... */
    int                    *data;
    int                     nbins;
    int                     fitstart;
    double                 *binwalls;
    int                     nphotons;
    int                     ndecays;
    int                     nparams;
    BayesUserFixedParams_t *paramfixing;
    double                  interval;
    double                  modulationperiod;
    BayesInstrRsp_t        *instr;
    double                  hyperparam;
    /* Probability 'normalization' parameters */
    double                  normalization;
};
typedef struct BayesContainerMultiExpMinusLogProbParams MultiExpMinusLogProbParams_t;

#if 0
int bayes_MultiExpDetermineMostProbParamValues(int                      *data,
                                               int                       nbins,
                                               int                       nphotons,
                                               ArrLikelihoodConstants_t *constants,
                                               int                       nparams,
                                               int                      *param_is_user_fixed,
                                               float                    *param_mp,//w0,w1,tau1,w2,tau2,...
                                               float                    *delay,
                                               float                    *width,
                                               float                     interval,
                                               float                     modulationperiod,
                                               int                       nhyperparams,
                                               float                    *hyperparams,
                                               float                    *val);
#else
int bayes_MultiExpDetermineMostProbParamValues(int                      *data,
                                               int                       nbins,
                                               double                   *binwalls,
                                               int                      *nphotons,
                                               int                       ndecays,
                                               double                   *weights,
                                               double                   *taus,
                                               BayesUserFixedParams_t   *paramfixing,
                                               double                    interval,
                                               double                    modulationperiod,
                                               BayesInstrRsp_t          *instr,
                                               double                    alpha,
                                               void                     *rapidgrid,
                                               double                   *val);
#endif

#if 0 // some code in preparation for doing integration (i.e. computing bayes factor) using transformed lifetimes
double bayes_MultiExpProbParamsTransformedTaus(double *x, int id, void *container);
#endif


#endif /* BAYES_MULTI_EXP_ANALYSIS_BINNED_H */

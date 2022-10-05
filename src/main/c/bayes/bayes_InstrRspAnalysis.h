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
#ifndef BAYES_INSTR_RSP_ANALYSIS_H
#define BAYES_INSTR_RSP_ANALYSIS_H

#include "bayes_Types.h"
#include "bayes_Interface.h"

#if 1
#define BAYES_MAX_FILENAME_LEN  256

struct BayesInstrRspEstimationStore
{
    /* Data in... */
    char            filepromptloaded[BAYES_MAX_FILENAME_LEN];
    char            filedataloaded[BAYES_MAX_FILENAME_LEN];
    int             allpixels;
    int             xcoordinate;
    int             ycoordinate;
    int             binningtype;
    int             binninglevel;
    float          *trans;
    int             transbins;
    float           transbinwidth;
    int             fitstart;
    int             fitend;
    int             nphotons;
    /* Decay model... */
    int             modeltype;
    /* Estimated instrument response... */
    BayesInstrRsp_t instr;
};
typedef struct BayesInstrRspEstimationStore BayesInstrRspEstimationStore_t;

BayesInstrRspEstimationStore_t *bayes_GetBayesInstrRspEstimationStorePtr(void);

/**
 * TODO: Thread Unsafe
 */
int bayes_UpdateEstimatedInstrRspDetailsStoreAdv(/* Image/prompt data used... */
                                              char           *fileprompt,
                                              char           *filedata,
                                              int             allpixels,
                                              int             xcoordinate,
                                              int             ycoordinate,
                                              int             binningtype,
                                              int             binninglevel,
                                              /* Estimation analysis settings... */
                                              float          *trans,
                                              int             transbins,
                                              float           transbinwidth,
                                              int             fitstart,
                                              int             fitend,
                                              int             nphotons,
                                              int             modeltype,
                                              /* Estimated instrument response... */
                                              BayesInstrRsp_t *instr);
#endif

int bayes_CopyInstrRspConfigParams(BayesInstrRsp_t *source,
                                   BayesInstrRsp_t *destination);

int bayes_CheckForDifferentInstrRspConfigParams(BayesInstrRsp_t *source,
                                                BayesInstrRsp_t *destination);
int bayes_SortInstrRspComponentsByWeight(BayesInstrRsp_t *instr);


#if 1
int bayes_FitPredictedApproxInstrRsp(float           *fitted,
                                     int              nbins,
                                     float            binwidth,
                                     BayesInstrRsp_t *instr);
#else
int bayes_FitPredictedApproxInstrRsp(float *fitted,
                                     int    nbins,
                                     float  binwidth,
                                     float  gamma1,
                                     float  delta1,
                                     float  sigma1,
                                     float  delay1,
                                     float  gamma2,
                                     float  delta2,
                                     float  sigma2,
                                     float  delay2,
                                     float  gamma3,
                                     float  delta3,
                                     float  sigma3,
                                     float  delay3);
#endif

/*
int bayes_DirectMonoExpMostProbW0W1WithInstrRspOptimizationBinLikelihood(int   *data,
                                                                         int    nbins,
					                                                     int    nphotons,
   				                                                         float *w0,
					                                                     float *w1,
					                                                     float *delay,
					                                                     float *width,
					                                                     float  interval,
                                                                         float  modulationperiod,
					                                                     float  alpha,
					                                                     float  precision);
*/

int bayes_DirectInstrRspAndMonoExpOptimization(/* Data in... */
                                               int                      *data,
                                               int                       nbins,
                                               int                       fitstart,
                                               double                   *binwalls,
                                               int                       nphotons,
                                               /* Estimates out... */
                                               double                   *w0,
                                               double                   *w1,
                                               BayesInstrRsp_t          *instr,
                                               double                   *minuslogprob,
                                               /* Instrument... */
                                               double                    interval,
                                               double                    modulationperiod,
                                               double                    alpha);


int bayes_DirectInstrRspAndMultiExpOptimization(/* Data in... */
                                                int                      *data,
                                                int                       nbins,
                                                int                       fitstart,
                                                double                   *binwalls,
                                                int                       nphotons,
                                                /* Estimates out... */
                                                int                       ndecays,
                                                double                   *weights,
                                                double                   *taus,
                                                BayesUserFixedParams_t   *paramfixing,
                                                BayesInstrRsp_t          *instr,
                                                double                   *minuslogprob,
                                                /* Instrument... */
                                                double                    interval,
                                                double                    modulationperiod,
                                                double                    alpha);


int bayes_InstrRspCoarseGuessValuesFromLoadedInstr(float *instr,
                                                   float  binwidth,
                                                   int    nbins,
                                                   float *delay,
                                                   float *width);

int bayes_InstrRspCoarseGuessWidthPriorValuesFromLoadedInstr(float *instr,
                                                             float  binwidth,
                                                             int    nbins,
                                                             float *width_ave,
                                                             float *width_sd);

#if 0 //configurable prior - code needs finishing
double MinusLogGaussianPrior(double  x, 
                int     nparams,
                double *params,
                void   *config);
#endif

#endif /* BAYES_INSTR_RSP_ANALYSIS_H */

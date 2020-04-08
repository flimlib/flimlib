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

#ifndef BAYES_MONO_EXP_RAPID_ANALYSIS_H
#define BAYES_MONO_EXP_RAPID_ANALYSIS_H


int bayes_RapidMonoExpDirectMostProbW0W1( int                           *data,
                                          int                            nbins,
                                          int                            fitstart,
                                          int                            nphotons,
                                          double                        *w0,
                                          double                        *w1,
                                          float                         *val,
                                          BayesInstrRsp_t               *instr,
                                          float                          interval,
                                          float                          modulationperiod,
                                          float                          alpha,
                                          BayesRapidMonoExpValueStore_t *grid,
                                          BayesAveErrDistn_t            *distr);



int bayes_RapidMonoExpAvgAndErrors( int                           *data,
                                    int                            nbins,
                                    int                            fitstart,
                                    int                            nphotons,
                                    double                        *w0_mp,
	                                double                        *w1_mp,
                                    double                        *w0_ave,
	                                double                        *w1_ave,
                                    double                        *dw0,
	                                double                        *dw1,
                                    BayesInstrRsp_t               *instr,
                                    float                          interval,
                                    float                          modulationperiod,
                                    float                          alpha,
	                                float                          precision,
	                                int                            p,/* hmm, this is actually nphotons */
	                                int                            quick,
	                                BayesAveErrDistn_t            *distr,
                                    BayesRapidMonoExpValueStore_t *grid,
                                    float                         *value);


#if 0
int  bayes_RapidMonoExpIndirectMostProbableW0W1(  int   *data,
                                                  int    nbins,
                                                  int    nphotons,
                                                  float *w0,
                                                  float *w1,
                                                  //float *val,
                                                  float                          delay,
                                                  float                          width,
                                                  float                          interval,
                                                  float                          modulationperiod,
                                                  float                          alpha,
                                                  BayesRapidMonoExpValueStore_t *grid);
#endif

#endif /* BAYES_MONO_EXP_RAPID_ANALYSIS_H */

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

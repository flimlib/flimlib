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
#ifndef BAYES_DISTR_FCT_LIKELIHOODS_H
#define BAYES_DISTR_FCT_LIKELIHOODS_H


#include "bayes_InstrRspAnalysis.h"

#define RESULT_VALID                     0
#define RESULT_INVALID_ARITHMETIC_ERROR -1


void bayes_CreateAndPopulateVectorInstrRspConstantGammaTilde(double **gammatilde, BayesInstrRsp_t *instr);


void bayes_ComputeArrBinLikelihoodConstantUpsilon1(double **upsilon1,/* array of length nbins+1 */
                                                   int     *data,
                                                   int      nbins,
                                                   double   interval,
                                                   double   width,
                                                   double   delay);

void bayes_CreateAndPopulateMatrixArrBinLikelihoodConstantUpsilon1(double        ***upsilon1,
                                                                   double          *binwalls,
                                                                   int              nbins,
                                                                   int             *data,
                                                                   int              ellmax,
                                                                   double           modperiod,
                                                                   BayesInstrRsp_t *instr);

double bayes_ComputeArrBinLikelihoodConstantUpsilon3(double repetitonperiod,
                                                     double tau);

#if 0
int  bayes_ArrBinLikelihoodsGivenTau(double *likelihoods,
                                     double *upsilon1_g,
                                     int    *data,
                                     int     nbins,
                                     double  interval,
                                     double  width,
                                     double  delay,
                                     double  tau,
                                     double  modperiod);
#else
int  bayes_ArrBinLikelihoodsGivenTau(double          *likelihoods,
                                     double          *binwalls_g,
                                     double          *upsilon1_g,
                                    // ArrLikelihoodConstants_t *constants_g,
                                     int             *data,
                                     int              nbins,
                                     double           interval,
                                     double           modperiod,
                                     BayesInstrRsp_t *instr,
                                     double           tau);
#endif


int  bayes_PopulateBinWallsVectorUniformIntervals(double *binwalls, /* Vector running from '0' to 'nbins'... */
                                                  int     nbins,
                                                  double  interval);


int bayes_ComputeBinLikelihoodsFromWeightsAndFluorescencePhotonLikelihoods(double           *binlikelihoods,
                                                                           int               nbins,
                                                                           double           *binwalls,
                                                                           int               ndecays,
                                                                           double          **fluorescencephotonlikelihoods,
                                                                           double           *weights,
                                                                           BayesInstrRsp_t  *instr,
                                                                           double            interval);


int bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(double          *fluorescencephotonlikelihoods,
                                                               int              nbins,
                                                               double          *binwalls,
                                                               int             *data,
                                                               double           interval,
                                                               double           modperiod,
                                                               BayesInstrRsp_t *instr,
                                                               double           tau,
                                                               int              ndecays,
                                                               double          *weights,
                                                               double          *taus);


int bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(double          *normalisation,
                                                              double           interval,
                                                              double           modperiod,
                                                              double           dithertime,
                                                              BayesInstrRsp_t *instr,
                                                              int              ndecays,
                                                              double          *weights,
                                                              double          *taus);


#endif /* BAYES_DISTR_FCT_LIKELIHOODS_H */

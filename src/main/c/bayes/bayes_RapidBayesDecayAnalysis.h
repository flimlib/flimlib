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
#ifndef BAYES_RAPID_BAYES_DECAY_ANALYSIS_H
#define BAYES_RAPID_BAYES_DECAY_ANALYSIS_H


#include "bayes_Types.h"


#define BAYES_RAPID_RESULT_NO_ERROR     0
#define BAYES_RAPID_RESULT_USER_CANCEL -99

#define BAYES_RAPID_GRID_MONO        1
#define BAYES_RAPID_GRID_BI          2
#define BAYES_RAPID_GRID_MONO_AND_BI 3


struct BayesRapidDiscreteGridSettings
{
    int     ntaus;         /* Number of lifetime values that span the grid...       */
    double *tau;           /* Lifetime values that span the grid...                 */
    int     nweights;      /* Number of weight values that span the grid...         */
    double *weight;        /* Weight values that span the grid...                   */
    double  backgroundmin; /* User defined minimum allowed background proportion... */
    double  backgroundmax; /* User defined maximum allowed background proportion... */
};
typedef struct BayesRapidDiscreteGridSettings BayesRapidGridSettings_t;


struct BayesRapidLikelihoodValuesContainer
{
    int     valid;
    double  tau;
    double *fluorescencedecayphotonlikelihoodsgiventau;
};    
typedef struct BayesRapidLikelihoodValuesContainer BayesRapidLikelihoodValues_t;


struct BayesRapidMonoExpDiscreteLikelihoodValueContainer
{
    double  tau;
    double  w0;
    double *logphotonlikelihoodgiventauandw0;
};    
typedef struct BayesRapidMonoExpDiscreteLikelihoodValueContainer BayesRapidMonoExpDiscreteValues_t;


struct BayesRapidMonoExpDiscreteValueStore
{
    int                                 validgrid;
    BayesRapidGridSettings_t           *settings;
    int                                 nbins;
    BayesInstrRsp_t                     instr;
    double                              interval;
    double                              modulationperiod;
    BayesRapidLikelihoodValues_t       *fluorescencelikelihoods;
    int                                 nlikelihoods;
    int                                 nvalidlikelihoods;
    int                                 ninvalidlikelihoods;
    BayesRapidMonoExpDiscreteValues_t **likelihoods;
    int                                 nstates;
    int                                 nvalidstates;
    int                                 ninvalidstates;
};
typedef struct BayesRapidMonoExpDiscreteValueStore BayesRapidMonoExpValueStore_t;


struct BayesRapidBiExpDiscreteLikelihoodValueContainer
{
    int     valid;
    double  weights[3];
    double  taus[3];
    double *logphotonlikelihoodgiventausandweights;
};    
typedef struct BayesRapidBiExpDiscreteLikelihoodValueContainer BayesRapidBiExpDiscreteValues_t;


struct BayesRapidBiExpDiscreteValueStoreMemoryPool
{
    int                               nvalid;
    BayesRapidBiExpDiscreteValues_t  *likelihoods;
    double                          **loglikelihoods;
};
typedef struct BayesRapidBiExpDiscreteValueStoreMemoryPool BayesRapidBiExpMemoryPool_t;


struct BayesRapidBiExpDiscreteValueStoreMemoryControl
{
    int                          npools;
    BayesRapidBiExpMemoryPool_t *pools;
};
typedef struct BayesRapidBiExpDiscreteValueStoreMemoryControl BayesRapidBiExpMemoryControl_t;


struct BayesRapidBiExpDiscreteValueStore
{
    int                                  validgrid;
    BayesRapidGridSettings_t            *settings;
    int                                 *low;
    int                                 *high;
    int                                  nbins;
    BayesRapidLikelihoodValues_t        *fluorescencelikelihoods;
    int                                  nlikelihoods;
    int                                  nvalidlikelihoods;
    int                                  ninvalidlikelihoods;
    BayesRapidBiExpDiscreteValues_t *****likelihoods;
    BayesInstrRsp_t                      instr;
    double                               interval;
    double                               modulationperiod;
    int                                  nstates;
    int                                  nvalidstates;
    int                                  ninvalidstates;
    double                               megabytes;
    BayesRapidBiExpMemoryControl_t      *memory;
};
typedef struct BayesRapidBiExpDiscreteValueStore BayesRapidBiExpValueStore_t;


struct BayesRapidLikelihoodsDiscreteValueStore
{
    int                           valid;
    BayesRapidGridSettings_t     *settings;
    BayesRapidLikelihoodValues_t *fluorescencelikelihoods;
    int                           nbins;
    BayesInstrRsp_t               instr;
    double                        interval;
    double                        modulationperiod;
    int                           nlikelihoods;
    int                           nvalidlikelihoods;
    int                           ninvalidlikelihoods;
};
typedef struct BayesRapidLikelihoodsDiscreteValueStore BayesRapidLikelihoodsValueStore_t;


struct BayesRapidDiscreteValueStore
{
    int                                validlikelihoodsgrid;
    BayesRapidLikelihoodsValueStore_t *likelihoodsvaluestore;
    int                                validmonoexpgrid;
    BayesRapidMonoExpValueStore_t     *monoexpvaluestore;
    int                                validbiexpgrid;
    BayesRapidBiExpValueStore_t       *biexpvaluestore;
};
typedef struct BayesRapidDiscreteValueStore BayesRapidValueStore_t;


BayesRapidValueStore_t *bayes_GetRapidValueStorePtr(void);
int bayes_InitializeRapidValueStore(BayesRapidValueStore_t *store);
int bayes_DestroyRapidValueStore(BayesRapidValueStore_t *store, int gridtype);
void bayes_InvalidateRapidValueStore(BayesRapidValueStore_t *store);
void bayes_InvalidateRapidBiExpValueStore(BayesRapidValueStore_t *store);
int bayes_CheckForValidRapidValueStore(BayesRapidValueStore_t *store, int gridtype);

int bayes_CreateRapidValueStore(BayesRapidValueStore_t *store,
                                int                     updatetype,
                                /* Mono-exponential settings... */
                                int                     ntaus,
                                double                 *taus,
                                int                     nweights,
                                double                 *weights,
                                double                  backgroundmax,
                                double                  backgroundmin,
                                /* Bi-exponential settings... */
                                int                     ntaus_bi,
                                double                 *taus_bi,
                                int                     nweights_bi,
                                double                 *weights_bi,
                                double                  backgroundmin_bi,
                                double                  backgroundmax_bi,
                                int                    *low,
                                int                    *high,
                                /* General (common) settings... */
                                int                     nbins,
                                double                 *binwalls,
                                BayesInstrRsp_t        *instr,
                                double                  interval,
                                double                  modulationperiod);


int bayes_DetermineIfBayesGridUpdateReqd(BayesRapidValueStore_t *store,
                                         int                     updatetype,
                                         /* Mono-exp... */
                                         int                     ntaus,
                                         double                 *taus,
                                         int                     nweights,
                                         double                 *weights,
                                         double                  backgroundmin,
                                         double                  backgroundmax,
                                         /* Bi-exp... */
                                         int                     ntausbi,
                                         double                 *tausbi,
                                         int                     nweightsbi,
                                         double                 *weightsbi,
                                         double                  backgroundminbi,
                                         double                  backgroundmaxbi,
                                         int                    *lowbi,
                                         int                    *highbi,
                                         /* Instrumentation and data... */
                                         int                     nbins,
                                         BayesInstrRsp_t        *instr,
                                         double                  interval,
                                         double                  modulationperiod);


int bayes_MapWeightValueToClosestRapidGridPoint(double value, int nweights, double *weights);
int bayes_MapLifetimeValueToClosestRapidGridPoint(double value, int ntaus, double *taus);

int bayes_RapidBiExpDetermineGridSizeAdv(/* Input... */
                                      int     ntaus,
                                      double *taus,
                                      int     nweights,
                                      double *weights,
                                      double  backgroundmin,
                                      double  backgroundmax,
                                      int    *low,
                                      int    *high,
                                      int     nbins,
                                      /* Output... */
                                      int    *npts,
                                      int    *nptsvalid,
                                      double *megabytesreqd);

int bayes_OutputGridToFile(char                   *filename,
                           char                   *comment,
                           BayesRapidValueStore_t *grid);


#endif /* BAYES_RAPID_BAYES_DECAY_ANALYSIS_H */







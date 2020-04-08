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
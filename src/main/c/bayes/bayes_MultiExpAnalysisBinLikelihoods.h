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
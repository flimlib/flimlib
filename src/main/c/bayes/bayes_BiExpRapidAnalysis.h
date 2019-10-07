#ifndef BAYES_BI_EXP_RAPID_ANALYSIS_H
#define BAYES_BI_EXP_RAPID_ANALYSIS_H


int bayes_RapidBiExpMostProbWeightsAndTaus(int                          *data,
                                           int                           nbins,
                                           int                           fitstart,
                                           double                       *binwalls,
                                           int                          *nphotons,
                                           int                           ndecays,
                                           double                       *weights_mp,
                                           double                       *taus_mp,
                                           double                       *weights_ave,
                                           double                       *taus_ave,
                                           double                       *weights_err,
                                           double                       *taus_err,
                                           BayesUserFixedParams_t       *paramfixing,
                                           double                        interval,
                                           double                        modulationperiod,
                                           BayesInstrRsp_t              *instr,
                                           double                        alpha,
                                           BayesRapidValueStore_t       *grid,
                                           double                       *val,
                                           BayesProbDistn_t             *distribution);


int bayes_RapidBiExpHyperParamOptimization(int                    *data,
                                           int                     nbins,
                                           int                     fitstart,
                                           int                     nphotons,
                                           double                 *binwalls,
                                           BayesInstrRsp_t        *instr,
                                           float                   interval,
                                           float                   modulationperiod,
                                           float                  *alphastar,
                                           float                   alphamin,
	                                       float                   precision,
                                           float                  *value,
                                           BayesRapidValueStore_t *grid);


#endif /* BAYES_BI_EXP_RAPID_ANALYSIS_H */
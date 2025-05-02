/*-
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2025 University of Oxford and Board of Regents of the
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

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
#ifndef BAYES_ANNALYSIS
#define BAYES_ANNALYSIS

#include "bayes/bayes_Interface.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * This function configure the search grid for Bayes_fitting_engine() by
 * specifying the lower and upper bound for each parameters (z, A, tau, ...).
 * The configuration takes effect only if nparam is 3 or 5 while the mode of
 * fitting is single or double exponential. The behavior is undefined if any of
 * the bounds are less than 0 or if parammax < parrammin.
 * 
 * \param[in] parammin The minimum of Bayesian search grid.
 * \param[in] parammax The maximum of Bayesian search grid.
 * \param[in] nparam The number of parameters (length of parammin and parammax).
 */
void Bayes_set_search_grid(float parammin[], float parammax[], int nparam);

/*=================================================================================*/
/*                                                                                 */
/*                 MAIN FITTING FUNCTIONS                                          */
/*                                                                                 */
/*=================================================================================*/

/**
* The main/convinient entry point for Bayesian parameter estimation. Wraps bayes_DoBayesFitting().
* This is designed to look similar to the entry point to other alogorithms (RLD, LMA), using similar var names and types.
* The Bayes algorithm does not accept the IRF (instr) here as othe alogorithms do, it must be loaded seperately.
*
* \param[in] xincr The time increment inbetween the values in the trans array.
* \param[in] laser_period Laser repetition or modualtion period.
* \param[in] trans The transient (time resolved) signal to be analysed, the 'data'.
* \param[in] ndata The number of data points.
* \param[in] fit_start The index into the trans array marking the start to the data to be used in the fit.
* \param[in] fit_end The index into the trans array marking the end of the data to be used in the fit.
* \param[in,out] param An array of parameters, the order of which must match fitfunc. Provide parameter estimates, these are overridden with the fitted values.
* \param[in] paramfree An array indicating which parameters are free (1), fixed (0)
* \param[in] nparam The number of parameters.
* \param[in] modeltype The fit model type, e.g. FIT_MONOEXP.
* \param[out] fitted An array containing values fitted to the data, the 'fit'. Fit points are coincident in time with the data points.
* \param[out] residuals An array containing the difference between the data and the fit. Can be NULL if not required.
* \param[out] error The estimated error in the parameters
* \param[out] minuslogprob The resulting negated log probability (-log(P)) of the estimate.
* \param[out] nphotons The total number of photons included in the fit.
* \param[out] chisq The resulting raw chi squared value of the fit. To get the reduced chisq, divide by the degrees of freedom (fit_end - fit_start - nparam). Requires residuals array. Can be NULL if not required.
* \return A negative error code on failure, non-negative on success.
*/
int Bayes_fitting_engine(/* Data in... */
                        float xincr,
                        float *trans,
                        int ndata,
                        int fit_start,
                        int fit_end,
                        float laser_period,
                        float instr[],
                        int ninstr,
                        /* Model... */
                        float param[],
                        int paramfree[],
                        int nparam,
                        /* Data out... */
                        float *fitted,
                        float *residuals,
                        float *error,
                        /* Metadata output */
                        float *minuslogprob,
                        int *nphotons,
                        float *chisq);

#ifdef __cplusplus
}
#endif

#endif /* BAYES_ANNALYSIS */

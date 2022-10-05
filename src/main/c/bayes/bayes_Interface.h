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
#ifndef BAYES__INTERFACE_H
#define BAYES__INTERFACE_H

#include "bayes_Types.h"
#include "bayes_Sizes.h"

#include "bayes_DataManagement.h"
#include "matrices.h"
#include "safe_globals.h"

#include "bayes_ModelTransformTools.h"
#include "bayes_DistributionFctsBinLikelihoods.h"
#include "bayes_MonoExpAnalysisBinLikelihoods.h"
#include "bayes_MultiExpAnalysisBinLikelihoods.h"
#include "bayes_InstrRspAnalysis.h"
#include "bayes_RapidBayesDecayAnalysis.h"
#include "bayes_MonoExpRapidAnalysis.h"
#include "bayes_BiExpRapidAnalysis.h"
#include "bayes_MultiExpRapidAnalysis.h"
#include "bayes_ModelSelection.h"

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
* \return An error code, 0 = success.
*/

int bayes_fitting_engine(   /* Data in... */
							float xincr,
							float laser_period,
							float *trans, 
							int ndata, 
							int fit_start, 
							int fit_end,
							/* Model... */
							float param[],
							int paramfree[], 
							int nparam, 
							int modeltype, 
							/* Data out... */
							float *fitted, 
							float *residuals, 
							float error[],
							/* Metadata output */
							float *minuslogprob, 
							int *nphotons,
							float *chisq);


/**
* The function to perform IRF estimation. Bayes requires a parameterised prompt in terms of gaussian components.
* These components are defined by a width, weight, delay (time position) and cut off (since a Gaussian will extend into negative time if not cut off)
* The estimation requires the prompt signal to estimate and a representative transient
*
* \param[in] trans The transient (time resolved) signal to be analysed, the 'data'.
* \param[in] ndata The number of data points in the transient.
* \param[in] xincr The time increment inbetween the values in the trans array.
* \param[in] fit_start The index into the trans array marking the start to the data to be used in the fit.
* \param[in] fit_end The index into the trans array marking the end of the data to be used in the fit.
* \param[out] nphotons The total number of photons included in the fit.
* \param[in] nprompt The number of data points in the prompt signal.
* \param[in] promptbinwidth The time increment inbetween the values in the prompt array.
* \param[in] modeltype The fit model type, e.g. FIT_IRFANDMONOEXP.
* \param[out] instr Fitted IRF values
* \param[in] laser_period Laser repetition or modualtion period.
* \param[out] param_mp An array of parameters that fit the trans dataset.
* \param[out] fitted An array containing values fitted to the data, the 'fit'. Fit points are coincident in time with the trans data points.
* \return An error code, 0 = success.
*/

int bayes_DoBayesInstrRspEstimation(/* Data in... */
                                        float                    *trans,
                                        int                       ndata,
                                        float                     xincr,
                                        int                       fitstart,
                                        int                       fitend,
                                        int                      *nphotons,
                                        /* Loaded prompt in... */
                                        float                    *prompt,
                                        int                       nprompt,
                                        float                     promptbinwidth,
                                        /* Decay model... */
                                        int                       modeltype,
                                        /* Instrument... */
                                        BayesInstrRsp_t          *instr,
                                        BayesIrEstConfig_t*       BayesIrEstConfig,
                                        float                     laser_period,
                                        /* Data out... */
										float                    *param_mp, /* Parameters input using conventional model and values (i.e. Z,A1,tau1,A2,tau2...) */
										float                    *fitted);


/*=================================================================================*/
/*                                                                                 */
/*                 INTERNAL FITTING FUNCTION                                       */
/*                                                                                 */
/*=================================================================================*/

int bayes_DoBayesFitting(    /* Data in... */
	float                    *trans,
	int                       transbins,
	float                     transbinwidth,
	int                       fitstart,
	int                       fitend,
	int                      *nphotons,
	/* Model... */
	int                       modeltype,
	int                       nparams,
	/* Instrument response... */
	BayesInstrRsp_t          *instr,
	float                     modulationperiod,
	/* Estimates... */
	int                      *param_free,
	float                    *param_mp,     //Z,A1,tau1,A2,tau2...
	float                    *param_err,
	/* Data out... */
	float                    *fitted,
	BayesAveErrDistn_t       *distr,
	int                       distr_xparam,
	int                       distr_yparam,
	float                    *val,
	/* Settings... */
	float                     precision,
	int                       quick,
	int                       rapidanalysis,
	BayesRapidValueStore_t   *rapid);

/*=================================================================================*/
/*                                                                                 */
/*                 ACCESS FUNCTIONS FOR ADVANCED CONFIGF VARIABLE VALUES           */
/*                                                                                 */
/*=================================================================================*/

void bayes_GetInstrRspParamValues(BayesInstrRsp_t *instr, BayesIrEstConfig_t* BayesIrEstConfig);
BayesIrEstConfig_t bayes_GetIrEstConfig(void);

float bayes_GetConfigParameterValueModulationPeriod(void);
void bayes_SetConfigParameterValueModulationPeriod(float val);
int bayes_UseRepetitionEffectsInAnalysis(void);
int bayes_GetIncludeRepetitionEffectsInAnalysisFlag(void);

/* Transient re-binning... */
void bayes_TransientRebinning(float *transin, int binsin, float *transout, int binsout);
int bayes_GetBayesTransientRebinningActiveFlag(void);
int bayes_GetBayesTransientRebinningFactor(void);

/* General minimization algorithm settings... */
int bayes_BiExpConfigGetNumberOfRestarts(void);
int bayes_ConfigGetMinimizationAlgorithm(void);

/* Mono-exp and downhill simplex configuration... */
float bayes_MonoExpConfigGetDownhillSimplexPrecision(void);

/* Bi-exp and simulated annealing configuration... */
int bayes_BiExpConfigGetRapidGridSearchInitialisationsFactor(void);
int bayes_BiExpConfigRapidGetGridSearchLocalisedExhaustiveSearchDelta(void);
float bayes_BiExpConfigGetSimAnnealTempReductionSchedule(void);
int bayes_BiExpConfigGetSimAnnealTempReductionNoOfSteps(void);
int bayes_BiExpConfigGetSimAnnealItersAtEachTemp(void);
float bayes_BiExpConfigGetSimAnnealFracConvergenceTolerance(void);

/* Instrument response determination configuration... */
float bayes_InstrConfigGetSimAnnealStartingTemp(void);
float bayes_InstrConfigGetSimAnnealTempReductionSchedule(void);
int bayes_InstrConfigGetSimAnnealTempReductionNoOfSteps(void);
int bayes_InstrConfigGetSimAnnealItersAtEachTemp(void);
float bayes_InstrConfigGetSimAnnealFracConvergenceTolerance(void);
int bayes_InstrConfigGetNumberOfRestarts(void);
float bayes_InstrConfigGetDownhillSimplexPrecision(void);
int bayes_InstrConfigGetMinimizationAlgorithm(void);

/* Rapid Bayesian analysis configuration... */
int bayes_RapidGetUseRapidBayesFlag(void);
int bayes_InitDiscreteGridForRapidBayes(void);
int bayes_ConfigureBayesianRapidGrid(int updatetype, float xincr, int fitend, BayesIrEstConfig_t* BayesIrEstConfig);
int bayes_DestroyDiscreteGridForRapidBayes(void);

void bayes_InvalidateBayesDiscreteGridsDueToParameterValueChange(void);
void bayes_RapidBiExpInvalidateGridDueToParameterValueChange(void);
int bayes_CheckForValidBayesDiscreteGrids(int gridtype);

// int bayes_RapidBiExpDetermineGridSize(int fitend);

int bayes_FitTypeToRapidGridUpdateType(int fittype);

/* General decay model configuration... */
int bayes_ConfigDecayModelMaxNumberOfModulationPeriodsConsidered(void);

/* Hyperparameter */
int bayes_ConfigUseFullBayesianHyperParamDetermination(void);

/* Model transform relations... */
float bayes_ToBayesModelWeightFromParamAAndTau( float A_i,
                                                float tau_i,
                                                float binsize,
                                                float interval,
                                                float delay,
                                                int   nphotons);

int bayes_ConvertConventionalToBayesModelParamValues(/* Model... */
                                                         int    modeltype,
                                                         int    nparams,
                                                         float *params_bayes,
                                                         float *params_conventional,
                                                         /* Instrument response... */
                                                         //float  delay,
                                                         //float  width,
                                                         BayesInstrRsp_t *instr,
                                                         /* Transient... */
                                                         float *trans,
                                                         int    transbins,
                                                         float  transbinwidth,
                                                         int    fitstart,
                                                         int    fitend);

int bayes_ConvertBayesModelToConventionalParamValues(/* Model... */
                                                         int    modeltype,
                                                         int    nparams,
                                                         float *params_bayes,
                                                         float *params_conventional,
                                                         /* Instrument response... */
                                                         //float  delay,
                                                         //float  width,
                                                         BayesInstrRsp_t *instr,
                                                         /* Transient... */
                                                         float *trans,
                                                         int    transbins,
                                                         float  transbinwidth,
                                                         int    fitstart,
                                                         int    fitend);

int bayes_ComputeBayesHyperParamsFromData(/* Data in... */
                                              float *trans,
                                              int    transbins,
                                              float  transbinwidth,
                                              int    fitstart,
                                              int    fitend,
                                              /* Hyperparameter values out... */
                                              int    nhyperparams,
                                              float *hyperparams);

/**
 * TODO: Thread Unsafe
 */
int bayes_UpdateEstimatedInstrRspDetailsStore(/* Image/prompt data used... */
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

/** 
	Calcualate the equivalent raw chisq for Bayes fit. Fills residuals array as a side effect. 

	\param[in] y The data
	\param[in] fitted The fit
	\param[in/out] residuals The difference between data and fit
	\param[in] fit_start Defines the start of the data in y
	\param[in] fit_end Defines the end of the data in y
*/

float bayes_CalculateResidualsAndEquivalentChisq(float y[], float fitted[], float *residuals, int fit_start, int fit_end);

/* Error management */
const char* bayes_GetBayesErrorDescription(int error);


/* Model selection */

int bayes_DoDecayModelSelection( /* Data in... */
                                     float                    *trans,
                                     int                       transbins,
                                     float                     transbinwidth,
                                     int                       fitstart,
                                     int                       fitend,
                                     int                      *nphotons,
                                     /* Instrumentation... */
                                     BayesInstrRsp_t          *instr,
                                     float                     modulationperiod,
                                     /* Model selection... */
                                     int                       modeltype,
                                     int                       nparams,
                                     /*  RLD estimates, used for initialising search... */
                                     float                    *param_rld,
                                     /* Data out... */
                                     float                    *decaymodellikelihoods,
                                     BayesParamValsAndFit_t   *paramvalsandfits,
                                     /* Configuration... */
                                     int                       rapidanalysis,
                                     BayesRapidValueStore_t   *rapidgrid);


int bayes_PerformBayesParameterEstimation(/* Data... */
                                              int                      *data,
                                              int                       nbins,
                                              int                       fitstart,
                                              double                   *binwalls,
                                              int                       nphotons,
                                              /* Instrument... */
                                              BayesInstrRsp_t          *instr,
                                              float                     modulationperiod,
                                              float                     interval,
                                              /* Decay model... */
                                              int                       modeltype,
                                              int                       ndecays,
                                              BayesUserFixedParams_t   *paramfixing,
                                              double                    alpha,
                                              /* Estimates... */
                                              double                   *weights_mp,
                                              double                   *taus_mp,
                                              double                   *weights_ave,
                                              double                   *taus_ave,
                                              double                   *weights_err,
                                              double                   *taus_err,
                                              float                    *minuslogprob,
                                              BayesAveErrDistn_t       *probdistr,
                                              /* Bayesian analysis settings... */
                                              int                       rapidanalysis,
                                              BayesRapidValueStore_t   *rapid);


int bayes_CheckParameterValueFixingForBayesFitting(BayesUserFixedParams_t *paramfixing,
                                                       int                     nparams,
                                                       int                    *paramsfree,
                                                       float                  *params,
                                                       int                     nbins,
                                                       int                     fitstart,
                                                       double                 *binwalls,
													   int					   nphotons,
                                                       double                  interval,
                                                       double                  modulationperiod,
                                                       BayesInstrRsp_t        *instr);

int bayes_FreeParameterValueFixingForBayesFitting(BayesUserFixedParams_t *paramfixing,
                                                      int                     ndecays,
													  int					  nbins);


#endif /* BAYES__INTERFACE_H */

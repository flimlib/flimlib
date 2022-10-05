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
#ifndef BAYES_TYPES_H
#define BAYES_TYPES_H

#define TRUE 1
#define FALSE 0

#ifndef NULL
#define NULL 0
#endif

// Fit types
#define FIT_UNKNOWN        0
#define FIT_MONOEXP        1
#define FIT_BIEXP          2
#define FIT_TRIEXP         3
#define FIT_STRETCHED      4
#define FIT_MULTIEXP       5
#define FIT_MODELSELECTION 6
#define FIT_IRFANDMONOEXP  7
#define FIT_IRFANDBIEXP    8
#define FIT_GLOBAL         10  // To indicate we've done a global fit

// Fit algorithms
#define ALG_3INTEGRAL  0
#define ALG_MARQUARDT  1
#define ALG_BAYES 	   2
#define ALG_PHASOR 	   3

/* Instrument response... */
typedef struct
{
	int bayesirnumcomponents;
	double bayesirweight1;
	double bayesircutoff1;
	double bayesirsigma1;
	double bayesiruc1;
	double bayesirweight2;
	double bayesircutoff2;
	double bayesirsigma2;
	double bayesiruc2;
	double bayesirweight3;
	double bayesircutoff3;
	double bayesirsigma3;
	double bayesiruc3;
} BayesIrEstConfig_t;

/* Rapid grids */
typedef struct
{
	int bayesrapidtaupts;
	float bayesrapidtaulow;
	float bayesrapidtauhigh;
	int bayesrapidwpts;
	float bayesrapidwlow;
	float bayesrapidwhigh;
	float bayesrapidbghigh;
} BayesMonoRapidGridConfig_t;

typedef struct
{
	int bayesrapidbitaupts;
	float bayesrapidbitaulow;
	float bayesrapidbitauhigh;
	int bayesrapidbiweightpts;
	float bayesrapidbiweightlow;
	float bayesrapidbiweighthigh;
	float bayesrapidbibgmin;
	float bayesrapidbibgmax;

	float bayesrapidbiw0low;
	float bayesrapidbiw1low;
	float bayesrapidbiw2low;
	float bayesrapidbitau1low;
	float bayesrapidbitau2low;

	float bayesrapidbiw0high;
	float bayesrapidbiw1high;
	float bayesrapidbiw2high;
	float bayesrapidbitau1high;
	float bayesrapidbitau2high;
} BayesBiRapidGridConfig_t;

/*====================================================================*/
/*    Generic multi-exponential probability distribution container    */
/*====================================================================*/

struct BayesProbDistnStateAndVal
{
    int     ndim;
    int    *state;
    double  val;
};    
typedef struct BayesProbDistnStateAndVal BayesProbDistnStateAndVal_t;

struct BayesProbDistributionMarginal
{
    int     indexlow;
    int     indexhigh;
    double *marginal;
};    
typedef struct BayesProbDistributionMarginal BayesProbMarginal_t;

struct BayesProbDistribution
{
    int                          ndim;
    int                          nweights;
    int                          ntaus;
    double                       dweight;
    double                       dtau;
    int                         *statesmin;
    int                         *statesmax;
    int                          nstates;
    BayesProbDistnStateAndVal_t *statesandvals;
    BayesProbMarginal_t         *marginals;
};    
typedef struct BayesProbDistribution BayesProbDistn_t;





struct BayesFunctionType
{
    double (*Function)(double *vars, void *params);
    void   *Params;
};
typedef struct BayesFunctionType BayesFct_t;



struct BayesContainerArrLikelihoodConstants
{
    /* Pre-calculated fixed instrument values, dependent only the instrument response */
    /* parameter values and the number of bins for the measured transient...          */
    double *upsilon1;
    double *gammatilde;
};
typedef struct BayesContainerArrLikelihoodConstants ArrLikelihoodConstants_t;



/* Approximated instrument response container types... */
#define BAYES_INSTR_RSP_MAX_COMPONENTS 3
#define BAYES_INSTR_RSP_PRIOR_MAX_PARAMS 5

struct BayesContainerInstrRspCutGaussianParams
{
    double weight;
    double width;
    double delay;
    double cutoff;
};
typedef struct BayesContainerInstrRspCutGaussianParams BayesInstrRspParams_t;

struct BayesContainerInstrRspParamFixing
{
    int weightfixed;
    int widthfixed;
    int delayfixed;
    int cutofffixed;
};
typedef struct BayesContainerInstrRspParamFixing BayesInstrRspParamFixing_t;

struct BayesContainerInstrRspParamPriorEvaluator
{
    int    nparams;                                  /* Number of parameters required to define the prior function... */
    double params[BAYES_INSTR_RSP_PRIOR_MAX_PARAMS]; /* Parameters that define the prior function, e.g. an average and width for a Gaussian... */
    double (*funk)(double, int, double *, void *);   /* Generic function, could be a Gaussian for example...  */
};
typedef struct BayesContainerInstrRspParamPriorEvaluator BayesInstrRspParamPriorEvaluator_t;
struct BayesContainerInstrRspParamPrior
{
    int weightprior;
    int widthprior;
    int delayprior;
    int cutoffprior;
    BayesInstrRspParamPriorEvaluator_t weightpriorevaluator;
    BayesInstrRspParamPriorEvaluator_t widthpriorevaluator;
    BayesInstrRspParamPriorEvaluator_t delaypriorevaluator;
    BayesInstrRspParamPriorEvaluator_t cutoffpriorevaluator;
};
typedef struct BayesContainerInstrRspParamPrior BayesInstrRspParamPrior_t;
struct BayesContainerApproximatedInstrRsp
{
    int                   ninstr;
    BayesInstrRspParams_t params[BAYES_INSTR_RSP_MAX_COMPONENTS];
    BayesInstrRspParamFixing_t paramsfixed[BAYES_INSTR_RSP_MAX_COMPONENTS];
    BayesInstrRspParamPrior_t  paramsprior[BAYES_INSTR_RSP_MAX_COMPONENTS];
};
typedef struct BayesContainerApproximatedInstrRsp BayesInstrRsp_t;


struct BayesPsuedoRapidDiscreteLikelihoodValueContainer
{
    int     valid;
    double  tau;
    double *fluorescencedecayphotonlikelihoodsgiventau;
};    
typedef struct BayesPsuedoRapidDiscreteLikelihoodValueContainer BayesPsuedoRapidDiscreteValues_t;


struct BayesContainerUserFixedParamsSettingsAndValues
{
    int                               nparams;
    int                               nparamsuserfixed;
    int                               nweightsuserfixed;
    int                               ntaususerfixed;
    int                              *weightuserfixed;
    int                              *tauuserfixed;
    double                           *weights;
    double                           *taus;
    BayesPsuedoRapidDiscreteValues_t *fluorescencelikelihoods; 
};
typedef struct BayesContainerUserFixedParamsSettingsAndValues BayesUserFixedParams_t;


/*====================================================================*/
/*            Generic decay model parameter values container          */
/*====================================================================*/

struct BayesDecayModelSelParamValuesAndFit
{
    int     ndecays;
    double  weights[3];
    double  taus[3];
    float  *fitted;
    float  *residuals;
};    
typedef struct BayesDecayModelSelParamValuesAndFit BayesParamValsAndFit_t;



#endif /* BAYES_TYPES_H */

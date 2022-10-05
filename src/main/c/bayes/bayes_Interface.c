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
/* File:       bayes_Interface.c                                                   */
/*---------------------------------------------------------------------------------*/
/* Purpose:    Interface for simple integration of Bayesian routines.              */
/*---------------------------------------------------------------------------------*/
/* References: None.                                                               */
/*---------------------------------------------------------------------------------*/
/* Notes:      None.                                                               */
/*=================================================================================*/
/* Revision history:                                                               */
/*---------------------------------------------------------------------------------*/
/* Date   | Modification                                                           */
/*---------------------------------------------------------------------------------*/
/* 091108 | Creation, mrowley.                                                     */
/*---------------------------------------------------------------------------------*/
/* 060719 | Modification for SLIM Curve, P Barber                                  */
/*---------------------------------------------------------------------------------*/
/*        |                                                                        */
/*=================================================================================*/


#include "bayes_Interface.h"

#include "math.h"
#include "stdlib.h"

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#ifndef NULL
#define NULL 0
#endif


struct BayesErrorString
{
    char *OptStr;
};
    
typedef struct BayesErrorString BayesError_t;

// Local variable store for advanced config
static const int bayesrapiduse = 1;
static float bayesrepperiod = 24.752108f;
static const int bayesrep = 1;
static const int bayesadvrebinning = 0;
static const int bayesadvrebinfactor = 1;
// static float bayesadvtempstart = 25.0f;
static const float bayesadvtempred = 0.990000f;
static const int bayesadvtempsteps = 200;
static const int bayesadvtempiters = 50;
static const float bayesadvtolerance = 0.0f;
static const int bayesadvrestarts = 0;
static const float bayesadvmonoprec = 5.0E-5f;
static const int bayesadvalgsel = 1;
static const int bayesadvbistarts = 30;
static const int bayesadvbidelta = 3;
static const float bayesadvirfsatempstart = 25.0f;
static const float bayesadvirfsatempred = 0.990000f;
static const int bayesadvirfsatempsteps = 200;
static const int bayesadvirfsatempiters = 50;
static const float bayesadvirfsatolerance = 0.0f;
static const float bayesadvirfsaprec = 0.000010f;
static const int bayesadvirfalgsel = 1;
// static int bayesadvirfsearchres = 0;

/*	bayes_fitting_engine

	The main entry point to the Bayes code 
	See bayes_Interface.h file for more detailed documentation.
*/

int bayes_fitting_engine(float xincr, float laser_period, float *trans, int ndata, int fit_start, int fit_end,
	float param[], int paramfree[], int nparam, int modeltype, float *fitted, float *residuals, float error[],
	float *minuslogprob, int *nphotons, float *chisq)
{
//	float         *param_ave;
	BayesInstrRsp_t          instr;
	int           ret, i;

	int quick = 0;

	float precision = bayes_MonoExpConfigGetDownhillSimplexPrecision();
	float modulationperiod = laser_period;
	bayes_SetConfigParameterValueModulationPeriod(modulationperiod);
	BayesIrEstConfig_t defaultIREstConfig = bayes_GetIrEstConfig();
	bayes_GetInstrRspParamValues(&instr, &defaultIREstConfig);

	// Setup the Bayesian rapid grid to search
	// TRI2 code does a quick check (now commented) the configure fn checks more things, always do that.
//	needgrid = bayes_RapidGetUseRapidBayesFlag();
//	validgrid = bayes_CheckForValidBayesDiscreteGrids(modeltype);
//	if (needgrid && !validgrid)
//	{
		bayes_ConfigureBayesianRapidGrid(modeltype, xincr, fit_end, &defaultIREstConfig);
//	}

	// Initialise parameters, and which are free to be optimised
	for (i = 0; i<nparam; i++)
	{
		param[i] = 0.0;
		error[i] = 0.0;
		paramfree[i] = 1;
	}

	// Perform the fit
	ret = bayes_DoBayesFitting(trans, ndata, xincr, fit_start, fit_end, nphotons,
		/* Model... */
		modeltype, nparam,
		/* Instrument response... */
		&instr, modulationperiod,
		/* Estimates... */
		paramfree,
		param,     //Z,A1,tau1,A2,tau2...
		error,
		/* Data out... */
		fitted,
		NULL, 0, 0,
		minuslogprob,
		/* Settings... */
		precision,
		quick,
		bayes_RapidGetUseRapidBayesFlag(),
		bayes_GetRapidValueStorePtr());

	// Calculate residuals and equivalent chisq
	if (residuals != NULL)
	{
		float chisq_local = bayes_CalculateResidualsAndEquivalentChisq(trans, fitted, residuals, fit_start, fit_end);
		if (chisq != NULL)
			*chisq = chisq_local;
	}

	return(ret);
}

static const char * const bayes_ErrorDescription[] =
{  // corresponds to the error code defines in bayes_Sizes.h, number of entries is used in the check in bayes_GetBayesErrorDescription below.
    "Bayes: OK",
    "Bayes: Invalid data",
    "Bayes: Invalid data window",
    "Bayes: Invalid model",
    "Bayes: Functionality not supported",
    "Bayes: Invalid fixed parameter value",
    "Bayes: All parameter values are fixed",
	"Bayes: Parameter estimation failure",
    "Bayes: No rapid grid for parameter estimation",
    "Bayes: Model selection parameter estimation failure",
    "Bayes: Model selection Hessian error",
	"Bayes: w max not found, pdf too sharp, too many counts?",
	"Bayes: Error in Ave & Errs (MP Vals only)",
	"Bayes: Error in Ave & Errs",
	"Bayes: Insufficient grid",
	""
};


const char* bayes_GetBayesErrorDescription(int error)
{
	if (error >= 0)
		return ("");
	else if (error == BAYES__RESULT_USER_CANCEL)  // -99
		return ("User Cancelled");
	else if (error < -14)
		return ("Unknown Error");
	else
		return (bayes_ErrorDescription[-error]);
}


            //this function wraps almost everything
            //bayes operates in a transformed parameter space
            //this routine takes and return estimates in the gci form (Z,A,tau),
            //transforming internally for bayes here
            //outside of here -  is pretty much as it was
            //inside - things are different, i.e 'transformed' then 'untransformed'
            //of course, certain bayes stuff needs to be defined outide too, distns etc...




int bayes_FitPredictedDecay(float           *fitted,
                            int              transbins,
                            int              fitstart,
                            int              fitend,
                            double          *binwalls,
                            BayesInstrRsp_t *instr,
                            float            interval,
                            float            modulationperiod,
                            int              ndecays,
                            double          *taus,
                            double          *weights,
                            int              nphotons)
{
    double timeunit, **likelihoods, norm;
    double dw0, dp;
    int    k, bin, nbins; 

    if ((!fitted) || (transbins <=0) || (fitstart<0) || (fitend>transbins) || (!binwalls))
    {
        return (-1);
    }

	nbins = fitend;

    for (k=1,dw0=1.0; k<=ndecays; k++)
        dw0 -= (double)weights[k];

    if (dw0<0.0)
        dw0=0.0;

    if(dw0>1.0)
        dw0=1.0;

    dp = (double)nphotons;

    likelihoods = Bayes_dmatrix(1,ndecays,0,nbins-1);
    timeunit    = interval/(float)nbins;

    for (k=1; k<=ndecays; k++)
    {
        //bayes_ArrBinLikelihoodsGivenTau(likelihoods[k],NULL,NULL,NULL,nbins,interval,modulationperiod,instr,taus[k]);

        bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(likelihoods[k],
                                                                   nbins,binwalls,NULL,
                                                                   interval,modulationperiod,instr,taus[k],
                                                                   0,NULL,NULL);
    }

    bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(&norm,interval,modulationperiod,binwalls[fitstart],instr,ndecays,weights,taus);

    //value = dw0/interval;

    for (bin=0; bin<fitstart; bin++)
        fitted[bin] = 0.0;

    for (bin=fitstart; bin<fitend; bin++)
    {
        fitted[bin] = (float)(dw0*(binwalls[bin+1]-binwalls[bin])/(interval-binwalls[fitstart]));

        for (k=1; k<=ndecays; k++)
            fitted[bin] += (float)(((float)(weights[k]*likelihoods[k][bin]))/norm);

        fitted[bin] *= (float)(dp);
    }

    for (bin=fitend; bin<nbins; bin++)
        fitted[bin] = 0.0;

    free_Bayes_dmatrix(likelihoods,1,ndecays,0,nbins-1);

    return (0);
}



//aim is to rebin the original transient in order to provide quicker (though possibly less accurate) transient analysis
void bayes_TransientRebinning(float *transin, int binsin, float *transout, int binsout)
{
    int   ratio, i, j, k;
    float temp;

    ratio = binsin/binsout;

    for (i=0, j=0, k=0; i<binsout; i++)
    {
        for (k=0, temp=0.0; k<ratio; k++)
            temp += transin[j+k];

        j           +=k;
        transout[i]  = temp;
    }

//	if (binsin%binsout)
//	{
//		for (j=0; j<ratio; j++)
//	}
}


int bayes_CheckAndTranformTransientDataForBayesFitting(/* Transient data in... */
                                                           float   *trans,
                                                           int      transbins,
                                                           float    transbinwidth,
                                                           int      fitstart,
                                                           int      fitend,
                                                           /* Re-binning settings */
                                                           int      rebinning,
                                                           int      rebinfactor,
                                                           /* Bayesian data out... */
                                                           int    **data,
                                                           int     *nbins,
                                                           double **binwalls,
                                                           int     *nphotons,
                                                           double  *interval)
{
    int    i;
    float *transdata;
    
    /* Check that the supplied data is valid... */
    if ((!trans) || (transbins<=0) || (transbinwidth<=0.0))
        return (BAYES_ERR_INVALID_DATA);

    if ((fitstart<0) || ((fitend-fitstart)>transbins))
        return (BAYES_ERR_INVALID_WINDOW);

    if (rebinning)
    {
        *nbins    = fitend/rebinfactor;
        transdata = Bayes_vector(0,*nbins-1);
        *interval = transbinwidth*(float)*nbins*(float)rebinfactor;
        bayes_TransientRebinning(trans,fitend,transdata,*nbins);
    }
    else
    {
        transdata = trans;
        *nbins    = fitend;
        *interval = transbinwidth*(float)*nbins;
    }

    /* Transform data for Bayesian algorithms... */
    *data = Bayes_ivector(0,*nbins-1);

    for (i=0; i<*nbins; i++)
        (*data)[i] = (int)(transdata[i]);

    for (i=fitstart, *nphotons=0; i<*nbins; i++)
        *nphotons += (*data)[i];

    if (rebinning)
        free_Bayes_vector(transdata,0,*nbins-1);

    if (*nphotons<=0)
    {
        free_Bayes_ivector(*data,0,*nbins-1);
		*data = NULL;
        return (BAYES_ERR_INVALID_DATA);
    }

    /* Binwalls... */
    *binwalls = Bayes_dvector(0,*nbins);
    bayes_PopulateBinWallsVectorUniformIntervals(*binwalls,*nbins,*interval);

    return (BAYES_ERR_NO_ERROR);
}


int bayes_PerformBayesHyperParameterOptimisation( /* Data... */
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
                                                      double                    alphamin,
                                                      double                    alphamax,
                                                      /* Estimates... */
                                                      float                    *alphastar,
                                                      float                    *minuslogprob,
                                                      /* Bayesian analysis settings... */
                                                      int                       rapidanalysis,
                                                      BayesRapidValueStore_t   *rapidstore)
{
    int     ret;

    if (rapidanalysis)
    {
        if (!rapidstore)
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);

        if ((modeltype==FIT_MONOEXP) && (!rapidstore->monoexpvaluestore))
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);
        
        if ((modeltype==FIT_MONOEXP) && (0==rapidstore->validmonoexpgrid))
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);

        if ((modeltype==FIT_BIEXP) && (!rapidstore->biexpvaluestore))
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);
        
        if ((modeltype==FIT_BIEXP) && (0==rapidstore->validbiexpgrid))
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);
    }

    switch (modeltype)
    {
        case FIT_MONOEXP:
        {
            if (rapidanalysis)
            {
                ret = bayes_RapidMonoExpHyperParamOptimization(data,nbins,fitstart,nphotons,binwalls,
                                                               instr,interval,modulationperiod,
                                                               alphastar,(float)alphamin,0.0f,minuslogprob,
                                                               rapidstore->monoexpvaluestore);
            }
            else
            {
                ret = BAYES_ERR_FUNCTIONALITY_NOT_SUPPORTED;
            }

            break;
        }

        case FIT_BIEXP:
        {
            if (rapidanalysis)
            {
                ret = bayes_RapidBiExpHyperParamOptimization(data,nbins,fitstart,nphotons,binwalls,
                                                             instr,interval,modulationperiod,
                                                             alphastar,(float)alphamin,0.0f,minuslogprob,
                                                             rapidstore);            
            }
            else
            {
                ret = BAYES_ERR_FUNCTIONALITY_NOT_SUPPORTED;
            }

            break;
        }

        default:
        {
            return (BAYES_ERR_INVALID_MODEL);
        }
    }

    return (ret);
}


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
                                              BayesRapidValueStore_t   *rapid)
{
    int     ret;
    int     i;
    double  dval;
    double  *w0_ave, *w1_ave, *dw0, *dw1;

    if (rapidanalysis)
    {
        if (!rapid)
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);

        if ((modeltype==FIT_MONOEXP) && (!rapid->monoexpvaluestore))
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);
        
        if ((modeltype==FIT_MONOEXP) && (0==rapid->validmonoexpgrid))
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);

        if ((modeltype==FIT_BIEXP) && (!rapid->biexpvaluestore))
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);
        
        if ((modeltype==FIT_BIEXP) && (0==rapid->validbiexpgrid))
            return (BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION);
    }

    switch (modeltype)
    {
        case FIT_MONOEXP:
        {
            //no parameter fixing in mono-exp
            //really want common interface and fct ptrs instead here

            if ((weights_ave) && (weights_err))
            {
                w0_ave = &(weights_ave[0]);
                dw0    = &(weights_err[0]);
            }
            else
            {
                w0_ave = NULL;
                dw0    = NULL;
            }

            if ((taus_ave) && (taus_err))
            {
                w1_ave = &(taus_ave[1]);
                dw1    = &(taus_err[1]);
            }
            else
            {
                w1_ave = NULL;
                dw1    = NULL;
            }
            
            if (!rapidanalysis)
            {
                ret = bayes_AveragesAndErrorBarsBinLikelihood(data,nbins,fitstart,binwalls,nphotons,
                                                              &(weights_mp[0]),&(taus_mp[1]),
                                                              w0_ave,w1_ave,
                                                              dw0,dw1,
                                                              instr,interval,modulationperiod,(float)alpha,bayes_MonoExpConfigGetDownhillSimplexPrecision(),
                                                              probdistr,minuslogprob);
            }
            else
            {
                ret = bayes_RapidMonoExpAvgAndErrors(data,nbins,fitstart,nphotons,
                                                     &(weights_mp[0]),&(taus_mp[1]),
                                                     w0_ave,w1_ave,
                                                     dw0,dw1,
                                                     instr,interval,modulationperiod,(float)alpha,/*precision*/bayes_MonoExpConfigGetDownhillSimplexPrecision(),nphotons,
                                                     /*quick*/0,probdistr,rapid->monoexpvaluestore,minuslogprob);
            }

            if (weights_ave)
            {
                weights_ave[1] = 1.0-weights_ave[0]; 
            }

            weights_mp[1]  = 1.0-weights_mp[0];

            if ((ret>=BAYES_ERR_NO_ERROR) && ((weights_mp[0]<0.0) || (weights_mp[0]>1.0) || (taus_mp[1]<0.0)))
                return (BAYES_ERR_PARAM_ESTIMATION_FAILURE);

            break;
        }

        case FIT_BIEXP:
        {
            if (!rapidanalysis)
            {
                ret = bayes_MultiExpDetermineMostProbParamValues(data,nbins,binwalls,&nphotons,
                                                                 ndecays,weights_mp,taus_mp,paramfixing,
                                                                 interval,modulationperiod,instr,
                                                                 alpha,
                                                                 NULL,&dval);            
            }
            else
            {
                ret = bayes_RapidBiExpMostProbWeightsAndTaus(data,nbins,fitstart,binwalls,&nphotons,
                                                             ndecays,
                                                             weights_mp,taus_mp,
                                                             weights_ave,taus_ave,
                                                             weights_err,taus_err,
                                                             paramfixing,
                                                             interval,modulationperiod,instr,
                                                             alpha,
                                                             bayes_GetRapidValueStorePtr(),&dval,NULL);
            }

            *minuslogprob = (float)dval;

            if ((ret>=0) && (!rapidanalysis)) //transform bayes model estimates for  use...
            {
                if (weights_err)
                    for (i=0; i<=ndecays; i++) //errors not currently determined...
                        weights_err[i] = -1.0;

                if (taus_err)
                    for (i=1; i<=ndecays; i++)
                        taus_err[i]    = -1.0;
            }

            break;
        }

        default:
        {
            return (BAYES_ERR_INVALID_MODEL);
        }
    }

    return (ret);
}






/* Input vector paramsfree concerns the  UI and therefore corresponds to the conventional model (Z,A1,tau1,A2,tau2,...) */
/* Assumption with weights in Bayes model is that they sum to unity, by design choose the weight with lowest index (i.e. usually background contribution 'w0') is fixed as a consequence of the other values, and is termed here 'model fixed'... */
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
                                                       BayesInstrRsp_t        *instr)
{

    // could also check that if more than one weight is fixed that the chosen values sum to allowed value...
    int i, k, nfixed, nweightsuserfixed, ntaususerfixed, nfree, ndecays, ret, error=0;

	// in case of error before allocation
	paramfixing->weightuserfixed=NULL;
	paramfixing->tauuserfixed=NULL;
	paramfixing->weights=NULL;
	paramfixing->taus=NULL;

    if ((nparams<1) || (!paramsfree) || (!paramfixing))
        return (BAYES_ERR_INVALID_MODEL);

    for (i=0,nfree=0,nfixed=0; i<nparams; i++)
    {
        if (paramsfree[i])
            nfree++;
        else
            nfixed++;
    }

    if (nfree<=0)
        return (BAYES_ERR_ALL_BAYES_PARAM_VALUES_FIXED);

    paramfixing->nparams           = nparams;
    paramfixing->nparamsuserfixed  = nfixed;
    ndecays                        = nparams/2;
    paramfixing->weightuserfixed   = Bayes_ivector(0,ndecays);
    paramfixing->tauuserfixed      = Bayes_ivector(1,ndecays);
    paramfixing->weights           = Bayes_dvector(0,ndecays);
    paramfixing->taus              = Bayes_dvector(1,ndecays);

    if (nfree<nparams)
    {
        nweightsuserfixed = 0;
        ntaususerfixed    = 0;

        /* User has fixed at least one parameter value... */
        /* ...extract fixed values and check for validity */
        paramfixing->weightuserfixed[0] = (paramsfree[0]) ? BAYES_PARAM_VALUE_FREE: BAYES_PARAM_VALUE_USER_FIXED;

        if (paramfixing->weightuserfixed[0] == BAYES_PARAM_VALUE_USER_FIXED)
        {
            nweightsuserfixed++;
            paramfixing->weights[0] = bayes_ToBayesModelTransformFromParamZ(params[0], nbins, nphotons);

            if ((paramfixing->weights[0]<0.0) || (paramfixing->weights[0]>1.0))
            {
				bayes_FreeParameterValueFixingForBayesFitting(paramfixing, ndecays, nbins);
                return (BAYES_ERR_INVALID_FIXED_PARAM_VALUE);
            }
        }

        for (k=1,i=1; k<=ndecays; k++)
        {
            /* Strip out any fixed weights... */
            paramfixing->weightuserfixed[k] = (paramsfree[i]) ? BAYES_PARAM_VALUE_FREE: BAYES_PARAM_VALUE_USER_FIXED;
            
            if (paramfixing->weightuserfixed[k] == BAYES_PARAM_VALUE_USER_FIXED)
            {
                nweightsuserfixed++;
#if 1 /* CURRENT: Decay component weight fixing is not supported, only background level fixing... */
				bayes_FreeParameterValueFixingForBayesFitting(paramfixing, ndecays, nbins);
                return (BAYES_ERR_FUNCTIONALITY_NOT_SUPPORTED);
#else /* FUTURE: Decay component weight fixing implementation... */
                paramfixing->weights[k] = params[i];

                if ((paramfixing->weights[k]<0.0) || (paramfixing->weights[k]>1.0))
                {
                    error = 1;
                    break;
                }
#endif
            }
            else
            {
                paramfixing->weights[k] = -1.0;
            }

            i++;

            /* Strip out any fixed lifetimes... */
            paramfixing->tauuserfixed[k] = (paramsfree[i]) ? BAYES_PARAM_VALUE_FREE: BAYES_PARAM_VALUE_USER_FIXED;

            if (paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_USER_FIXED)
            {
                ntaususerfixed++;
                paramfixing->taus[k] = params[i];

                if (paramfixing->taus[k]<=0.0)
                {
                    error = 1;
                    break;
                }
            }
            else
            {
                paramfixing->taus[k] = -1.0;

            }

            i++;
        }

        if (error)
        {
			bayes_FreeParameterValueFixingForBayesFitting(paramfixing, ndecays, nbins);
            return (BAYES_ERR_INVALID_FIXED_PARAM_VALUE);
        }

        paramfixing->nweightsuserfixed = nweightsuserfixed;
        paramfixing->ntaususerfixed    = ntaususerfixed;

        /* Pre-compute fluorescence likelihoods for use in subsequent algorithms... */
        paramfixing->fluorescencelikelihoods = (BayesPsuedoRapidDiscreteValues_t*)malloc(sizeof(BayesPsuedoRapidDiscreteValues_t)*(1+ndecays));
        
        for (k=1,i=1; k<=ndecays; k++)
        {
            if (paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_USER_FIXED)
            {
                paramfixing->fluorescencelikelihoods[k].fluorescencedecayphotonlikelihoodsgiventau = Bayes_dvector(0,nbins-1);

                ret = bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(paramfixing->fluorescencelikelihoods[k].fluorescencedecayphotonlikelihoodsgiventau,
                                                                                 nbins,binwalls,NULL,
                                                                                 interval,modulationperiod,instr,
                                                                                 paramfixing->taus[k],ndecays,NULL,NULL);

                if (ret>=0)
                {
                    paramfixing->fluorescencelikelihoods[k].valid = 1;
                    paramfixing->fluorescencelikelihoods[k].tau   = paramfixing->taus[k];
                }
                else
                {
					bayes_FreeParameterValueFixingForBayesFitting(paramfixing, ndecays, nbins);
                    free_Bayes_dvector(paramfixing->fluorescencelikelihoods[k].fluorescencedecayphotonlikelihoodsgiventau,0,nbins-1);
                    return (BAYES_ERR_INVALID_MODEL);
                }
            }
            else
            {
                paramfixing->fluorescencelikelihoods[k].valid                                      = 0;
                paramfixing->fluorescencelikelihoods[k].tau                                        = -1.0;
                paramfixing->fluorescencelikelihoods[k].fluorescencedecayphotonlikelihoodsgiventau = NULL;        
            }        
        }
    }
    else
    {
        paramfixing->nweightsuserfixed = 0;
        paramfixing->ntaususerfixed    = 0;

        for (k=0; k<=ndecays; k++)
        {
            paramfixing->weightuserfixed[k] = BAYES_PARAM_VALUE_FREE;
            paramfixing->weights[k]         = -1.0;
        }

        for (k=1; k<=ndecays; k++)
        {
            paramfixing->tauuserfixed[k]    = BAYES_PARAM_VALUE_FREE;
            paramfixing->taus[k]            = -1.0;
        }

        paramfixing->fluorescencelikelihoods = NULL;
    }

    return (BAYES_ERR_NO_ERROR);
}


int bayes_FreeParameterValueFixingForBayesFitting(BayesUserFixedParams_t *paramfixing,
                                                      int                     ndecays,
													  int					  nbins)
{
	int i;
	
	if (paramfixing->weightuserfixed){
		free_Bayes_ivector(paramfixing->weightuserfixed,0,ndecays);
		paramfixing->weightuserfixed=NULL;
	}

	if (paramfixing->tauuserfixed){
		free_Bayes_ivector(paramfixing->tauuserfixed,1,ndecays);
		paramfixing->tauuserfixed=NULL;
	}

	if (paramfixing->weights){
		free_Bayes_dvector(paramfixing->weights,0,ndecays);
		paramfixing->weights=NULL;
	}

	if (paramfixing->taus){
		free_Bayes_dvector(paramfixing->taus,1,ndecays);
		paramfixing->taus=NULL;
	}

	if (paramfixing->fluorescencelikelihoods){
	 for (i=1; i<=ndecays; i++){
		 if (paramfixing->fluorescencelikelihoods[i].fluorescencedecayphotonlikelihoodsgiventau){
			free_Bayes_dvector(paramfixing->fluorescencelikelihoods[i].fluorescencedecayphotonlikelihoodsgiventau,0,nbins-1);
			paramfixing->fluorescencelikelihoods[i].fluorescencedecayphotonlikelihoodsgiventau = NULL;
		 }
	 }

	 free(paramfixing->fluorescencelikelihoods);
	 paramfixing->fluorescencelikelihoods = NULL;
	}

    return (BAYES_ERR_NO_ERROR);
}

int bayes_DoBayesFitting(/* Data in... */
                             float                    *trans,
                             int                       transbins,
                             float                     transbinwidth,
                             int                       fitstart,
                             int                       fitend,
                             int                      *nphotons,
                             /* Model... */
                             int                       modeltype,
                             int                       nparams,
                             /* Instrument... */
                             BayesInstrRsp_t          *instr,
                             float                     modulationperiod,
                             /* Estimates... */
                             int                      *param_free,
                             float                    *param_mp, /* Parameters input using conventional model and values (i.e. Z,A1,tau1,A2,tau2...) */
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
                             BayesRapidValueStore_t   *rapid)
{
    int     ret = BAYES_ERR_NO_ERROR;
    int    *data=NULL, nbins, rebinning, rebinfactor, ndecays, b_nparams, i;
    float   avg, alpha, *b_param_mp=NULL, *b_param_ave=NULL, *b_param_err=NULL;
    double  interval, sum, *binwalls=NULL;
    double *weights_mp=NULL, *taus_mp=NULL, *weights_ave=NULL, *taus_ave=NULL, *weights_err=NULL, *taus_err=NULL;

    float   alphastar, alphamin=0.005f, alphamax=1.5f;

    BayesUserFixedParams_t paramfixing;

    /* Check input transient for validity and convert and compute   */
    /* as required for subsequent calling of Bayesian algorithms... */
    rebinning   = bayes_GetBayesTransientRebinningActiveFlag();
    rebinfactor = bayes_GetBayesTransientRebinningFactor();
    ret         = bayes_CheckAndTranformTransientDataForBayesFitting(trans,transbins,transbinwidth,fitstart,fitend,
                                                                         rebinning,rebinfactor,
                                                                         &data,&nbins,&binwalls,nphotons,&interval);

    if (ret<BAYES_ERR_NO_ERROR) 
	{
		if (binwalls) free_Bayes_dvector(binwalls,0,nbins);
		if (data)     free_Bayes_ivector(data,0,nbins-1);
		return (ret); /* Data related error has occured... */
	}

    /* Data okay, so estimate hyperparameter value from the data... */
    avg   = data_ComputeBinnedDataAverageArrTime(data,nbins,fitstart,*nphotons,(float)interval);
    alpha = 1.0f/avg;

    /* Convert the  parameter estimates to Bayesian model counterparts (used as starting point in search)... */
    b_nparams   = nparams /*1+2*ndecays*/;
    b_param_mp  = Bayes_vector(0,nparams-1);
    b_param_ave = Bayes_vector(0,nparams-1);
    b_param_err = Bayes_vector(0,nparams-1);

    ndecays     = nparams/2;
    weights_mp  = Bayes_dvector(0,ndecays);
    weights_ave = Bayes_dvector(0,ndecays);
    weights_err = Bayes_dvector(0,ndecays);
    taus_mp     = Bayes_dvector(1,ndecays);
    taus_ave    = Bayes_dvector(1,ndecays);
    taus_err    = Bayes_dvector(1,ndecays);

    bayes_ToBayesModelParamValuesFromConventionalModelParamValues(nparams,b_param_mp,param_mp,
                                                                  /**delay*/(float)instr->params[0].delay,nbins,*nphotons,(float)interval);

    ret = bayes_CheckParameterValueFixingForBayesFitting(&paramfixing,nparams,param_free,param_mp,
                                                             nbins,fitstart,binwalls,*nphotons,
                                                             interval,modulationperiod,instr);

    if (ret<BAYES_ERR_NO_ERROR)
	{
		if (binwalls)	free_Bayes_dvector(binwalls,0,nbins);
		if (data)		free_Bayes_ivector(data,0,nbins-1);

    	if (b_param_mp) free_Bayes_vector(b_param_mp,0,nparams-1);
    	if (b_param_ave)free_Bayes_vector(b_param_ave,0,nparams-1);
    	if (b_param_err)free_Bayes_vector(b_param_err,0,nparams-1);
    	if (weights_mp) free_Bayes_dvector(weights_mp,0,ndecays);
    	if (weights_ave)free_Bayes_dvector(weights_ave,0,ndecays);
    	if (weights_err)free_Bayes_dvector(weights_err,0,ndecays);
    	if (taus_mp)    free_Bayes_dvector(taus_mp,1,ndecays);
    	if (taus_ave)   free_Bayes_dvector(taus_ave,1,ndecays);
    	if (taus_err)   free_Bayes_dvector(taus_err,1,ndecays);
		
		bayes_FreeParameterValueFixingForBayesFitting(&paramfixing, ndecays, nbins);

		return (ret); /* Parameter fixing related error has occured... */
	}

    bayes_PopulateWeightsAndTausVectorsFromParamVector(2*ndecays,param_mp,weights_mp,taus_mp);

    for (i=0,sum=0.0; i<=ndecays; i++)
        sum += weights_mp[i];

    for (i=0; i<=ndecays; i++)
        weights_mp[i] /= sum;

    /* No problems yet, perform full hyperparameter optimisation if required... */
    if (bayes_ConfigUseFullBayesianHyperParamDetermination())
    {
        ret = bayes_PerformBayesHyperParameterOptimisation(data,nbins,fitstart,binwalls,*nphotons,
                                                               instr,modulationperiod,(float)interval,
                                                               modeltype,ndecays,&paramfixing,alphamin,alphamax,
                                                               &alphastar,val,
                                                               rapidanalysis,rapid);

        if (ret>=0)
            alpha = alphastar; 
    }

    /* Everything okay so far, now perform Bayesian parameter estimation... */
	ret = bayes_PerformBayesParameterEstimation(data,nbins,fitstart,binwalls,*nphotons,
                                                    instr,modulationperiod,(float)interval,
                                                    modeltype,ndecays,&paramfixing,alpha,
                                                    weights_mp,taus_mp,weights_ave,taus_ave,weights_err,taus_err,
                                                    val,
                                                    distr,rapidanalysis,rapid);

    if ((ret<BAYES_ERR_NO_ERROR) && (ret!=BAYES_AVE_ERR_RAPID_INSUFFICIENT_GRID))
	{
		if (binwalls) free_Bayes_dvector(binwalls,0,nbins);
		if (data)     free_Bayes_ivector(data,0,nbins-1);
    	if (b_param_mp) free_Bayes_vector(b_param_mp,0,nparams-1);
    	if (b_param_ave)free_Bayes_vector(b_param_ave,0,nparams-1);
    	if (b_param_err)free_Bayes_vector(b_param_err,0,nparams-1);
    	if (weights_mp) free_Bayes_dvector(weights_mp,0,ndecays);
    	if (weights_ave)free_Bayes_dvector(weights_ave,0,ndecays);
    	if (weights_err)free_Bayes_dvector(weights_err,0,ndecays);
    	if (taus_mp)    free_Bayes_dvector(taus_mp,1,ndecays);
    	if (taus_ave)   free_Bayes_dvector(taus_ave,1,ndecays);
    	if (taus_err)   free_Bayes_dvector(taus_err,1,ndecays);
		
		bayes_FreeParameterValueFixingForBayesFitting(&paramfixing, ndecays, nbins);
		return (ret); /* Parameter estimation error has occured... */
	}

    if (FIT_MONOEXP==modeltype)
    {
        /* Reporting average values to user; avoid coarse discretisation due to grid values... */
        weights_ave[1] = 1.0-weights_ave[0]; 
        weights_mp[0]  = weights_ave[0];
        weights_mp[1]  = weights_ave[1];
        taus_mp[1]     = taus_ave[1];

        /* Populate  model vector using converted Bayes weights and taus... */
        param_mp[0]  = bayes_FromBayesModelTransformToParamZ((float)weights_mp[0],nbins,*nphotons);
        
        if (rebinning)
            param_mp[0] /= (float)rebinfactor;

        param_mp[1]  = bayes_FromBayesModelWeightAndTauToParamA((float)(1.0-weights_mp[0]), (float)taus_mp[1],transbinwidth, (float)interval,/**delay*/(float)instr->params[0].delay,*nphotons);
        param_mp[2]  = (float)taus_mp[1];
        
        param_err[0] = bayes_FromBayesModelTransformToParamZ((float)weights_err[0],nbins,*nphotons);
        param_err[1] = -1.0f /* Need to determine this uncertainty */;
        param_err[2] = (float)taus_err[1];
    }
    else
    {
        if (ndecays==2)
        {
            if (rapidanalysis)
            {
                param_mp[0] = (float)weights_ave[0];
                param_mp[1] = (float)weights_ave[1];
                param_mp[2] = (float)taus_ave[1];
                param_mp[3] = (float)weights_ave[2];
                param_mp[4] = (float)taus_ave[2];

                param_err[0] = -1.0f;
                param_err[1] = (float)weights_err[1];
                param_err[2] = (float)taus_err[1];
                param_err[3] = (float)weights_err[2];
                param_err[4] = (float)taus_err[2];
            }
            else
            {
                param_mp[0] = (float)weights_mp[0];
                param_mp[1] = (float)weights_mp[1];
                param_mp[2] = (float)taus_mp[1];
                param_mp[3] = (float)weights_mp[2];
                param_mp[4] = (float)taus_mp[2];

                param_err[0] = -1.0f;
                param_err[1] = -1.0f;
                param_err[2] = -1.0f;
                param_err[3] = -1.0f;
                param_err[4] = -1.0f;            
            }

            param_mp[0] = bayes_FromBayesModelTransformToParamZ((float)weights_mp[0],nbins,*nphotons);

            if (rebinning)
                param_mp[0] /= (float)rebinfactor;

            param_mp[1] = bayes_FromBayesModelWeightAndTauToParamA((float)(1.0-weights_mp[0]-weights_mp[2]), (float)taus_mp[1],transbinwidth, (float)interval,/**delay*/(float)instr->params[0].delay,*nphotons);
            param_mp[3] = bayes_FromBayesModelWeightAndTauToParamA((float)(1.0-weights_mp[0]-weights_mp[1]), (float)taus_mp[2],transbinwidth, (float)interval,/**delay*/(float)instr->params[0].delay,*nphotons);
        }
    }

    if ((fitted) && ((ret>=BAYES_ERR_NO_ERROR) || (ret == BAYES_AVE_ERR_RAPID_INSUFFICIENT_GRID)))
    {
        if (rebinning) /* Need to recompute binwalls... */
        {
            interval = transbinwidth*(float)transbins;
            free_Bayes_dvector(binwalls,0,nbins);
            binwalls = Bayes_dvector(0,transbins);
            bayes_PopulateBinWallsVectorUniformIntervals(binwalls,transbins,interval);
        }

        bayes_FitPredictedDecay(fitted,transbins,fitstart,fitend,binwalls,
                                instr, (float)interval,modulationperiod,
                                ndecays,taus_mp,weights_mp,*nphotons);
    }

	if (binwalls) free_Bayes_dvector(binwalls,0,nbins);
	if (data)     free_Bayes_ivector(data,0,nbins-1);

	if (b_param_mp) free_Bayes_vector(b_param_mp,0,nparams-1);
	if (b_param_ave)free_Bayes_vector(b_param_ave,0,nparams-1);
	if (b_param_err)free_Bayes_vector(b_param_err,0,nparams-1);
	if (weights_mp) free_Bayes_dvector(weights_mp,0,ndecays);
	if (weights_ave)free_Bayes_dvector(weights_ave,0,ndecays);
	if (weights_err)free_Bayes_dvector(weights_err,0,ndecays);
	if (taus_mp)    free_Bayes_dvector(taus_mp,1,ndecays);
	if (taus_ave)   free_Bayes_dvector(taus_ave,1,ndecays);
	if (taus_err)   free_Bayes_dvector(taus_err,1,ndecays);
		
	bayes_FreeParameterValueFixingForBayesFitting(&paramfixing, ndecays, nbins);

    return (ret);
}

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
                                     BayesRapidValueStore_t   *rapidgrid)
{
    int     ret = BAYES_ERR_NO_ERROR;
    int    *data, nbins, i;
    int     rebinning, rebinfactor;
    double *binwalls, interval;
    BayesParamValsAndFit_t decayestimates[3];

    double weights[3]; //hack to use A's rather than w's...

    rebinning   = bayes_GetBayesTransientRebinningActiveFlag();
    rebinfactor = bayes_GetBayesTransientRebinningFactor();
    ret         = bayes_CheckAndTranformTransientDataForBayesFitting(trans,transbins,transbinwidth,fitstart,fitend,
                                                                         rebinning,rebinfactor,
                                                                         &data,&nbins,&binwalls,nphotons,&interval);

    if (ret<0) return (ret); /* Data related error... */

    decayestimates[1].weights[0] = bayes_ToBayesModelTransformFromParamZ(param_rld[0],nbins,*nphotons);
    decayestimates[1].weights[1] = 1.0-decayestimates[1].weights[0];
    decayestimates[1].taus[1]    = param_rld[2];
    decayestimates[1].fitted     = NULL;
    decayestimates[1].residuals  = NULL;

    decayestimates[2].weights[0] = bayes_ToBayesModelTransformFromParamZ(param_rld[0],nbins,*nphotons);
    decayestimates[2].weights[1] = 0.5*(1.0-decayestimates[2].weights[0]);
    decayestimates[2].weights[2] = 1.0-decayestimates[2].weights[0]-decayestimates[2].weights[1];
    decayestimates[2].taus[1]    = 1.5*param_rld[2];
    decayestimates[2].taus[2]    = 0.5*param_rld[2];
    decayestimates[2].fitted     = NULL;
    decayestimates[2].residuals  = NULL;

    ret = bayes_DetemineDecayModelRelativeLikelihoods(nbins,fitstart,binwalls,data,*nphotons,
                                                      interval,modulationperiod,instr,
                                                      decayestimates,
                                                      decaymodellikelihoods,paramvalsandfits,
                                                      rapidanalysis,rapidgrid);

    if (ret>=BAYES_ERR_NO_ERROR)
    {
        if (paramvalsandfits)
        {
            if (rebinning) /* Need to recompute binwalls... */
            {
                interval = transbinwidth*(float)transbins;
                free_Bayes_dvector(binwalls,0,nbins);
                binwalls = Bayes_dvector(0,transbins);
                bayes_PopulateBinWallsVectorUniformIntervals(binwalls,transbins,interval);
            }

            for (i=1; i<=2; i++) //need to make generic
            {
                if (paramvalsandfits[i].fitted)
                {
                    bayes_FitPredictedDecay(paramvalsandfits[i].fitted,transbins,fitstart,fitend,binwalls,
                                            instr, (float)interval,modulationperiod,
                                            i,paramvalsandfits[i].taus,paramvalsandfits[i].weights,*nphotons);                
                }

                if (i==1)
                {
                    weights[0] = paramvalsandfits[1].weights[0];
                    weights[1] = paramvalsandfits[1].weights[1];

                    paramvalsandfits[1].weights[0] = bayes_FromBayesModelTransformToParamZ((float)weights[0],nbins,*nphotons);

                    if (rebinning)
                        paramvalsandfits[1].weights[0] /= (float)rebinfactor;

                    paramvalsandfits[1].weights[1] = bayes_FromBayesModelWeightAndTauToParamA((float)(1.0-weights[0]), (float)paramvalsandfits[1].taus[1],
						(float)interval/ (float)nbins, (float)interval,/**delay*/(float)instr->params[0].delay,*nphotons);
                }
                else
                {
                    weights[0] = paramvalsandfits[2].weights[0];
                    weights[1] = paramvalsandfits[2].weights[1];
                    weights[2] = paramvalsandfits[2].weights[2];

                    paramvalsandfits[2].weights[0] = bayes_FromBayesModelTransformToParamZ((float)weights[0],nbins,*nphotons);

                    if (rebinning)
                        paramvalsandfits[2].weights[0] /= (float)rebinfactor;
                    
                    paramvalsandfits[2].weights[1] = bayes_FromBayesModelWeightAndTauToParamA((float)(1.0-weights[0]-weights[2]), (float)paramvalsandfits[2].taus[1],
						(float)interval/(float)nbins, (float)interval,/**delay*/(float)instr->params[0].delay,*nphotons);
                    paramvalsandfits[2].weights[2] = bayes_FromBayesModelWeightAndTauToParamA((float)(1.0-weights[0]-weights[1]), (float)paramvalsandfits[2].taus[2],
						(float)interval/(float)nbins, (float)interval,/**delay*/(float)instr->params[0].delay,*nphotons);
                }
            }
        }    
    }

    free_Bayes_ivector(data,0,nbins-1);
    free_Bayes_dvector(binwalls,0,nbins);

    return (ret);
}

/*=================================================================================*/
/*                                                                                 */
/*                    BAYESIAN (MINUS LOG) PROBABILITY CALCULATOR                  */
/*                                                                                 */
/*=================================================================================*/

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
                                                         int    fitend)
{
    int     ndecays, i, nphotons;
    double *weights_bayes, *taus_bayes, *weights_conventional, *taus_conventional;
    int    *data, nbins;
    float   interval;
    int     rebinning, rebinfactor;
    float  *transdata;

    if ((!trans) || (transbins<=0) || (transbinwidth<=0.0))
        return (BAYES_ERR_INVALID_DATA);

    if ((fitstart<0) || ((fitend-fitstart)>transbins))
        return (BAYES_ERR_INVALID_WINDOW);

    /* Transient data and settings valid, do rebinning, windowing and determine window characteristics... */
    rebinning = bayes_GetBayesTransientRebinningActiveFlag();
	
    if (rebinning)
    {
        rebinfactor = bayes_GetBayesTransientRebinningFactor();
        nbins       = (fitend-fitstart)/rebinfactor;
        transdata   = Bayes_vector(0,nbins-1);
        bayes_TransientRebinning(trans,fitend-fitstart,transdata,nbins);

        interval    = transbinwidth*(float)nbins*(float)rebinfactor;
    }
    else
    {
        nbins     = fitend-fitstart;
        transdata = trans;

        interval  = transbinwidth*(float)nbins;
    }

    data  = Bayes_ivector(0,nbins-1);

    for (i=0; i<nbins; i++)
        data[i] = (int)(transdata[fitstart+i]);

    for (i=0, nphotons=0; i<nbins; i++)
        nphotons += data[i];

    if (nphotons<=0)
    {
        free_Bayes_ivector(data,0,nbins-1);

        if (rebinning)
            free_Bayes_vector(transdata,0,nbins-1);

        return (BAYES_ERR_INVALID_DATA);
    }

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_conventional,&taus_conventional);
    bayes_PopulateWeightsAndTausVectorsFromParamVector(nparams,params_conventional,weights_conventional,taus_conventional);

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_bayes,&taus_bayes);

    for (i=1; i<=ndecays; i++)
    {
//        weights_bayes[i] = (double)bayes_ToBayesModelWeightFromParamAAndTau(weights_conventional[i],taus_conventional[i],
  //                                                                          interval/(float)nbins,interval,delay,nphotons);
        weights_bayes[i] = (double)bayes_ToBayesModelWeightFromParamAAndTau((float)weights_conventional[i], (float)taus_conventional[i],
            interval/(float)nbins, (float)interval, (float)instr->params[0].delay,nphotons);

        taus_bayes[i]    = taus_conventional[i];
    }

    for (i=1,weights_bayes[0]=1.0; i<=ndecays; i++)
        weights_bayes[0] -= weights_bayes[i];

    bayes_PopulateParamVectorFromWeightsAndTausVectors(params_bayes,nparams,weights_bayes,taus_bayes);

//    bayes_FreeWeightsAndTausVectors(ndecays,weights_conventional,

    return (0);
}


/* Calcualate the equivalent raw chisq for Bayes fit. Fills residuals array as a side effect. */
float bayes_CalculateResidualsAndEquivalentChisq(float y[], float fitted[], float *residuals, int fit_start, int fit_end)
{
	float chisq_local = 0.0f, res, sigma2;
	int i;

	// ignore the first part of the fit for fair comparison of chisq with with LMA, after all Bayes does not use the chisq to fit!
	for (i = 0; i < fit_start; i++) {
		res = y[i] - fitted[i];
		if (residuals != NULL)
			residuals[i] = res;
	}

	// Chisq for noise model POISSON_FIT
	for (; i < fit_end; i++) {
		res = y[i] - fitted[i];
		if (residuals != NULL)
			residuals[i] = res;
		// don't let variance drop below 1 
		sigma2 = (fitted[i] > 1 ? 1.0f / fitted[i] : 1.0f);
		chisq_local += res * res * sigma2;
	}

	return (chisq_local);
}



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
                                                         int    fitend)
{
    int     ndecays, i, nphotons;
    double *weights_bayes, *taus_bayes, *weights_conventional, *taus_conventional;
    int    *data, nbins;
    float   interval;
    int     rebinning, rebinfactor;
    float  *transdata;

    if ((!trans) || (transbins<=0) || (transbinwidth<=0.0))
        return (BAYES_ERR_INVALID_DATA);

    if ((fitstart<0) || ((fitend-fitstart)>transbins))
        return (BAYES_ERR_INVALID_WINDOW);

    /* Transient data and settings valid, do rebinning, windowing and determine window characteristics... */
    rebinning = bayes_GetBayesTransientRebinningActiveFlag();
	
    if (rebinning)
    {
        rebinfactor = bayes_GetBayesTransientRebinningFactor();
        nbins       = (fitend-fitstart)/rebinfactor;
        transdata   = Bayes_vector(0,nbins-1);
        bayes_TransientRebinning(trans,fitend-fitstart,transdata,nbins);

        interval    = transbinwidth*(float)nbins*(float)rebinfactor;
    }
    else
    {
        nbins     = fitend-fitstart;
        transdata = trans;

        interval  = transbinwidth*(float)nbins;
    }

    data  = Bayes_ivector(0,nbins-1);

    for (i=0; i<nbins; i++)
        data[i] = (int)(transdata[fitstart+i]);

    for (i=0, nphotons=0; i<nbins; i++)
        nphotons += data[i];

    if (nphotons<=0)
    {
        free_Bayes_ivector(data,0,nbins-1);

        if (rebinning)
            free_Bayes_vector(transdata,0,nbins-1);

        return (BAYES_ERR_INVALID_DATA);
    }

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_conventional,&taus_conventional);

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_bayes,&taus_bayes);
    bayes_PopulateWeightsAndTausVectorsFromParamVector(nparams,params_bayes,weights_bayes,taus_bayes);

    for (i=1; i<=ndecays; i++)
    {
//        weights_conventional[i] = (double)bayes_FromBayesModelWeightAndTauToParamA(weights_bayes[i],taus_bayes[i],
  //                                                                          interval/(float)nbins,interval,delay,nphotons);
        weights_conventional[i] = (double)bayes_FromBayesModelWeightAndTauToParamA((float)weights_bayes[i], (float)taus_bayes[i],
            interval/(float)nbins,interval, (float)instr->params[0].delay,nphotons);

        taus_conventional[i]    = taus_bayes[i];        
    }

    for (i=1,weights_bayes[0]=1.0; i<=ndecays; i++)
        weights_bayes[0] -= weights_bayes[i];

    bayes_PopulateParamVectorFromWeightsAndTausVectors(params_conventional,nparams,weights_conventional,taus_conventional);

//    bayes_FreeWeightsAndTausVectors(ndecays,weights_conventional,

    return (0);
}

int bayes_ComputeBayesHyperParamsFromData(/* Data in... */
                                              float *trans,
                                              int    transbins,
                                              float  transbinwidth,
                                              int    fitstart,
                                              int    fitend,
                                              /* Hyperparameter values out... */
                                              int    nhyperparams,
                                              float *hyperparams)
{
    int   ret = BAYES_ERR_NO_ERROR;

    int   *data, rebinning, rebinfactor, nbins, nphotons, i;
    float *transdata, interval, avg, alpha;

    if ((!hyperparams) || (nhyperparams<=0))
        return (BAYES_ERR_INVALID_MODEL);
    
    if ((!trans) || (transbins<=0) || (transbinwidth<=0.0))
        return (BAYES_ERR_INVALID_DATA);

    if ((fitstart<0) || ((fitend-fitstart)>transbins))
        return (BAYES_ERR_INVALID_WINDOW);

    /* Transient data and settings valid, do rebinning, windowing and determine window characteristics... */
    rebinning = bayes_GetBayesTransientRebinningActiveFlag();
	
    if (rebinning)
    {
        rebinfactor = bayes_GetBayesTransientRebinningFactor();
        nbins       = (fitend-fitstart)/rebinfactor;
        transdata   = Bayes_vector(0,nbins-1);
        bayes_TransientRebinning(trans,fitend-fitstart,transdata,nbins);

        interval    = transbinwidth*(float)nbins*(float)rebinfactor;
    }
    else
    {
        nbins     = fitend-fitstart;
        transdata = trans;

        interval  = transbinwidth*(float)nbins;
    }

    data  = Bayes_ivector(0,nbins-1);

    for (i=0; i<nbins; i++)
        data[i] = (int)(transdata[fitstart+i]);

    for (i=0, nphotons=0; i<nbins; i++)
        nphotons += data[i];

    if (nphotons<=0)
    {
        free_Bayes_ivector(data,0,nbins-1);

        if (rebinning)
            free_Bayes_vector(transdata,0,nbins-1);

        return (BAYES_ERR_INVALID_DATA);
    }

   	avg   = data_ComputeBinnedDataAverageArrTime(data,nbins,fitstart,nphotons,interval);
    alpha = 1.0f/avg;

    for (i=0; i<=nhyperparams; i++)
        hyperparams[i] = alpha;

    free_Bayes_ivector(data,0,nbins-1);

    if (rebinning)
        free_Bayes_vector(transdata,0,nbins-1);

    return (ret);
}


/*=================================================================================*/
/*                                                                                 */
/*                 ACCESS FUNCTIONS FOR VARIABLE VALUES STORED BY UI               */
/*                                                                                 */
/*=================================================================================*/

float bayes_GetConfigParameterValueModulationPeriod(void)
{
	return (bayesrepperiod);
}
void bayes_SetConfigParameterValueModulationPeriod(float val)
{
	bayesrepperiod = val;
}

int bayes_UseRepetitionEffectsInAnalysis(void)
{
	return (bayes_GetIncludeRepetitionEffectsInAnalysisFlag());
}
int bayes_GetIncludeRepetitionEffectsInAnalysisFlag(void)
{
	return (bayesrep);
}

int bayes_GetBayesTransientRebinningActiveFlag(void)
{
    return (bayesadvrebinning);
}

int bayes_GetBayesTransientRebinningFactor(void)
{
	return (bayesadvrebinfactor);
}

float bayes_BiExpConfigGetSimAnnealTempReductionSchedule(void)
{
	return (bayesadvtempred);
}

int bayes_BiExpConfigGetSimAnnealTempReductionNoOfSteps(void)
{
	return (bayesadvtempsteps);
}

int bayes_BiExpConfigGetSimAnnealItersAtEachTemp(void)
{
    return (bayesadvtempiters);
}

float bayes_BiExpConfigGetSimAnnealFracConvergenceTolerance(void)
{
    return (bayesadvtolerance);
}

int bayes_BiExpConfigGetNumberOfRestarts(void)
{
    return (bayesadvrestarts);
}

float bayes_MonoExpConfigGetDownhillSimplexPrecision(void)
{
	return (bayesadvmonoprec);
}

int bayes_ConfigGetMinimizationAlgorithm(void)
{
    return (bayesadvalgsel);
}

/*===========================================================================*/
/*                                                                           */
/*             Multi-exp rapid search settings access functions              */
/*                                                                           */
/*===========================================================================*/

int bayes_BiExpConfigGetRapidGridSearchInitialisationsFactor(void)
{
	return (bayesadvbistarts);
}

int bayes_BiExpConfigRapidGetGridSearchLocalisedExhaustiveSearchDelta(void)
{
    return (bayesadvbidelta);
}


/*===========================================================================*/
/*                                                                           */
/*          Instrument response determination setting access functions       */
/*                                                                           */
/*===========================================================================*/
float bayes_InstrConfigGetSimAnnealStartingTemp(void)
{
    return (bayesadvirfsatempstart);
}

float bayes_InstrConfigGetSimAnnealTempReductionSchedule(void)
{
    return (bayesadvirfsatempred);
}

int bayes_InstrConfigGetSimAnnealTempReductionNoOfSteps(void)
{
    return (bayesadvirfsatempsteps);
}

int bayes_InstrConfigGetSimAnnealItersAtEachTemp(void)
{
    return (bayesadvirfsatempiters);
}

float bayes_InstrConfigGetSimAnnealFracConvergenceTolerance(void)
{
    return (bayesadvirfsatolerance);
}

int bayes_InstrConfigGetNumberOfRestarts(void)
{
    return (bayesadvrestarts);
}

float bayes_InstrConfigGetDownhillSimplexPrecision(void)
{
    return (bayesadvirfsaprec);
}

int bayes_InstrConfigGetMinimizationAlgorithm(void)
{
    return (bayesadvirfalgsel);
}


/*===========================================================================*/
/*                 RAPID BAYES (DISCRETE GRID IMPLEMENTATION)                */
/*===========================================================================*/

int bayes_RapidGetUseRapidBayesFlag(void)
{
    return (bayesrapiduse);
}


int bayes_FitTypeToRapidGridUpdateType(int fittype)
{
    int gridupdatetype;

    if (fittype==FIT_MONOEXP)
    {
        gridupdatetype = BAYES_RAPID_GRID_MONO;
    }
    else if (fittype==FIT_BIEXP)
    {
        gridupdatetype = BAYES_RAPID_GRID_BI;
    }
    else if (fittype==FIT_MODELSELECTION)
    {
        gridupdatetype = BAYES_RAPID_GRID_MONO_AND_BI;
    }
    else
    {
        gridupdatetype = BAYES_RAPID_GRID_MONO_AND_BI;
    }

    return (gridupdatetype);
}

int bayes_ConfigureBayesianRapidGrid(int updatetype, float xincr, int fitend, BayesIrEstConfig_t *BayesIrEstConfig)
{
    int     i, nbins, fitstart, rebinfactor, update, ret=0;
    int     ntaus_mono=100, nweights_mono=200, ntaus_bi=0, nweights_bi=0;
    float   taulow_mono= 0.1f, tauhigh_mono=4.0f, weightlow_mono=0.0f, weighthigh_mono=1.0f, bglow_mono=0.0;
    float   taulow_bi=0.0, tauhigh_bi=0.0, weightlow_bi=0.0, weighthigh_bi=0.0, bghigh_mono=0.3f, bglow_bi=0.0, bghigh_bi=0.0, val;
    float   interval, modulationperiod;
    double *weights_mono=NULL, *taus_ono=NULL, *weights_bi=NULL, *taus_bi=NULL, *binwalls=NULL, dw, dt;
    int     low[5], high[5];

    BayesInstrRsp_t instr;

    if (bayes_RapidGetUseRapidBayesFlag())
    {
        //* Mono-exponential grid settings... 
        if ((updatetype==BAYES_RAPID_GRID_MONO) ||
            (updatetype==BAYES_RAPID_GRID_MONO_AND_BI))
        {
            BayesMonoRapidGridConfig_t *BayesMonoRapidGridConfig = bayes_GetMonoRapidGridConfigPtrSafe();
			ntaus_mono = BayesMonoRapidGridConfig->bayesrapidtaupts;
			taulow_mono = BayesMonoRapidGridConfig->bayesrapidtaulow;
			tauhigh_mono = BayesMonoRapidGridConfig->bayesrapidtauhigh;
			nweights_mono = BayesMonoRapidGridConfig->bayesrapidwpts;
			weightlow_mono = BayesMonoRapidGridConfig->bayesrapidwlow;
			weighthigh_mono = BayesMonoRapidGridConfig->bayesrapidwhigh;
			bghigh_mono = BayesMonoRapidGridConfig->bayesrapidbghigh;

            weights_mono   = Bayes_dvector(0,nweights_mono-1);
            taus_ono  = Bayes_dvector(0,ntaus_mono-1);

            for (i=0,dw=(weighthigh_mono-weightlow_mono)/(double)(nweights_mono-1); i<nweights_mono; i++)
                weights_mono[i] = weightlow_mono + dw*(double)i;

            for (i=0,dt=(tauhigh_mono-taulow_mono)/(double)(ntaus_mono-1); i<ntaus_mono; i++)
                taus_ono[i] = taulow_mono + dt*(double)i;
        }

        //* Bi-exponential grid settings... 
        if ((updatetype==BAYES_RAPID_GRID_BI) ||
            (updatetype==BAYES_RAPID_GRID_MONO_AND_BI))
        {
            BayesBiRapidGridConfig_t *BayesBiRapidGridConfig = bayes_GetBiRapidGridConfigPtrSafe();
			ntaus_bi = BayesBiRapidGridConfig->bayesrapidbitaupts;
			taulow_bi = BayesBiRapidGridConfig->bayesrapidbitaulow;
			tauhigh_bi = BayesBiRapidGridConfig->bayesrapidbitauhigh;
			nweights_bi = BayesBiRapidGridConfig->bayesrapidbiweightpts;
			weightlow_bi = BayesBiRapidGridConfig->bayesrapidbiweightlow;
			weighthigh_bi = BayesBiRapidGridConfig->bayesrapidbiweighthigh;
			bglow_bi = BayesBiRapidGridConfig->bayesrapidbibgmin;
			bghigh_bi = BayesBiRapidGridConfig->bayesrapidbibgmax;

            weights_bi = Bayes_dvector(0,nweights_bi-1);
            taus_bi    = Bayes_dvector(0,ntaus_bi-1);

            for (i=0,dw=(weighthigh_bi-weightlow_bi)/(double)(nweights_bi-1); i<nweights_bi; i++)
                weights_bi[i] = weightlow_bi + dw*(double)i;

            for (i=0,dt=(tauhigh_bi-taulow_bi)/(double)(ntaus_bi-1); i<ntaus_bi; i++)
                taus_bi[i]    = taulow_bi + dt*(double)i;

            //GetCtrlVal(BayesConfigPanel,BAYES_RAPID_BI_W0_PRE_L,&val);
			val = BayesBiRapidGridConfig->bayesrapidbiw0low;
            if (val<bglow_bi)  val = bglow_bi;
            if (val>bghigh_bi) val = bghigh_bi;
            low[0] = bayes_MapWeightValueToClosestRapidGridPoint(val,nweights_bi,weights_bi);

            low[1] = bayes_MapWeightValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbiw1low,nweights_bi,weights_bi);
			low[2] = bayes_MapWeightValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbiw2low,nweights_bi,weights_bi);
			low[3] = bayes_MapLifetimeValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbitau1low,ntaus_bi,taus_bi);
			low[4] = bayes_MapLifetimeValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbitau2low,ntaus_bi,taus_bi);

            //GetCtrlVal(BayesConfigPanel,BAYES_RAPID_BI_W0_PRE_H,&val);
			val = BayesBiRapidGridConfig->bayesrapidbiw0high;
            if (val<bglow_bi)  val = bglow_bi;
            if (val>bghigh_bi) val = bghigh_bi;
            high[0] = bayes_MapWeightValueToClosestRapidGridPoint(val,nweights_bi,weights_bi);
			high[1] = bayes_MapWeightValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbiw1high,nweights_bi,weights_bi);
			high[2] = bayes_MapWeightValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbiw2high,nweights_bi,weights_bi);
			high[3] = bayes_MapLifetimeValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbitau1high,ntaus_bi,taus_bi);
			high[4] = bayes_MapLifetimeValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbitau2high,ntaus_bi,taus_bi);

            high[0] = (high[0]>=nweights_bi) ? (nweights_bi) : (high[0]);
            high[1] = (high[1]>=nweights_bi) ? (nweights_bi) : (high[1]);
            high[2] = (high[2]>=nweights_bi) ? (nweights_bi) : (high[2]);
            high[3] = (high[3]>=ntaus_bi) ? (ntaus_bi) : (high[3]);
            high[4] = (high[4]>=ntaus_bi) ? (ntaus_bi) : (high[4]);
        }

        //* Common settings... 
        bayes_GetInstrRspParamValues(&instr, BayesIrEstConfig);
        modulationperiod = bayes_GetConfigParameterValueModulationPeriod();

        fitstart = 0;
        nbins    = fitend-fitstart;

        if (bayes_GetBayesTransientRebinningActiveFlag())
            rebinfactor = bayes_GetBayesTransientRebinningFactor();
        else
            rebinfactor = 1;

        nbins    = (fitend-fitstart)/rebinfactor;
        interval = xincr *(float)nbins*(float)rebinfactor;

        update   = bayes_DetermineIfBayesGridUpdateReqd(bayes_GetRapidValueStorePtr(),updatetype,
                                                        //* Mono-exp... 
                                                        ntaus_mono,taus_ono,
                                                        nweights_mono,weights_mono,
                                                        bglow_mono,bghigh_mono,
                                                        //* Bi-exp... 
                                                        ntaus_bi,taus_bi,
                                                        nweights_bi,weights_bi,bglow_bi,bghigh_bi,
                                                        low,high,
                                                        //* Instrumentation and data... 
                                                        nbins,&instr,interval,modulationperiod);

        if (update)
        {
            bayes_DestroyRapidValueStore(bayes_GetRapidValueStorePtr(),update);

            binwalls = Bayes_dvector(0,nbins);
            bayes_PopulateBinWallsVectorUniformIntervals(binwalls,nbins,interval);

            ret = bayes_CreateRapidValueStore(bayes_GetRapidValueStorePtr(),update,
                                              ntaus_mono,taus_ono,
                                              nweights_mono,weights_mono,
                                              bglow_mono,bghigh_mono,
                                              ntaus_bi,taus_bi,
                                              nweights_bi,weights_bi,bglow_bi,bghigh_bi,low,high,
                                              nbins,binwalls,&instr,interval,modulationperiod);
        }
    }
    
    if (weights_mono)        free_Bayes_dvector(weights_mono,0,nweights_mono-1);
    if (taus_ono)       free_Bayes_dvector(taus_ono,0,ntaus_mono-1);
    if (weights_bi) free_Bayes_dvector(weights_bi,0,nweights_bi-1);
    if (taus_bi)    free_Bayes_dvector(taus_bi,0,ntaus_bi-1);
    if (binwalls)   free_Bayes_dvector(binwalls,0,nbins);

    return (ret);
}

/*
int bayes_RapidBiExpDetermineGridSize(int fitend)
{
    int     i, ntaus, nweights, nbins, fitstart, rebinfactor, npts, nptsvalid;
    float   taulow, tauhigh, weightlow, weighthigh, bglow, bghigh, val;
    double *weights, *taus, dw, dt, memoryreqd;
    int     low[5], high[5];

    //* Bi-exponential grid settings... 
	ntaus = BayesBiRapidGridConfig->bayesrapidbitaupts;
	taulow = BayesBiRapidGridConfig->bayesrapidbitaulow;
	tauhigh = BayesBiRapidGridConfig->bayesrapidbitauhigh;
	nweights = BayesBiRapidGridConfig->bayesrapidbiweightpts;
	weightlow = BayesBiRapidGridConfig->bayesrapidbiweightlow;
	weighthigh = BayesBiRapidGridConfig->bayesrapidbiweighthigh;

    weights = Bayes_dvector(0,nweights-1);
    taus    = Bayes_dvector(0,ntaus-1);

    for (i=0,dw=(weighthigh-weightlow)/(double)(nweights-1); i<nweights; i++)
        weights[i] = weightlow + dw*(double)i;

    for (i=0,dt=(tauhigh-taulow)/(double)(ntaus-1); i<ntaus; i++)
        taus[i]    = taulow + dt*(double)i;

    //* 'Snap' the user defined background minimum and maximum to values on the weight vector... 
	bglow = BayesBiRapidGridConfig->bayesrapidbibgmin;
	if (bglow<weightlow)  bglow = weightlow;
    if (bglow>weighthigh) bglow = weighthigh;
    bglow = (float)weights[bayes_MapWeightValueToClosestRapidGridPoint(bglow,nweights,weights)];
	BayesBiRapidGridConfig->bayesrapidbibgmin = bglow;

	bghigh = BayesBiRapidGridConfig->bayesrapidbibgmax;
	if (bghigh<weightlow)  bghigh = weightlow;
    if (bghigh>weighthigh) bghigh = weighthigh;
    bghigh = (float)weights[bayes_MapWeightValueToClosestRapidGridPoint(bghigh,nweights,weights)];
	BayesBiRapidGridConfig->bayesrapidbibgmax = bghigh;

    //* Pre-computed bi-exponential log likelihoods... 
	val = BayesBiRapidGridConfig->bayesrapidbiw0low;
    if (val<bglow)  val = bglow;
    if (val>bghigh) val = bghigh;
    low[0] = bayes_MapWeightValueToClosestRapidGridPoint(val,nweights,weights);

	low[1] = bayes_MapWeightValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbiw1low,nweights,weights);
	low[2] = bayes_MapWeightValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbiw2low,nweights,weights);
	low[3] = bayes_MapLifetimeValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbitau1low,ntaus,taus);
	low[4] = bayes_MapLifetimeValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbitau2low,ntaus,taus);

    //GetCtrlVal(BayesConfigPanel,BAYES_RAPID_BI_W0_PRE_H,&val);
	val = BayesBiRapidGridConfig->bayesrapidbiw0high;
	if (val<bglow)  val = bglow;
    if (val>bghigh) val = bghigh;
    high[0] = bayes_MapWeightValueToClosestRapidGridPoint(val,nweights,weights);
	high[1] = bayes_MapWeightValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbiw1high,nweights,weights);
	high[2] = bayes_MapWeightValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbiw2high,nweights,weights);
	high[3] = bayes_MapLifetimeValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbitau1high,ntaus,taus);
	high[4] = bayes_MapLifetimeValueToClosestRapidGridPoint(BayesBiRapidGridConfig->bayesrapidbitau2high,ntaus,taus);

    high[1] = (high[1]>=nweights) ? (nweights) : (high[1]);
    high[0] = (high[0]>=nweights) ? (nweights) : (high[0]);
    high[2] = (high[2]>=nweights) ? (nweights) : (high[2]);
    high[3] = (high[3]>=ntaus) ? (ntaus) : (high[3]);
    high[4] = (high[4]>=ntaus) ? (ntaus) : (high[4]);

    //* Data and re-binning... 
    fitstart = 0;
    nbins    = fitend-fitstart;

    if (bayes_GetBayesTransientRebinningActiveFlag())
        rebinfactor = bayes_GetBayesTransientRebinningFactor();
    else
        rebinfactor = 1;

    nbins = (fitend-fitstart)/rebinfactor;

    //* Determine pre-computed grid size... 
    bayes_RapidBiExpDetermineGridSizeAdv(ntaus,taus,nweights,weights,bglow,bghigh,low,high,nbins,
                                      &npts,&nptsvalid,&memoryreqd);

    //* Update, snap the pre-computed user selected values to grid co-ordinates
	BayesBiRapidGridConfig->bayesrapidbiw0low = (float)weights[low[0]];
	BayesBiRapidGridConfig->bayesrapidbiw1low = (float)weights[low[1]];
	BayesBiRapidGridConfig->bayesrapidbiw2low = (float)weights[low[2]];
	BayesBiRapidGridConfig->bayesrapidbitau1low = (float)taus[low[3]];
	BayesBiRapidGridConfig->bayesrapidbitau2low = (float)taus[low[4]];
	BayesBiRapidGridConfig->bayesrapidbiw0high = (float)weights[high[0]];
	BayesBiRapidGridConfig->bayesrapidbiw1high = (float)weights[high[1]];
	BayesBiRapidGridConfig->bayesrapidbiw2high = (float)weights[high[2]];
	BayesBiRapidGridConfig->bayesrapidbitau1high = (float)taus[high[3]];
	BayesBiRapidGridConfig->bayesrapidbitau2high = (float)taus[high[4]];

    free_Bayes_dvector(weights,0,nweights-1);
    free_Bayes_dvector(taus,0,ntaus-1);

    return (0);
}
*/

void bayes_InvalidateBayesDiscreteGridsDueToParameterValueChange(void)
{
    bayes_InvalidateRapidValueStore(bayes_GetRapidValueStorePtr());
}


void bayes_RapidBiExpInvalidateGridDueToParameterValueChange(void)
{
    bayes_InvalidateRapidBiExpValueStore(bayes_GetRapidValueStorePtr());
}


int bayes_CheckForValidBayesDiscreteGrids(int fittype)
{
    return (bayes_CheckForValidRapidValueStore(bayes_GetRapidValueStorePtr(),
                                               bayes_FitTypeToRapidGridUpdateType(fittype)));
}


int bayes_InitDiscreteGridForRapidBayes(void)
{
    return (bayes_InitializeRapidValueStore(bayes_GetRapidValueStorePtr()));
}

int bayes_DestroyDiscreteGridForRapidBayes(void)
{
    /* Destroy all grids as stopping time resolved analysis... */
    return (bayes_DestroyRapidValueStore(bayes_GetRapidValueStorePtr(),BAYES_RAPID_GRID_MONO_AND_BI));
}

void bayes_GetInstrRspParamValues(BayesInstrRsp_t *instr, BayesIrEstConfig_t *BayesIrEstConfig)
{
    BayesInstrRsp_t temp;
    int             i, j, ninstr;
    double          norm;

    //GetCtrlVal(BayesConfigPanel,BAYES_IR_COMPONENTS_MAX,&temp.ninstr);
	temp.ninstr = BayesIrEstConfig->bayesirnumcomponents;

	temp.params[0].weight = BayesIrEstConfig->bayesirweight1;
	temp.params[0].cutoff = BayesIrEstConfig->bayesircutoff1;
	temp.params[0].width = BayesIrEstConfig->bayesirsigma1;
	temp.params[0].delay = BayesIrEstConfig->bayesiruc1;

	temp.params[1].weight = BayesIrEstConfig->bayesirweight2;
	temp.params[1].cutoff = BayesIrEstConfig->bayesircutoff2;
	temp.params[1].width = BayesIrEstConfig->bayesirsigma2;
	temp.params[1].delay = BayesIrEstConfig->bayesiruc2;

	temp.params[2].weight = BayesIrEstConfig->bayesirweight3;
	temp.params[2].cutoff = BayesIrEstConfig->bayesircutoff3;
	temp.params[2].width = BayesIrEstConfig->bayesirsigma3;
	temp.params[2].delay = BayesIrEstConfig->bayesiruc3;

    //* Check for valid values... 
    for (i=0,ninstr=0,norm=0.0; i<temp.ninstr; i++)
    {
        if ((temp.params[i].weight>0.0) && (temp.params[i].width>0.0))
        {
            ninstr++;
            norm += temp.params[i].weight;
        }
        else
        {
            temp.params[i].weight = 0.0;
            temp.params[i].delay  = 0.0;
            temp.params[i].cutoff = 0.0;
            temp.params[i].width  = 0.0;
        }
    }

    for (; i<BAYES_INSTR_RSP_MAX_COMPONENTS; i++)
    {
        temp.params[i].weight = 0.0;
        temp.params[i].delay  = 0.0;
        temp.params[i].cutoff = 0.0;
        temp.params[i].width  = 0.0;
    }

    //* Only copy the valid components into the container, normalise weights... 
    instr->ninstr = ninstr;

    for (i=0, j=0; i<ninstr; i++)
    {
        if (temp.params[i].weight>0.0)
        {
            instr->params[j].weight = temp.params[i].weight/norm;
            instr->params[j].delay  = temp.params[i].delay;
            instr->params[j].cutoff = temp.params[i].cutoff;
            instr->params[j].width  = temp.params[i].width;
            j++;
        }
    }

    //* Neutralize any empty component spaces... 
    for (; j<BAYES_INSTR_RSP_MAX_COMPONENTS; j++)
    {
        instr->params[j].weight = 0.0;
        instr->params[j].delay  = 0.0;
        instr->params[j].cutoff = 0.0;
        instr->params[j].width  = 0.0;
    }
}

/*===============================================================================*/
/*                                                                               */
/*                GENERAL DECAY MODEL CONFIGURATION PARAMETERS                   */
/*                                                                               */
/*===============================================================================*/

int bayes_ConfigDecayModelMaxNumberOfModulationPeriodsConsidered(void)
{
    return (2);
}

/*===============================================================================*/
/*                                                                               */
/*                        HYPERPARAMETER OPTIMISATION                            */
/*                                                                               */
/*===============================================================================*/

int bayes_ConfigUseFullBayesianHyperParamDetermination(void)
{
    return (0);
}


/*===============================================================================*/
/*                                                                               */
/*                        INSTRUMENT RESPONSE ESTIMATION                         */
/*                                                                               */
/*===============================================================================*/

int bayes_InstrRspCoarseGuessValuesFromLoadedTransient(float *trans,
                                                       float  binwidth,
                                                       int    nbins,
                                                       float *delay,
                                                       float *width)
{
    double p=1.0, max;
    int    i, maxi;

    if ((!trans) || (binwidth<=0.0) || (nbins<=0))
        return (-1);

    /* Initial delay parameter estimate... */
    for (i=0,max=0.0,maxi=0; i<nbins; i++)
    {
        if (trans[i]>max)
        {
            max  = trans[i];
            maxi = i;
        }
    }

    *delay = (float)(((double)maxi+0.5)*binwidth);

    /* Initial width parameter estimate... */
	for (i=maxi; i>=0; i--)
	{
		if (trans[i]<max/2.0) 
		{
			p = (double)i;
			break;
		}
	}

    *width = (float)(2.0*((double)maxi-p)*binwidth/BAYES_INSTR_FWHM_TO_SIGMA_FACTOR);
/*
    for (i=0; i<nbins; i++)
	{
		if (trans[i]>max/2.0) 
		{
			p = (double)i;
			break;
		}
	}
*/
//	*width = (*delay-(p*binwidth))/(sqrt(-log(0.5)));

    return (0);
}


int bayes_InstrRspCoarseGuessValuesFromSmoothedTransient(/* Transient data in... */
                                                         float   *trans,
                                                         int      transbins,
                                                         float    transbinwidth,
                                                         int      fitstart,
                                                         int      fitend,
                                                         /* Guesses out... */
                                                         float   *delay1,
                                                         float   *delay2,
                                                         float   *width1,
                                                         float   *width2)
{
    int    ret;
    int    i, maxi, j, max1, max2;
    int    *data, nbins, rebinning, rebinfactor, nphotons;         
    double interval, *binwalls, temp;
    int    *delta, *convsum, *deltaconvsum, convwidth=3, pos2neg, pos2negconvsum;

    /* Coarse rebinning to smooth the transient... */
    rebinning   = 1;
    rebinfactor = max(transbins/16, 1);
    ret         = bayes_CheckAndTranformTransientDataForBayesFitting(trans,transbins,transbinwidth,fitstart,fitend,
                                                                         rebinning,rebinfactor,
                                                                         &data,&nbins,&binwalls,&nphotons,&interval);

    delta        = Bayes_ivector(0,nbins-1);
    convsum      = Bayes_ivector(0,nbins-1);
    deltaconvsum = Bayes_ivector(0,nbins-1);

    /* Data gradient... */
    for (i=fitstart+1,delta[fitstart]=0; i<nbins; i++)
    {
        delta[i] = data[i]-data[i-1];
    }

    /* Determine how many times data gradient turns +ve to -ve... */
    for (i=fitstart+1,pos2neg=0; i<nbins; i++)
    {
        if ((delta[i]<0) && (delta[i-1]>=0))
            pos2neg++;
    }

    /* Convolution with a box, further smoothes the noise... */
    for (i=fitstart; i<nbins; i++)
    {
        for (j=i,convsum[i]=0; j<i+convwidth; j++)
            convsum[i] += data[j];
    }

    /* Convoluted data gradient... */
    for (i=fitstart+1,deltaconvsum[fitstart]=0; i<nbins; i++)
    {
        deltaconvsum[i] = convsum[i]-convsum[i-1];
    }

    /* Determine how many times data gradient turns +ve to -ve... */
    for (i=fitstart+1,pos2negconvsum=0; i<nbins-convwidth; i++)
    {
        if ((deltaconvsum[i]<0) && (deltaconvsum[i-1]>=0))
            pos2negconvsum++;
    }

    /* Locate the 'pos2negconvsum' largest +ve gradients in transient... */
    for (i=fitstart,j=0,max1=0,maxi=fitstart; i<nbins; i++)
    {
        if (delta[i]>max1)
        {
            max1  = delta[i];
            maxi = i;
        }
    }

	*delay1 = (float)(((float)maxi + 0.5)*(float)rebinfactor*transbinwidth);

    for (i=fitstart,j=0,max2=0; i<nbins; i++)
    {
        if ((delta[i]>max2) && (delta[i]<max1))
        {
            max2  = delta[i];
            maxi = i;
        }
    }

	*delay2 = (float)(((float)maxi + 0.5)*(float)rebinfactor*transbinwidth);

    if (*delay2<*delay1)
    {
        temp    = *delay1;
        *delay1 = *delay2;
        *delay2 = (float)temp;
    }

    *width1 = 0.1f;
    *width2 = 0.2f;


#if 0
    /* Initial delay parameter estimate... */
    *delay = ((float)maxi+0.5)*(float)rebinfactor*transbinwidth;

    /* Initial width parameter estimate... */
	for (i=maxi,max=data[maxi]; i>=0; i--)
	{
		if (data[i]<max/2.0) 
		{
			p = (double)i;
			break;
		}
	}
#endif
	//*width = *delay-(p*(float)rebinfactor*transbinwidth)/(sqrt(-log(0.5)));

    free_Bayes_ivector(delta,0,nbins-1);
    free_Bayes_ivector(convsum,0,nbins-1);
    free_Bayes_ivector(deltaconvsum,0,nbins-1);

    return (0);
}


float bayes_InstrRspCoarseGuessAvgDecayTimeRelativeToInstrRsp(/* Transient data in... */
                                                              float   *trans,
                                                              int      transbins,
                                                              float    transbinwidth,
                                                              int      fitstart,
                                                              int      fitend)
{
    int    ret;
    int    i, maxi, j, max;


    int    *data, nbins, rebinning, rebinfactor, nphotons, *leftsum, leftwidth=3;         
    double interval, avg, ttl, *binwalls;

    /* Coarse rebinning to smooth the transient... */
    rebinning   = 1;
	rebinfactor = max(transbins / 16, 1);
	ret         = bayes_CheckAndTranformTransientDataForBayesFitting(trans,transbins,transbinwidth,fitstart,fitend,
                                                                         rebinning,rebinfactor,
                                                                         &data,&nbins,&binwalls,&nphotons,&interval);

    

    /* Not perfect, but search from the right for maximum then the drop... */
    leftsum = Bayes_ivector(0,nbins-1);

    for (i=nbins-1; i>=fitstart; i--)
    {
        for (j=i,leftsum[i]=0; j>=i-leftwidth; j--)
        {
            leftsum[i] += data[j];
        }
    }

    for (i=nbins-1,max=0,maxi=0; i>=fitstart; i--)
    {
        if (leftsum[i]>max)
        {
            max  = leftsum[i];
            maxi = i;
        }
    }

    /* Compute average arrival time following the peak... */
    for (i=maxi,avg=0.0,ttl=0; i<nbins; i++)
    {
        if (leftsum[i])
        {
            avg += (0.5+(float)(i-maxi))*leftsum[i];
            ttl += leftsum[i];
        }
    }

    free_Bayes_ivector(leftsum,0,nbins-1);


    return ((float)rebinfactor*transbinwidth*(float)avg/(float)ttl);  // weighted average bin = avg/nPhotons, x by interval/nbins to get average time
}


int bayes_DoBayesInstrRspEstimation(/* Data in... */
                                        float                    *trans,
                                        int                       transbins,
                                        float                     transbinwidth,
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
                                        BayesIrEstConfig_t       *BayesIrEstConfig,
                                        float                     modulationperiod,
                                        /* Data out... */
										float                    *param_mp, /* Parameters input using conventional model and values (i.e. Z,A1,tau1,A2,tau2...) */
										float                    *fitted)
{
    int   ret = BAYES_ERR_NO_ERROR;
    
    float uc, sigma, uc1, sigma1, uc2, sigma2, tau, alpha;
    double interval, *binwalls=NULL, minuslogprob, weights_mp[3], taus_mp[3];
    int    ncomponents, *data, nbins, rebinning, rebinfactor, ndecays, paramsfree[] = {1,1,1,1,1,1,1};
    BayesUserFixedParams_t paramfixing;

    /* Check input transient for validity and convert and compute   */
    /* as required for subsequent calling of Bayesian algorithms... */
    rebinning   = bayes_GetBayesTransientRebinningActiveFlag();
    rebinfactor = bayes_GetBayesTransientRebinningFactor();
    ret         = bayes_CheckAndTranformTransientDataForBayesFitting(trans,transbins,transbinwidth,fitstart,transbins,
                                                                         rebinning,rebinfactor,
                                                                         &data,&nbins,&binwalls,nphotons,&interval);

    if (ret<BAYES_ERR_NO_ERROR) 
    {
	    if (binwalls) free_Bayes_dvector(binwalls,0,nbins);
	    if (data)     free_Bayes_ivector(data,0,nbins-1);
	    return (ret); /* Data related error has occured... */
    }

    /* Initial decay and IRF parameter estimates are required for search algorithm... */
    /* Data certainly exists, use for initial estimates... */
    /* NEEDS TO BE UPDATED TO USE THE INTEGER DATA, POTENTIALLY REBINNED!!! */
    /* Coarse guess of centre and width of instr. rsp. from the decay data... */
    ret = bayes_InstrRspCoarseGuessValuesFromLoadedTransient(trans,transbinwidth,transbins,&uc,&sigma);

    bayes_InstrRspCoarseGuessValuesFromSmoothedTransient(trans,transbins,transbinwidth,fitstart,transbins,
                                                         &uc1,&uc2,&sigma1,&sigma2);

    /* A prompt is loaded, improve initial estimates if possible and establish prior if required... */
    if (prompt)
    {
	    /* Coarse guess of centre and width of the instr rsp pulse... */
	    ret = bayes_InstrRspCoarseGuessValuesFromLoadedInstr(prompt,promptbinwidth,nprompt,&uc1,&sigma1);
        uc2    = 1.5f*uc1;
        sigma2 = 2.0f*sigma1;

#if 0 //configurable prior - code needs finishing
        instr->paramsprior[0].widthprior                    = 1;
        instr->paramsprior[0].widthpriorevaluator.nparams   = 2;
        instr->paramsprior[0].widthpriorevaluator.params[0] = sigma1;
        instr->paramsprior[0].widthpriorevaluator.params[1] = sigma1/2.0;
        instr->paramsprior[0].widthpriorevaluator.funk      = MinusLogGaussianPrior;
#endif
    }

    /* Initialise the instr. rsp. container for subsequent searching... */
    ncomponents = BayesIrEstConfig->bayesirnumcomponents;
    instr->paramsfixed[0].cutofffixed = 1;
    instr->params[0].cutoff           = 0.0;
    instr->paramsfixed[1].cutofffixed = 1;
    instr->params[1].cutoff           = 0.0;
    
    instr->ninstr = ncomponents;
    if (!instr->paramsfixed[0].weightfixed) instr->params[0].weight = (1==ncomponents)?(1.00):(0.80);
    if (!instr->paramsfixed[0].delayfixed)  instr->params[0].delay  = uc;
    if (!instr->paramsfixed[0].widthfixed)  instr->params[0].width  = sigma1;
    if (!instr->paramsfixed[0].cutofffixed) instr->params[0].cutoff = 0.0;
    
    if (ncomponents>1)
    {
        if (!instr->paramsfixed[1].weightfixed) instr->params[1].weight = (2==ncomponents)?(0.20):(0.10);
        if (!instr->paramsfixed[1].delayfixed)  instr->params[1].delay  = uc2;
        if (!instr->paramsfixed[1].widthfixed)  instr->params[1].width  = sigma2;
        if (!instr->paramsfixed[1].cutofffixed) instr->params[1].cutoff = 0.0/*uc1*/;

        if (ncomponents>2)
        {
            if (!instr->paramsfixed[2].weightfixed) instr->params[2].weight = 0.10;
            if (!instr->paramsfixed[2].delayfixed)  instr->params[2].delay  = uc1*0.90;
            if (!instr->paramsfixed[2].widthfixed)  instr->params[2].width  = sigma1;
            if (!instr->paramsfixed[2].cutofffixed) instr->params[2].cutoff = 0.0;        
        }
        else
        {
            instr->params[2].weight = 0.0;
            instr->params[2].delay  = 0.0;
            instr->params[2].width  = 0.0;
            instr->params[2].cutoff = 0.0;
        }
    }
    else
    {
        instr->params[1].weight = 0.0;
        instr->params[1].delay  = 0.0;
        instr->params[1].width  = 0.0;
        instr->params[1].cutoff = 0.0;
        instr->params[2].weight = 0.0;
        instr->params[2].delay  = 0.0;
        instr->params[2].width  = 0.0;
        instr->params[2].cutoff = 0.0;
    }

    /* Coarse estimate of decay parameters as initialisation for search... */
    //avg = data_ComputeBinnedDataAverageArrTime(data,nbins,fitstart,*nphotons,interval);

    tau   = bayes_InstrRspCoarseGuessAvgDecayTimeRelativeToInstrRsp(trans,transbins,transbinwidth,fitstart,transbins);
    alpha = 1.0f/tau;

    weights_mp[0] = 0.05; //this could (and should really!!!) be guessed...
    weights_mp[1] = 0.70;
    weights_mp[2] = 0.25;
    taus_mp[1]    = tau;
    taus_mp[2]    = tau/2.0;

    /* Estimate the instrument response and the decay parameters... */
    if (modeltype==FIT_IRFANDMONOEXP)
    //if (bayes_InstrConfigUseMonoExpDecayModelForInstrEstimation())
    {
        ndecays = 1;
        //Determine optimal values using a mono-exponential fit
        ret = bayes_DirectInstrRspAndMonoExpOptimization(/* Data in... */ 
                                                         data,nbins,fitstart,binwalls,*nphotons,
                                                         /* Estimates out... */
                                                         &weights_mp[0],&taus_mp[1],instr,&minuslogprob,
                                                         /* Instrument... */
                                                         interval,modulationperiod,alpha);

        if (ret>=0)
        {
            weights_mp[1] = 1.0-weights_mp[0];
        
            param_mp[0] = bayes_FromBayesModelTransformToParamZ((float)weights_mp[0],nbins,*nphotons);
            param_mp[1] = bayes_FromBayesModelWeightAndTauToParamA((float)weights_mp[1], (float)taus_mp[1], (float)interval/(float)nbins, (float)interval,/**delay*/(float)instr->params[0].delay,*nphotons);
            param_mp[2] = (float)taus_mp[1];
        }
    }
    else
    {
        //Determine optimal values using a multi-exponential fit
        ndecays = 2;
        int nparams = 4;

        ret = bayes_CheckParameterValueFixingForBayesFitting(&paramfixing,nparams,paramsfree,NULL,
                                                                 nbins,fitstart,binwalls,*nphotons,
                                                                 interval,modulationperiod,instr);

        if (ret<BAYES_ERR_NO_ERROR) 
        {
	        if (binwalls) free_Bayes_dvector(binwalls,0,nbins);
	        if (data)     free_Bayes_ivector(data,0,nbins-1);
	        return (ret); /* Parameter fixing related error has occured... */
        }

        ret = bayes_DirectInstrRspAndMultiExpOptimization(/* Data in... */ 
                                                          data,nbins,fitstart,binwalls,*nphotons,
                                                          /* Estimates out... */
                                                          ndecays,weights_mp,taus_mp,&paramfixing,instr,&minuslogprob,
                                                          /* Instrument... */
                                                          interval,modulationperiod,alpha);
        if (ret>=0)
        {
            param_mp[0] = bayes_FromBayesModelTransformToParamZ((float)weights_mp[0],nbins,*nphotons);
            param_mp[1] = bayes_FromBayesModelWeightAndTauToParamA((float)weights_mp[1], (float)taus_mp[1], (float)interval/(float)nbins, (float)interval,/**delay*/(float)instr->params[0].delay,*nphotons);
            param_mp[2] = (float)taus_mp[1];
            param_mp[3] = bayes_FromBayesModelWeightAndTauToParamA((float)weights_mp[2], (float)taus_mp[2], (float)interval/(float)nbins, (float)interval,/**delay*/(float)instr->params[0].delay,*nphotons);
            param_mp[4] = (float)taus_mp[2];
        }
    }

    if (ret>=BAYES_ERR_NO_ERROR)
        bayes_SortInstrRspComponentsByWeight(instr);

    if ((fitted) && (ret>=BAYES_ERR_NO_ERROR))
    {
        if (rebinning) /* Need to recompute binwalls... */
        {
            interval = transbinwidth*(float)transbins;
            free_Bayes_dvector(binwalls,0,nbins);
            binwalls = Bayes_dvector(0,transbins);
            bayes_PopulateBinWallsVectorUniformIntervals(binwalls,transbins,interval);
        }

        bayes_FitPredictedDecay(fitted,transbins,fitstart,fitend,binwalls,
                                instr, (float)interval,modulationperiod,
                                ndecays,taus_mp,weights_mp,*nphotons);
    }
	
	// Record the values
	BayesIrEstConfig->bayesirweight1 = instr->params[0].weight;
	BayesIrEstConfig->bayesirsigma1 = instr->params[0].width;
	BayesIrEstConfig->bayesiruc1 = instr->params[0].delay;
	BayesIrEstConfig->bayesircutoff1 = instr->params[0].cutoff;
	BayesIrEstConfig->bayesirweight2 = instr->params[1].weight;
	BayesIrEstConfig->bayesirsigma2 = instr->params[1].width;
	BayesIrEstConfig->bayesiruc2 = instr->params[1].delay;
	BayesIrEstConfig->bayesircutoff2 = instr->params[1].cutoff;
	BayesIrEstConfig->bayesirweight3 = instr->params[2].weight;
	BayesIrEstConfig->bayesirsigma3 = instr->params[2].width;
	BayesIrEstConfig->bayesiruc3 = instr->params[2].delay;
	BayesIrEstConfig->bayesircutoff3 = instr->params[2].cutoff;

	return (ret);
}

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
                                                  BayesInstrRsp_t *instr)
{
    return (bayes_UpdateEstimatedInstrRspDetailsStoreAdv(fileprompt,filedata,
                                                      allpixels,xcoordinate,ycoordinate,binningtype,binninglevel,
                                                      trans,transbins,transbinwidth,fitstart,fitend,nphotons,modeltype,
                                                      instr));
}
/* Rapid grids setup */

/* IRF estimation setup */

BayesIrEstConfig_t bayes_GetIrEstConfig(void)
{
	return (BayesIrEstConfig_t) {
        .bayesirnumcomponents = 1,
        .bayesirweight1 = 1.0,
        .bayesircutoff1 = 0.0,
        .bayesirsigma1 = 0.001,
        .bayesiruc1 = 0.0,
        .bayesirweight2 = 0.0,
        .bayesircutoff2 = 0.0,
        .bayesirsigma2 = 0.0,
        .bayesiruc2 = 0.0,
        .bayesirweight3 = 0.0,
        .bayesircutoff3 = 0.0,
        .bayesirsigma3 = 0.0,
        .bayesiruc3 = 0.0
    };
}

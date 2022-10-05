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
/***********************************************************************************/
/*                                                                                 */
/*                             DEFINES / SWITCHES                                  */
/*                                                                                 */
/***********************************************************************************/

/* None. */


/***********************************************************************************/
/*                                                                                 */
/*                             HEADER FILE INCLUSION                               */
/*                                                                                 */
/***********************************************************************************/

#include "stdlib.h"
#include "bayes_Interface.h"

#include "extmath.h"
#include "random.h"
#include "matrices.h"
#include "DTYPE.h"

#include "bayes_Sizes.h"
#include "bayes_Types.h"
#include "bayes_DataManagement.h"
#include "bayes_DistributionFctsBinLikelihoods.h"
#include "bayes_ModelTransformTools.h"
#include "bayes_MultiExpAnalysisBinLikelihoods.h"


/***********************************************************************************/
/*                                                                                 */
/*                         STATIC FUNCTION PROTOTYPES                              */
/*                                                                                 */
/***********************************************************************************/

/* None. */


/***********************************************************************************/
/*                                                                                 */
/*        (BINNED DATA) MULTI-EXPONENTIAL ANALYSIS SPECIFIC MODEL FUNCTIONS        */
/*                                                                                 */
/***********************************************************************************/



#define	FREE_DECAYPHOTONLIKELIHOODS	{\
	int ii;\
	if (decayphotonlikelihoods[0]) {free_Bayes_dvector(decayphotonlikelihoods[0],0,nbins-1);}\
	for (ii=1;ii<=ndecays;ii++) {\
	if (decayphotonlikelihoods[ii])\
	if (params->paramfixing->tauuserfixed[ii] == BAYES_PARAM_VALUE_FREE)\
	free_Bayes_dvector(decayphotonlikelihoods[ii],0,nbins-1);}\
	free(decayphotonlikelihoods);\
}

double bayes_MultiExpMinusLogProbParams(double *x, int id, void *container)
{
    int    *data, bin, nbins, nphotons, nphotonsbin, ndecays, i, k, modelfixedweightindex, type, ret;
    double  interval, period, *binwalls, bH, bL, bjoverT, norm;
    double  hyperparam, val, temp, *weights, *taus, **decayphotonlikelihoods;

    MultiExpMinusLogProbParams_t *params;
    BayesInstrRsp_t              *instr;

    params  = (MultiExpMinusLogProbParams_t *)(container);
    ndecays = params->ndecays;

    /* Determine which weight value is fixed by the model (i.e. weights summing to '1')... */
    for (k=0; k<=ndecays; k++)
    {
        if (params->paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            modelfixedweightindex = k;
            break;
        }
    }

    weights = Bayes_dvector(0,ndecays);
    taus    = Bayes_dvector(1,ndecays);

    /* Extract and validate the weights from 'x' first... */
    for (k=1,i=1; k<=ndecays; k++)
    {
        if (params->paramfixing->weightuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            weights[k] = x[i];
            i++;
        }
        else
        {
            weights[k] = params->paramfixing->weights[k];
        }

        if ((weights[k]<0.0) || (weights[k]>1.0)) /* Check that weights are valid... */
        {
            free_Bayes_dvector(weights,0,ndecays);
            free_Bayes_dvector(taus,1,ndecays);
            return (BAYES_SIZE_DOUBLE_HUGE);
        }
    }

    for (k=0,val=1.0; k<=ndecays; k++)
    {
        if (k != modelfixedweightindex)
            val -= weights[k];
    }

    if ((val<0.0) || (val>1.0))
    {
        free_Bayes_dvector(weights,0,ndecays);
        free_Bayes_dvector(taus,1,ndecays);
        return (BAYES_SIZE_DOUBLE_HUGE);
    }

    weights[modelfixedweightindex] = val;

    /* Now extract and validate the lifetime indices from 'x'... */
    for (k=1; k<=ndecays; k++)
    {
        if (params->paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            taus[k] = x[i];
            i++;
        }
        else
        {
            taus[k] = params->paramfixing->taus[k]; /* User fixed value... */
        }

        if ((taus[k]<SMALLTIME) || (taus[k]>params->interval)) /* Check that lifetimes are valid... */
        {
            free_Bayes_dvector(weights,0,ndecays);
            free_Bayes_dvector(taus,1,ndecays);
            return (BAYES_SIZE_DOUBLE_HUGE);
        }
    }

    for (k=1; k<ndecays; k++)
    {
        if (taus[k]<taus[k+1]) /* Reduce the search space by requiring decreasing lifetime with increasing index... */
        {
            free_Bayes_dvector(weights,0,ndecays);
            free_Bayes_dvector(taus,1,ndecays);
            return (BAYES_SIZE_DOUBLE_HUGE);
        }
    }

    /* Weight and lifetime values okay... */
    data        = params->data;
    nbins       = params->nbins;
    //fitstart    = params->fitstart;
    binwalls    = params->binwalls;
    nphotons    = params->nphotons;
    interval    = params->interval;
    period      = params->modulationperiod;
    instr       = params->instr;
    hyperparam  = params->hyperparam;

    /* Compute the minus log likelihood using pre-computed bin likelihoods... */
//    decayphotonlikelihoods = (double **) malloc((1+ndecays)*sizeof(double*));
    decayphotonlikelihoods = (double **) calloc((1+ndecays), sizeof(double*));  // calloc to init everything to NULL
    
    for (k=1; k<=ndecays; k++)
    {
        if (params->paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            decayphotonlikelihoods[k] = Bayes_dvector(0,nbins-1);

            ret = bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(decayphotonlikelihoods[k],
                                                                             nbins,binwalls,data,
                                                                             interval,period,instr,
                                                                             taus[k],
                                                                             0,NULL,NULL);

            if (ret<0)
            {
                free_Bayes_dvector(weights,0,ndecays);
                free_Bayes_dvector(taus,1,ndecays);
				FREE_DECAYPHOTONLIKELIHOODS;
	            return (BAYES_SIZE_DOUBLE_HUGE);
            }
        }
        else
        {
            decayphotonlikelihoods[k] = params->paramfixing->
                fluorescencelikelihoods[k].fluorescencedecayphotonlikelihoodsgiventau;
        }
    }
    

#if 0
    likelihoods = Bayes_dmatrix(1,ndecays,0,nbins);

//    for (i=1; i<=ndecays; i++)
  //      bayes_ArrBinLikelihoodsGivenTau(likelihoods[i],NULL,NULL/*upsilon1*/,data,nbins,interval,period,instr,taus[i]);

    for (i=1; i<=ndecays; i++)
        bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(likelihoods[i],
                                                                   nbins,binwalls,data,
                                                                   interval,period,instr,
                                                                   taus[i],
                                                                   0,NULL,NULL);
#endif

    ret = bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(&norm,interval,params->modulationperiod,/*binwalls[fitstart]*/0.0,params->instr,ndecays,weights,taus);

    if (ret<0)
    {
        free_Bayes_dvector(weights,0,ndecays);
        free_Bayes_dvector(taus,1,ndecays);
		FREE_DECAYPHOTONLIKELIHOODS;
        return (BAYES_SIZE_DOUBLE_HUGE);
    }

    for (i=1, val=0.0; i<=ndecays; i++)
        val += (hyperparam*taus[i]); 

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            bL      = binwalls[bin];
            bH      = binwalls[bin+1];
            bjoverT = (bH-bL)/interval;

            for (i=1, temp=weights[0]*bjoverT; i<=ndecays; i++)
                temp += (weights[i]*decayphotonlikelihoods[i][bin])/norm;

            val -= (double)(nphotonsbin)*(log(temp));

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(val,&type))
            {
                free_Bayes_dvector(weights,0,ndecays);
                free_Bayes_dvector(taus,1,ndecays);
				FREE_DECAYPHOTONLIKELIHOODS;
	            return (BAYES_SIZE_DOUBLE_HUGE);
            }
        }
    }

    free_Bayes_dvector(weights,0,ndecays);
    free_Bayes_dvector(taus,1,ndecays);
	FREE_DECAYPHOTONLIKELIHOODS;

	return (val);
}

#undef	FREE_DECAYPHOTONLIKELIHOODS

double bayes_MultiExpProbParams(double *x, int id, void *container)
{
    double                        val, norm;
    MultiExpMinusLogProbParams_t *params;

    params = (MultiExpMinusLogProbParams_t *)(container);
    norm   = params->normalization;

    val = bayes_MultiExpMinusLogProbParams(x,id,container);
    val = exp(val+norm);

    return (val);
}
/*
int bayes_RapidMultiExpMostProbWeightsAndTausSimulatedAnnealing(int                          *data,
                                                                int                           nbins,
                                                                double                       *binwalls,
                                                                int                          *nphotons,
                                                                int                           ndecays,
                                                                double                       *weights,
                                                                double                       *taus,
                                                                BayesUserFixedParams_t       *paramfixing,
                                                                double                        interval,
                                                                double                        modulationperiod,
                                                                BayesInstrRsp_t              *instr,
                                                                double                        alpha,
                                                                BayesPsuedoRapidValueStore_t *grid,
                                                                double                       *val)*/
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
                                               double                   *val)

{
	int     ret, id=0;
	double  value, bestvalue, *x, *xbest, temp;
    int     i, k, modelfixedweightindex, ndim, restarts, nrestarts, nweightsfree, ntausfree;
    int   (*minimizer)(double (*)(double *, int, void *), int, void *, int, double *, double *, void *);
    void   *config=NULL;
    double deltas[] = {-1.0,0.05,0.05,0.25,0.25,0.05,0.25,0.05,0.25,0.05,0.25};

    MultiExpMinusLogProbParams_t container;

    container.data             =  data;
    container.nbins            =  nbins;
    container.binwalls         =  binwalls;
    container.nphotons         =  *nphotons;
    container.nparams          =  2*ndecays+1;
    container.ndecays          =  ndecays;
    container.paramfixing      =  paramfixing;
    container.interval         =  interval;
    container.modulationperiod =  modulationperiod;
    container.instr            =  instr;
    container.hyperparam       =  alpha;
    container.normalization    =  0.0;

    /* Initialise input starting location 'paramsfree' first... */
    ndim                       = 2*ndecays;
    ndim                      -= paramfixing->nparamsuserfixed;
    x                          = Bayes_dvector(1,ndim);

    for (k=1,i=1; k<=ndecays; k++)
    {
        if (paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            x[i] = weights[k];
            i++;
        }
    }

    for (k=1; k<=ndecays; k++)
    {
        if (paramfixing->tauuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            x[i] = taus[k];
            i++;
        }
    }

    value = BAYES_SIZE_DOUBLE_HUGE;

    minimizer = &math_MinimiseFctDoubleWithGenericContainer;
    config = malloc(sizeof(AmoebaConfigParams_t));
    ((AmoebaConfigParams_t*)config)->monitor   = 0;
    ((AmoebaConfigParams_t*)config)->tolerance = bayes_MonoExpConfigGetDownhillSimplexPrecision();

	ret = minimizer(bayes_MultiExpMinusLogProbParams,id,(void*)(&container),ndim,x,&value,(void*)(config));

	if (ret >= MATH_MINIMISATION_RESULT_SUCCESS)
	{
        nrestarts = bayes_BiExpConfigGetNumberOfRestarts();

		if (nrestarts)
	    {
	        xbest = Bayes_dvector(1,ndim);

	        for (i=1; i<=ndim; i++)
	            xbest[i] = x[i];

	        bestvalue = value;

            for (k=0, nweightsfree=ndecays; k<=ndim; k++)
                if (paramfixing->weightuserfixed[k] == BAYES_PARAM_VALUE_USER_FIXED)
                    nweightsfree--;

            for (k=1, ntausfree=ndecays; k<=ndim; k++)
                if (paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_USER_FIXED)
                    ntausfree--;

            rand_InitializeRandomSeed();
		    
            for (restarts=0; restarts<nrestarts; restarts++)
		    {
	            /* Randomly initialize the search starting location in parameter space... */
                for (i=1; i<=nweightsfree; i++)
                    x[i] = rand_RandomDouble();

                for (; i<=nweightsfree; i++)
                    x[i] = rand_RandomDouble();

	            /* Perform the search... */
                ret = minimizer(bayes_MultiExpMinusLogProbParams,id,(void*)(&container),ndim,x,&value,(void*)(config));

	            /* Update the best result if required... */
		        if ((value < bestvalue) && (ret >= MATH_MINIMISATION_RESULT_SUCCESS))
		        {
	                bestvalue = value;

	                for (i=1; i<=ndim; i++)
				        xbest[i] = x[i];
	            }
	        }

	        /* Put our best results back into the free parameter array... */
	        for (i=1; i<=ndim; i++)
	            x[i] = xbest[i];

	        value = bestvalue;

	        free_Bayes_dvector(xbest,1,ndim);
	    }	
    }

    /* Recombine the fixed and optimal free parameter values to yield parameter vector... */
    if (ret >= MATH_MINIMISATION_RESULT_SUCCESS)
    {
        for (k=1,i=1; k<=ndecays; k++)
        {
            if (paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
            {
                weights[k] = x[i];
                i++;
            }
            else
            {
                weights[k] = paramfixing->weights[k];
            }
        }

        for (k=0; k<=ndecays; k++)
        {
            if (paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
            {
                modelfixedweightindex = k;
                break;
            }
        }

        for (k=0,temp=1.0; k<=ndecays; k++)
        {
            if (k != modelfixedweightindex)
                temp -= weights[k];
        }

        weights[modelfixedweightindex] = temp;

        for (k=1; k<=ndecays; k++)
        {
            if (paramfixing->tauuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
            {
                taus[k] = x[i];
                i++;
            }
            else
            {
                taus[k] = paramfixing->taus[k];
            }
        }

	    *val = (float)(value);
	}
	
    free_Bayes_dvector(x,1,ndim);
	if (config) free(config);

    return (ret);
}


#if 0  // some code in preparation for doing integration (i.e. computing bayes factor) using transformed lifetimes
double bayes_ComputeTauFromTransformedLifetimeValue(double y)
{
    return ((1.0-y)/y);
}


double bayes_MultiExpProbParamsTransformedTaus(double *x, int id, void *container)
{
    double                        norm, minuslogprob, prob, temp;
    int                          *isparamfixed, nparams, nparamsfixed;
    int                          *isweightfixed, *istaufixed, ndecays, i;
    double                       *paramsfixed, *weights, *taus, *taus_transformed;
    float                        *y;
    MultiExpMinusLogProbParams_t *params;

    int    type;

    params        = (MultiExpMinusLogProbParams_t *)(container);
    nparams       = params->nparams;
    isparamfixed  = params->isparamfixed;
    nparamsfixed  = params->nparamsfixed;
    paramsfixed   = params->paramsfixed;
    isweightfixed = params->isweightfixed;
    istaufixed    = params->istaufixed;

    y = Bayes_vector(0,nparams-1);

    bayes_PopulateParamVectorFromFreeAndFixedVectors //reconstruct the parameter vector from free and fixed
        (y,nparams,x,nparams-nparamsfixed,paramsfixed,nparamsfixed,isparamfixed);

    bayes_AllocateWeightsAndTausVectors //separate the weights and lifetimes
        (nparams,&ndecays,&weights,&taus_transformed);

    bayes_PopulateWeightsAndTausVectorsFromParamVector(nparams,y,weights,taus_transformed);

    taus = Bayes_dvector(1,ndecays);

    for (i=1; i<=ndecays; i++) // transform the lifetimes to actual taus
        taus[i] = bayes_ComputeTauFromTransformedLifetimeValue(taus_transformed[i]);

    bayes_UpdateWeightsVectorModelDefinedValue(weights,ndecays+1,isweightfixed);

    // now need to recombine so that we can pass things to the minus log prob fct
    bayes_PopulateParamVectorFromWeightsAndTausVectors
        (y,nparams,weights,taus);
    
    bayes_PopulateFreeAndFixedVectorsFromParamVector
        (y,nparams,x,nparams-nparamsfixed,paramsfixed,nparamsfixed,isparamfixed);


////////////////////////
    // extract transformed lifetimes
    // transform the lifetimes to actual taus
    // compute minus log prob using taus
    // divide result by each transformed lifetime squared
    
    params = (MultiExpMinusLogProbParams_t *)(container);
    norm   = params->normalization;

    minuslogprob = bayes_MultiExpMinusLogProbParams(x,id,container);

    for (i=1, temp=0.0; i<=ndecays; i++)
        temp -= 2.0*log(taus_transformed[i]);

    prob = exp(-minuslogprob+norm+temp);

    if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(prob, &type))
        prob = 1.0e-30; //return (BAYES_SIZE_DOUBLE_HUGE);

    if (prob <= 1.0e-30)
        prob = 1.0e-30;

    if (prob >= 1.0e+30)
        prob = 1.0e+30;

    free_Bayes_vector(y,0,nparams-1);
    free_Bayes_dvector(taus,1,ndecays);
    bayes_FreeWeightsAndTausVectors(ndecays,weights,taus_transformed);

    return (prob);
}
#endif

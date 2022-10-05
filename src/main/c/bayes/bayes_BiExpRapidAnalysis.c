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

#include "stdio.h"
#include "stdlib.h"
#include "bayes_Interface.h"

#include "extmath.h"
#include "random.h"
#include "matrices.h"
#include "DTYPE.h"

#include "bayes_Sizes.h"
#include "bayes_DataManagement.h"
#include "bayes_DistributionFctsBinLikelihoods.h"
#include "bayes_ModelTransformTools.h"
#include "bayes_ModelSelection.h"
#include "bayes_MultiExpRapidAnalysis.h"
#include "bayes_RapidBayesDecayAnalysis.h"



struct BayesContainerRapidBiExpMinusLogProbParams
{
//    int                           ndecays;
    int                     *data;
    int                      nbins;
    int                      fitstart;
    int                      nphotons;
    double                  *binwalls;
    double                   interval;
    double                   modulationperiod;
    BayesInstrRsp_t         *instr;
    double                   hyperparam;
    BayesRapidValueStore_t  *rapidparamsandlikelihoods;
    BayesUserFixedParams_t  *paramfixing;
    int                      fast;
    /* Hyperparameter optimisation... */
    double                   alphamin;
    double                   alphamax;
    double               ****datalikelihoods;
};
typedef struct BayesContainerRapidBiExpMinusLogProbParams RapidBiExpMinusLogProbParams_t;


struct BayesContainerRapidBiExpSearchReporting
{
    int fast;
};
typedef struct BayesContainerRapidBiExpSearchReporting RapidBiExpSearchReporting_t;


/***********************************************************************************/
/*                                                                                 */
/*                      GLOBALS FOR ENHANCED DEBUG OUTPUT                          */
/*                                                                                 */
/***********************************************************************************/

int                             gggExtendedDebug=0;
int                             gggExtendedDebugNumOfDecays;
int                             gggExtendedDebugNumOfWeightsAsParams;
int                             gggExtendedDebugNumOfTausAsParams;
RapidBiExpMinusLogProbParams_t *gggExtendedDebugParamsContainer;






int bayes_BiExpRapidChechStateWithinPrecomputedGrid(int *x,
                                                    int *low,
                                                    int *high)
{
    int i;

    for (i=1; i<=4; i++)
    {
        if ((x[i]<low[i]) ||
            (x[i]>high[i]))
        {
            return (0);
        }
    }

    return (1);
}



double bayes_BiExpRapidMinusLogProbGivenParamsFast(int *x, int id, void *container, void *report)
{
    int     i, k, modelfixedweightindex;
    int     *data, bin, nbins, fitstart, nphotonsbin;
    double  val, hyperparam;
    int     type;
    int     weightindexes[3], tauindexes[3];
    double  weights[3], taus[3];
    double  *logphotonlikelihoodgivenparams;

    RapidBiExpMinusLogProbParams_t *params;
    RapidBiExpSearchReporting_t    *reporting;

    reporting = (RapidBiExpSearchReporting_t*) report;
    params  = (RapidBiExpMinusLogProbParams_t *)(container);

    /* A null pointer in pre-computed grid means that the weight/tau combination is not valid (or an error occured)... */
    if (!params->rapidparamsandlikelihoods->biexpvaluestore->likelihoods[x[1]][x[2]][x[3]][x[4]])
        return (BAYES_SIZE_DOUBLE_HUGE);

    logphotonlikelihoodgivenparams =
        params->rapidparamsandlikelihoods->biexpvaluestore->likelihoods[x[1]][x[2]][x[3]][x[4]]->logphotonlikelihoodgiventausandweights;

    /* No user weight fixing for pre-computed grid... */
    modelfixedweightindex = 0;

    /* Extract and validate the weight indices from 'x' first... */
    for (k=1,i=1; k<=/*ndecays*/2; k++)
    {
        weightindexes[k] = x[i];
        i++;
    }

    /* Now extract and validate the lifetime indices from 'x'... */
    for (k=1; k<=/*ndecays*/2; k++)
    {
        tauindexes[k] = x[i];
        i++;
    }

    /* Populate weights vector with free, user fixed, and model fixed values... */
    for (k=0; k<=/*ndecays*/2; k++)
    {
        if (k != modelfixedweightindex)
        {
            weights[k] = params->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[weightindexes[k]];
        }
    }

    for (k=0,val=1.0; k<=/*ndecays*/2; k++)
    {
        if (k != modelfixedweightindex)
            val -= weights[k];
    }

    weights[modelfixedweightindex] = val;

    /* Populate lifetimes vector with free and user fixed values... */
    for (k=1; k<=/*ndecays*/2; k++)
    {
        taus[k] = params->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[tauindexes[k]];
    }

    /* Now compute the minus log likelihood... */
    hyperparam = params->hyperparam;
    data       = params->data;
    nbins      = params->nbins;
    fitstart   = params->fitstart;

    for (k=1, val=0.0; k<=/*ndecays*/2; k++)
        val += (hyperparam*taus[k]);

    for (bin=fitstart; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            /* Log bin likelihoods have been precomputed... */
            val -= (double)(nphotonsbin)*logphotonlikelihoodgivenparams[bin];

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(val,&type))
            {
                val = BAYES_SIZE_DOUBLE_HUGE;
                break;
            }
        }
    }

    return (val);
}


double bayes_BiExpRapidMinusLogProbGivenParams(int *x, int id, void *container, void *report)
{
    int    i, k, modelfixedweightindex, ndecays;
    int    *data, bin, nbins, fitstart, nphotonsbin;
    double val, temp, hyperparam, *binwalls, **decayphotonlikelihoods, interval, norm, backgroundmax;
    int    type;
    int    weightindexes[3], tauindexes[3], fast;
    double weights[3], taus[3];

    RapidBiExpMinusLogProbParams_t *params;
    RapidBiExpSearchReporting_t    *reporting;

    reporting = (RapidBiExpSearchReporting_t*) report;

    params  = (RapidBiExpMinusLogProbParams_t *)(container);

    fast = bayes_BiExpRapidChechStateWithinPrecomputedGrid(x,
                                                           params->rapidparamsandlikelihoods->biexpvaluestore->low,
                                                           params->rapidparamsandlikelihoods->biexpvaluestore->high);

    reporting->fast = fast;

    if (fast) /* Use pre-computed log likelihoods if they exist... */
    {
        return (bayes_BiExpRapidMinusLogProbGivenParamsFast(x,id,container,report));
    }

    ndecays = 2;

    /* Determine which weight value is fixed by the model (i.e. weights summing to '1')... */
    for (k=0; k<=ndecays; k++)
    {
        if (params->paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            modelfixedweightindex = k;
            break;
        }
    }

    /* Extract and validate the weight indices from 'x' first... */
    for (k=1,i=1; k<=ndecays; k++)
    {
        if (params->paramfixing->weightuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            weightindexes[k] = x[i];
            i++;

            if ((weightindexes[k]<0) ||
                (weightindexes[k]>params->rapidparamsandlikelihoods->biexpvaluestore->settings->nweights))
            {
                return (BAYES_SIZE_DOUBLE_HUGE); /* Grid location does not exist... */
            }
        }
    }

    /* Now extract and validate the lifetime indices from 'x'... */
    for (k=1; k<=ndecays; k++)
    {
        if (params->paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            tauindexes[k] = x[i];
            i++;

            if ((tauindexes[k]<0) ||
                (tauindexes[k]>params->rapidparamsandlikelihoods->biexpvaluestore->settings->ntaus))
            {
                return (BAYES_SIZE_DOUBLE_HUGE); /* Grid location does not exist... */
            }
        }
    }

    /* Populate weights vector with free, user fixed, and model fixed values... */
    for (k=0; k<=ndecays; k++)
    {
        if (k != modelfixedweightindex)
        {
            if (params->paramfixing->weightuserfixed[k] == BAYES_PARAM_VALUE_FREE)
                weights[k] = params->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[weightindexes[k]];
            else
                weights[k] = params->paramfixing->weights[k];
        }
    }

    for (k=0,val=1.0; k<=ndecays; k++)
    {
        if (k != modelfixedweightindex)
            val -= weights[k];
    }

    weights[modelfixedweightindex] = val;

    /* Populate lifetimes vector with free and user fixed values... */
    for (k=1; k<=ndecays; k++)
    {
        if (params->paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_FREE)
            taus[k] = params->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[tauindexes[k]];
        else
            taus[k] = params->paramfixing->taus[k];
    }

    /* Check that weight values are valid, background can be zero but decay */
    /* component weights should be greater than zero really...              */
    backgroundmax = params->rapidparamsandlikelihoods->biexpvaluestore->settings->backgroundmax;

    if ((weights[0]<0.0) || (weights[0]>backgroundmax) || (weights[0]>1.0))
    {
        return (BAYES_SIZE_DOUBLE_HUGE);
    }

    for (k=1; k<=ndecays; k++)
    {
        if ((weights[k]<=0.0) || (weights[k]>1.0))
        {
            return (BAYES_SIZE_DOUBLE_HUGE);
        }
    }

    /* Check that lifetimes are valid */
    for (k=1; k<=ndecays; k++)
    {
        if (taus[k]<SMALLTIME)
        {
            return (BAYES_SIZE_DOUBLE_HUGE);
        }
    }

    /* Check that lifetime values are in descending order... */
    if (params->paramfixing->nparamsuserfixed==0)
    {
        for (k=1; k<ndecays; k++)
        {
            if ((taus[k]<taus[k+1]))
            {
                return (BAYES_SIZE_DOUBLE_HUGE);
            }
        }    
    }

    /* Now compute the minus log likelihood... */
    hyperparam = params->hyperparam;
    data       = params->data;
    nbins      = params->nbins;
    fitstart   = params->fitstart;
    binwalls   = params->binwalls;
    interval   = params->interval;
    
    /* Weight and lifetime values okay, compute the minus log likelihood using pre-computed bin likelihoods... */
    decayphotonlikelihoods = (double **) malloc((1+ndecays)*sizeof(double*));

    for (k=1; k<=ndecays; k++)
    {
        if (params->paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            decayphotonlikelihoods[k] = params->rapidparamsandlikelihoods->biexpvaluestore->
                fluorescencelikelihoods[tauindexes[k]].fluorescencedecayphotonlikelihoodsgiventau;
        }
        else
        {
            decayphotonlikelihoods[k] = params->paramfixing->
                fluorescencelikelihoods[k].fluorescencedecayphotonlikelihoodsgiventau;
        }
    }

    bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(&norm,interval,params->modulationperiod,binwalls[fitstart],params->instr,ndecays,weights,taus);

    for (k=1, val=0.0; k<=ndecays; k++)
        val += (hyperparam*taus[k]);

    for (bin=fitstart; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            /* Fluorescence decay photon likelihoods have been precomputed... */
            temp = weights[0]*(binwalls[bin+1]-binwalls[bin])/(interval-binwalls[fitstart]);

            for (k=1; k<=ndecays; k++)
                temp += (weights[k]*((decayphotonlikelihoods[k])[bin]))/norm;

            val -= (double)(nphotonsbin)*(log(temp));

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(val,&type))
            {
                val = BAYES_SIZE_DOUBLE_HUGE;
                break;
            }
        }
    }
            
    free(decayphotonlikelihoods);

    return (val);
}



int bayes_BiExpDiscreteSpaceMinimisationExhaustiveSearch(double (*funk)(int *, int, void *, void *),
                                                         int      id,
                                                         void    *container,
                                                         int      ndim,
                                                         int     *where,
                                                         double  *value,
                                                         void    *config,
                                                         int     *nfast,
                                                         int     *nslow)
{
    int    ret = MATH_MINIMISATION_RESULT_SUCCESS;

    double y, yb;
    int    x[5], xb[5], *xmin, *xmax, i, j, k, ell, m, nit;
    BayesProbDistn_t *distn=NULL;

    RapidBiExpSearchReporting_t report;
    
    int    DebugTrace=0, d;
    FILE   *fp = NULL;

    xmin  = ((MultiExpDiscreteGridSearchConfigParams_t *)config)->gridmins;
    xmax  = ((MultiExpDiscreteGridSearchConfigParams_t *)config)->gridmaxs;
    distn = ((MultiExpDiscreteGridSearchConfigParams_t *)config)->distn;

    if ((ndim<=0) || (!where) || (!funk) || (!xmin) || (!xmax))
        return (MATH_MINIMISATION_RESULT_ERROR_INVALID_INPUT);

    /* Use 'where' location as provided by calling routine for the optimal initial values... */
    for (k=1; k<=ndim; k++)
    {
        xb[k] = x[k] = where[k];
    }

    yb = y = funk(x,id,container,(void*)(&report));

    if (DebugTrace)
    {
        if ((fp = fopen("BiExpDiscreteSpaceMinimisationExhaustiveSearchDebugOutput.txt", "w")))
        {
            DebugTrace = 0;
        }
        else
        {
            fprintf(fp,"bayes_BiExpDiscreteSpaceMinimisationExhaustiveSearch ==>\n");

            fprintf(fp,"ndim: %d\n", ndim);

            if (gggExtendedDebug)
            {
                gggExtendedDebugParamsContainer  = (RapidBiExpMinusLogProbParams_t *)(container);
                gggExtendedDebugNumOfDecays      = gggExtendedDebugParamsContainer->paramfixing->nparams/2;

                for (d=1, gggExtendedDebugNumOfTausAsParams=0; d<=gggExtendedDebugNumOfDecays; d++)
                {
                    if (gggExtendedDebugParamsContainer->paramfixing->tauuserfixed[d] != BAYES_PARAM_VALUE_USER_FIXED)
                        gggExtendedDebugNumOfTausAsParams++;
                }

                gggExtendedDebugNumOfWeightsAsParams = gggExtendedDebugNumOfDecays-(gggExtendedDebugParamsContainer->paramfixing->nparamsuserfixed-(gggExtendedDebugNumOfDecays-gggExtendedDebugNumOfTausAsParams));

                fprintf(fp,"\nExtended debugging trace output ==>\n");
                fprintf(fp,"gExtendedDebugNumOfWeightsAsParams: %d, gExtendedDebugNumOfTausAsParams: %d\n", gggExtendedDebugNumOfWeightsAsParams, gggExtendedDebugNumOfTausAsParams);
            }

            fprintf(fp,"Starting location ==> \n");

            for (d=1; d<=ndim; d++)
                fprintf(fp,"%d\t", where[d]);

            if (gggExtendedDebug)
            {
                fprintf(fp,"\t(");

                for (d=1; d<=gggExtendedDebugNumOfWeightsAsParams; d++)
                {
                    fprintf(fp,"%g", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[where[d]]);

                    if (d!=ndim)
                        fprintf(fp,"\t");
                }

                for (; d<=ndim; d++)
                {
                    fprintf(fp,"%g", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[where[d]]);

                    if (d!=ndim)
                        fprintf(fp,"\t");
                }

                fprintf(fp,")\t");
            }

            fprintf(fp,"y: %g\n\n",y);

            fprintf(fp,"Performing exhaustive search ==>\n");
        }
    }

    for (i=xmin[1],nit=0; i<=xmax[1]; i++)
    {
        for (j=xmin[2]; j<=xmax[2]; j++)
        {
            for (k=xmin[3]; k<=xmax[3]; k++)
            {
                for (ell=xmin[4]; ell<=xmax[4]; ell++,nit++)
                {
                    x[1] = i;
                    x[2] = j;
                    x[3] = k;
                    x[4] = ell;

                    y = funk(x,id,container,(void*)(&report));

                    if (report.fast)
                        (*nfast)++;
                    else
                        (*nslow)++;

                    if (distn)
                    {
                        distn->statesandvals[nit].ndim = ndim;
                        distn->statesandvals[nit].val  = y;

                        for (m=1; m<=ndim; m++)
                            distn->statesandvals[nit].state[m] = x[m];
                    }

                    if (y<yb)
                    {
                        yb = y;

                        for (m=1; m<=ndim; m++)
                            xb[m] = x[m];
                    }
                    
                    if (DebugTrace)
                    {
                        for (d=1; d<=ndim; d++)
                            fprintf(fp,"%d ",x[d]);

                        if (gggExtendedDebug)
                        {
                            fprintf(fp,"\t(");

                            for (d=1; d<=gggExtendedDebugNumOfWeightsAsParams; d++)
                            {
                                fprintf(fp,"%g", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[x[d]]);

                                if (d!=ndim)
                                    fprintf(fp,"\t");
                            }

                            for (; d<=ndim; d++)
                            {
                                fprintf(fp,"%g", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[x[d]]);

                                if (d!=ndim)
                                    fprintf(fp,"\t");
                            }

                            fprintf(fp,")\t");
                        }

                        fprintf(fp,"y: %g\t",y);

                        if (report.fast)
                            fprintf(fp,"fast\n");
                        else
                            fprintf(fp,"slow\n");
                    }
                }
            }
        }
    }
/////////////////////////
#if 0
    /* Start at the lowest desired parameter values... */
    for (k=1; k<=ndim; k++)
        x[k] = xmin[k];

    p1 = ndim;

    y = yb = *value;

    terminate = 0;
    nit       = 0;

    while (!terminate)
    {
        y = funk(x,id,container,(void*)(&report));

        if (report.fast)
            (*nfast)++;
        else
            (*nslow)++;

        if (distn)
        {
            distn->statesandvals[nit].ndim = ndim;
            distn->statesandvals[nit].val  = y;

            for (k=1; k<=ndim; k++)
                distn->statesandvals[nit].state[k] = x[k];
        }

        if (y<yb)
        {
            yb = y;

            for (k=1; k<=ndim; k++)
                xb[k] = x[k];
        }
        
        if (DebugTrace)
        {
            for (d=1; d<=ndim; d++)
                fprintf(fp,"%d ",x[d]);

            if (gggExtendedDebug)
            {
                fprintf(fp,"\t(");

                for (i=1; i<=gggExtendedDebugNumOfWeightsAsParams; i++)
                {
                    fprintf(fp,"%g", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->settings->weight[x[i]]);

                    if (i!=ndim)
                    fprintf(fp,"\t");
                }

                for (; i<=ndim; i++)
                {
                    fprintf(fp,"%g", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->settings->tau[x[i]]);

                    if (i!=ndim)
                    fprintf(fp,"\t");
                }

                fprintf(fp,")\t");
            }

            fprintf(fp,"y: %g\n",y);
        }

        if (x[ndim]<xmax[ndim])
        {
            x[ndim]++;
        }
        else
        {
            for (k=ndim; k>=p1; k--)
            {
                if (x[k]<xmax[k])
                {
                    x[k]++;

                    for (i=ndim; i>k; i--)
                        x[i] = xmin[i];
                    
                    break;
                }
            }

            if (k==p1-1)
            {
                p1--;

                if (p1>=1)
                {
                    x[p1]++;

                    for (i=ndim; i>p1; i--)
                        x[i] = xmin[i];

                    terminate = 0;
                }
                else
                {
                    terminate = 1;
                }   
            }        
        }

        nit++;
    }
#endif

    for (k=1; k<=ndim; k++) /* Start at the lowest desired parameter values... */
        where[k] = xb[k];

    *value = yb;

    if (DebugTrace)
    {
        fprintf(fp,"\n\n\n=============== Overall Results =================\n");
        fprintf(fp,"Minimum (best) location ==>");

        for (i=1; i<=ndim; i++)
            fprintf(fp,"%d ", where[i]);

        if (gggExtendedDebug)
        {
            fprintf(fp,"\t(");

            for (i=1; i<=gggExtendedDebugNumOfWeightsAsParams; i++)
            {
                fprintf(fp,"%g", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[where[i]]);

                if (i!=ndim)
                    fprintf(fp,"\t");
            }

            for (; i<=ndim; i++)
            {
                fprintf(fp,"%g", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[where[i]]);

                if (i!=ndim)
                    fprintf(fp,"\t");
            }

            fprintf(fp,")\t");
        }

        fprintf(fp,"\nValue at minimum: %g\n", *value);
        fprintf(fp,"=================================================");
        fprintf(fp,"\nTotal iterations: %d\n", nit);
        fprintf(fp,"=================================================");
        fclose(fp);
    }

    return (0);
}




int bayes_BiExpDiscreteSpaceSearchForPreferableNeighbourState(double (*funk)(int *, int, void *, void *),
                                                              int      id,
                                                              void    *container,
                                                              int      ndim,
                                                              int     *where,
                                                              double  *value,
                                                              void    *config,
                                                              int     *nfast,
                                                              int     *nslow)
{
    int    i, j, x[5], xb[5];
    double y, yb;

    RapidBiExpSearchReporting_t report;

    yb = *value;

    /* Positive changes... */
    for (i=1; i<=ndim; i++)
    {
        for (j=1; j<=ndim; j++)
            x[j] = where[j];

        x[i] += 1;
        y     = funk(x,id,container,(void*)(&report));

        if (report.fast)
            (*nfast)++;
        else
            (*nslow)++;

        if (y<yb)
        {
            for (j=1; j<=ndim; j++)
                xb[j] = x[j];

            yb = y;
        }
    }

    /* Negative changes (followed by the 'minimum' at centre)... */
    for (i=1; i<=ndim; i++)
    {
        for (j=1; j<=ndim; j++)
            x[j] = where[j];

        x[i] -= 1;
        y     = funk(x,id,container,(void*)(&report));

        if (report.fast)
            (*nfast)++;
        else
            (*nslow)++;

        if (y<yb)
        {
            for (j=1; j<=ndim; j++)
                xb[j] = x[j];

            yb = y;
        }
    }
    
    if (yb<*value)
    {
        for (j=1; j<=ndim; j++)
            where[j] = xb[j];

        *value = yb;

        return (1);
    }
    else
    {
        return (0);
    }
}


int bayes_BiExpDiscreteSpaceRapidGenerateRandomValidState(int  *x,
                                                          void *container)
{
    int    nweights, ntaus, temp;
    double *weights, *taus, bg, bgmin, bgmax;

    RapidBiExpMinusLogProbParams_t *params;

    params = (RapidBiExpMinusLogProbParams_t *)(container);

    nweights = params->rapidparamsandlikelihoods->biexpvaluestore->settings->nweights;
    weights  = params->rapidparamsandlikelihoods->biexpvaluestore->settings->weight;
    ntaus    = params->rapidparamsandlikelihoods->biexpvaluestore->settings->ntaus;
    taus     = params->rapidparamsandlikelihoods->biexpvaluestore->settings->tau;
    bgmin    = params->rapidparamsandlikelihoods->biexpvaluestore->settings->backgroundmin;
    bgmax    = params->rapidparamsandlikelihoods->biexpvaluestore->settings->backgroundmax;

    bg = bgmin+rand_RandomDouble()*(bgmax-bgmin);

	x[1] = bayes_MapWeightValueToClosestRapidGridPoint(rand_RandomDouble()*(1.0-bg),nweights,weights); // w1
    x[2] = bayes_MapWeightValueToClosestRapidGridPoint((1.0-bg-weights[x[1]]),nweights,weights);          // w2
    x[3] = (int)(rand_RandomFloat()*(float)ntaus);                                                    // tau1
    x[4] = (int)(rand_RandomFloat()*(float)ntaus);                                                    // tau2

    if (x[3]<x[4]) //decreasing tau with increasing index
    {
        temp = x[4];
        x[4] = x[3];
        x[3] = temp;
    }

    return (0);
}


int bayes_BiExpDiscreteSpaceRapidGenerateRandomValidStateOutsidePreComputedGrid(int  *x,
                                                                                void *container)
{
    int    nweights, ntaus, temp, *low, *high;
    double *weights, *taus, bg, bgmin, bgmax;

    RapidBiExpMinusLogProbParams_t *params;

    params = (RapidBiExpMinusLogProbParams_t *)(container);

    nweights = params->rapidparamsandlikelihoods->biexpvaluestore->settings->nweights;
    weights  = params->rapidparamsandlikelihoods->biexpvaluestore->settings->weight;
    ntaus    = params->rapidparamsandlikelihoods->biexpvaluestore->settings->ntaus;
    taus     = params->rapidparamsandlikelihoods->biexpvaluestore->settings->tau;
    bgmin    = params->rapidparamsandlikelihoods->biexpvaluestore->settings->backgroundmin;
    bgmax    = params->rapidparamsandlikelihoods->biexpvaluestore->settings->backgroundmax;
    low      = params->rapidparamsandlikelihoods->biexpvaluestore->low;
    high     = params->rapidparamsandlikelihoods->biexpvaluestore->high;

	bg = bgmin+rand_RandomDouble()*(bgmax-bgmin); // w0

	if (rand_RandomDouble()<0.5)
    {
    }
    else
    {
    }

    if (rand_RandomDouble()<((weights[low[1]]-weights[0])/(weights[nweights-1]-weights[high[1]])))
        x[1] = bayes_MapWeightValueToClosestRapidGridPoint(rand_RandomDouble()*(weights[low[1]]-bg),nweights,weights); // w1
    else
        x[1] = bayes_MapWeightValueToClosestRapidGridPoint(rand_RandomDouble()*(1.0-weights[high[1]]-bg),nweights,weights); // w2

    x[2] = bayes_MapWeightValueToClosestRapidGridPoint((1.0-bg-weights[x[1]]),nweights,weights);          // w2
    
    if (rand_RandomDouble()<((taus[low[3]]-taus[0])/(taus[ntaus-1]-taus[high[3]]))) // tau1
        x[3] = bayes_MapLifetimeValueToClosestRapidGridPoint(rand_RandomDouble()*(taus[low[3]]-taus[0]),ntaus,taus);
    else
        x[3] = bayes_MapLifetimeValueToClosestRapidGridPoint(rand_RandomDouble()*(taus[ntaus-1]-taus[high[3]]),ntaus,taus);

    if (rand_RandomDouble()<((taus[low[4]]-taus[0])/(taus[ntaus-1]-taus[high[4]]))) // tau2
        x[4] = bayes_MapLifetimeValueToClosestRapidGridPoint(rand_RandomDouble()*(taus[low[4]]-taus[0]),ntaus,taus);
    else
        x[4] = bayes_MapLifetimeValueToClosestRapidGridPoint(rand_RandomDouble()*(taus[ntaus-1]-taus[high[4]]),ntaus,taus);

    if (x[3]<x[4]) //decreasing tau with increasing index
    {
        temp = x[4];
        x[4] = x[3];
        x[3] = temp;
    }

    return (0);
}


int bayes_BiExpDiscreteSpaceMinimisationStochasticSearch(double (*funk)(int *, int, void *, void *),
                                                         int      id,
                                                         void    *container,
                                                         int      ndim,
                                                         int     *where,
                                                         double  *value,
                                                         void    *config,
                                                         int     *nfast,
                                                         int     *nslow)
{
    int    ret = MATH_MINIMISATION_RESULT_MAX_FCT_CALLS_REACHED;
	
    int    i, j, nit;
    double y, ynew, yb, ybb, t, dt;
    int    x[5], xnew[5], xb[5], xbb[5], *grid_extents;
    int    initialisation, ninitialisations;

    RapidBiExpSearchReporting_t report;

    int    DebugTrace=0, d;
    FILE   *fp = NULL;

    grid_extents  = ((MultiExpDiscreteGridSearchConfigParams_t *)config)->gridextents;
	
    if ((ndim<=0) || (!where) || (!funk) || (!grid_extents))
        return (MATH_MINIMISATION_RESULT_ERROR_INVALID_INPUT);

	if (DebugTrace)
    {
        if ((fp = fopen("BiExpDiscreteSpaceMinimisationStochasticSearchOutput.txt", "w")))
        {
            DebugTrace = 0;
        }
        else
        {
            fprintf(fp,"\nbayes_BiExpDiscreteSpaceMinimisationStochasticSearch ==>\n");

            fprintf(fp,"ndim: %d\n", ndim);

            if (gggExtendedDebug)
            {
                gggExtendedDebugParamsContainer  = (RapidBiExpMinusLogProbParams_t *)(container);
                gggExtendedDebugNumOfDecays      = gggExtendedDebugParamsContainer->paramfixing->nparams/2;

                for (i=1, gggExtendedDebugNumOfTausAsParams=0; i<=gggExtendedDebugNumOfDecays; i++)
                {
                    if (gggExtendedDebugParamsContainer->paramfixing->tauuserfixed[i] != BAYES_PARAM_VALUE_USER_FIXED)
                        gggExtendedDebugNumOfTausAsParams++;
                }

                gggExtendedDebugNumOfWeightsAsParams = gggExtendedDebugNumOfDecays-(gggExtendedDebugParamsContainer->paramfixing->nparamsuserfixed-(gggExtendedDebugNumOfDecays-gggExtendedDebugNumOfTausAsParams));

                fprintf(fp,"\nExtended debugging trace output ==>\n");
                fprintf(fp,"gExtendedDebugNumOfWeightsAsParams: %d, gExtendedDebugNumOfTausAsParams: %d\n", gggExtendedDebugNumOfWeightsAsParams, gggExtendedDebugNumOfTausAsParams);
            }

            fprintf(fp,"\nStarting location ==> \n");

            for (i=1; i<=ndim; i++)
                fprintf(fp,"%d\t", where[i]);

            if (gggExtendedDebug)
            {
                fprintf(fp,"\t(");

                for (i=1; i<=gggExtendedDebugNumOfWeightsAsParams; i++)
                    fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[where[i]]);

                for (; i<=ndim; i++)
                    fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[where[i]]);

                fprintf(fp,")");
            }

            fprintf(fp,"\ty: %g\n\n",funk(where,id,container,(void*)(&report)));
        }
    }

    /* Initialize routine...                               */	
    /* Starting location as provided by calling routine... */
    for (j=1; j<=ndim; j++)
    {
        xbb[j] = xb[j] = x[j] = where[j];
    }

    ybb = yb = y = funk(x,id,container,(void*)(&report));

    nit = 0;
    rand_InitializeRandomSeed();

    /* Iterations show grow with increasing 'ndim'... */
    for (j=1,ninitialisations=2; j<=ndim; j++)
        ninitialisations *= 2;

    ninitialisations *= bayes_BiExpConfigGetRapidGridSearchInitialisationsFactor();

    t  = 15.0;
    dt = t/(double)ninitialisations;

    if (DebugTrace)
    {
        fprintf(fp,"ninitialisations: %d\n",ninitialisations);
        fprintf(fp,"t: %g\n",t);
        fprintf(fp,"dt: %g\n\n",dt);
    }

    for (initialisation=0; initialisation<ninitialisations; initialisation++)
    {
        t -= dt;

        bayes_BiExpDiscreteSpaceRapidGenerateRandomValidState(x,container);

        for (j=1; j<=ndim; j++)
            xnew[j] = x[j];

        ynew = funk(xnew,id,container,(void*)(&report));

        if (DebugTrace)
        {
            fprintf(fp,"\nRandom starting point (%d, t: %g, ybb: %g) ==>\n", initialisation, t, ybb);

            for (i=1; i<=ndim; i++)
                fprintf(fp,"%d\t", xnew[i]);

            if (gggExtendedDebug)
            {
                fprintf(fp,"\t(");

                for (i=1; i<=gggExtendedDebugNumOfWeightsAsParams; i++)
                    fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[xnew[i]]);

                for (; i<=ndim; i++)
                    fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[xnew[i]]);

                fprintf(fp,")");
            }

            fprintf(fp,"\ty: %g\n", ynew);
        }

        if (report.fast)
            (*nfast)++;
        else
            (*nslow)++;

        if ((ynew<=ybb) ||
            (ynew-15.0<ybb))
//            ((ynew>ybb) && (exp(-(ynew-ybb)/t)>rand_NrModifiedRan1(0))))
        {
            /* While one exists...                        */
            /* ...search for the optimal neighbour state. */
            while (bayes_BiExpDiscreteSpaceSearchForPreferableNeighbourState(funk,id,container,ndim,xnew,&ynew,config,nfast,nslow))
            {
                /* Save the neighbouring state... */
                for (j=1; j<=ndim; j++)
                    x[j] = xnew[j];

                y = ynew;

                if (DebugTrace)
                {
                    fprintf(fp,"Neighbour (better) state saved...\n");                

                    for (i=1; i<=ndim; i++)
                        fprintf(fp,"%d\t", xnew[i]);

                    if (gggExtendedDebug)
                    {
                        fprintf(fp,"\t(");

                        for (i=1; i<=gggExtendedDebugNumOfWeightsAsParams; i++)
                            fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[xnew[i]]);

                        for (; i<=ndim; i++)
                            fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[xnew[i]]);

                        fprintf(fp,")\t");
                    }

                    fprintf(fp,"ynew: %g\n", ynew);
                }

                if (y<yb)
                {
                    /* New best point found in this inner loop... */
                    for (j=1; j<=ndim; j++)
                    {
                        xb[j] = x[j];
                    }

                    yb = y;

                    if (DebugTrace)
                    {
                        fprintf(fp,"New (inner) best state saved...\n");
                    }
                }

                nit++;
            }
        }
        else
        {
            if (DebugTrace)
            {
                fprintf(fp,"Ignored\n");
            }
        }

        if (yb < ybb)
        {
            for (j=1; j<=ndim; j++)
                xbb[j] = xb[j];

            ybb = yb;

            if (DebugTrace)
            {
                fprintf(fp,"New (outer) best state saved\n");

                for (d=1; d<=ndim; d++)
                    fprintf(fp,"%d ", xbb[d]);

                if (gggExtendedDebug)
                {
                    fprintf(fp,"\t(");

                    for (i=1; i<=gggExtendedDebugNumOfWeightsAsParams; i++)
                        fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[xbb[i]]);

                    for (; i<=ndim; i++)
                        fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[xbb[i]]);

                    fprintf(fp,")\t");
                }

                fprintf(fp,"ybb: %g\n", ybb);        
            }
        }
    }

    *value = ybb; 

    for (j=1; j<=ndim; j++)
        where[j] = xbb[j];
    
	if (DebugTrace)
    {
        fprintf(fp,"\n\n\n=============== Overall Results =================\n");
        fprintf(fp,"Minimum (best) location ==>");

        for (i=1; i<=ndim; i++)
            fprintf(fp,"%d ", where[i]);

        if (gggExtendedDebug)
        {
            fprintf(fp,"\t(");

            for (i=1; i<=gggExtendedDebugNumOfWeightsAsParams; i++)
                fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->weight[where[i]]);

            for (; i<=ndim; i++)
                fprintf(fp,"%g\t", gggExtendedDebugParamsContainer->rapidparamsandlikelihoods->biexpvaluestore->settings->tau[where[i]]);

            fprintf(fp,")\t");
        }

        fprintf(fp,"\nValue at minimum: %g\n", *value);

        fprintf(fp,"=================================================");
        fclose(fp);
    }

    return (ret);
}





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
                                           BayesProbDistn_t             *distribution)
{
    double value, valuebest, temp, dweight, dtau;
    int    x[5], xbest[5], xmin[5], xmax[5], id=0, i, j, k, ndim, modelfixedweightindex, extents[5];
    int    ret, nfast=0, nslow=0;

    RapidBiExpMinusLogProbParams_t            container;
    MultiExpDiscreteGridSearchConfigParams_t  config;
    BayesProbDistn_t                         *distn=NULL;

    ndim     = 2*ndecays;
    ndim    -= paramfixing->nparamsuserfixed;

    container.data                      = data;
    container.nbins                     = nbins;
    container.fitstart					= fitstart;
    container.binwalls                  = binwalls;
    container.instr                     = instr;
    container.interval                  = interval;
    container.modulationperiod          = modulationperiod;
    container.rapidparamsandlikelihoods = grid;
    container.hyperparam                = alpha;
    container.paramfixing               = paramfixing;

    /* Initialise input starting location 'x' first... */
    for (k=1,i=1; k<=ndecays; k++)
    {
        if (paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            x[i]       = bayes_MapWeightValueToClosestRapidGridPoint(weights_mp[k],
                                                                     grid->biexpvaluestore->settings->nweights,
                                                                     grid->biexpvaluestore->settings->weight);

            extents[i] = grid->biexpvaluestore->settings->nweights;
            i++;
        }
    }

    for (k=1; k<=ndecays; k++)
    {
        if (paramfixing->tauuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            x[i]       = bayes_MapLifetimeValueToClosestRapidGridPoint(taus_mp[k],
                                                                       grid->biexpvaluestore->settings->ntaus,
                                                                       grid->biexpvaluestore->settings->tau);

            extents[i] = grid->biexpvaluestore->settings->ntaus;
            i++;
        }
    }

    config.gridextents = extents;
    config.gridmins    = grid->biexpvaluestore->low;
    config.gridmaxs    = grid->biexpvaluestore->high;

    value = BAYES_SIZE_DOUBLE_HUGE;

    /* Exhaustive search of the pre-computed region initially, should be fast... */
    config.gridmins = grid->biexpvaluestore->low;
    config.gridmaxs = grid->biexpvaluestore->high;
    config.distn    = NULL/*distn*/; //don't want the exhaustive search to be slow, may need attention if called from bpda

    ret = bayes_BiExpDiscreteSpaceMinimisationExhaustiveSearch(bayes_BiExpRapidMinusLogProbGivenParams,
                                                               id,(void*)(&container),
                                                               ndim,x,&value,(void*)(&config),&nfast,&nslow);

    /* Search space using stochastic algorithm... */
    ret = bayes_BiExpDiscreteSpaceMinimisationStochasticSearch(bayes_BiExpRapidMinusLogProbGivenParams,
                                                               id,(void*)(&container),
                                                               ndim,x,&value,(void*)(&config),&nfast,&nslow);

    
	if (ret == MATH_MINIMISATION_RESULT_USERCANCEL)
		goto Cancel;


    if ((ret >= MATH_MINIMISATION_RESULT_SUCCESS) && 
        (value < BAYES_SIZE_DOUBLE_HUGE) &&
        (bayes_BiExpConfigRapidGetGridSearchLocalisedExhaustiveSearchDelta()))
    {
        /* Record current best state and look for a better one... */
        valuebest = value;
        
        for (i=1; i<=ndim; i++)
            xbest[i] = x[i];

        /* Setup the limits of the region around 'xbest' to be searched exhaustively... */
        for (i=1; i<=ndim; i++)
        {
	        xmin[i] = x[i]-bayes_BiExpConfigRapidGetGridSearchLocalisedExhaustiveSearchDelta();
            xmax[i] = x[i]+bayes_BiExpConfigRapidGetGridSearchLocalisedExhaustiveSearchDelta();
        }

        /* Constrain the limits of the region around 'xbest' to defined grid points... */
        for (i=1; i<=ndecays-paramfixing->nweightsuserfixed; i++)
        {
            xmin[i] = (xmin[i]<0)?(0):(xmin[i]);
            xmax[i] = (xmax[i]>grid->biexpvaluestore->settings->nweights-1)?(grid->biexpvaluestore->settings->nweights-1):(xmax[i]);
        }

        for (; i<=ndim; i++)
        {
            xmin[i] = (xmin[i]<0)?(0):(xmin[i]);
            xmax[i] = (xmax[i]>grid->biexpvaluestore->settings->ntaus-1)?(grid->biexpvaluestore->settings->ntaus-1):(xmax[i]);
        }

        if (distribution)
        {
            distn = distribution; /* Container has been included by calling function for use later (e.g. BPDA)... */
        }
        else
        {
            /* Only allocate if parameter average/error values are required... */
            if (((weights_ave) && (weights_err)) || ((taus_ave) && (taus_err)))
            {
                distn   = (BayesProbDistn_t*)malloc(sizeof(BayesProbDistn_t));
                dweight = (grid->biexpvaluestore->settings->weight[grid->biexpvaluestore->settings->nweights-1]-grid->biexpvaluestore->settings->weight[0])/(double)grid->biexpvaluestore->settings->nweights;
                dtau    = (grid->biexpvaluestore->settings->tau[grid->biexpvaluestore->settings->ntaus-1]-grid->biexpvaluestore->settings->tau[0])/(double)grid->biexpvaluestore->settings->ntaus;

                ret = bayes_AllocateForMultiExpDiscreteProbDistn(distn,ndim,
                                                                 ndecays-paramfixing->nweightsuserfixed,
                                                                 ndecays-paramfixing->ntaususerfixed,
                                                                 dweight,dtau,xmin,xmax);

                if (ret<0)
                {
                    *val = BAYES_SIZE_DOUBLE_HUGE;

                    //free_Bayes_ivector(x,1,ndim);
                    //free_Bayes_ivector(xmin,1,ndim);
                    //free_Bayes_ivector(xmax,1,ndim);
                    //free_Bayes_ivector(extents,1,ndim);
					free (distn);

                    return (BAYES_ERR_PARAM_ESTIMATION_FAILURE);
                }            
            }
        }

        /* Exhaustive search of region around the claimed minimum for a better solution... */
        config.gridmins = xmin;
        config.gridmaxs = xmax;
        config.distn    = distn;

        ret = bayes_BiExpDiscreteSpaceMinimisationExhaustiveSearch(bayes_BiExpRapidMinusLogProbGivenParams,
                                                                   id,(void*)(&container),
                                                                   ndim,x,&value,(void*)(&config),&nfast,&nslow);

        if ((ret>=0) && (value<valuebest))
        {
            for (i=1; i<=ndim; i++) //if at edge of exhaustive search region allow one more greedy better neighbour search...
            {
                if ((x[i]==config.gridmins[i]) || (x[i]==config.gridmaxs[i]))
                {
                    /* While one exists...                        */
                    /* ...search for the optimal neighbour state. */
                    while (bayes_BiExpDiscreteSpaceSearchForPreferableNeighbourState(bayes_BiExpRapidMinusLogProbGivenParams,
                                                                                     id,(void*)(&container),
                                                                                     ndim,x,&value,(void*)(&config),&nfast,&nslow))
                    {
                        /* Save the neighbouring state... */
                        for (j=1; j<=ndim; j++)
                            xbest[j] = x[j];

                        valuebest = value;
                    }
                }
            }

            valuebest = value;
            
            for (i=1; i<=ndim; i++)
	            xbest[i] = x[i];
        }

        /* Put our best results back into the free parameter array... */
        value = valuebest;
        
        for (i=1; i<=ndim; i++)
            x[i] = xbest[i];
    }

    /* Check that the optimal values are sensible, i.e. at least exist on the grid... */
    for (i=1; i<=ndecays-paramfixing->nweightsuserfixed; i++)
    {
        if ((x[i]<0) && (x[i]>=grid->biexpvaluestore->settings->nweights))
            ret = BAYES_ERR_PARAM_ESTIMATION_FAILURE;
    }

    for (; i<=ndim; i++)
    {
        if ((x[i]<0) && (x[i]>=grid->biexpvaluestore->settings->ntaus))
            ret = BAYES_ERR_PARAM_ESTIMATION_FAILURE;
    }

    if (ret==BAYES_ERR_PARAM_ESTIMATION_FAILURE)
    {
        *val = BAYES_SIZE_DOUBLE_HUGE;

        if ((!distribution) && (distn))
            bayes_FreeForMultiExpDiscreteProbDistn(distn);

        return (BAYES_ERR_PARAM_ESTIMATION_FAILURE);
    }

    /* Compute marginal distributions in order to return 'continuous' space parameter estimates... */
    if (distn)
    {
        bayes_NormaliseMultiExpDiscreteProbDistn(distn,value);
        bayes_DetermineMarginalsForMultiExpDiscreteProbDistn(distn);
    }

    /* Extract the (optimal) weight and lifetimes from 'x' for return to user... */
    for (k=1,i=1; k<=ndecays; k++)
    {
        if (paramfixing->weightuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            if ((distn) && (weights_ave) && (weights_err))
            {
                bayes_ComputeParamAveAndErrUsingMultiExpDiscreteProbDistnMarginal(
                    distn->marginals[i].marginal,
                    grid->biexpvaluestore->settings->weight,
                    dweight,
                    distn->marginals[i].indexlow,
                    distn->marginals[i].indexhigh,
                    &weights_ave[k],
                    &weights_err[k]);
            }

            weights_mp[k] = grid->biexpvaluestore->settings->weight[x[i]];
            i++;
        }
        else
        {
            weights_mp[k]  = paramfixing->weights[k];
            weights_ave[k] = paramfixing->weights[k];
            weights_err[k] = 0.0;
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

    /* Most probable model fixed weight value... */
    for (k=0,temp=1.0; k<=ndecays; k++)
    {
        if (k != modelfixedweightindex)
            temp -= weights_mp[k];
    }

    weights_mp[modelfixedweightindex] = temp;

    /* Average model fixed weight value (determined from averages of the remaining weights)... */
    if (weights_ave)
    {
        for (k=0,temp=1.0; k<=ndecays; k++)
        {
            if (k != modelfixedweightindex)
                temp -= weights_ave[k];
        }

        weights_ave[modelfixedweightindex] = temp;
    }

    /* Now extract and validate the lifetime indices from 'x'... */
    for (k=1; k<=ndecays; k++)
    {
        if (paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            if ((distn) && (taus_ave) && (taus_err))
            {
                bayes_ComputeParamAveAndErrUsingMultiExpDiscreteProbDistnMarginal(
                    distn->marginals[i].marginal,
                    grid->biexpvaluestore->settings->tau,
                    dtau,
                    distn->marginals[i].indexlow,
                    distn->marginals[i].indexhigh,
                    &taus_ave[k],
                    &taus_err[k]);
            }

            taus_mp[k] = grid->biexpvaluestore->settings->tau[x[i]];
            i++;
        }
        else
        {
            taus_mp[k]  = paramfixing->taus[k];
            taus_ave[k] = paramfixing->taus[k];
            taus_err[k] = 0.0;
        }
    }
 
    *val = value;

Cancel:

	if (!distribution && distn)
        bayes_FreeForMultiExpDiscreteProbDistn(distn);

    return (BAYES_ERR_NO_ERROR);
}







////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
#if 1 //hyperparam optimization


double bayes_RapidBiExpMinusLogProbDataLikelihood(int *x, int id, void *container)
{
    int     *data, bin, nbins, fitstart, nphotonsbin;
    double   value, *likelihoods;

    RapidBiExpMinusLogProbParams_t *params;

    params      = (RapidBiExpMinusLogProbParams_t*)(container);
    data        = params->data;
    nbins       = params->nbins;
    fitstart    = params->fitstart;

    if (!params->rapidparamsandlikelihoods->biexpvaluestore->likelihoods[x[1]][x[2]][x[3]][x[4]])
        return (BAYES_SIZE_DOUBLE_HUGE);

    likelihoods = params->rapidparamsandlikelihoods->biexpvaluestore->likelihoods[x[1]][x[2]][x[3]][x[4]]->logphotonlikelihoodgiventausandweights;

    for (bin=fitstart,value=0.0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            value -= (double)(nphotonsbin) * likelihoods[bin];
        }
    }

    return (value);
}


int bayes_RapidBiExpPopulateDataLikelihoodGrid(double                     ****datalikelihoods,
                                               int                           *xmin,
                                               int                           *xmax,
                                               double                        *weights,
                                               double                        *taus,
                                               double                         w0low,
                                               double                         w0high,
                                               int                           *data,
                                               int                            nbins,
                                               int                            fitstart,
                                               int                            nphotons,
                                               double                        *binwalls,
                                               BayesInstrRsp_t               *instr,
                                               float                          interval,
                                               float                          modulationperiod,
                                               BayesRapidValueStore_t        *grid)
{
    double value, sum, w0;
    int    x[5], i, j, k, ell, id=0;
    RapidBiExpMinusLogProbParams_t container;

    container.data                      = data;
    container.nbins                     = nbins;
    container.fitstart                  = fitstart;
    container.rapidparamsandlikelihoods = grid;

    for (i=xmin[1]; i<=xmax[1]; i++)
    {
        for (j=xmin[2]; j<=xmax[2]; j++)
        {
            sum = weights[i]+weights[j];
            w0  = 1.0-sum;

            if ((sum>=0.0) && (sum<=1.0) && (w0>=w0low) && (w0<=w0high))
            {
                for (k=xmin[3]; k<=xmax[3]; k++)
                {
                    for (ell=xmin[4]; ell<=xmax[4]; ell++)
                    {
                        if (taus[k]>taus[ell])
                        {
                            x[1] = i;
                            x[2] = j;
                            x[3] = k;
                            x[4] = ell;

                            value                         = bayes_RapidBiExpMinusLogProbDataLikelihood(x,id,(void*)(&container));;
                            datalikelihoods[i][j][k][ell] = value;                        
                        }
                    }
                }
            }

        }
    }

    return (0);
}


int bayes_RapidBiExpMostProbWeightsAndTausPreComputedDataLikelihood(int                            *xout,
                                                                    int                            *xmin,
                                                                    int                            *xmax,
                                                                    double                         *weights,
                                                                    double                         *taus,
                                                                    double                          w0low,
                                                                    double                          w0high,
                                                                    float                          *val,
                                                                    float                           alpha,
                                                                    BayesRapidBiExpValueStore_t    *grid,
                                                                    double                      ****datalikelihoods)
{
    double value, valuebest, sum, w0;
    int    i, j, k, ell, d, id=0, x[5], xbest[5];

    for (i=xmin[1],valuebest=BIG; i<=xmax[1]; i++)
    {
        for (j=xmin[2]; j<=xmax[2]; j++)
        {
            sum = weights[i]+weights[j];
            w0  = 1.0-sum;

            if ((sum>=0.0) && (sum<=1.0) && (w0>=w0low) && (w0<=w0high))
            {
                for (k=xmin[3]; k<=xmax[3]; k++)
                {
                    for (ell=xmin[4]; ell<=xmax[4]; ell++)
                    {
                        if (taus[k]>taus[ell])
                        {
                            x[1] = i;
                            x[2] = j;
                            x[3] = k;
                            x[4] = ell;

                            if (!grid->likelihoods[i][j][k][ell])
                            {
                                value  = BIG;
                            }
                            else
                            {
                                value  = alpha*(grid->likelihoods[i][j][k][ell]->taus[1]+grid->likelihoods[i][j][k][ell]->taus[2]);
                                value += datalikelihoods[i][j][k][ell];
                            }

                            if (value<valuebest)
                            {
                                valuebest = value;

                                for (d=1; d<=4; d++)
                                    xbest[d] = x[d];
                            }                        
                        }
                    }
                }
            }
        }
    }

    for (d=1; d<=4; d++)
        xout[d] = xbest[d];

    *val = (float)valuebest;

    return (0);
}



double bayes_RapidBiExpMinusLogProbAlphaTimesModelEvidence(double *x, int id, void *container)
{
    int     ret, *data, nbins, fitstart, nphotons, weightsandtaus[5], d;
    float   minuslogprob;
    double  alpha, alphamin, logevidence, weights[3], taus[3], interval, modperiod, *binwalls, ****datalikelihoods;

    RapidBiExpMinusLogProbParams_t     *paramscontainer;
    BayesRapidValueStore_t     *grid;
    BayesInstrRsp_t               *instr;

    paramscontainer = (RapidBiExpMinusLogProbParams_t *)(container);
    alphamin        = paramscontainer->alphamin;
    alpha           = (float)x[1];

    if (alpha<alphamin)
        return (BIG);
    
    data            = paramscontainer->data;
    nbins           = paramscontainer->nbins;
    fitstart        = paramscontainer->fitstart;
    nphotons        = paramscontainer->nphotons;
    instr           = paramscontainer->instr;
    interval        = paramscontainer->interval;
    modperiod       = paramscontainer->modulationperiod;
    binwalls        = paramscontainer->binwalls;
    grid            = paramscontainer->rapidparamsandlikelihoods;
    datalikelihoods = paramscontainer->datalikelihoods;

    /* Determine most-probable parameter weights and taus given the data and the hyperparameter alpha... */
    ret = bayes_RapidBiExpMostProbWeightsAndTausPreComputedDataLikelihood(weightsandtaus,
                                                                          grid->biexpvaluestore->low,grid->biexpvaluestore->high,
                                                                          grid->biexpvaluestore->settings->weight,grid->biexpvaluestore->settings->tau,
                                                                          grid->biexpvaluestore->settings->backgroundmin,grid->biexpvaluestore->settings->backgroundmax,
                                                                          &minuslogprob,(float)alpha,grid->biexpvaluestore,datalikelihoods);

    if (ret<0)
        return (BIG);

    /* Grab the actual weight and tau values from the optimal grid point... */
    for (d=1; d<=2; d++)
        weights[d] = grid->biexpvaluestore->settings->weight[weightsandtaus[d]];

    for (; d<=4; d++)
        taus[d-2]    = grid->biexpvaluestore->settings->tau[weightsandtaus[d]];

    /* Use the Gaussian approximation to determine the (log of) integral (i.e. log of model evidence)... */
    ret = bayes_DetemineDecayModelEvidence(/*ndecays*/2,weights,taus,NULL,minuslogprob,
                                           nbins,binwalls,data,
                                           interval,modperiod,instr,&logevidence);

    if (ret<0)
        return (BIG);
    else
        return (-(log(alpha)+logevidence));
}



#if 0
double **Bayes_dmatrix(int nrl,int nrh,int ncl,int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) matrices_error("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) matrices_error("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
#endif

#if 1
#if 0
double ****bayes_Allocate4dmatrix(int *low, int *high)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, j;
    int        dim[5];
	double ****m;

	/* allocate pointers to rows */
	m = (double****)malloc(dim[1]*sizeof(double***));

    if (!m){return (NULL);}
	
	m -= low[1];

    m[low[1]] = (double***)malloc((dim[1]*dim[2])*sizeof(double**));
	
    if (!m[low[1]]){return (NULL);}
	m[low[1]] -= low[2];

	for(i=low[1]+1; i<=high[1]; i++)
        m[i] = m[i-1]+dim[2];

    for(i=low[1]; i<=high[1]; i++)
        for (j=low[2]; j<=high[2]; j++)
            m[i][j] = Bayes_dmatrix(low[3],high[3],low[4],high[4]);

	return (m);
}
#endif

double ****bayes_Allocate4dDoubleMatrix(int *low, int *high)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    double ****all_x, ***all_y, **all_r, *all_c, ****result;
    int        x, y, r, dim[5];

    for (x=1; x<=4; x++)
        dim[x] = high[x]-low[x]+1;

    all_x  = malloc(dim[1]*sizeof(double***));
    all_y  = malloc(dim[1]*dim[2]*sizeof(double**));
    all_r  = malloc(dim[1]*dim[2]*dim[3]*sizeof(double*));
    all_c  = malloc(dim[1]*dim[2]*dim[3]*dim[4]*sizeof(double));

    if ((!all_x) || (!all_y) || (!all_r) || (!all_c))
        return (NULL);

    all_x -= low[1];
    all_y -= low[2];
    all_r -= low[3];
    all_c -= low[4];

    result = all_x;    

    for (x=low[1]; x<high[1]; x++,all_y+=dim[2])
    {
        result[x] = all_y;
        
        for (y=low[2]; y<high[2]; y++,all_r+=dim[3])
        {
            result[x][y] = all_r;
        
            for (r=low[3]; r<high[3]; r++,all_c+=dim[4])
            {
                result[x][y][r] = all_c;
            }
        }
    }

    return (result);
}

void bayes_Free4dDoubleMatrix(double ****m, int *low, int *high)
{
    free(m[low[1]][low[2]][low[3]]);
    free(m[low[1]][low[2]]);
    free(m[low[1]]);
    free(m);
}


double ****bayes_AllocateDataLikelihoodsMatrix(int *low, int *high/*int nrl,int nrh,int ncl,int nch*/)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, j;
	double ****m;

    int ndim[5];

    for (i=1; i<=4; i++)
        ndim[i] = high[i]-low[i]+1;

    m  = (double****)malloc(ndim[1]*sizeof(double***));
    m -= low[1]; /* Not zero-based... */

    for (i=low[1]; i<=high[1]; i++)
    {
        m[i] = (double***)malloc(ndim[2]*sizeof(double**));

        if (!m[i]){return (NULL);}

	    m[i] -= low[2];

        for (j=low[2]; j<=high[2]; j++)
        {
            m[i][j] = Bayes_dmatrix(low[3],high[3],low[4],high[4]);

         /*   if (!m[i][j])
                return (NULL);

	        m[i][j] -= low[3];

            for (k=low[3]; k<=high[3]; k++)
            {
                m[i][j][k] = (BayesRapidBiExpDiscreteValues_t**)malloc((size_t)(ndim[4]*sizeof(BayesRapidBiExpDiscreteValues_t*)));

                if (!m[i][j][k])
                    return (NULL);

	            m[i][j][k] -= low[4];
            }*/
        }
    }

	return (m);
}


void bayes_FreeDataLikelihoodsMatrix(double ****m, int *low, int *high/*int nrl,int nrh,int ncl,int nch*/)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, j;
    int ndim[5];

    for (i=1; i<=4; i++)
        ndim[i] = high[i]-low[i]+1;

    for (i=low[1]; i<=high[1]; i++)
    {
        for (j=low[2]; j<=high[2]; j++)
        {
            free_Bayes_dmatrix(m[i][j],low[3],high[3],low[4],high[4]);
        }

        if (m[i]+low[2])
            free(m[i]+low[2]);
    }

    if (m+low[1])
        free(m+low[1]);
}
#endif


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
                                           BayesRapidValueStore_t *grid)
{
    double x[2], val, valbest, xbest[2], ****datalikelihoods;
    int    i, id=0, ret;
    double alpha, dalpha=0.05;

    double tempval[100], tempalpha[100];

    RapidBiExpMinusLogProbParams_t container;

    container.data             = data;
    container.nbins            = nbins;
    container.fitstart         = fitstart;
    container.nphotons         = nphotons;
    container.binwalls         = binwalls;
    container.instr            = instr;
    container.interval         = interval;
    container.modulationperiod = modulationperiod;
    
    container.rapidparamsandlikelihoods = grid;
    container.alphamin                  = alphamin;
    

    datalikelihoods = bayes_AllocateDataLikelihoodsMatrix/*bayes_Allocate4dDoubleMatrix*/(grid->biexpvaluestore->low,grid->biexpvaluestore->high);

    ret = bayes_RapidBiExpPopulateDataLikelihoodGrid(datalikelihoods,grid->biexpvaluestore->low,grid->biexpvaluestore->high,
                                                     grid->biexpvaluestore->settings->weight,grid->biexpvaluestore->settings->tau,
                                                     grid->biexpvaluestore->settings->backgroundmin,grid->biexpvaluestore->settings->backgroundmax,
                                                     data,nbins,fitstart,nphotons,
                                                     binwalls,instr,interval,modulationperiod,
                                                     grid);



    container.datalikelihoods  = datalikelihoods;

    for (i=0,valbest=BIG; i<15; i++)
    {
        alpha = alphamin+dalpha*(float)i;
        x[1] = (float)alpha;

        val = bayes_RapidBiExpMinusLogProbAlphaTimesModelEvidence(x,id,(void*)(&container));

        if (val<valbest)
        {
            valbest  = val;
            xbest[1] = x[1];
        }

        tempalpha[i] = alpha;
        tempval[i]   = val;
    }

    *alphastar = (float)xbest[1];
    *value     = (float)valbest;

    bayes_FreeDataLikelihoodsMatrix(datalikelihoods,grid->biexpvaluestore->low,grid->biexpvaluestore->high);

    //free_Bayes_dmatrix(datalikelihoods,0,grid->nw0s,0,grid->ntaus);

    return (0);
}




#endif

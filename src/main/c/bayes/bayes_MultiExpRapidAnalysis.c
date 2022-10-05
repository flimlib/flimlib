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
#include "bayes_MultiExpRapidAnalysis.h"
#include "bayes_RapidBayesDecayAnalysis.h"



/***********************************************************************************/
/*                                                                                 */
/*                      GLOBALS FOR ENHANCED DEBUG OUTPUT                          */
/*                                                                                 */
/***********************************************************************************/

int                                      ggExtendedDebug=0;
int                                      ggExtendedDebugNumOfDecays;
int                                      ggExtendedDebugNumOfWeightsAsParams;
int                                      ggExtendedDebugNumOfTausAsParams;
PsuedoRapidMultiExpMinusLogProbParams_t *ggExtendedDebugParamsContainer;



/***********************************************************************************/
/*                                                                                 */
/*        MULTI-EXPONENTIAL PSUEDO-RAPID ANALYSIS SPECIFIC MODEL FUNCTIONS         */
/*                                                                                 */
/***********************************************************************************/

double bayes_PsuedoRapidMultiExpMinusLogProbGivenWeightsAndTaus(int *x, int id, void *container)
{
    int    i, k, modelfixedweightindex, ndecays, *weightindexes, *tauindexes;
    int    *data, bin, nbins, fitstart, nphotonsbin;
    double val, temp, hyperparam, *weights, *taus, *binwalls, **decayphotonlikelihoods, interval, norm, backgroundmax;
    int type;

    PsuedoRapidMultiExpMinusLogProbParams_t *params;

    params           = (PsuedoRapidMultiExpMinusLogProbParams_t *)(container);
    ndecays          = params->ndecays;
    weightindexes    = Bayes_ivector(1,ndecays);
    tauindexes       = Bayes_ivector(1,ndecays);

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

            if ((weightindexes[k]<0) || (weightindexes[k]>params->rapidparamsandlikelihoods->likelihoodsvaluestore->settings->nweights))
            {
                /* Grid location does not exist... */
                free_Bayes_ivector(weightindexes,1,ndecays);
                free_Bayes_ivector(tauindexes,1,ndecays);
                return (BAYES_SIZE_DOUBLE_HUGE);
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

            if ((tauindexes[k]<0) || (tauindexes[k]>params->rapidparamsandlikelihoods->likelihoodsvaluestore->settings->ntaus))
            {
                /* Grid location does not exist... */
                free_Bayes_ivector(weightindexes,1,ndecays);
                free_Bayes_ivector(tauindexes,1,ndecays);
                return (BAYES_SIZE_DOUBLE_HUGE);
            }
        }
    }
    
    weights = Bayes_dvector(0,ndecays);
    taus    = Bayes_dvector(1,ndecays);

    /* Populate weights vector with free, user fixed, and model fixed values... */
    for (k=0; k<=ndecays; k++)
    {
        if (k != modelfixedweightindex)
        {
            if (params->paramfixing->weightuserfixed[k] == BAYES_PARAM_VALUE_FREE)
                weights[k] = params->rapidparamsandlikelihoods->likelihoodsvaluestore->settings->weight[weightindexes[k]];
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
            taus[k] = params->rapidparamsandlikelihoods->likelihoodsvaluestore->settings->tau[tauindexes[k]];
        else
            taus[k] = params->paramfixing->taus[k];
    }

    /* Check that weight values are valid, background can be zero but decay */
    /* component weights should be greater than zero really...              */

    backgroundmax = params->rapidparamsandlikelihoods->likelihoodsvaluestore->settings->backgroundmax;

    if ((weights[0]<0.0) || (weights[0]>backgroundmax) || (weights[0]>1.0))
    {
        free_Bayes_dvector(weights,0,ndecays);
        free_Bayes_dvector(taus,1,ndecays);
        free_Bayes_ivector(weightindexes,1,ndecays);
        free_Bayes_ivector(tauindexes,1,ndecays);
        return (BAYES_SIZE_DOUBLE_HUGE);
    }

    for (k=1; k<=ndecays; k++)
    {
        if ((weights[k]<=0.0) || (weights[k]>1.0))
        {
            free_Bayes_dvector(weights,0,ndecays);
            free_Bayes_dvector(taus,1,ndecays);
            free_Bayes_ivector(weightindexes,1,ndecays);
            free_Bayes_ivector(tauindexes,1,ndecays);
            return (BAYES_SIZE_DOUBLE_HUGE);
        }
    }

    /* Check that lifetimes are valid */
    for (k=1; k<=ndecays; k++)
    {
        if (taus[k]<SMALLTIME)
        {
            free_Bayes_dvector(weights,0,ndecays);
            free_Bayes_dvector(taus,1,ndecays);
            free_Bayes_ivector(weightindexes,1,ndecays);
            free_Bayes_ivector(tauindexes,1,ndecays);
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
                free_Bayes_dvector(weights,0,ndecays);
                free_Bayes_dvector(taus,1,ndecays);
                free_Bayes_ivector(weightindexes,1,ndecays);
                free_Bayes_ivector(tauindexes,1,ndecays);
                return (BAYES_SIZE_DOUBLE_HUGE);
            }
        }    
    }

    /* Weight and lifetime values okay, compute the minus log likelihood using pre-computed bin likelihoods... */
    decayphotonlikelihoods = (double **) malloc((1+ndecays)*sizeof(double*));

    for (k=1; k<=ndecays; k++)
    {
        if (params->paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            decayphotonlikelihoods[k] = params->rapidparamsandlikelihoods->likelihoodsvaluestore->
                fluorescencelikelihoods[tauindexes[k]].fluorescencedecayphotonlikelihoodsgiventau;
        }
        else
        {
            decayphotonlikelihoods[k] = params->paramfixing->
                fluorescencelikelihoods[k].fluorescencedecayphotonlikelihoodsgiventau;
        }
    }

    /* Now compute the minus log likelihood... */
    hyperparam = params->hyperparam;
    data       = params->data;
    nbins      = params->nbins;
    fitstart   = params->fitstart;
    binwalls   = params->binwalls;
    interval   = params->interval;

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
            
    free_Bayes_dvector(weights,0,ndecays);
    free_Bayes_dvector(taus,1,ndecays);
    free_Bayes_ivector(weightindexes,1,ndecays);
    free_Bayes_ivector(tauindexes,1,ndecays);
    free(decayphotonlikelihoods);

    return (val);
}

/////////////////// save prob distn... start

int bayes_ComputeParamAveAndErrUsingMultiExpDiscreteProbDistnMarginal(double *probx,
                                                                      double *x,
                                                                      double  dx,
                                                                      int     indexlow,
                                                                      int     indexhigh,
                                                                      double *ave,
                                                                      double *err)
{
    int    i;
    double val, temp;

    if ((!probx) || (!x) || (dx<=0.0) ||
        (indexlow<0) || (indexhigh<=0) || (indexhigh<=indexlow) ||
        (!ave) || (!err))
    {
        return (-1);
    }

    for (i=indexlow,val=0.0; i<=indexhigh; i++)
    {
        val += x[i]*probx[i];
    }

    *ave = val;

    for (i=indexlow,val=0.0; i<=indexhigh; i++)
    {
        temp  = x[i]-*ave;
        val  += temp*temp*probx[i];
    }

    *err = sqrt(val);

    return (0);
}


int bayes_DetermineMarginalsForMultiExpDiscreteProbDistn(BayesProbDistn_t *distn)
{
    int    i, j, ndim, s, nstates, *state;
    double val, *marginal, norm;

    if (!distn)
        return (-1);

    ndim    = distn->ndim;
    nstates = distn->nstates;

    for (s=0; s<nstates; s++)
    {
        state = distn->statesandvals[s].state;
        val   = distn->statesandvals[s].val;

        for (i=1; i<=ndim; i++)
        {
            marginal = distn->marginals[i].marginal;

            for (j=1; j<=ndim; j++)
            {
                if ((state[j]>=distn->marginals[i].indexlow) && (state[j]<=distn->marginals[i].indexhigh))
                    marginal[state[j]] += val;
            }
        }
    }

    /* 'Normalise' the marginal... */
    for (i=1; i<=ndim; i++)
    {
        marginal = distn->marginals[i].marginal;

        for (j=distn->marginals[i].indexlow,norm=0.0; j<=distn->marginals[i].indexhigh; j++)
            norm += marginal[j];

        for (j=distn->marginals[i].indexlow; j<=distn->marginals[i].indexhigh; j++)
            marginal[j] /= norm;
    }

    return (0);
}


int bayes_NormaliseMultiExpDiscreteProbDistn(BayesProbDistn_t *distn,
                                             double            min)
{
    int    i;
    double temp, val, norm;

    if (!distn)
        return (-1);
    
    for (i=0, norm=0.0; i<distn->nstates; i++)
    {
        val   = distn->statesandvals[i].val;
        temp  = exp(min-val);
        norm += temp;

        distn->statesandvals[i].val = temp;
    }

    for (i=1; i<=distn->nweights; i++)
    {
        norm *= distn->dweight;
    }

    for (i=1; i<=distn->ntaus; i++)
    {
        norm *= distn->dtau;
    }

    for (i=0; i<distn->nstates; i++)
    {
        distn->statesandvals[i].val /= norm;
    }

    return (0);
}


int bayes_AllocateForMultiExpDiscreteProbDistn(BayesProbDistn_t *distn,
                                               int               ndim,
                                               int               nweights,
                                               int               ntaus,
                                               double            dweight,
                                               double            dtau,
                                               int              *statesmin,
                                               int              *statesmax)
{
    int nstates, i, j, min, max;

    if ((!distn) || 
        (ndim<=0) || (nweights<=0) || (ntaus<=0) || 
        (dweight<=0.0) || (dtau<=0.0) || 
        (!statesmin) || (!statesmax))
    {
        return (-1);
    }

    for (i=1,nstates=1; i<=ndim; i++)
        nstates *= (1+statesmax[i]-statesmin[i]);

    if (nstates<=0)
    {
        return (-2);
    }

    distn->statesandvals = (BayesProbDistnStateAndVal_t *)malloc(nstates*sizeof(BayesProbDistnStateAndVal_t));

    if (!distn->statesandvals)
    {
        return (-3);
    }

    distn->marginals     = (BayesProbMarginal_t *)malloc((1+ndim)*sizeof(BayesProbMarginal_t));

    if (!distn->marginals)
    {
        return (-4);
    }

    /* Parameters and allocation okay, proceed... */
    distn->ndim      = ndim;
    distn->nweights  = nweights;
    distn->ntaus     = ntaus;
    distn->dweight   = dweight;
    distn->dtau      = dtau;
    distn->statesmin = Bayes_ivector(1,ndim);
    distn->statesmax = Bayes_ivector(1,ndim);
    distn->nstates   = nstates;

    for (i=1; i<=ndim; i++)
    {
        distn->statesmin[i] = statesmin[i];
        distn->statesmax[i] = statesmax[i];
    }

    for (i=0; i<nstates; i++)
    {
        distn->statesandvals[i].state = Bayes_ivector(1,ndim);
    }

    for (i=0; i<=ndim; i++)
    {
        distn->marginals[i].marginal = Bayes_dvector(distn->statesmin[i],distn->statesmax[i]);
    }

    for (i=0; i<=ndim; i++)
    {
        min = distn->statesmin[i];
        max = distn->statesmax[i];

        distn->marginals[i].indexlow  = min;
        distn->marginals[i].indexhigh = max;

        for (j=min; j<=max; j++)
        {
            distn->marginals[i].marginal[j] = 0.0;
        }
    }

    return (0);
}


int bayes_FreeForMultiExpDiscreteProbDistn(BayesProbDistn_t *distn)
{
    int i;

    if ((!distn) || (!distn->statesandvals) || (!distn->marginals))
        return (-1);

    for (i=0; i<distn->nstates; i++)
        free_Bayes_ivector(distn->statesandvals[i].state,1,distn->ndim);

    for (i=0; i<=distn->ndim; i++)
        free_Bayes_dvector(distn->marginals[i].marginal,distn->statesmin[i],distn->statesmax[i]);

    free_Bayes_ivector(distn->statesmin,1,distn->ndim);
    free_Bayes_ivector(distn->statesmax,1,distn->ndim);

    free(distn->statesandvals);
    free(distn->marginals);
    free(distn);
    
    return (0);
}

/////////////////// save prob distn... end

int bayes_MultiExpDiscreteSpaceMinimisationExhaustiveSearch(double (*funk)(int *, int, void *),
                                                            int      id,
                                                            void    *container,
                                                            int      ndim,
                                                            int     *where,
                                                            double  *value,
                                                            void    *config)
{
    int    ret = MATH_MINIMISATION_RESULT_SUCCESS;

    double y, yb;
    int    *x, *xb, *xmin, *xmax, i, k, nit;
    int    p1, terminate;
    BayesProbDistn_t *distn;
    
    xmin  = ((MultiExpDiscreteGridSearchConfigParams_t *)config)->gridmins;
    xmax  = ((MultiExpDiscreteGridSearchConfigParams_t *)config)->gridmaxs;
    distn = ((MultiExpDiscreteGridSearchConfigParams_t *)config)->distn;

    if ((ndim<=0) || (!where) || (!funk) || (!xmin) || (!xmax))
        return (MATH_MINIMISATION_RESULT_ERROR_INVALID_INPUT);

    x  = Bayes_ivector(1,ndim);
    xb = Bayes_ivector(1,ndim);

    /* Use 'where' location as provided by calling routine for the optimal initial values... */
    for (k=1; k<=ndim; k++)
    {
        xb[k] = x[k] = where[k];
    }

    yb = y = funk(x,id,container);

    /* Start at the lowest desired parameter values... */
    for (k=1; k<=ndim; k++)
        x[k] = xmin[k];

    p1 = ndim;

    y = yb = *value;

    terminate = 0;
    nit       = 0;

    while (!terminate)
    {
        y = funk(x,id,container);

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

    for (k=1; k<=ndim; k++) /* Start at the lowest desired parameter values... */
        where[k] = xb[k];

    *value = yb;

    free_Bayes_ivector(x,1,ndim);
    free_Bayes_ivector(xb,1,ndim);

    return (0);
}


int bayes_MultiExpDiscreteSpaceSearchForPreferableNeighbourState(double (*funk)(int *, int, void *),
                                                                 int      id,
                                                                 void    *container,
                                                                 int      ndim,
                                                                 int     *where,
                                                                 double  *value,
                                                                 void    *config)
{
    int    i, j, *x, *xb;
    double y, yb;

    x    = Bayes_ivector(1,ndim);
    xb   = Bayes_ivector(1,ndim);

    yb = *value;

    /* Positive changes... */
    for (i=1; i<=ndim; i++)
    {
        for (j=1; j<=ndim; j++)
            x[j] = where[j];

        x[i] += 1;
        y     = funk(x,id,container);

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
        y     = funk(x,id,container);

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

        free_Bayes_ivector(x,1,ndim);
        free_Bayes_ivector(xb,1,ndim);
        return (1);
    }
    else
    {
        free_Bayes_ivector(x,1,ndim);
        free_Bayes_ivector(xb,1,ndim);
        return (0);
    }
}


/* Custom routine: stochastic in performs 'best neighbour' search from multiple random starting states... */
int bayes_MultiExpDiscreteSpaceMinimisationStochasticSearch(double (*funk)(int *, int, void *),
                                                            int      id,
                                                            void    *container,
                                                            int      ndim,
                                                            int     *where,
                                                            double  *value,
                                                            void    *config)
{
    int    ret = MATH_MINIMISATION_RESULT_MAX_FCT_CALLS_REACHED;

    int    j, nit;
    double y, ynew, yb, ybb, t, dt;
    int    *x, *xnew, *xb, *xbb, *grid_extents;
    int    initialisation, ninitialisations;

    grid_extents  = ((MultiExpDiscreteGridSearchConfigParams_t *)config)->gridextents;
	
    if ((ndim<=0) || (!where) || (!funk) || (!grid_extents))
        return (MATH_MINIMISATION_RESULT_ERROR_INVALID_INPUT);
	
    /* Initialize routine... */	
    x    = Bayes_ivector(1,ndim);
    xb   = Bayes_ivector(1,ndim);
    xbb  = Bayes_ivector(1,ndim);
    xnew = Bayes_ivector(1,ndim);

    /* Starting location as provided by calling routine... */
    for (j=1; j<=ndim; j++)
    {
        xbb[j] = xb[j] = x[j] = where[j];
    }

    ybb = yb = y = funk(x,id,container);

    nit = 0;
    rand_InitializeRandomSeed();

    /* Iterations show grow with increasing 'ndim'... */
    for (j=1,ninitialisations=2; j<=ndim; j++)
        ninitialisations *= 2;

    ninitialisations *= bayes_BiExpConfigGetRapidGridSearchInitialisationsFactor();

    t  = 150.0;
    dt = t/(double)ninitialisations;

    for (initialisation=0; initialisation<ninitialisations; initialisation++)
    {
        for (j=1; j<=ndim; j++)//new starting point in outer iteration, easier to explore more space here
            x[j] = (int)(rand_RandomFloat()*(float)grid_extents[j]);

        if (x[3]<x[4])
            x[4]=x[3];

        for (j=1; j<=ndim; j++)
            xnew[j] = x[j];

        ynew = funk(xnew,id,container);

        /*if ((ynew<=ybb) ||
            ((ynew>ybb) && (exp(-(ynew-ybb)/t)>rand_NrModifiedRan1(0))))
        {
            t -= dt;*/

            /* While one exists...                        */
            /* ...search for the optimal neighbour state. */
            while (bayes_MultiExpDiscreteSpaceSearchForPreferableNeighbourState(funk,id,container,ndim,xnew,&ynew,config))
            {
                /* Save the neighbouring state... */
                for (j=1; j<=ndim; j++)
                    x[j] = xnew[j];

                y = ynew;

                if (y<yb)
                {
                    /* New best point found in this inner loop... */
                    for (j=1; j<=ndim; j++)
                    {
                        xb[j] = x[j];
                    }

                    yb = y;

                    /*if (DebugTrace)
                    {
                        fprintf(fp,"New (inner) best state saved...\n");
                    }*/
                }

                nit++;
            }

            if (yb < ybb)
            {
                for (j=1; j<=ndim; j++)
                    xbb[j] = xb[j];

                ybb = yb;

            }
        /*}*/
    }

    *value = ybb; 

    for (j=1; j<=ndim; j++)
        where[j] = xbb[j];

    free_Bayes_ivector(x,1,ndim);
    free_Bayes_ivector(xb,1,ndim);
    free_Bayes_ivector(xbb,1,ndim);
    free_Bayes_ivector(xnew,1,ndim);

    return (ret);
}




#define BAYES_ERR_NO_ERROR                         0
#define BAYES_ERR_PARAM_ESTIMATION_FAILURE        -7



int bayes_RapidMultiExpMostProbWeightsAndTaus(int                          *data,
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
    int    *x, *xbest, *xmin, *xmax, id=0, i, k, ndim, modelfixedweightindex, *extents;
    int    ret;

    PsuedoRapidMultiExpMinusLogProbParams_t   container;
    MultiExpDiscreteGridSearchConfigParams_t  config;
    BayesProbDistn_t                         *distn=NULL;

    ndim     = 2*ndecays;
    ndim    -= paramfixing->nparamsuserfixed;
    x        = Bayes_ivector(1,ndim);
    extents  = Bayes_ivector(1,ndim);

    container.data                      = data;
    container.nbins                     = nbins;
    container.fitstart					= fitstart;
    container.binwalls                  = binwalls;
    container.instr                     = instr;
    container.interval                  = interval;
    container.modulationperiod          = modulationperiod;
    container.ndecays                   = ndecays;
    container.rapidparamsandlikelihoods = grid;
    container.hyperparam                = alpha;
    container.paramfixing               = paramfixing;

    /* Initialise input starting location 'x' first... */
    for (k=1,i=1; k<=ndecays; k++)
    {
        if (paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            x[i]       = bayes_MapWeightValueToClosestRapidGridPoint(weights_mp[k],grid->likelihoodsvaluestore->settings->nweights,grid->likelihoodsvaluestore->settings->weight);
            extents[i] = grid->likelihoodsvaluestore->settings->nweights;
            i++;
        }
    }

    for (k=1; k<=ndecays; k++)
    {
        if (paramfixing->tauuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            x[i]       = bayes_MapLifetimeValueToClosestRapidGridPoint(taus_mp[k],grid->likelihoodsvaluestore->settings->ntaus,grid->likelihoodsvaluestore->settings->tau);
            extents[i] = grid->likelihoodsvaluestore->settings->ntaus;
            i++;
        }
    }

    config.gridextents = extents;

    value = BAYES_SIZE_DOUBLE_HUGE;

    /* Search space using stochastic algorithm initially... */
    ret = bayes_MultiExpDiscreteSpaceMinimisationStochasticSearch(bayes_PsuedoRapidMultiExpMinusLogProbGivenWeightsAndTaus,
                                                                  id,(void*)(&container),
                                                                  ndim,x,&value,(void*)(&config));

    
	if (ret == MATH_MINIMISATION_RESULT_USERCANCEL)
		goto Cancel;

    if ((ret >= MATH_MINIMISATION_RESULT_SUCCESS) && 
        (value < BAYES_SIZE_DOUBLE_HUGE) &&
        (bayes_BiExpConfigRapidGetGridSearchLocalisedExhaustiveSearchDelta()))
    {
        /* Record current best state and look for a better one... */
        valuebest = value;
        xbest     = Bayes_ivector(1,ndim);
        
        for (i=1; i<=ndim; i++)
            xbest[i] = x[i];

        /* Setup the limits of the region around 'xbest' to be searched exhaustively... */
        xmin  = Bayes_ivector(1,ndim);
        xmax  = Bayes_ivector(1,ndim);

        for (i=1; i<=ndim; i++)
        {
	        xmin[i] = x[i]-bayes_BiExpConfigRapidGetGridSearchLocalisedExhaustiveSearchDelta();
            xmax[i] = x[i]+bayes_BiExpConfigRapidGetGridSearchLocalisedExhaustiveSearchDelta();
        }

        /* Constrain the limits of the region around 'xbest' to defined grid points... */
        for (i=1; i<=ndecays-paramfixing->nweightsuserfixed; i++)
        {
            xmin[i] = (xmin[i]<0)?(0):(xmin[i]);
            xmax[i] = (xmax[i]>grid->likelihoodsvaluestore->settings->nweights-1)?(grid->likelihoodsvaluestore->settings->nweights-1):(xmax[i]);
        }

        for (; i<=ndim; i++)
        {
            xmin[i] = (xmin[i]<0)?(0):(xmin[i]);
            xmax[i] = (xmax[i]>grid->likelihoodsvaluestore->settings->ntaus-1)?(grid->likelihoodsvaluestore->settings->ntaus-1):(xmax[i]);
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
                dweight = (grid->likelihoodsvaluestore->settings->weight[grid->likelihoodsvaluestore->settings->nweights-1]-grid->likelihoodsvaluestore->settings->weight[0])/(double)grid->likelihoodsvaluestore->settings->nweights;
                dtau    = (grid->likelihoodsvaluestore->settings->tau[grid->likelihoodsvaluestore->settings->ntaus-1]-grid->likelihoodsvaluestore->settings->tau[0])/(double)grid->likelihoodsvaluestore->settings->ntaus;

                ret = bayes_AllocateForMultiExpDiscreteProbDistn(distn,ndim,
                                                                 ndecays-paramfixing->nweightsuserfixed,
                                                                 ndecays-paramfixing->ntaususerfixed,
                                                                 dweight,dtau,xmin,xmax);

                if (ret<0)
                {
                    *val = BAYES_SIZE_DOUBLE_HUGE;

                    free_Bayes_ivector(x,1,ndim);
                    free_Bayes_ivector(xmin,1,ndim);
                    free_Bayes_ivector(xmax,1,ndim);
                    free_Bayes_ivector(extents,1,ndim);

                    return (BAYES_ERR_PARAM_ESTIMATION_FAILURE);
                }            
            }
        }

        /* Exhaustive search of region around the claimed minimum for a better solution... */
        config.gridmins = xmin;
        config.gridmaxs = xmax;
        config.distn    = distn;

        ret = bayes_MultiExpDiscreteSpaceMinimisationExhaustiveSearch(bayes_PsuedoRapidMultiExpMinusLogProbGivenWeightsAndTaus,
                                                                      id,(void*)(&container),
                                                                      ndim,x,&value,(void*)(&config));

        if ((ret>=0) && (value<valuebest))
        {
            valuebest = value;

            for (i=1; i<=ndim; i++)
	            xbest[i] = x[i];
        }

        /* Put our best results back into the free parameter array... */
        value = valuebest;
        
        for (i=1; i<=ndim; i++)
            x[i] = xbest[i];

        free_Bayes_ivector(xbest,1,ndim);
        free_Bayes_ivector(xmin,1,ndim);
        free_Bayes_ivector(xmax,1,ndim);
    }

    /* Check that the optimal values are sensible, i.e. at least exist on the grid... */
    for (i=1; i<=ndecays-paramfixing->nweightsuserfixed; i++)
    {
        if ((x[i]<0) && (x[i]>=grid->likelihoodsvaluestore->settings->nweights))
            ret = BAYES_ERR_PARAM_ESTIMATION_FAILURE;
    }

    for (; i<=ndim; i++)
    {
        if ((x[i]<0) && (x[i]>=grid->likelihoodsvaluestore->settings->ntaus))
            ret = BAYES_ERR_PARAM_ESTIMATION_FAILURE;
    }

    if (ret==BAYES_ERR_PARAM_ESTIMATION_FAILURE)
    {
        *val = BAYES_SIZE_DOUBLE_HUGE;

        free_Bayes_ivector(x,1,ndim);
        free_Bayes_ivector(extents,1,ndim);

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
                    grid->likelihoodsvaluestore->settings->weight,
                    dweight,
                    distn->marginals[i].indexlow,
                    distn->marginals[i].indexhigh,
                    &weights_ave[k],
                    &weights_err[k]);
            }

            weights_mp[k] = grid->likelihoodsvaluestore->settings->weight[x[i]];
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
                    grid->likelihoodsvaluestore->settings->tau,
                    dtau,
                    distn->marginals[i].indexlow,
                    distn->marginals[i].indexhigh,
                    &taus_ave[k],
                    &taus_err[k]);
            }

            taus_mp[k] = grid->likelihoodsvaluestore->settings->tau[x[i]];
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

    free_Bayes_ivector(x,1,ndim);
    free_Bayes_ivector(extents,1,ndim);

	if (!distribution && distn)
        bayes_FreeForMultiExpDiscreteProbDistn(distn);

    return (BAYES_ERR_NO_ERROR);
}

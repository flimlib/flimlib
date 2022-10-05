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
#include "stdlib.h"
#include "string.h"

#include "math.h"
#include "extmath.h"
#include "random.h"
#include "matrices.h"

#include "bayes_Types.h"
#include "bayes_Sizes.h"
#include "bayes_DataManagement.h"
#include "bayes_InstrRspAnalysis.h"
#include "bayes_DistributionFctsBinLikelihoods.h"
#include "bayes_MonoExpAnalysisBinLikelihoods.h"
#include "bayes_MultiExpAnalysisBinLikelihoods.h"
#include "bayes_ModelTransformTools.h"
#include "bayes_Interface.h"



BayesInstrRspEstimationStore_t gBayesEstimatedInstrRsp;


BayesInstrRspEstimationStore_t *bayes_GetBayesInstrRspEstimationStorePtr(void)
{
    return (&gBayesEstimatedInstrRsp);
}

/**
 * TODO: Thread Unsafe
 */
int bayes_UpdateEstimatedInstrRspDetailsStoreAdv(/* Image/prompt data used... */
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
    /* Image/prompt data used... */
    strcpy(gBayesEstimatedInstrRsp.filepromptloaded,fileprompt);
    strcpy(gBayesEstimatedInstrRsp.filedataloaded,filedata);
    gBayesEstimatedInstrRsp.allpixels    = allpixels;
    gBayesEstimatedInstrRsp.xcoordinate  = xcoordinate;
    gBayesEstimatedInstrRsp.ycoordinate  = ycoordinate;
    gBayesEstimatedInstrRsp.binningtype  = binningtype;
    gBayesEstimatedInstrRsp.binninglevel = binninglevel;

    /* Estimation analysis settings... */
    gBayesEstimatedInstrRsp.transbins     = transbins;
    gBayesEstimatedInstrRsp.transbinwidth = transbinwidth;
    gBayesEstimatedInstrRsp.fitstart      = fitstart;
    gBayesEstimatedInstrRsp.fitend        = fitend;
    gBayesEstimatedInstrRsp.nphotons      = nphotons;
    gBayesEstimatedInstrRsp.modeltype     = modeltype;

    /* Estimated instrument response... */
    bayes_CopyInstrRspConfigParams(instr,&gBayesEstimatedInstrRsp.instr);

    return (0);
}

#define BAYES_IR_SMALL_DELAY  0.10
#define BAYES_IR_SMALL_WIDTH  0.01
#define BAYES_IR_SMALL_WEIGHT 0.00
#define BAYES_IR_SMALL_CUTOFF 0.00


int bayes_CopyInstrRspConfigParams(BayesInstrRsp_t *source,
                                   BayesInstrRsp_t *destination)
{
    int i;
    
    if ((!source) || (!destination))
        return (-1);

    destination->ninstr = source->ninstr;

    for (i=0; i<source->ninstr; i++)
    {
        destination->params[i].weight = source->params[i].weight;
        destination->params[i].width  = source->params[i].width;
        destination->params[i].delay  = source->params[i].delay;
        destination->params[i].cutoff = source->params[i].cutoff;
    }

    //could nullify all other values, but not really reqd...
    return (0);
}


int bayes_SortInstrRspComponentsByWeight(BayesInstrRsp_t *instr)
{
    int ninstr, sortkeys[]={0,1,2};
    double          temp;

    if (!instr)
        return (-1);

    ninstr = instr->ninstr;

    if (ninstr<=1)
        return (0);

    if (instr->params[0].weight<instr->params[1].weight)
    {
        temp                    = instr->params[0].weight;
        instr->params[0].weight = instr->params[1].weight;
        instr->params[1].weight = temp;

        temp                    = instr->params[0].width;
        instr->params[0].width  = instr->params[1].width;
        instr->params[1].width  = temp;

        temp                    = instr->params[0].delay;
        instr->params[0].delay  = instr->params[1].delay;
        instr->params[1].delay  = temp;

        temp                    = instr->params[0].cutoff;
        instr->params[0].cutoff = instr->params[1].cutoff;
        instr->params[1].cutoff = temp;
    }

    return (0);
}


//return '1' if instrument rsp containers hold a different configuration...
int bayes_CheckForDifferentInstrRspConfigParams(BayesInstrRsp_t *source,
                                                BayesInstrRsp_t *destination)
{
    int i;

    if (destination->ninstr != source->ninstr)
        return (1);

    for (i=0; i<source->ninstr; i++)
    {
        if ((destination->params[i].weight != source->params[i].weight) ||
            (destination->params[i].width  != source->params[i].width) ||
            (destination->params[i].delay  != source->params[i].delay) ||
            (destination->params[i].cutoff != source->params[i].cutoff))
        {
            return (1);
        }
    }

    //could nullify all other values, but not really reqd...

    return (0);
}


double bayes_LogProbInstrRspLikelihoodGivenDelayAndWidth(double u,
			                                             double uc,
			                                             double sigma)
{
	double val;
	
	val  = (u-uc)/(ROOTTWO*sigma);
	val *= -val;
	
	val += log(2.0);
	
	val -= log(sigma*ROOTTWO*ROOTPI);
	
	val -= log(1.0+erf(uc/(sigma*ROOTTWO)));
	
	return (val);
}


double bayes_LogProbInstrRspLikelihoodGivenDelayAndWidthAndCutoff(double u,
			                                                      double uc,
			                                                      double sigma,
                                                                  double delta)
{
	double val;
	
	val  = (u-uc)/(ROOTTWO*sigma);
	val *= -val;
	
	val += log(2.0);
	
	val -= log(sigma*ROOTTWO*ROOTPI);
	
	val -= log(1.0+erf((uc-delta)/(sigma*ROOTTWO)));
	
	return (val);
}


double bayes_InstrRspLikelihoodGivenDelayAndWidth(double u,
                                                  double uc,
                                                  double sigma)
{
    double val;

    if ((u<0.0) || (uc<0.0) || (sigma<0.0))
        return (0.0);
	
	val = bayes_LogProbInstrRspLikelihoodGivenDelayAndWidth(u,uc,sigma);
	val = exp(val);

    return (val);
}


double bayes_InstrRspLikelihoodGivenDelayAndWidthAndCutoff(double u,
                                                           double uc,
                                                           double sigma,
                                                           double delta)
{
    double val;

    if ((u<0.0) || (uc<0.0) || (sigma<0.0) || (u<delta))
        return (0.0);
	
	val = bayes_LogProbInstrRspLikelihoodGivenDelayAndWidthAndCutoff(u,uc,sigma,delta);
	val = exp(val);

    return (val);
}



#if 1
int bayes_FitPredictedApproxInstrRsp(float           *fitted,
                                     int              nbins,
                                     float            binwidth,
                                     BayesInstrRsp_t *instr)
{
    int    i, j;
    double u, val, weight, width, delay, cutoff;

	for (i=0; i<nbins; i++)
	{
        u = binwidth*(0.5+(float)i);

        for (j=0, val=0.0; j<BAYES_INSTR_RSP_MAX_COMPONENTS; j++)
        {
            weight = instr->params[j].weight;
            width  = instr->params[j].width;

            if ((weight>0.0) && (width>0.0))
            {  
                delay  = instr->params[j].delay;
                cutoff = instr->params[j].cutoff;

                val   += weight*bayes_InstrRspLikelihoodGivenDelayAndWidthAndCutoff(u,delay,width,cutoff);
            }
        }
		
        fitted[i] = (float)val;
	}
	
    return (0);
}
#else
int bayes_FitPredictedApproxInstrRsp(float *fitted,
                                     int    nbins,
                                     float  binwidth,
                                     float  gamma1,
                                     float  delta1,
                                     float  sigma1,
                                     float  delay1,
                                     float  gamma2,
                                     float  delta2,
                                     float  sigma2,
                                     float  delay2,
                                     float  gamma3,
                                     float  delta3,
                                     float  sigma3,
                                     float  delay3)
{
    int   i;
	float u, val;

    if (((delay1>0.0) && (sigma1>0.0)) || ((delay2>0.0) && (sigma2>0.0)) || ((delay3>0.0) && (sigma3>0.0)))
	{
		for (i=0; i<nbins; i++)
		{
            u          = binwidth*(0.5+(float)i);
        	val        = gamma1*bayes_InstrRspLikelihoodGivenDelayAndWidthAndCutoff(u,delay1,sigma1,delta1);
            
            if (gamma2>0.0)
                val += gamma2*bayes_InstrRspLikelihoodGivenDelayAndWidthAndCutoff(u,delay2,sigma2,delta2);

            if (gamma3>0.0)
                val += gamma3*bayes_InstrRspLikelihoodGivenDelayAndWidthAndCutoff(u,delay3,sigma3,delta3);
			
            fitted[i]  = val;
		}
	}
	else if ((delay1>0.0) || (delay2>0.0))
	{
		for (i=0; i<nbins; i++)
			fitted[i] = 0.0;
		
		fitted[(int)(delay1/binwidth)] = 1.0;
	}
	else
	{
		for (i=0; i<nbins; i++)
			fitted[i] = 0.0;
	}
	
    return (0);
}
#endif




/////////////////////////// NEW INSTR RSP //////////////////
/* Instrument response and mono-exponential parameter estimation */



double bayes_InstrRspAndMonoExpParamsMinusLogProb(double *x, int id, void *container)
{
    int     *data, bin, nbins, nphotons, nphotonsbin, i, j, type, ret;
    double   interval, modperiod, temp;
    double   alpha, value, w0, oneminusw0, w1, bjoverT, bL, bH;
    double  *likelihoods, *binwalls;
    int      ndecays;
    double   weights[2], taus[2];

    MonoExpMinusLogProbW0W1_t *params1;
    BayesInstrRsp_t           *instr;

    w0 = x[1]; i=2;
    w1 = x[2]; i++;

    if ((w0 < 0.0) || (w0 > 1.0))
        return (BAYES_SIZE_DOUBLE_HUGE);

    oneminusw0 = 1.0-w0;
    
    if (w1 < SMALLTIME)
        return (BAYES_SIZE_DOUBLE_HUGE);
	
	params1  = (MonoExpMinusLogProbW0W1_t *)(container);
	interval = params1->interval;
	
    instr = params1->instr;

	if ((instr->ninstr<1) || (instr->ninstr>BAYES_INSTR_RSP_MAX_COMPONENTS))
		return (BAYES_SIZE_DOUBLE_HUGE);

    if (!instr->paramsfixed[0].delayfixed)  {instr->params[0].delay  = x[i]; i++;}
    if (!instr->paramsfixed[0].cutofffixed) {instr->params[0].cutoff = x[i]; i++;}
    if (!instr->paramsfixed[0].widthfixed)  {instr->params[0].width  = x[i]; i++;}

    if ((instr->params[0].delay<BAYES_IR_SMALL_DELAY)   || (instr->params[0].delay>=interval) ||
        (instr->params[0].cutoff<BAYES_IR_SMALL_CUTOFF) || (instr->params[0].cutoff>=interval) ||
        (instr->params[0].width<BAYES_IR_SMALL_WIDTH)   || (instr->params[0].width>=interval))
    {
        return (BAYES_SIZE_DOUBLE_HUGE);
    }

    if (instr->ninstr>1)
    {
        for (j=1,temp=1.0; j<instr->ninstr; j++)
        {
            if (!instr->paramsfixed[j].weightfixed) {instr->params[j].weight = x[i]; i++;}
            if (!instr->paramsfixed[j].delayfixed)  {instr->params[j].delay  = x[i]; i++;}
            if (!instr->paramsfixed[j].cutofffixed) {instr->params[j].cutoff = x[i]; i++;}
            if (!instr->paramsfixed[j].widthfixed)  {instr->params[j].width  = x[i]; i++;}

            if ((instr->params[j].weight<BAYES_IR_SMALL_WEIGHT) || (instr->params[j].weight>1.0) ||
                (instr->params[j].delay<BAYES_IR_SMALL_DELAY)   || (instr->params[j].delay>=interval) ||
                (instr->params[j].cutoff<BAYES_IR_SMALL_CUTOFF) || (instr->params[j].cutoff>=interval) ||
                (instr->params[j].width<BAYES_IR_SMALL_WIDTH)   || (instr->params[j].width>=interval))
            {
                return (BAYES_SIZE_DOUBLE_HUGE);
            }

            temp -= instr->params[j].weight;
        }

        instr->params[0].weight = temp;

        if ((temp<0.0) || (temp>1.0))
            return (BAYES_SIZE_DOUBLE_HUGE);
    }
    else
    {
        instr->params[0].weight = 1.0;
        j = 1;
    }

    for (; j<BAYES_INSTR_RSP_MAX_COMPONENTS; j++)
    {
        instr->params[j].weight = 0.0;
        instr->params[j].delay  = 0.0;
        instr->params[j].cutoff = 0.0;
        instr->params[j].width  = 0.0;        
    }

    modperiod   = params1->modulationperiod;
    data        = params1->data;
    nbins       = params1->nbins;
    binwalls    = params1->binwalls;
    nphotons    = params1->nphotons;
    alpha       = params1->hyperparam;

    likelihoods = Bayes_dvector(0,nbins);

    ndecays     = 1;
    weights[0]  = w0;
    weights[1]  = oneminusw0;
    taus[1]     = w1;

    ret = bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(likelihoods,
                                                                     nbins,binwalls,data,
                                                                     interval,modperiod,instr,
                                                                     w1,ndecays,weights,taus);
//    bayes_ArrBinLikelihoodsGivenTau(likelihoods,NULL,NULL,data,nbins,interval,modperiod,&instr,w1);
    
    if (ret<0)
    {
        free_Bayes_dvector(likelihoods,0,nbins);
        return (BAYES_SIZE_DOUBLE_HUGE);
    }

    value = 0.0;

#if 0 //configurable prior - code needs finishing
    /* Compute prior likelihhod of instr parameter if required */
    for (j=0,temp=0.0; j<instr->ninstr; j++)
    {
        if (instr->paramsprior[j].weightprior)
        {
            temp += (*instr->paramsprior[j].weightpriorevaluator.funk)(instr->params[j].weight,
                                                                      instr->paramsprior[j].weightpriorevaluator.nparams,
                                                                      instr->paramsprior[j].weightpriorevaluator.params,
                                                                      NULL);
        }

        if (instr->paramsprior[j].delayprior)
        {
            temp += (*instr->paramsprior[j].delaypriorevaluator.funk)(instr->params[j].delay,
                                                                     instr->paramsprior[j].delaypriorevaluator.nparams,
                                                                     instr->paramsprior[j].delaypriorevaluator.params,
                                                                     NULL);
        }

        if (instr->paramsprior[j].cutoffprior)
        {
            temp += (*instr->paramsprior[j].cutoffpriorevaluator.funk)(instr->params[j].cutoff,
                                                                      instr->paramsprior[j].cutoffpriorevaluator.nparams,
                                                                      instr->paramsprior[j].cutoffpriorevaluator.params,
                                                                      NULL);
        }

        if (instr->paramsprior[j].widthprior)
        {
            temp += (*instr->paramsprior[j].widthpriorevaluator.funk)(instr->params[j].width,
                                                                    instr->paramsprior[j].widthpriorevaluator.nparams,
                                                                    instr->paramsprior[j].widthpriorevaluator.params,
                                                                    NULL);
        }
    }
#else
    temp = 0.0;
#endif

    value += temp;

    for (bin=0,value+=alpha*w1; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            bL       = binwalls[bin];
            bH       = binwalls[bin+1];
            bjoverT  = (bH-bL)/interval;
            value   -= (double)(nphotonsbin) * (log(w0*bjoverT + oneminusw0*likelihoods[bin]));

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value,&type))
            {
                free_Bayes_dvector(likelihoods,0,nbins);
	            return (BAYES_SIZE_DOUBLE_HUGE);
            }
        }
    }

    free_Bayes_dvector(likelihoods,0,nbins);

    return (value);
}



#define	FREE_DECAYPHOTONLIKELIHOODS	{int ii; for(ii=0;ii<(1+ndecays);ii++) if(decayphotonlikelihoods[ii]) free_Bayes_dvector(decayphotonlikelihoods[ii],0,nbins-1); free(decayphotonlikelihoods);}

double bayes_InstrRspAndMultiExpParamsMinusLogProb(double *x, int id, void *container)
{
    int     *data, bin, nbins, nphotons, nphotonsbin, i, j, k, modelfixedweightindex, fitstart, type, ndecays, ret;
    double   interval, modperiod, *binwalls, temp;
    double   value, bjoverT, bL, bH;
    double  *weights, *taus, **decayphotonlikelihoods, hyperparam;

    MultiExpMinusLogProbParams_t *params;
    BayesUserFixedParams_t       *paramfixing;
    BayesInstrRsp_t              *instr;

	params = (MultiExpMinusLogProbParams_t *)(container);
    ndecays = params->ndecays;
    weights     = Bayes_dvector(0,ndecays);
    taus        = Bayes_dvector(1,ndecays);
    paramfixing = params->paramfixing;

    /* Determine which weight value is fixed by the model (i.e. weights summing to '1')... */
    for (k=0; k<=ndecays; k++)
    {
        if (params->paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            modelfixedweightindex = k;
            break;
        }
    }

    /* Extract and validate the weights from 'x' first... */
    for (k=1,i=1; k<=ndecays; k++)
    {
        if (paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            weights[k] = x[i];
            i++;
        }
        else
        {
            weights[k] = params->paramfixing->weights[k];
        }
    }

    for (k=0,value=1.0; k<=ndecays; k++)
    {
        if (k != modelfixedweightindex)
            value -= weights[k];
    }

    weights[modelfixedweightindex] = value;
    for (k=0; k<=ndecays; k++)
    {
        if ((weights[k]<0.0) || (weights[k]>1.0)) /* Check that weights are valid... */
        {
            free_Bayes_dvector(weights,0,ndecays);
            free_Bayes_dvector(taus,1,ndecays);
            return (BAYES_SIZE_DOUBLE_HUGE);
        }
    }

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

    /* Now extract the instrument response parameters... */
    instr        = params->instr;
    interval     = params->interval;

	if ((instr->ninstr<1) || (instr->ninstr>BAYES_INSTR_RSP_MAX_COMPONENTS))
		return (BAYES_SIZE_DOUBLE_HUGE);

    if (!instr->paramsfixed[0].delayfixed)  {instr->params[0].delay  = x[i]; i++;}
    if (!instr->paramsfixed[0].cutofffixed) {instr->params[0].cutoff = x[i]; i++;}
    if (!instr->paramsfixed[0].widthfixed)  {instr->params[0].width  = x[i]; i++;}
    for (j=1; j<instr->ninstr; j++)
    {
        if (!instr->paramsfixed[j].weightfixed) {instr->params[j].weight = x[i]; i++;}
        if (!instr->paramsfixed[j].delayfixed)  {instr->params[j].delay  = x[i]; i++;}
        if (!instr->paramsfixed[j].cutofffixed) {instr->params[j].cutoff = x[i]; i++;}
        if (!instr->paramsfixed[j].widthfixed)  {instr->params[j].width  = x[i]; i++;}
    }

    if (instr->ninstr>1)
    {
        for (j=1,temp=1.0; j<instr->ninstr; j++)
            temp -= instr->params[j].weight;

        instr->params[0].weight = temp;
    }
    else
    {
        instr->params[0].weight = 1.0;
    }

    for (j=0; j<instr->ninstr; j++)
    {
        if ((instr->params[j].weight<BAYES_IR_SMALL_WEIGHT) || (instr->params[j].weight>1.0) ||
            (instr->params[j].delay<BAYES_IR_SMALL_DELAY)   || (instr->params[j].delay>=interval) ||
            (instr->params[j].cutoff<BAYES_IR_SMALL_CUTOFF) || (instr->params[j].cutoff>=interval) ||
            (instr->params[j].width<BAYES_IR_SMALL_WIDTH)   || (instr->params[j].width>=interval))
        {
            free_Bayes_dvector(weights,0,ndecays);
            free_Bayes_dvector(taus,1,ndecays);
            return (BAYES_SIZE_DOUBLE_HUGE);
        }    
    }

    for (; j<BAYES_INSTR_RSP_MAX_COMPONENTS; j++)
    {
        instr->params[j].weight = 0.0;
        instr->params[j].delay  = 0.0;
        instr->params[j].cutoff = 0.0;
        instr->params[j].width  = 0.0;        
    }



    /* Weights, lifetimes, instrument values okay... */
    data        = params->data;
    nbins       = params->nbins;
    fitstart    = params->fitstart;
    binwalls    = params->binwalls;
    nphotons    = params->nphotons;
    modperiod   = params->modulationperiod;
    hyperparam  = params->hyperparam;

    /* Compute the minus log likelihood using pre-computed bin likelihoods... */
    decayphotonlikelihoods = (double **) calloc((1+ndecays), sizeof(double*));

    for (k=1; k<=ndecays; k++)
    {
        if (params->paramfixing->tauuserfixed[k] == BAYES_PARAM_VALUE_FREE)
        {
            decayphotonlikelihoods[k] = Bayes_dvector(0,nbins-1);

            ret = bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(decayphotonlikelihoods[k],
                                                                             nbins,binwalls,data,
                                                                             interval,modperiod,instr,
                                                                             taus[k],
                                                                             ndecays,weights,taus);

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

    if (ret<0)
    {
        free_Bayes_dvector(weights,0,ndecays);
        free_Bayes_dvector(taus,1,ndecays);
		FREE_DECAYPHOTONLIKELIHOODS;
        return (BAYES_SIZE_DOUBLE_HUGE);
    }

    for (i=1, value=0.0; i<=ndecays; i++)
        value += (hyperparam*taus[i]); 

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            bL      = binwalls[bin];
            bH      = binwalls[bin+1];
            bjoverT = (bH-bL)/interval;

            for (i=1, temp=weights[0]*bjoverT; i<=ndecays; i++)
                temp += (weights[i]*decayphotonlikelihoods[i][bin]);

            value -= (double)(nphotonsbin)*(log(temp));

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value,&type))
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

    return (value);
}


double bayes_GetRandomLifetimeValue(double lifetimemax)
{
    return (lifetimemax*rand_RandomDouble());
}


double bayes_GetRandomWeightValue(double weightmax)
{
    return (weightmax*rand_RandomDouble());
}


double bayes_GetRandomInstrRspDelayValue(double delaymax)
{
    return (delaymax*rand_RandomDouble());
}


double bayes_GetRandomInstrRspWeightValue(double weightmax)
{
    return (weightmax*rand_RandomDouble());
}


double bayes_GetRandomInstrRspWidthValue(double widthmax)
{
    return (widthmax*rand_RandomDouble());
}


#define BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_WIDTH    0.02
#define BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_WEIGHT   0.05
#define BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_DELAY    0.05
#define BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_CUTOFF   0.05
#define BAYES_INSTR_SEARCH_LENGTH_SCALE_DECAY_WEIGHT 0.01
#define BAYES_INSTR_SEARCH_LENGTH_SCALE_DECAY_TAU    0.05

int bayes_DirectInstrRspAndMonoExpOptimization(/* Data in... */
                                               int                      *data,
                                               int                       nbins,
                                               int                       fitstart,
                                               double                   *binwalls,
                                               int                       nphotons,
                                               /* Estimates out... */
                                               double                   *w0,
                                               double                   *w1,
                                               BayesInstrRsp_t          *instr,
                                               double                   *minuslogprob,
                                               /* Instrument... */
                                               double                    interval,
                                               double                    modulationperiod,
                                               double                    alpha)
{
    double bestvalue, value, *bestx, *x, *deltas, temp;
    int    ret, id=0, i, j, ninstr, ninstrfreeparams, ninstrfixedparams, ndim, restarts, nrestarts;
    void   *config;

    MonoExpMinusLogProbW0W1_t container;

    int   (*minimizer)(double (*)(double *, int, void *), int, void *, int, double *, double *, void *);

    container.data             = data;
    container.nbins            = nbins;
    container.binwalls         = binwalls;
    container.nphotons         = nphotons;
    container.interval         = interval;
    container.modulationperiod = modulationperiod;
    container.hyperparam       = alpha;
	container.upsilon1         = NULL;
    container.instr            = instr;

    ninstr                     = instr->ninstr;
    for (i=0,ninstrfreeparams=0,ninstrfixedparams=0; i<ninstr; i++)
    {
        if (instr->paramsfixed[i].weightfixed) ninstrfixedparams++; else ninstrfreeparams++;
        if (instr->paramsfixed[i].delayfixed)  ninstrfixedparams++; else ninstrfreeparams++;
        if (instr->paramsfixed[i].cutofffixed) ninstrfixedparams++; else ninstrfreeparams++;
        if (instr->paramsfixed[i].widthfixed)  ninstrfixedparams++; else ninstrfreeparams++;
    }
    
    ninstrfreeparams -= 1; //-1 as weights must sum to '1'...
    ndim                       = 2+ninstrfreeparams;

    bestx     = Bayes_dvector(1,ndim);
    deltas    = Bayes_dvector(1,ndim); //required for the minimization algorithms
    bestx[1]  = (double)(*w0);
    deltas[1] = BAYES_INSTR_SEARCH_LENGTH_SCALE_DECAY_WEIGHT;
    bestx[2]  = (double)(*w1);
    deltas[2] = BAYES_INSTR_SEARCH_LENGTH_SCALE_DECAY_TAU;

    for (i=0,j=3; i<ninstr; i++)
    {
        if ((i>0) && /* Will always have the weight of component '0' being determined by the other weights... */
            (!instr->paramsfixed[i].weightfixed))
        {
            bestx[j]  = instr->params[i].weight;
            deltas[j] = BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_WEIGHT;
            j++;
        }

        if (!instr->paramsfixed[i].delayfixed)
        {
            bestx[j]  = instr->params[i].delay;
            deltas[j] = BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_DELAY;
            j++;
        }

        if (!instr->paramsfixed[i].cutofffixed)
        {
            bestx[j]  = instr->params[i].cutoff;
            deltas[j] = BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_CUTOFF;
            j++;
        }

        if (!instr->paramsfixed[i].widthfixed)
        {
            bestx[j]  = instr->params[i].width;
            deltas[j] = BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_WIDTH;
            j++;
        }    
    }

    minimizer = &math_MinimiseFctDoubleWithGenericContainer;
    config = malloc(sizeof(AmoebaConfigParams_t));
    ((AmoebaConfigParams_t*)config)->monitor   = 0;
    ((AmoebaConfigParams_t*)config)->tolerance = bayes_InstrConfigGetDownhillSimplexPrecision();
    ((AmoebaConfigParams_t*)config)->deltas    = deltas;
 
    bestvalue = BAYES_SIZE_DOUBLE_HUGE;

    ret = minimizer(bayes_InstrRspAndMonoExpParamsMinusLogProb,id,(void*)(&container),ndim,bestx,&bestvalue,(void*)(config));

	if (ret >= MATH_MINIMISATION_RESULT_SUCCESS)
	{
		nrestarts = bayes_InstrConfigGetNumberOfRestarts();

		if (nrestarts)
	    {
            x = Bayes_dvector(1,ndim);

		    for (restarts=0; restarts<nrestarts; restarts++)
		    {
	            /* Randomly initialize the search starting location in parameter space... */
                x[1] = bayes_GetRandomWeightValue(0.20);
                x[2] = bayes_GetRandomLifetimeValue(interval/2.0);
                x[3] = bayes_GetRandomInstrRspDelayValue(interval/3.0); //delay
                x[4] = bayes_GetRandomInstrRspDelayValue(interval/3.0); //cutoff
                x[5] = bayes_GetRandomInstrRspWidthValue(interval/3.0); //width

                if (ninstr>1)
                {
                    for (i=0,j=1; j<ninstr; i++,j++)
                    {
                        x[6+i*4] = bayes_GetRandomInstrRspWeightValue(1.0);
                        x[7+i*4] = bayes_GetRandomInstrRspDelayValue(interval/3.0); //delay
                        x[8+i*4] = bayes_GetRandomInstrRspDelayValue(interval/3.0); //cutoff
                        x[9+i*4] = bayes_GetRandomInstrRspWidthValue(interval/3.0); //width   
                    }    
                }
                
	            /* Perform the search... */
                value = BAYES_SIZE_DOUBLE_HUGE;
                minimizer(bayes_InstrRspAndMonoExpParamsMinusLogProb,id,(void*)(&container),ndim,x,&value,(void*)(config));

	            /* Update the best result if required... */
                if ((ret >= MATH_MINIMISATION_RESULT_SUCCESS) && (value<BAYES_SIZE_DOUBLE_HUGE))
                {
		            if (value < bestvalue)
		            {
	                    bestvalue = value;

	                    for (i=1; i<=ndim; i++)
				            bestx[i] = x[i];
	                }
                }
	        }

	        free_Bayes_dvector(x,1,ndim);
	    }
    }

    *w0    = (float)bestx[1]; j=2;
    *w1    = (float)bestx[2]; j++;

    if (!instr->paramsfixed[0].delayfixed)  {instr->params[0].delay  = bestx[j]; j++;}
    if (!instr->paramsfixed[0].cutofffixed) {instr->params[0].cutoff = bestx[j]; j++;}
    if (!instr->paramsfixed[0].widthfixed)  {instr->params[0].width  = bestx[j]; j++;}

    for (i=1; i<ninstr; i++)
    {
        if (!instr->paramsfixed[i].weightfixed) {instr->params[i].weight = bestx[j]; j++;}
        if (!instr->paramsfixed[i].delayfixed)  {instr->params[i].delay  = bestx[j]; j++;}
        if (!instr->paramsfixed[i].cutofffixed) {instr->params[i].cutoff = bestx[j]; j++;}
        if (!instr->paramsfixed[i].widthfixed)  {instr->params[i].width  = bestx[j]; j++;}
    }

    if (ninstr>1)
    {
        for (i=1,temp=1.0; i<ninstr; i++)
            temp -= instr->params[i].weight;

        instr->params[0].weight = temp;
    }
    else
    {
        instr->params[0].weight = 1.0;
    }

    for (; i<BAYES_INSTR_RSP_MAX_COMPONENTS; i++)
    {
        instr->params[i].weight = 0.0;
        instr->params[i].delay  = 0.0;
        instr->params[i].cutoff = 0.0;
        instr->params[i].width  = 0.0;        
    }

    *minuslogprob = bestvalue;

    free_Bayes_dvector(bestx,1,ndim);
    free_Bayes_dvector(deltas,1,ndim);
	
	return (0);
}


int bayes_DirectInstrRspAndMultiExpOptimization(/* Data in... */
                                                int                      *data,
                                                int                       nbins,
                                                int                       fitstart,
                                                double                   *binwalls,
                                                int                       nphotons,
                                                /* Estimates out... */
                                                int                       ndecays,
                                                double                   *weights,
                                                double                   *taus,
                                                BayesUserFixedParams_t   *paramfixing,
                                                BayesInstrRsp_t          *instr,
                                                double                   *minuslogprob,
                                                /* Instrument... */
                                                double                    interval,
                                                double                    modulationperiod,
                                                double                    alpha)
{
    double bestvalue, *bestx, temp, *deltas;
    int    ret, id=0, i, j, k, modelfixedweightindex, ninstr, ninstrfreeparams, ninstrfixedparams, ndim;
    void   *config;

    MultiExpMinusLogProbParams_t container;

    int   (*minimizer)(double (*)(double *, int, void *), int, void *, int, double *, double *, void *);

    container.data             =  data;
    container.nbins            =  nbins;
    container.binwalls         =  binwalls;
    container.nphotons         =  nphotons;
    container.nparams          =  2*ndecays+1;
    container.ndecays          =  ndecays;
    container.paramfixing      =  paramfixing;
    container.interval         =  interval;
    container.modulationperiod =  modulationperiod;
    container.instr            =  instr;
    container.hyperparam       =  alpha;
    container.normalization    =  0.0;

    /* Initialise input starting location (weight and lifetime) 'paramsfree' first... */
    ninstr = instr->ninstr;

    for (i=0,ninstrfreeparams=0,ninstrfixedparams=0; i<ninstr; i++)
    {
        if (instr->paramsfixed[i].weightfixed) ninstrfixedparams++; else ninstrfreeparams++;
        if (instr->paramsfixed[i].delayfixed)  ninstrfixedparams++; else ninstrfreeparams++;
        if (instr->paramsfixed[i].cutofffixed) ninstrfixedparams++; else ninstrfreeparams++;
        if (instr->paramsfixed[i].widthfixed)  ninstrfixedparams++; else ninstrfreeparams++;
    }
    
    ninstrfreeparams -= 1; //-1 as weights must sum to '1'...
    ndim              = (2*ndecays)-paramfixing->nparamsuserfixed+ninstrfreeparams; //instrument and decay parameters...

    bestx  = Bayes_dvector(1,ndim);
    deltas = Bayes_dvector(1,ndim); //required for the minimization algorithms

    for (k=1,i=1; k<=ndecays; k++)
    {
        if (paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            bestx[i] = weights[k];
            deltas[i] = BAYES_INSTR_SEARCH_LENGTH_SCALE_DECAY_WEIGHT;
            i++;
        }
    }

    for (k=1; k<=ndecays; k++)
    {
        if (paramfixing->tauuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
        {
            bestx[i] = taus[k];
            deltas[i] = BAYES_INSTR_SEARCH_LENGTH_SCALE_DECAY_TAU;
            i++;
        }
    }

    /* ...now initialise input starting location instrument response values */
    for (j=0; j<ninstr; j++)
    {
        if ((j>0) && /* Will always have the weight of component '0' being determined by the other weights... */
            (!instr->paramsfixed[j].weightfixed))
        {
            bestx[i]  = instr->params[j].weight;
            deltas[i] = BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_WEIGHT;
            i++;
        }

        if (!instr->paramsfixed[j].delayfixed)
        {
            bestx[i]  = instr->params[j].delay;
            deltas[i] = BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_DELAY;
            i++;
        }

        if (!instr->paramsfixed[j].cutofffixed)
        {
            bestx[i]  = instr->params[j].cutoff;
            deltas[i] = BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_CUTOFF;
            i++;
        }

        if (!instr->paramsfixed[j].widthfixed)
        {
            bestx[i]  = instr->params[j].width;
            deltas[i] = BAYES_INSTR_SEARCH_LENGTH_SCALE_IRF_WIDTH;
            i++;
        }
    }

    minimizer = &math_MinimiseFctDoubleWithGenericContainer;
    config = malloc(sizeof(AmoebaConfigParams_t));
    ((AmoebaConfigParams_t*)config)->monitor   = 0;
    ((AmoebaConfigParams_t*)config)->tolerance = bayes_InstrConfigGetDownhillSimplexPrecision();
    ((AmoebaConfigParams_t*)config)->deltas    = deltas;

    bestvalue = BAYES_SIZE_DOUBLE_HUGE;

    ret = minimizer(bayes_InstrRspAndMultiExpParamsMinusLogProb,id,(void*)(&container),ndim,bestx,&bestvalue,(void*)(config));
#if 0
	if (ret >= MATH_MINIMISATION_RESULT_SUCCESS)
	{
		nrestarts = bayes_InstrConfigGetNumberOfRestarts();

		if (nrestarts)
	    {
            x = Bayes_dvector(1,ndim);

		    for (restarts=0; restarts<nrestarts; restarts++)
		    {
	            /* Randomly initialize the search starting location in parameter space... */
                x[1] = bayes_GetRandomWeightValue(0.20);
                x[2] = bayes_GetRandomLifetimeValue(interval/2.0);
                x[3] = bayes_GetRandomInstrRspDelayValue(interval/3.0);
                x[4] = bayes_GetRandomInstrRspWidthValue(interval/3.0);

                if (ninstr>1)
                {
                    for (i=0,j=2; j<=ninstr; i++,j++)
                    {
                        x[5+i*3] = bayes_GetRandomInstrRspWeightValue(1.0);
                        x[6+i*3] = bayes_GetRandomInstrRspDelayValue(interval/3.0);
                        //x[7+i*4] = instr->params[j].cutoff;
                        x[7+i*3] = bayes_GetRandomInstrRspWidthValue(interval/3.0);    
                    }    
                }
                
	            /* Perform the search... */
                value = BAYES_SIZE_DOUBLE_HUGE;
                minimizer(bayes_InstrRspAndMonoExpParamsMinusLogProb,id,(void*)(&container),ndim,x,&value,(void*)(config));

	            /* Update the best result if required... */
                if ((ret >= MATH_MINIMISATION_RESULT_SUCCESS) && (value<BAYES_SIZE_DOUBLE_HUGE))
                {
		            if (value < bestvalue)
		            {
	                    bestvalue = value;

	                    for (i=1; i<=ndim; i++)
				            bestx[i] = x[i];
	                }
                }
	        }

	        free_Bayes_dvector(x,1,ndim);
	    }
    }
#endif

    /* Recombine the fixed and optimal free parameter values to yield parameter vector... */
    if (ret >= MATH_MINIMISATION_RESULT_SUCCESS)
    {
        for (k=1,i=1; k<=ndecays; k++)
        {
            if (paramfixing->weightuserfixed[k] != BAYES_PARAM_VALUE_USER_FIXED)
            {
                weights[k] = bestx[i];
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
                taus[k] = bestx[i];
                i++;
            }
            else
            {
                taus[k] = paramfixing->taus[k];
            }
        }

        if (!instr->paramsfixed[0].delayfixed)  {instr->params[0].delay  = bestx[i]; i++;}
        if (!instr->paramsfixed[0].cutofffixed) {instr->params[0].cutoff = bestx[i]; i++;}
        if (!instr->paramsfixed[0].widthfixed)  {instr->params[0].width  = bestx[i]; i++;}

        for (j=1; j<ninstr; j++)
        {
            if (!instr->paramsfixed[j].weightfixed) {instr->params[j].weight = bestx[i]; i++;}
            if (!instr->paramsfixed[j].delayfixed)  {instr->params[j].delay  = bestx[i]; i++;}
            if (!instr->paramsfixed[j].cutofffixed) {instr->params[j].cutoff = bestx[i]; i++;}
            if (!instr->paramsfixed[j].widthfixed)  {instr->params[j].width  = bestx[i]; i++;}
        }

        if (ninstr>1)
        {
            for (j=1,temp=1.0; j<ninstr; j++)
                temp -= instr->params[j].weight;

            instr->params[0].weight = temp;
        }
        else
        {
            instr->params[0].weight = 1.0;
        }

        for (; j<BAYES_INSTR_RSP_MAX_COMPONENTS; j++)
        {
            instr->params[j].weight = 0.0;
            instr->params[j].delay  = 0.0;
            instr->params[j].cutoff = 0.0;
            instr->params[j].width  = 0.0;        
        }

	    *minuslogprob = (float)bestvalue;
	}

    free_Bayes_dvector(bestx,1,ndim);
	
	return (0);
}

//renamed from DoBayesInstrRspCoarseGuessValuesFromLoadedInstr
#if 0
int bayes_InstrRspCoarseGuessValuesFromLoadedInstr(float *instr,
                                                   float  binwidth,
                                                   int    nbins,
                                                   float *delay1,
                                                   float *width1,
                                                   float *delay2,
                                                   float *width2)
{
    double p1, p2, max;
    int    i, maxi;

    if ((!instr) || (binwidth<=0.0) || (nbins<=0))
        return (-1);

    /* Search for the dominant delay and spreading around it... */
    for (i=0, max=0.0,maxi=0; i<nbins; i++)
    {
        if (instr[i]>max)
        {
            max  = instr[i];
            maxi = i;
        }
    }

    *delay = ((double)maxi+0.5)*binwidth;

    for (i=maxi; i>=0; i--) /* Look for drop below FWHM to left of dominant delay... */
	{
		if (instr[i]<max/2.0) 
		{
			p1 = i;
			break;
		}
	}

	for (i=maxi; i<nbins; i++) /* ...and now to the right of dominant delay */
	{
		if (instr[i]<max/2.0) 
		{
			p2 = i;
			break;
		}
	}

	*width = (p2-p1)*binwidth/(2.0*sqrt(-log(0.5)));




    /* Initial delay parameter estimate... */
    for (i=0, max=0.0,maxi=0; i<nbins; i++)
    {
        if (instr[i]>max)
        {
            max  = instr[i];
            maxi = i;
        }
    }

    *delay = ((double)maxi+0.5)*binwidth;

    /* Initial width parameter estimate... */
	for (i=0; i<nbins; i++)
	{
		if (instr[i] > max/2.0) 
		{
			p1 = i;
			break;
		}
	}

	for (i=nbins-1; i>=0; i--)
	{
		if (instr[i] > max/2.0) 
		{
			p2 = i;
			break;
		}
	}

	*width = (p2-p1)*binwidth / (2.0 * sqrt(-log(0.5)));

    return (0);
}
#else
int bayes_InstrRspCoarseGuessValuesFromLoadedInstr(float *instr,
                                                   float  binwidth,
                                                   int    nbins,
                                                   float *delay,
                                                   float *width)
{
    double p1, p2, max;
    int    i, maxi, quick=1;

    if ((!instr) || (binwidth<=0.0) || (nbins<=0))
        return (-1);

    /* Initial delay parameter estimate... */
    for (i=0, max=0.0,maxi=0; i<nbins; i++)
    {
        if (instr[i] > max)
        {
            max  = instr[i];
            maxi = i;
        }
    }

    *delay = ((float)maxi+0.5f)*binwidth;

    /* Initial width parameter estimate... */
	p1 = 0; // for safety, set to first val
	for (i=0; i<nbins; i++)
	{
		if (instr[i] > max/2.0) 
		{
			p1 = (double)i;
			break;
		}
	}

	p2 = nbins-1; // for safety, set to last val
	for (i=maxi; i<nbins; i++)
	{
		if (instr[i]<max/2.0) 
		{
			p2 = (double)i;
			break;
		}
	}

//	*width = (p2-p1)*binwidth / (2.0 * sqrt(-log(0.5)));
    *width = (float)((p2-p1)*binwidth/BAYES_INSTR_FWHM_TO_SIGMA_FACTOR);

    return (0);
}
#endif


#if 0 //configurable prior - code needs finishing
double MinusLogGaussianPrior(double  x, 
                int     nparams,
                double *params,
                void   *config)
{
    double ave, width, val;

    ave   = params[0];
    width = params[1];

    val = (x-ave)/width;

    return (val*val);
}
#endif

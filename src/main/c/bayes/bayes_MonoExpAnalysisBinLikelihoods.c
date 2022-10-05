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
/* File:       bayes_MonoExpAnalysisBinned.c                                       */
/*---------------------------------------------------------------------------------*/
/* Purpose:    Routines for analysis of binned mono-exponential data.              */
/*---------------------------------------------------------------------------------*/
/* References: [1] - "Bayesian data analysis and parameter inference and           */
/*                   applications to imaging based on photon counting",            */
/*                   ACC Coolen, KCL, December 8th 2006.                           */
/*             --------------------------------------------------------------------*/
/*             [2] - "Numerical Recipes in C++",                                   */
/*                   Second Edition, W. H. Press et al.                            */
/*---------------------------------------------------------------------------------*/
/* Notes:      None.                                                               */
/*=================================================================================*/
/* Revision history:                                                               */
/*---------------------------------------------------------------------------------*/
/* Date   | Modification                                                           */
/*---------------------------------------------------------------------------------*/
/* 090424 | Creation, mrowley.                                                     */
/*---------------------------------------------------------------------------------*/
/*        |                                                                        */
/*=================================================================================*/


/***********************************************************************************/
/*                                                                                 */
/*                             DEFINES / SWITCHES                                  */
/*                                                                                 */
/***********************************************************************************/



/***********************************************************************************/
/*                                                                                 */
/*                             HEADER FILE INCLUSION                               */
/*                                                                                 */
/***********************************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"

#include "extmath.h"
#include "matrices.h"
#include "DTYPE.h"
#include "bayes_Types.h"
#include "bayes_Sizes.h"
#include "bayes_DataManagement.h"
#include "bayes_DistributionFctsBinLikelihoods.h"
#include "bayes_MonoExpAnalysisBinLikelihoods.h"
#include "bayes_ModelTransformTools.h"
#include "bayes_Interface.h"

#ifndef MAX
#define MAX(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)            (((a) < (b)) ? (a) : (b))
#endif

/***********************************************************************************/
/*                                                                                 */
/*                          (LOCAL) GLOBAL VARIABLES                               */
/*                                                                                 */
/***********************************************************************************/

/* None. */


/***********************************************************************************/
/*                                                                                 */
/*                      EXTERNAL GLOBAL VARIABLES ACCESS                           */
/*                                                                                 */
/***********************************************************************************/

/* None. */


/***********************************************************************************/
/*                                                                                 */
/*                         STATIC FUNCTION PROTOTYPES                              */
/*                                                                                 */
/***********************************************************************************/

/* None. */


/***********************************************************************************/
/*                                                                                 */
/*        (BINNED DATA) MONO-EXPONENTIAL ANALYSIS SPECIFIC MODEL FUNCTIONS         */
/*                                                                                 */
/***********************************************************************************/
       
double bayes_MonoExpMinusLogProbW0W1BinLikelihood(double *x, int id, void *container)
{
    int    *data, bin, nbins, fitstart, nphotons, nphotonsbin, ret;
    double  interval, modperiod, *binwalls, bL, bH, bjoverT;
    double  alpha, value, w0, w1, oneminusw0;
    double *likelihoods, w[2], t[2];

    MonoExpMinusLogProbW0W1_t *params1;
    BayesInstrRsp_t           *instr;

    w0 = x[1];
    w1 = x[2];

    if ((w0 < 0.0) || (w0 > 1.0))
        return (BAYES_SIZE_DOUBLE_HUGE);
    
    if (w1 <= 0.0)
        return (BAYES_SIZE_DOUBLE_HUGE);

    params1     = (MonoExpMinusLogProbW0W1_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    fitstart    = params1->fitstart;
    binwalls    = params1->binwalls;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    modperiod   = params1->modulationperiod;
    instr       = params1->instr;
    alpha       = params1->hyperparam;

    value       = alpha*w1;
    oneminusw0  = (1.0-w0);

    likelihoods = Bayes_dvector(0,nbins);

    w[0] = w0; //normalisation
    w[1] = oneminusw0;
    t[1] = w1;
    //bayes_ArrBinLikelihoodsGivenTau(likelihoods,NULL,NULL/*upsilon1*/,data,nbins,interval,modperiod,instr,w1);
    //bayes_ArrTimeBinLikelihoodsGivenTauFixedInstrFixedLifetime(likelihoods,data,nbins,interval,width,delay,w1);
    ret = bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(likelihoods,
                                                                     nbins,binwalls,data,
                                                                     interval,modperiod,instr,
                                                                     w1,
                                                                     //0,NULL,NULL);
                                                                     1,w,t);
    if (ret<0)
    {
        free_Bayes_dvector(likelihoods,0,nbins);
        return (BAYES_SIZE_DOUBLE_HUGE);
    }

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            bL       = binwalls[bin];
            bH       = binwalls[bin+1];
            bjoverT  = (bH-bL)/interval;
            value   -= (double)(nphotonsbin) * (log((w0*bjoverT)+(oneminusw0*likelihoods[bin])));
        }
    }

    free_Bayes_dvector(likelihoods,0,nbins);

    return (value);
}


void bayes_DirectMonoExpMostProbW0W1BinLikelihood(int             *data,
                                                  int              nbins,
                                                  int              fitstart,
                                                  double          *binwalls,
                                                  int              nphotons,
                                                  double          *w0,
                                                  double          *w1,
                                                  double          *val,
                                                  BayesInstrRsp_t *instr,
                                                  double           interval,
                                                  double           modulationperiod,
                                                  double           alpha,
                                                  double           precision)
{
    double value, *x, deltas[3]={0.0,0.05,0.05};
    int    id=0;

    int   (*minimizer)(double (*)(double *, int, void *), int, void *, int, double *, double *, void *);
    void   *config;

    MonoExpMinusLogProbW0W1_t container;

    container.data             = data;
    container.nbins            = nbins;
    container.fitstart         = fitstart;
    container.binwalls         = binwalls;
    container.nphotons         = nphotons;
    container.interval         = interval;
    container.modulationperiod = modulationperiod;
    container.instr            = instr;
    container.hyperparam       = alpha;

    x    = Bayes_dvector(1,2);
    x[1] = (double)(*w0);
    x[2] = (double)(*w1);

	minimizer = &math_MinimiseFctDoubleWithGenericContainer;
    config = malloc(sizeof(AmoebaConfigParams_t));
    ((AmoebaConfigParams_t*)config)->monitor   = 0;
    ((AmoebaConfigParams_t*)config)->tolerance = bayes_MonoExpConfigGetDownhillSimplexPrecision();
    ((AmoebaConfigParams_t*)config)->deltas    = deltas;

	minimizer(bayes_MonoExpMinusLogProbW0W1BinLikelihood,id,(void*)(&container),2,x,&value,(void*)(config));
 
    *w0  = (float)x[1]; 
    *w1  = (float)x[2]; 
    *val = (float)value;

    free_Bayes_dvector(x,1,2);
}


double bayes_MapW0Fast(double  w0,
                       double  w1,
                       double  interval,
                       int    *data,
                       double *likelihood,
                       int     nbins,
                       int     p)
{
    double dum, value, factor;
    int    bin, nphotonsbin;
    int    type;
   
    if (w1 <= 0.0)
        return (1.0);
    
    if ((w0 <= 0.0) || (w0 >= 1.0))
        return (0.0);
   
    factor = interval*(1.0/w0 - 1.0);
    value  = 0.0;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            dum    = 1.0 + factor*likelihood[bin];
            value += nphotonsbin/dum;
        }

        if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value,&type))
	    {
		    bayes_dm_CorrectInvalidDoubleValue(&value,type);
        }
    }

    return (value/(double)p);
}


double bayes_MonoExpMinusLogProbW1BinLikelihood(double *x, int id, void *container)
{
    double  value, alternative, diff, neww0, w0, w1, oldw0, oneminusw0;
    double  dp, interval, modperiod, alpha, *likelihood, *binwalls;
    double *upsilon1;
    int    *data, nbins, nphotons, bin, nphotonsbin, type, ret;

    MonoExpMinusLogProbW0W1_t *params1;
    BayesInstrRsp_t           *instr;

    w1 = x[1];

    //if (w1 <= 0.0)
    if (w1 < SMALLTIME/2.0)
        return (BIG);

    params1     = (MonoExpMinusLogProbW0W1_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    binwalls    = params1->binwalls;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    modperiod   = params1->modulationperiod;
    instr       = params1->instr;
    alpha       = params1->hyperparam;
    upsilon1    = params1->upsilon1;


    likelihood = Bayes_dvector(0,nbins-1);
    //ret = bayes_ArrBinLikelihoodsGivenTau(likelihood,NULL,upsilon1,data,nbins,interval,modperiod,instr,w1);
    ret = bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(likelihood,
                                                                     nbins,binwalls,data,
                                                                     interval,modperiod,instr,
                                                                     w1,
                                                                     0,NULL,NULL);
    
    if (ret<0)
    {
        free_Bayes_dvector(likelihood,0,nbins-1);
        return (BIG); //an error has occured, but give minimizer a chance to recover with the next tested param values...
    } 

    dp       = (double)(nphotons);
    diff     = 1.0;
    w0       = 0.001;
    oldw0    = 0.001;
    diff     = 0.1;

    while (diff > 0.000000001)
    {
        w0 = bayes_MapW0Fast(oldw0,w1,interval,data,likelihood,nbins,nphotons);

        /* Check that w0 is a valid number... */
        if ((w0<0.0) || (w0>1.0) || (BAYES_DM_DOUBLE_TYPE_INVALID==bayes_dm_CheckDoubleValueValid(w0,&type)))
        {
            free_Bayes_dvector(likelihood,0,nbins-1);
            return (BIG);
        }

        diff  = w0-oldw0;

        if (diff<0.0)
            diff = -diff;

        oldw0 = w0;
    }

    /* Computation of [1], Eqn. (88), assuming noisy environment (i.e. 0.0<w0<1.0)... */
    value      = alpha * w1/dp;
    oneminusw0 = 1.0-w0;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            value -= (double)(nphotonsbin) * (log(oneminusw0*interval*likelihood[bin]+w0)/dp);

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value,&type))
		    {
			    bayes_dm_CorrectInvalidDoubleValue(&value,type);
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    /* Computation of [1], Eqn. (88), assuming perfect noiseless environment (i.e. w0 = 0.0)... */
    neww0       = 0.0;
    alternative = alpha * w1/dp;
    
    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            alternative -= (double)(nphotonsbin) * (log(interval*likelihood[bin])/dp);

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value,&type))
		    {
			    alternative = value + 100.0;
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    /* Computation of [1], Eqn. (88), assuming overwhelming noisy environment (i.e. w0 = 1.0)... */
    neww0       = 1.0;
    alternative = alpha*w1/dp;

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    free_Bayes_dvector(likelihood,0,nbins-1);

    params1->w0 = w0;
    
    return (value);
}


int  bayes_IndirectMostProbableW0W1BinLikelihood(int             *data,
                                                 int              nbins,
                                                 double          *binwalls,
                                                 int              nphotons,
                                                 double          *w0,
                                                 double          *w1,
                                                 double          *val,
                                                 BayesInstrRsp_t *instr,
                                                 double           interval,
                                                 double           modulationperiod,
                                                 double           alpha,
                                                 double           precision)
{
    double  value, *x;
    int     id=0, error;
    void    *config;

    MonoExpMinusLogProbW0W1_t container;

    container.data             = data;
    container.nbins            = nbins;
    container.binwalls         = binwalls;
    container.nphotons         = nphotons;
    container.interval         = interval;
    container.modulationperiod = modulationperiod;
    container.instr            = instr;
    container.hyperparam       = alpha;
    container.error            = 0;

    x    = Bayes_dvector(1,1);
    x[1] = (double)(*w1);

    config = malloc(sizeof(AmoebaConfigParams_t));
    ((AmoebaConfigParams_t*)config)->monitor   = 0;
    ((AmoebaConfigParams_t*)config)->tolerance = precision;

    math_MinimiseFctDoubleWithGenericContainer
        (bayes_MonoExpMinusLogProbW1BinLikelihood,id,(void*)(&container),1,x,&value,(void*)(config));

    error = container.error;

    if (error >= 0)
    {
        math_MinimiseFctDoubleWithGenericContainer
            (bayes_MonoExpMinusLogProbW1BinLikelihood,id,(void*)(&container),1,x,&value,(void*)(config));
    }

    error = container.error;

    if (error < 0)
    {
        free_Bayes_dvector(x,1,1);
        free(config);

        return (-1);
    }

    *w1  = (float)x[1];
    *w0  = (float)container.w0;
    *val = (float)value;
   
    free_Bayes_dvector(x,1,1);
    free(config);

    return (0);
}


double bayes_MonoExpProbW0W1BinLikelihoodFixedInstrFixedLifetime(double   w0,
	     				 		                                 double   w1,
                                                                 double  *likelihoods,
							                                     int     *data,
							                                     int      nbins,
                                                                 double  *binwalls,
							                                     double   interval,
							                                     double   alpha,
							                                     double   norm,
							                                     int     *result)
{
    int     bin, nphotonsbin;
    double  value, oneminusw0, bL, bH, bjoverT;
    int     type;

    if ((w0 < 0.0) || (w0 > 1.0) || (w1 <= 0.0))
        return (TINY);

    oneminusw0 = (1.0-w0);
    value      = -alpha*w1;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            bL       = binwalls[bin];
            bH       = binwalls[bin+1];
            bjoverT  = (bH-bL)/interval;
            value   += (double)(nphotonsbin) * log((w0*bjoverT)+(oneminusw0*likelihoods[bin]));      
        }
    }

    value = exp(value+norm);

    if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value,&type))
    {
	    bayes_dm_CorrectInvalidDoubleValue(&value,type);
        *result = RESULT_INVALID_ARITHMETIC_ERROR;
    }
    else
	{
        *result = RESULT_VALID;
	}

    return (value);
}


double bayes_MonoExpProbW0W1BinLikelihood(double           w0,
	 				 		              double           w1,
							              int             *data,
							              int              nbins,
                                          double          *binwalls,
							              double           interval,
							              BayesInstrRsp_t *instr,
                                          double           alpha,
							              double           norm,
							              int             *result,
                                          double          *upsilon1_g,
                                          double           modperiod)
{
    int     bin, nphotonsbin;
    double *likelihoods, bL, bH, bjoverT;
    double  value, oneminusw0;
    int     type;

    if ((w0 < 0.0) || (w0 > 1.0) || (w1 <= 0.0))
        return (TINY);

    oneminusw0 = (1.0-w0);
    value      = -alpha*w1;

    likelihoods = Bayes_dvector(0,nbins-1);
    //bayes_ArrBinLikelihoodsGivenTau(likelihoods,upsilon1,data,nbins,interval,width,delay,w1,modperiod);
    //bayes_ArrBinLikelihoodsGivenTau(likelihoods,NULL,NULL,data,nbins,interval,modperiod,instr,w1);

    bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(likelihoods,
                                                               nbins,binwalls,data,
                                                               interval,modperiod,instr,
                                                               w1,
                                                               0,NULL,NULL);
   
    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            bL       = binwalls[bin];
            bH       = binwalls[bin+1];
            bjoverT  = (bH-bL)/interval;
            value   += (double)(nphotonsbin) * (log((w0*bjoverT)+(oneminusw0*likelihoods[bin])));      
        }
    }

    value = exp(value+norm);

    if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
    {
	    bayes_dm_CorrectInvalidDoubleValue(&value, type);
        *result = RESULT_INVALID_ARITHMETIC_ERROR;
    }
    else
	{
        *result = RESULT_VALID;
	}

    free_Bayes_dvector(likelihoods,0,nbins-1);

    return (value);
}




int bayes_AveragesAndErrorBarsBinLikelihood(int                *data,
                                            int                 nbins,
                                            int                 fitstart,
                                            double             *binwalls,
                                            int                 nphotons,
		                                    double             *w0_mp,
			                                double             *w1_mp,
		                                    double             *w0_ave,
			                                double             *w1_ave,
                                            double             *dw0,
			                                double             *dw1,
		                                    BayesInstrRsp_t    *instr,
                                            float               interval,
                                            float               modulationperiod,
                                            float               alpha,
			                                float               precision,
			                                BayesAveErrDistn_t *distr,
                                            float              *minuslogprob)
{
	// Calculates the average of the pdf and 'error bars' from it's width
	// The starting point is the most probable (mp) values
	// Then it roughly finds the range fo the data
	// Then calculates the average and width from a 2D grid of values
	
	
	int    ret = BAYES_AVE_ERRS_ROUTINE_NO_ERRORS;
	double average, error, value, temp, min, dxdy, val;
	double *likelihoods;
    double norm, w0min, w0max, w1min, w1max, X[3], **Prob;
	double *w0values, *w1values, marginal, *marginals, dx, dy, range, maxprob, newprob, neww0, neww1;
    int    noise, hi, lo;
	int    i, j, xpts, ypts, manualrange, maxCount=1000;
    int    id=0;
    int    result, quick=0;
    int    type;

    MonoExpMinusLogProbW0W1_t container1;
    ArrLikelihoodConstants_t  constants;

#if 0
    if (constants_g)
    {
        upsilon1 = constants_g->upsilon1;
        constants.upsilon1 = upsilon1;
    }
    else
    {
        /* Compute constant quantites depending only on instrument response parameters and number of bins... */
        upsilon1 = Bayes_dvector(0,nbins);
        bayes_ComputeArrBinLikelihoodConstantUpsilon1(&upsilon1,data,nbins,interval,width,delay);
        constants.upsilon1 = upsilon1;
    }
#else
    constants.upsilon1 = NULL;
#endif


    /* Determine most-probable values as the starting point               */
    /* for computation of all distribution parameter statistics...        */
    //bayes_DirectMonoExpMostProbW0W1BinLikelihood - could make direct/inderect selectable from UI...
    //bayes_IndirectMostProbableW0W1BinLikelihood

    bayes_DirectMonoExpMostProbW0W1BinLikelihood    
        (data,nbins,fitstart,binwalls,nphotons,w0_mp,w1_mp,&val,instr,interval,modulationperiod,alpha,precision);

    if ((BAYES_DM_DOUBLE_TYPE_VALID != bayes_dm_CheckDoubleValueValid(*w0_mp,&type)) ||
        (BAYES_DM_DOUBLE_TYPE_VALID != bayes_dm_CheckDoubleValueValid(*w1_mp,&type)) ||
        (*w0_mp<0.0) || (*w0_mp>1.0) || (*w1_mp<0.0))
    {
        *w0_mp  = -1.0;
        *w1_mp  = -1.0;

        if (w0_ave) *w0_ave = -1.0;
        if (w1_ave) *w1_ave = -1.0;
        if (dw0)    *dw0    = -1.0;
        if (dw1)    *dw1    = -1.0;
        
        *minuslogprob = (float)BIG;

        return (BAYES_AVE_ERRS_ROUTINE_ERROR);
    }

    X[1] = *w0_mp;
	X[2] = *w1_mp;
    *minuslogprob = (float)val;

    /* ...now determine the 'normalization' constant for subsequent probability computations... */
    container1.data             = data;
    container1.nbins            = nbins;
    container1.binwalls         = binwalls;
    container1.nphotons         = nphotons;
    container1.interval         = interval;
    container1.modulationperiod = modulationperiod;
    container1.instr            = instr;
    container1.hyperparam       = alpha;
    container1.upsilon1         = /*upsilon1*/NULL;

    min     = bayes_MonoExpMinusLogProbW0W1BinLikelihood(X,id,(void *)(&container1));
	maxprob = bayes_MonoExpProbW0W1BinLikelihood(X[1],X[2],data,nbins,binwalls,interval,instr,alpha,min,&result,NULL,modulationperiod);

    /* Now start analysis of the distribution and generate its statistics... */
	if (distr)
	{
		xpts = distr->m;
		ypts = distr->n;
	}
	else
	{
    	xpts = MAX((int)(0.5+0.001/precision),1);
		ypts = xpts;
	}
	
    Prob     = Bayes_dmatrix(0,xpts-1,0,ypts-1);
    w0values = Bayes_dvector(0,xpts-1);
	w1values = Bayes_dvector(0,ypts-1);

    neww0    = X[1];
	neww1    = X[2];
    
    /* Estimate the w0-width around the most likely point (w0,w1) */
    range = 0.1;

    likelihoods = Bayes_dvector(0,nbins-1);
    //bayes_ArrBinLikelihoodsGivenTau(likelihoods,NULL,/*upsilon1*/NULL,data,nbins,interval,modulationperiod,instr,neww1);
    bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(likelihoods,
                                                               nbins,binwalls,data,
                                                               interval,modulationperiod,instr,
                                                               neww1,
                                                               0,NULL,NULL);

    manualrange = (distr) ? (distr->manualrange) : (0);

    if (manualrange)
    {
        w0min = distr->xmin;
        w0max = distr->xmax;
        w1min = distr->ymin;
        w1max = distr->ymax;        
    }
    else
    {
        lo=0;
        hi=1;

        noise = (X[1]<0.5)?(lo):(hi);

        while (range > 0.00005) /* Determine the width of the background parameter marginal distribution... */
        {
            neww0 += (lo==noise)?(range):(-range);
    		
	        if (((neww0 < 0.0)) || (neww0 > 1.0))
            {
                quick = 1;
                ret   = BAYES_AVE_ERRS_ROUTINE_OUTOFRANGE;
                break;
            }

            newprob = bayes_MonoExpProbW0W1BinLikelihoodFixedInstrFixedLifetime(neww0,neww1,likelihoods,data,nbins,binwalls,interval,alpha,min,&result);
          
            if (newprob < (0.4*maxprob))
	        {
		        neww0 -= (lo==noise)?(range):(-range);
		        range *= 0.1;
	        }
			else if (newprob >= maxprob)
			{
				range *= 1.5;
			}

			if (range > DBL_MAX)
				break;

        }

        neww0 += (lo==noise)?(5.0*range):(-5.0*range);
        range  = (lo==noise)?(neww0-X[1]):(X[1]-neww0);

        w0min  = X[1]-6.0*range;
	    w0max  = X[1]+6.0*range;

        if (w0min < 0.0)
	        w0min = 0.0;
    	
        if (w0max > 1.0)
	        w0max = 1.0;

        if (dw0) *dw0 = (float)range;    

        /* Estimate the w1-width around the most likely point (w0,w1) */
        neww0 = X[1];
        range = 0.1;

	    while (range > 0.00005) /* Detemine 'w1' distribution width limits... */
	    {
            neww1   += range;

            newprob = bayes_MonoExpProbW0W1BinLikelihood
                          (neww0,neww1,data,nbins,binwalls,interval,instr,alpha,min,&result,/*upsilon1*/NULL,modulationperiod);

		    if (newprob < (0.4*maxprob))
		    {
			    neww1 -= range;
			    range *= 0.1;
		    }
			else if (newprob >= maxprob)
			{
				range *= 10.0;
			}

			if (range > DBL_MAX)
				break;
        }

	    neww1 += (5.0*range);
	    range  = neww1-X[2];
        w1min  = X[2]-6.0*range;
	    w1max  = X[2]+6.0*range;
       
	    if (w1min < 0.0)
		    w1min = SMALLTIME;
       
	    if (dw1) *dw1 = (float)range;
    }

    if (quick)
		return (ret);
   
    /* Now do the proper calculation of average and error bars */
	/* Fill 2D array with p(w|D) values... */
    dx   = (w0max - w0min) / (xpts - 1.0);
	dy   = (w1max - w1min) / (ypts - 1.0);
	dxdy = dx*dy;
	norm = 0.0;

	for (i=0; i<xpts; i++)
		w0values[i] = w0min+(double)i*dx;
	
	for (j=0; j<ypts; j++)
		w1values[j] = w1min+(double)j*dy;


	for (j=0; j<ypts; j++)   /* Run through the lifetime values in outer loop...   */
    {                         /* ...so that inner loop is over a constant lifetime. */
        bayes_ArrBinLikelihoodsGivenTau(likelihoods,NULL,NULL/*upsilon1*/,data,nbins,interval,modulationperiod,instr,w1values[j]);

        for (i=0; i<xpts; i++)
		{
            value = bayes_MonoExpProbW0W1BinLikelihoodFixedInstrFixedLifetime(
                        w0values[i],w1values[j],likelihoods,data,nbins,binwalls,interval,alpha,min,&result);

            Prob[i][j] = value;
			norm      += value;
        }
    }

    norm *= dxdy;

	for (i=0; i<xpts; i++){
//		printf("\n");
		for (j=0; j<ypts; j++){
			Prob[i][j] = Prob[i][j]/norm;
//			printf("%.1f ", Prob[i][j]);
		}
	}

    if (distr)//could be done in reverse, and then use distr values later
	{
        distr->x = Bayes_vector(0,xpts-1);
        distr->y = Bayes_vector(0,ypts-1);
        distr->z = Bayes_matrix(0,xpts-1,0,ypts-1);
        
        for (i=0; i<xpts; i++)
		{
            (distr->x)[i] = (float)(w0values[i]);
            
            for (j=0; j<ypts; j++)
			{
				(distr->y)[j]    = (float)(w1values[j]);
                (distr->z)[i][j] = (float)(Prob[i][j]);
			}
        }
    }

   // Get the averages and errors
    marginals = Bayes_dvector(0,xpts-1);

    average = 0.0;

	for (i=0; i<xpts; i++)
	{
        marginal = 0.0;
		
		for (j=0; j<ypts; j++)
			marginal += Prob[i][j];
      
        marginal     *= dy;
        average      += w0values[i]*marginal;
        marginals[i]  = marginal;
    }
    
	average *= dx;
	error    = 0.0;

    for (i=0; i<xpts; i++)
	{
        temp   = w0values[i]-average;
        error += temp*temp*marginals[i];
    }
    
	error  *= dx;
	error   = sqrt(error);
	if (w0_ave) *w0_ave = (float)average;
	if (dw0)    *dw0    = (float)error;
	
	free_Bayes_dvector(marginals,0,xpts-1);
	marginals = Bayes_dvector(0,ypts-1);

    average = 0.0;

    for (j=0; j<ypts; j++)
	{
        marginal = 0.0;
		
		for (i=0; i<xpts; i++)
			marginal += Prob[i][j];
      
		marginal     *= dx;
		average      += w1values[j]*marginal;
        marginals[j]  = marginal;
    }
    
	average *= dy;
	error    = 0.0;

    for (j=0; j<ypts; j++)
    {
        temp   = w1values[j]-average;
        error += temp*temp*marginals[j];
    }

	error  *= dy;
	error   = sqrt(error);
	if (w1_ave) *w1_ave = (float)average;
	if (dw1)    *dw1    = (float)error;

    free_Bayes_dvector(marginals,0,ypts-1);
    free_Bayes_dvector(likelihoods,0,nbins);

	free_Bayes_dmatrix(Prob,0,xpts-1,0,xpts-1);
    free_Bayes_dvector(w0values,0,xpts-1);
	free_Bayes_dvector(w1values,0,ypts-1);

	return (ret);
}





//free instrument





#if 0
int bayes_DirectMonoExpMostProbW0W1WithInstrRspOptimizationBinLikelihood(int   *data,
                                                                         int    nbins,
					                                                     int    nphotons,
   				                                                         float *w0,
					                                                     float *w1,
					                                                     //float *delay,
					                                                     //float *width,
					                                                     BayesInstrRsp_t *instr,
                                                                         float  interval,
                                                                         float  modulationperiod,
					                                                     float  alpha,
					                                                     float  precision)
{
    double value, *x;
    int    id=0;
    void   *config;

    MonoExpMinusLogProbW0W1_t container;

    container.data             = data;
    container.nbins            = nbins;
    container.nphotons         = nphotons;
    container.interval         = interval;
    container.modulationperiod = modulationperiod;
    container.hyperparam       = alpha;
	container.upsilon1         = NULL;

    x    = Bayes_dvector(1,4);
    x[1] = (double)(*w0);
    x[2] = (double)(*w1);
	//x[3] = (double)(*delay);
	//x[4] = (double)(*width);
    x[3] = instr->params[0].delay;
    x[4] = instr->params[0].width;

    (AmoebaConfigParams_t*)config = (AmoebaConfigParams_t*)malloc(sizeof(AmoebaConfigParams_t));
    ((AmoebaConfigParams_t*)config)->monitor   = 0;
    ((AmoebaConfigParams_t*)config)->tolerance = precision;

    math_MinimiseFctDoubleWithGenericContainer
        (bayes_MonoExpMinusLogProbW0W1WithInstrRspOptimizationBinLikelihood,id,(void*)(&container),4,x,&value,(void*)(config));
    
    math_MinimiseFctDoubleWithGenericContainer
        (bayes_MonoExpMinusLogProbW0W1WithInstrRspOptimizationBinLikelihood,id,(void*)(&container),4,x,&value,(void*)(config));

    *w0    = (float)x[1]; 
    *w1    = (float)x[2];
	*delay = (float)x[3];
	*width = (float)x[4];

    free_Bayes_dvector(x,1,4);
	
	return (0);
}
#endif
///very optimistic!!!!!

















#if 0
int bayes_IndirectMostProbableW0W1WithInstrRspOptimizationBinLikelihood(int   *data,
						                                                int    nbins,
						                                                int    nphotons,
						                                                float *w0,
						                                                float *w1,
						                                                float *delay,
						                                                float *width,
						                                                float  interval,
                                                                        float  modulationperiod,
						                                                float  alpha,
						                                                float  precision)
{
    double  value, *x;
    int     id=0;
    void    *config;

    MonoExpMinusLogProbW0W1_t container;

    container.data             = data;
    container.nbins            = nbins;
    container.nphotons         = nphotons;
    container.interval         = interval;
    container.modulationperiod = modulationperiod;
    container.width            = -1.0;
    container.delay            = -1.0;
    container.hyperparam       = alpha;
	container.upsilon1         = NULL;

    x    = Bayes_dvector(1,3);
    x[1] = (double)(*w1);
	x[2] = (double)(*delay);
	x[3] = (double)(*width);

    (AmoebaConfigParams_t*)config = (AmoebaConfigParams_t*)malloc(sizeof(AmoebaConfigParams_t));
    ((AmoebaConfigParams_t*)config)->monitor   = 0;
    ((AmoebaConfigParams_t*)config)->tolerance = precision;

    math_MinimiseFctDoubleWithGenericContainer
        (bayes_MonoExpMinusLogProbW1WithInstrRspOptimizationBinLikelihood,id,(void*)(&container),3,x,&value,(void*)(config));

    math_MinimiseFctDoubleWithGenericContainer
        (bayes_MonoExpMinusLogProbW1WithInstrRspOptimizationBinLikelihood,id,(void*)(&container),3,x,&value,(void*)(config));

    *w1    = (float)x[1];
    *w0    = (float)container.w0;
	*delay = (float)x[2];
	*width = (float)x[3];
   
    free_Bayes_dvector(x,1,3);
	
	return (0);
}
#endif



/***********************************************************************************/
/*                                                                                 */
/*                         WRAPPERS FOR A COMMON INTERFACE                         */
/*                                                                                 */
/***********************************************************************************/


#if 0
int direct_most_probable_wrapped_bin_likelihood (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, BayesAveErrDistn_t *distr, float *val, float *params, float *fitted, float *residuals)
{
    float timeunit;
    int   i;

    if (w0_ave) *w0_ave = 0.0;
    if (w1_ave) *w1_ave = 0.0;
	if (dw0)    *dw0    = 0.0;
	if (dw1)    *dw1    = 0.0;
	
    bayes_DirectMonoExpMostProbW0W1BinLikelihood
        (data,nbins,nphotons,w0_mp,w1_mp,val,p,*delay,*width,interval,alpha,precision,NULL);

    if (params)
    {
		params[0] = bayes_ComputeRawBackgroundZ(*w0_mp,interval/(float)(nbins),interval,nphotons);            /* Z   */
		params[1] = bayes_ComputeRawAmplitudeA(*w0_mp,*w1_mp,interval/(float)nbins,interval,*delay,nphotons); /* A   */
        params[2] = *w1_mp;                                                                                   /* tau */
    }
#if 0
    if (fitted)
    {
        timeunit = interval/(float)nbins;

        for (i=0; i<nbins; i++)
        {
            fitted[i] = bayes_ModelTransformComputePhotonCount
                            (*delay,*width,interval,timeunit*(float)i,*w0_mp,*w1_mp,timeunit,nphotons);
        }
    }

    if (residuals)
    {
        for (i=0; i<nbins; i++)
        {
            residuals[i] = (float)(data[i])-fitted[i];
        }

        scaleDataAccordingToSignalNoise(residuals,nbins,fitted);    
    }
#endif
    return (0);
}




int indirect_most_probable_wrapped_bin_likelihood (float *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, BayesAveErrDistn_t *distr, float *val/*, float *params, float *fitted, float *residuals*/)
{
    float timeunit;
    int   i, *idata;

    if (w0_ave) *w0_ave = 0.0;
    if (w1_ave) *w1_ave = 0.0;
	if (dw0)    *dw0    = 0.0;
	if (dw1)    *dw1    = 0.0;
	
	idata = Bayes_ivector(0,nbins-1);
	
	for (i=0; i<nbins; i++)
		idata[i] = (int)(data[i]);

    bayes_IndirectMostProbableW0W1BinLikelihood
        (idata,nbins,nphotons,w0_mp,w1_mp,val,p,*delay,*width,interval,alpha,precision,NULL);
#if 0
    if (params)
    {
		params[0] = bayes_ComputeRawBackgroundZ(*w0_mp,interval/(float)(nbins),interval,nphotons);            /* Z   */
		params[1] = bayes_ComputeRawAmplitudeA(*w0_mp,*w1_mp,interval/(float)nbins,interval,*delay,nphotons); /* A   */
        params[2] = *w1_mp;                                                                                   /* tau */
    }

    if (fitted)
    {
        timeunit = interval/(float)nbins;

        for (i=0; i<nbins; i++)
        {
            fitted[i] = bayes_ModelTransformComputePhotonCount
                            (*delay,*width,interval,timeunit*(float)i,*w0_mp,*w1_mp,timeunit,nphotons);    
        }
    }

    if (residuals)
    {
        for (i=0; i<nbins; i++)
        {
            residuals[i] = (float)(data[i])-fitted[i];
        }

        scaleDataAccordingToSignalNoise(residuals,nbins,fitted);    
    }
#endif
    return 0;
}




#if 0
double bayes_MonoExpMinusLogProbW1BinnedDataFastBinLikelihood(double *x, int id, void *container)
{
    double  value, alternative, diff, neww0, w0, w1, oldw0;
    double  dp, interval, width, delay, alpha, *likelihood;
    int    *data, nbins, nphotons, bin, nphotonsbin, type;

    MonoExpMinusLogProbW0W1BinnedData_t *params1;

    w1 = x[1];

    if (w1 <= 0.0)
        return (BIG);

    params1     = (MonoExpMinusLogProbW0W1BinnedData_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    width       = params1->width;
    delay       = params1->delay;
    alpha       = params1->hyperparam;

    likelihood = dvector(0,nbins-1);

    bayes_ArrTimeBinLikelihoodsGivenTauFixedInstrFixedLifetime(likelihood,data,nbins,interval,width,delay,w1);

    dp       = (double)(nphotons);
    diff     = 1.0;
    w0       = 0.001;
    oldw0    = 0.001;
    diff     = 0.1;

    while (diff > 0.000000001)
    {
        w0 = bayes_MapW0BinnedDataFast(oldw0, w1, interval, data, likelihood, nbins, nphotons);

		/* Check that w0 is a valid number... */
		bayes_dm_CheckDoubleValueValid(w0, &type);

	    if (type == BAYES_DM_DTYPE_NOT_A_NUMBER)
            w0 = 1.0;

        diff  = sqrt((w0-oldw0)*(w0-oldw0));
        oldw0 = w0;
    }

    /* Computation of [1], Eqn. (88), assuming noisy environment (i.e. 0.0<w0<1.0)... */
    value = alpha * w1/dp;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            value       -= (double)(nphotonsbin) * (log((1.0-w0) * interval * likelihood[bin] + w0) / dp);

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    bayes_dm_CorrectInvalidDoubleValue(&value, type);
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    /* Computation of [1], Eqn. (88), assuming perfect noiseless environment (i.e. w0 = 0.0)... */
    neww0       = 0.0;
    alternative = alpha * w1/dp;
    
    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            alternative -= (double)(nphotonsbin) * (log(interval * likelihood[bin]) / dp);

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    alternative = value + 100.0;
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    /* Computation of [1], Eqn. (88), assuming overwhelming noisy environment (i.e. w0 = 1.0)... */
    neww0       = 1.0;
    alternative = alpha*w1/dp;

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    free_dvector(likelihood,0,nbins-1);

    params1->w0 = w0;
    
    return (value);
}
#endif




#if 0
int most_probable_freeinstrument_wrapped_binned (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, int distr, float *val, float *params, float *fitted, float *residuals)
{
    float timeunit;
    int   i;
    int   trial, trials;
    float bestw0, bestw1, w0start, w1start, bestwidth, bestdelay, bestval;

    if (w0_ave) *w0_ave = 0.0;
    if (w1_ave) *w1_ave = 0.0;
	if (dw0)    *dw0    = 0.0;
	if (dw1)    *dw1    = 0.0;
#if 1
//original	bayes_MostProbableFreeInstrument(data, nbins, nphotons, w0_mp, w1_mp, val, p, delay, width, interval, alpha, precision);
	
    *width = 0.0;
    *delay = 1.8;
    bayes_MostProbableFreeInstrDelayNoSpread(data, nbins, nphotons, w0_mp, w1_mp, val, p, delay, width, interval, alpha, precision);
#else

    /* Don't really know where to start from for width and delay estimates...          */
    /* ...and don't want to find a local minimum due to poor initial search parameters */
    /* So, after the first attempt starting at calling routines suggested location...  */
    /* ...introduce a stochastic element.                                              */

    w0start = *w0_mp;
    w1start = *w1_mp;

    bayes_MostProbableFreeInstrument(data, nbins, nphotons, w0_mp, w1_mp, val, p, delay, width, interval, alpha, precision);

    bestval   = *val;
    bestwidth = *width;
    bestdelay = *delay;
    bestw0    = *w0_mp;
    bestw1    = *w1_mp;

    trials = 100;

    if (trials > 1)
        SETRANDOM();
    
    trial = 0;

    *delay = 0.0;
    *width = 0.0;
    bayes_MostProbableFreeInstrument(data, nbins, nphotons, w0_mp, w1_mp, val, p, delay, width, interval, alpha, precision);
    printf("\nNo response: val=%g, width=%g, delay=%g", *val, *width, *delay);

    while (trial<trials)
    {
        *delay = (interval/2.0)*RANDOM();
        *width = 2.0*RANDOM();
        *w0_mp  = w0start;
        *w1_mp  = w1start;

        bayes_MostProbableFreeInstrument(data, nbins, nphotons, w0_mp, w1_mp, val, p, delay, width, interval, alpha, precision);

        if (*val<bestval)
        {
            printf("\nBetter (%d): val=%g, width=%g, delay=%g", trial, *val, *width, *delay);
            bestval   = *val;
            bestwidth = *width;
            bestdelay = *delay;
            bestw0    = *w0_mp;
            bestw1    = *w1_mp;
        }

        trial++;
    }

    *width = bestwidth;
    *delay = bestdelay;
    *val   = bestval;
    *w0_mp = bestw0;
    *w1_mp = bestw1;

#endif
    if (params)
    {
		params[0] = bayes_ComputeRawBackgroundZ(*w0_mp, interval/(float)(nbins), interval, nphotons);                 /* Z   */
		params[1] = bayes_ComputeRawAmplitudeA(*w0_mp, *w1_mp, interval/(float)nbins, interval, *delay, nphotons);    /* A   */
        params[2] = *w1_mp;                                                                                           /* tau */
    }

    if (fitted)
    {
        timeunit = interval/(float)nbins;

        for (i=0; i<nbins; i++)
        {
            fitted[i] = bayes_ModelTransformComputePhotonCount
                            (*delay,*width,interval,timeunit*(float)i,*w0_mp,*w1_mp,timeunit,nphotons);    
        }
    }

    if (residuals)
    {
        for (i=0; i<nbins; i++)
            residuals[i] = (float)(data[i])-fitted[i];

        scaleDataAccordingToSignalNoise(residuals,nbins,fitted);    
    }
	
    return (0);
}
#endif



int averages_and_errorbars_wrapped_bin_likelihood (float *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, BayesAveErrDistn_t *distr, float *val/*, float *params, float *fitted, float *residuals*/)
{
    float  timeunit, avg;
    double *likelihoods;
    int    *datadummy;
    int    i, ret, *idata;

	if (val) *val=0;
	
	idata = Bayes_ivector(0,nbins-1);
	
	for (i=0; i<nbins; i++)
		idata[i] = (int)(data[i]);
	
	if (1)
	{
		FILE *fp;

		fp = fopen("BayesFormatHisto", "w");
		if (fp == NULL) return -6;

		fprintf(fp, "%d\n\n", nbins);

		for (i=0; i<nbins; i++)
		{
			fprintf(fp, "%d\n", idata[i]);
		}

		fclose(fp);	
	}

	avg = data_ComputeBinnedDataAverageArrTime(idata,nbins,nphotons,interval);
	
	alpha = 1.0/avg;

    ret = bayes_AveragesAndErrorBarsBinLikelihood(idata, nbins, nphotons,
                                                  w0_mp, w1_mp, w0_ave, w1_ave, dw0, dw1,
                                                  *delay, *width, interval, alpha, precision,
                                                  p, quick, distr, NULL);

    if (-1 == ret) /* Too many points to determine error bars... */
    {
        *w0_ave = *w0_mp;
        *w1_ave = *w1_mp;
        *dw0    = 0.0;
        *dw1    = 0.0;
    }
#if 0
    if (params)
    {
		params[0] = bayes_ComputeRawBackgroundZ(*w0_ave,interval/(float)(nbins),interval,nphotons);             /* Z   */
		params[1] = bayes_ComputeRawAmplitudeA(*w0_ave,*w1_ave,interval/(float)nbins,interval,*delay,nphotons); /* A   */
        params[2] = *w1_ave;                                                                                    /* tau */

		params[3] = bayes_ComputeRawBackgroundZ(*w0_mp,interval/(float)(nbins),interval,nphotons);            /* Z (MP)   */
		params[4] = bayes_ComputeRawAmplitudeA(*w0_mp,*w1_mp,interval/(float)nbins,interval,*delay,nphotons); /* A (MP)   */
        params[5] = *w1_mp;                                                                                   /* tau (MP) */
    }

    if (fitted)
    {
        timeunit    = interval/(float)nbins;
        likelihoods = Bayes_dvector(0,nbins);
        datadummy   = Bayes_ivector(0,nbins);

        for (i=0; i<nbins; i++)
            datadummy[i] = 1;

        bayes_ArrBinLikelihoodsGivenTau(likelihoods,NULL,NULL,datadummy,nbins,interval,*width,*delay,*w1_mp);

        for (i=0; i<nbins; i++)
            fitted[i] = bayes_ModelTransformComputePhotonCount
                            (interval,*w0_mp,(float)(likelihoods[i]),timeunit,nphotons);    
    }

    if (residuals)
    {
        for (i=0; i<nbins; i++)
        {
            residuals[i] = (float)(data[i])-fitted[i];
        }

        scaleDataAccordingToSignalNoise(residuals,nbins,fitted);    
    }
#endif
    return (ret);
}
#endif


#if 0
int bayes_MonoExpDetermineMostProbAveErrParamValues(float                    *data,
                                                    int                       nbins,
                                                    int                       fitstart,
                                                    int                       fitend,
                                                    int                       nphotons,
                                                    ArrLikelihoodConstants_t *constants,
                                                    int                       nparams,
                                                    int                      *param_free,
                                                    float                    *param_mp,//w1,tau1,w2,tau2
                                                    float                    *param_ave,
                                                    float                    *param_err,
                                                    //float                    *delay,
                                                    //float                    *width,
                                                    BayesInstrRsp_t           *instr,
                                                    float                     interval,
                                                    float                     modulationperiod,
                                                    float                     alpha,
                                                    float                     precision,
                                                    int                       quick,
                                                    BayesAveErrDistn_t       *distr,
                                                    float                    *val)
{
    int    i, *idata, ret;

	if (val)
        *val=0;

    idata = Bayes_ivector(0,fitend-fitstart-1);

    for (i=0; i<fitend-fitstart; i++)
        idata[i] = (int)(data[fitstart+i]);

    ret = bayes_AveragesAndErrorBarsBinLikelihood(idata,fitend-fitstart,nphotons,
                                                  &param_mp[1],&param_mp[2],
                                                  &param_ave[1],&param_ave[2],
                                                  &param_err[1],&param_err[2],
                                                  /**delay,*width*/instr,
                                                  interval*(((float)fitend-(float)fitstart+1.0)/(float)nbins),
                                                  modulationperiod,
                                                  alpha,precision,
                                                  nphotons,quick,distr,NULL);

    if (-1 == ret) /* Too many points to determine error bars... */
    {
        param_ave[1] = param_mp[1];
        param_ave[2] = param_mp[2];
        param_err[1] = 0.0;
        param_err[2] = 0.0;
    }

    return (ret);
}
#endif





////////////////////////////////////////////////////////////////////
///////////  old code from here onwards ////////////////////////////
////////////////////////////////////////////////////////////////////


//fast indirect binned likelihood implementation, with an init fct

//container now includes pre-calculated erf values where lifetime is not an argument
//egn 6, term 3, an array
//eqn 7, the erf sum, just one constant value

#if 0
void bayes_MonoExpMinusLogProbW1FixedInstrBinLikelihoodPreCalcInit(double **const1,/* array of length nbins+1 */
                                                            double  *const2,
                                                            int     *data,
                                                            int      nbins,
                                                            double   interval,
                                                            double   width,
                                                            double   delay)
{
    double binwidth, oneoverwidthroottwo;
    int    bin, valid;

    if (width>0.0)
    {
        binwidth            = interval/(double)nbins;
        oneoverwidthroottwo = 1.0/(width*ROOTTWO);

        /* Lifetime independent bin endpoint values (occupied bins only)... */
        for (bin=0,valid=0; bin<nbins; bin++)
        {
            if (data[bin])
            {
                if (valid)
                {
                    (*const1)[bin+1] = erf((delay-(bin+1)*binwidth)*oneoverwidthroottwo);
                    valid            = 1;
                }
                else
                {
                    (*const1)[bin]   = erf((delay-bin*binwidth)*oneoverwidthroottwo);
                    (*const1)[bin+1] = erf((delay-(bin+1)*binwidth)*oneoverwidthroottwo);
                    valid            = 1;
                }
            }
            else
            {
                valid = 0;
            }
        }

        /* Lifetime independent term in the denominator (i.e. normalisation constant)... */
        *const2 = erf((interval-delay)*oneoverwidthroottwo)+erf(delay*oneoverwidthroottwo);    
    }
}
#endif



#if 0
double bayes_MonoExpMinusLogProbW0W1BinnedData(double *x, int id, void *container)
{
    int     *data, bin, nbins, nphotons, nphotonsbin;
    double  arrivaltime, halfbinsize, interval, width, delay;
    double  alpha, value, w0, w1, dp, oneminusw0T;
    int     type;

    double  binmax, binmin, binsize;

    MonoExpMinusLogProbW0W1BinnedData_t *params1;

    w0 = x[1];
    w1 = x[2];

    if ((w0 < 0.0) || (w0 > 1.0))
        return (BIG);
    
    if (w1 <= 0.0)
        return (BIG);

    params1     = (MonoExpMinusLogProbW0W1BinnedData_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    width       = params1->width;
    delay       = params1->delay;
    alpha       = params1->hyperparam;
    halfbinsize = interval/(double)(nbins<<1);

    binsize     = interval/(double)nbins;

    dp          = (double)nphotons;
    value       = alpha*w1/dp;
    oneminusw0T = (1.0-w0)*interval;

    bayes_ArrTimeLikelihoodGivenTauFixedInstrumentFixedLifetimeUpdate(w1,interval,width,delay);

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            //value       -= (double)(nphotonsbin) * (log(oneminusw0T*bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay)+w0)/dp);
            value       -= (double)(nphotonsbin) * (log(oneminusw0T*bayes_ArrTimeLikelihoodGivenTauFixedInstrumentFixedLifetime(arrivaltime,w1,interval,width,delay)+w0)/dp);

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    bayes_dm_CorrectInvalidDoubleValue(&value, type);
			    return (value);
		    }
        }
    }
  
    return (value);
}
#endif


#if 0
double bayes_MonoExpMinusLogProbW0W1BinnedDataOverflowCorrected(double *x, int id, void *container)
{
    int     *data, bin, nbins, nphotons, nphotonsbin;
    double  arrivaltime, halfbinsize, interval, width, delay;
    double  alpha, value, w0, w1, oneminusw0T;
    int     type;

    MonoExpMinusLogProbW0W1BinnedData_t *params1;

    w0 = x[1];
    w1 = x[2];

    if ((w0 < 0.0) || (w0 > 1.0))
        return (BIG);
    
    if (w1 <= 0.0)
        return (BIG);

    params1     = (MonoExpMinusLogProbW0W1BinnedData_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    width       = params1->width;
    delay       = params1->delay;
    alpha       = params1->hyperparam;
    halfbinsize = interval/(double)(nbins<<1);

    value       = alpha*w1;
    oneminusw0T = (1.0-w0)*interval;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            value       -= (double)(nphotonsbin) * (log(oneminusw0T*bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay)+w0));

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    bayes_dm_CorrectInvalidDoubleValue(&value, type);
			    return (value);
		    }
        }
    }
       
    return (value);
}
#endif



#if 0
void bayes_MonoExpDirectMostProbW0W1BinnedData( int   *data,
                                                int    nbins,
                                                int    nphotons,
                                                float *w0,
                                                float *w1,
                                                float *val,
                                                int    p,
                                                float  delay,
                                                float  width,
                                                float  interval,
                                                float  alpha,
                                                float  precision)
{
    double value, *x;
    int    id=0;

    MonoExpMinusLogProbW0W1BinnedData_t container;

    container.data       = data;
    container.nbins      = nbins;
    container.nphotons   = nphotons;
    container.interval   = interval;
    container.width      = width;
    container.delay      = delay;
    container.hyperparam = alpha;

    x    = Bayes_dvector(1,2);
    x[1] = (double)(*w0);
    x[2] = (double)(*w1);

    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW0W1BinnedData, 2, id, (void *)(&container), &value, x, (double)precision);
    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW0W1BinnedData, 2, id, (void *)(&container), &value, x, (double)precision);

    *w0  = (float)x[1]; 
    *w1  = (float)x[2]; 
    *val = (float)value;

    free_Bayes_dvector(x,1,2);
}
#endif




#if 0
double bayes_MonoExpMinusLogProbW1BinnedData(double *x, int id, void *container)
{
    double  value, alternative, diff, neww0, w0, w1, oldw0;
    double  dp, interval, width, delay, alpha, halfbinsize, arrivaltime;
    int    *data, nbins, nphotons, bin, nphotonsbin, type;

    MonoExpMinusLogProbW0W1BinnedData_t *params1;

    w1 = x[1];

    if (w1 <= 0.0)
        return (BIG);

    params1     = (MonoExpMinusLogProbW0W1BinnedData_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    width       = params1->width;
    delay       = params1->delay;
    alpha       = params1->hyperparam;
    halfbinsize = interval/(double)(nbins<<1);

    dp       = (double)(nphotons);
    diff     = 1.0;
    w0       = 0.001;
    oldw0    = 0.001;
    diff     = 0.1;

    while (diff > 0.000000001)
    {
        w0 = bayes_MapW0BinnedData(oldw0, w1, halfbinsize, interval, width, delay, data, nbins, nphotons);

		/* Check that w0 is a valid number... */
		bayes_dm_CheckDoubleValueValid(w0, &type);

	    if (type == BAYES_DM_DTYPE_NOT_A_NUMBER)
            w0 = 1.0;

        diff  = sqrt((w0-oldw0)*(w0-oldw0));
        oldw0 = w0;
    }

    /* Computation of [1], Eqn. (88), assuming noisy environment (i.e. 0.0<w0<1.0)... */
    value = alpha * w1/dp;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            value       -= (double)(nphotonsbin) * (log((1.0-w0) * interval * bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay) + w0) / dp);

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    bayes_dm_CorrectInvalidDoubleValue(&value, type);
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    /* Computation of [1], Eqn. (88), assuming perfect noiseless environment (i.e. w0 = 0.0)... */
    neww0       = 0.0;
    alternative = alpha * w1/dp;
    
    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            alternative -= (double)(nphotonsbin) * (log(interval * bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay)) / dp);

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    alternative = value + 100.0;
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    /* Computation of [1], Eqn. (88), assuming overwhelming noisy environment (i.e. w0 = 1.0)... */
    neww0       = 1.0;
    alternative = alpha*w1/dp;

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    params1->w0 = w0;
    
    return (value);
}
#endif


#if 0
double bayes_MapW0BinnedDataOverflowCorrected(double  w0,
                                              double  w1,
                                              double  halfbinsize,
                                              double  interval,
                                              double  width,
                                              double  delay, 
                                              int    *data,
                                              int     nbins,
                                              int     p)
{
    double dum, value, factor, arrivaltime;
    int    bin, nphotonsbin;
   
    if (w1 <= 0.0)
        return (1.0);
    
    if ((w0 <= 0.0) || (w0 >= 1.0))
        return (0.0);
   
    factor = interval*(1.0/w0 - 1.0);
    value  = 0.0;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            dum          = 1.0 + factor*bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay);
            value       += nphotonsbin/dum;
        }
    }

    return (value);
}
#endif

#if 0
double bayes_MonoExpMinusLogProbW1BinnedDataOverflowCorrected(double *x, int id, void *container)
{
    double  value, alternative, diff, neww0, w0, w1, oldw0, alphatimesw1;
    double  dp, interval, width, delay, alpha, halfbinsize, arrivaltime;
    int    *data, nbins, nphotons, bin, nphotonsbin, type;

    MonoExpMinusLogProbW0W1BinnedData_t *params1;

    w1 = x[1];

    if (w1 <= 0.0)
        return (BIG);

    params1     = (MonoExpMinusLogProbW0W1BinnedData_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    width       = params1->width;
    delay       = params1->delay;
    alpha       = params1->hyperparam;
    halfbinsize = interval/(double)(nbins<<1);

    dp       = (double)(nphotons);
    diff     = 1.0;
    w0       = 0.001;
    oldw0    = 0.001;
    diff     = 0.1;

    while (diff > 0.000000001)
    {
        w0 = bayes_MapW0BinnedDataOverflowCorrected(oldw0, w1, halfbinsize, interval, width, delay, data, nbins, nphotons);

		/* Check that w0 is a valid number... */
		bayes_dm_CheckDoubleValueValid(w0, &type);

	    if (type == BAYES_DM_DTYPE_NOT_A_NUMBER)
            w0 = 1.0;

        diff  = sqrt((w0-oldw0)*(w0-oldw0));
        oldw0 = w0;
    }

    /* Computation of [1], Eqn. (88), assuming noisy environment (i.e. 0.0<w0<1.0)... */
    alphatimesw1 = alpha*w1;
    value = alphatimesw1;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            value       -= (double)(nphotonsbin) * (log((1.0-w0) * interval * bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay) + w0));

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    bayes_dm_CorrectInvalidDoubleValue(&value, type);
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    /* Computation of [1], Eqn. (88), assuming perfect noiseless environment (i.e. w0 = 0.0)... */
    neww0       = 0.0;
    alternative = alphatimesw1;
    
    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            alternative -= (double)(nphotonsbin) * (log(interval * bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay)));

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    alternative = value + 100.0;
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    /* Computation of [1], Eqn. (88), assuming overwhelming noisy environment (i.e. w0 = 1.0)... */
    neww0       = 1.0;
    alternative = alphatimesw1;

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    params1->w0 = w0;
    
    return (value);
}
#endif



#if 0
double bayes_MonoExpMinusLogProbW1BinnedDataFast(double *x, int id, void *container)
{
    double  value, alternative, diff, neww0, w0, w1, oldw0;
    double  dp, interval, width, delay, alpha, halfbinsize, arrivaltime, *likelihood;
    int    *data, nbins, nphotons, bin, nphotonsbin, type, update;

    MonoExpMinusLogProbW0W1BinnedData_t *params1;

    w1 = x[1];

    if (w1 <= 0.0)
        return (BIG);

    params1     = (MonoExpMinusLogProbW0W1BinnedData_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    width       = params1->width;
    delay       = params1->delay;
    alpha       = params1->hyperparam;
    halfbinsize = interval/(double)(nbins<<1);

    likelihood = dvector(0,nbins-1);

    bayes_ArrTimeLikelihoodGivenTauFixedInstrumentFixedLifetimeUpdate(w1,interval,width,delay);

    for (bin=0, update=1; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime     = halfbinsize * (double)(1+(bin<<1));
            //likelihood[bin] = bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay);
            likelihood[bin] = bayes_ArrTimeLikelihoodGivenTauFixedInstrumentFixedLifetime(arrivaltime,w1,interval,width,delay);
            //update          = 0;
        }
        else
        {
            likelihood[bin] = 0.0;        
        }
    }

    dp       = (double)(nphotons);
    diff     = 1.0;
    w0       = 0.001;
    oldw0    = 0.001;
    diff     = 0.1;

    while (diff > 0.000000001)
    {
        w0 = bayes_MapW0BinnedDataFast(oldw0, w1, interval, data, likelihood, nbins, nphotons);

		/* Check that w0 is a valid number... */
		bayes_dm_CheckDoubleValueValid(w0, &type);

	    if (type == BAYES_DM_DTYPE_NOT_A_NUMBER)
            w0 = 1.0;

        diff  = sqrt((w0-oldw0)*(w0-oldw0));
        oldw0 = w0;
    }

    /* Computation of [1], Eqn. (88), assuming noisy environment (i.e. 0.0<w0<1.0)... */
    value = alpha * w1/dp;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            value       -= (double)(nphotonsbin) * (log((1.0-w0) * interval * likelihood[bin] + w0) / dp);

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    bayes_dm_CorrectInvalidDoubleValue(&value, type);
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    /* Computation of [1], Eqn. (88), assuming perfect noiseless environment (i.e. w0 = 0.0)... */
    neww0       = 0.0;
    alternative = alpha * w1/dp;
    
    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            alternative -= (double)(nphotonsbin) * (log(interval * likelihood[bin]) / dp);

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    alternative = value + 100.0;
			    bin = nbins; /* i.e. exit the outer loop... */
		    }
        }
    }

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    /* Computation of [1], Eqn. (88), assuming overwhelming noisy environment (i.e. w0 = 1.0)... */
    neww0       = 1.0;
    alternative = alpha*w1/dp;

    if (alternative < value)
    {
        w0    = neww0;
        value = alternative;
    }

    free_dvector(likelihood,0,nbins-1);

    params1->w0 = w0;
    
    return (value);
}
#endif








#if 0
double bayes_MapW0BinnedData(double w0, double w1, double halfbinsize, double interval, double width, double delay, int *data, int nbins, int p)
{
    double dum, value, factor, dp, arrivaltime;
    int    bin, nphotonsbin;
   
    if (w1 <= 0.0)
        return (1.0);
    
    if ((w0 <= 0.0) || (w0 >= 1.0))
        return (0.0);
   
    factor = interval*(1.0/w0 - 1.0);
    dp     = (double)p;
    value  = 0.0;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            dum          = 1.0 + factor*bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay);
            value       += nphotonsbin/(dum*dp);
        }
    }

    return (value);
}
#endif







#if 0
void bayes_IndirectMostProbableBinnedData(  int   *data,
                                            int    nbins,
                                            int    nphotons,
                                            float *w0,
                                            float *w1,
                                            float *val,
                                            int    p,
                                            float  delay,
                                            float  width,
                                            float  interval,
                                            float  alpha,
                                            float  precision)
{
    double value, *x;
    int    id=0;

    MonoExpMinusLogProbW0W1BinnedData_t container;

    container.data       = data;
    container.nbins      = nbins;
    container.nphotons   = nphotons;
    container.interval   = interval;
    container.width      = width;
    container.delay      = delay;
    container.hyperparam = alpha;

    x    = Bayes_dvector(1,1);
    x[1] = (double)(*w1);

    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW1BinnedData, 1, id, (void *)(&container), &value, x, (double)precision);
    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW1BinnedData, 1, id, (void *)(&container), &value, x, (double)precision);

    *w1  = (float)x[1];
    *w0  = (float)container.w0;
    *val = (float)value;
   
    free_Bayes_dvector(x,1,1);
}
#endif

#if 0
void bayes_IndirectMostProbableBinnedDataFast(  int   *data,
                                                int    nbins,
                                                int    nphotons,
                                                float *w0,
                                                float *w1,
                                                float *val,
                                                int    p,
                                                float  delay,
                                                float  width,
                                                float  interval,
                                                float  alpha,
                                                float  precision)
{
    double value, *x;
    int    id=0;

    MonoExpMinusLogProbW0W1BinnedData_t container;

    container.data       = data;
    container.nbins      = nbins;
    container.nphotons   = nphotons;
    container.interval   = interval;
    container.width      = width;
    container.delay      = delay;
    container.hyperparam = alpha;

    x    = Bayes_dvector(1,1);
    x[1] = (double)(*w1);
#ifndef BAYES_BINNED_LIKELIHOOD_IMPLEMENTATION
    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW1BinnedDataFast, 1, id, (void *)(&container), &value, x, (double)precision);
    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW1BinnedDataFast, 1, id, (void *)(&container), &value, x, (double)precision);
#else
    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW1BinnedDataFastBinLikelihood, 1, id, (void *)(&container), &value, x, (double)precision);
    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW1BinnedDataFastBinLikelihood, 1, id, (void *)(&container), &value, x, (double)precision);
#endif
    *w1  = (float)x[1];
    *w0  = (float)container.w0;
    *val = (float)value;
   
    free_Bayes_dvector(x,1,1);
}
#endif





#if 0
double bayes_MonoExpMinusLogProbW0W1WidthDelayBinnedData(double *x, int id, void *container)
{
    double  oneminusw0T, value, w0, w1, dp;
    double  interval, width, delay, alpha, halfbinsize, arrivaltime;
    int    *data, nbins, nphotons, bin, nphotonsbin, type;

    MonoExpMinusLogProbW0W1WidthDelayBinnedData_t *params1;
    
	w0    = x[1];
	w1    = x[2];
	delay = x[3];
	width = x[4];

	if ((w0 < 0.0) || (w0 > 1.0) || (w1 <= 0.0))
		return (BIG);

    if ((delay < 0.0) || (width < 0.0))
		return (BIG);

    params1     = (MonoExpMinusLogProbW0W1WidthDelayBinnedData_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    alpha       = params1->hyperparam;
    halfbinsize = interval/(double)(nbins<<1);
    dp          = (double)nphotons;

    oneminusw0T = (1.0-w0)*interval;
    value       = alpha*(w1+width+delay)/dp;

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            value       -= (double)(nphotonsbin) * (log(oneminusw0T*bayes_ArrTimeLikelihoodGivenTau(arrivaltime,w1,interval,width,delay)+w0)/dp); 

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    bayes_dm_CorrectInvalidDoubleValue(&value, type);
			    return value;
		    }
        }
    }

    return (value);
}
#endif


#if 0
double bayes_MonoExpMinusLogProbW0W1DelayBinnedData(double *x, int id, void *container)
{
    double  oneminusw0T, value, w0, w1, dp;
    double  interval, width, delay, alpha, halfbinsize, arrivaltime, *likelihoods;
    int    *data, nbins, nphotons, bin, nphotonsbin, type;

    MonoExpMinusLogProbW0W1DelayBinnedData_t *params1;
    
	w0    = x[1];
	w1    = x[2];
	delay = x[3];

	if ((w0 < 0.0) || (w0 > 1.0) || (w1 <= 0.0))
		return (BIG);

    if (delay < 0.0)
		return (BIG);

    params1     = (MonoExpMinusLogProbW0W1DelayBinnedData_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    nphotons    = params1->nphotons;
    interval    = params1->interval;
    width       = params1->width;
    alpha       = params1->hyperparam;
    halfbinsize = interval/(double)(nbins<<1);
    dp          = (double)nphotons;

    oneminusw0T = (1.0-w0)*interval;
    value       = (alpha*w1/*+alpha*0.5*delay*/)/dp;

    likelihoods = dvector(0,nbins-1);
    bayes_ArrTimeBinLikelihoodsGivenTauFixedInstrFixedLifetimeUpdate(interval,width,delay,w1);
    bayes_ArrTimeBinLikelihoodsGivenTauFixedInstrFixedLifetime(likelihoods,data,nbins,interval,width,delay,w1);

    for (bin=0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            arrivaltime  = halfbinsize * (double)(1+(bin<<1)); /* i.e. faster implementation of bin midpoint value determination... */
            value       -= (double)(nphotonsbin) * (log(oneminusw0T*likelihoods[bin]+w0)/dp); 

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value, &type))
		    {
			    bayes_dm_CorrectInvalidDoubleValue(&value, type);
			    return value;
		    }
        }
    }

    free_dvector(likelihoods,0,nbins-1);

    return (value);
}
#endif


#if 0
void bayes_MostProbableFreeInstrument(  int   *data,
                                        int    nbins,
                                        int    nphotons,
                                        float *w0,
                                        float *w1,
                                        float *val,
                                        int    p,
                                        float *delay,
                                        float *width,
                                        float  interval,
                                        float  alpha,
                                        float  precision)
{
    double value, *x;
    int    id=0;

    MonoExpMinusLogProbW0W1WidthDelayBinnedData_t container;

    container.data       = data;
    container.nbins      = nbins;
    container.nphotons   = nphotons;
    container.interval   = interval;
    container.hyperparam = alpha;

    x    = Bayes_dvector(1,4);
    x[1] = (double)(*w0);
    x[2] = (double)(*w1);
    x[3] = (double)(*delay);
    x[4] = (double)(*width);

    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW0W1WidthDelayBinnedData, 4, id, (void *)(&container), &value, x, (double)precision);
    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW0W1WidthDelayBinnedData, 4, id, (void *)(&container), &value, x, (double)precision);
   
    *w0    = (float)x[1];
    *w1    = (float)x[2];
    *delay = (float)x[3];
    *width = (float)x[4];
    *val   = (float)value;
   
    free_Bayes_dvector(x,1,4);
}
#endif



#if 0
void bayes_MostProbableFreeInstrDelayNoSpread(int   *data,
                                             int    nbins,
                                             int    nphotons,
                                             float *w0,
                                             float *w1,
                                             float *val,
                                             int    p,
                                             float *delay,
                                             float *width,
                                             float  interval,
                                             float  alpha,
                                             float  precision)
{
    double value, *x;
    int    id=0;

    MonoExpMinusLogProbW0W1DelayBinnedData_t container;

    container.data       = data;
    container.nbins      = nbins;
    container.nphotons   = nphotons;
    container.interval   = interval;
    container.width      = *width;
    container.hyperparam = alpha;

    x    = Bayes_dvector(1,3);
    x[1] = (double)(*w0);
    x[2] = (double)(*w1);
    x[3] = (double)(*delay);

    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW0W1DelayBinnedData, 3, id, (void *)(&container), &value, x, (double)precision);
    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW0W1DelayBinnedData, 3, id, (void *)(&container), &value, x, (double)precision);
   
    *w0    = (float)x[1];
    *w1    = (float)x[2];
    *delay = (float)x[3];
    *val   = (float)value;
   
    free_Bayes_dvector(x,1,3);
}
#endif

#if 0
int bayes_AveragesAndErrorBarsBinnedData(   int   *data,
                                            int    nbins,
                                            int    nphotons,
						                    float *w0,
							                float *w1,
							                float *dw0,
							                float *dw1,
							                float  delay,
							                float  width,
						                    float  interval,
                                            float  alpha,
							                float  precision,
							                int    p,/* hmm, this is actually nphotons */
							                int    quick,
							                int    distr)
{
	double average, error, value, dxdy;
	double norm, w0min, w0max, w1min, w1max, X[3], **Prob;
	double *w0values, *w1values, marginal, dx, dy, range, maxprob, newprob, neww0, neww1;
	int    i, j, points, maxCount=1000;
    int    id=0;
    int    result;
	char   filename[100];
	FILE   *fp;

    MonoExpMinusLogProbW0W1BinnedData_t container1;

    /* Locate the most probable parameter values initially... */
    container1.data       = data;
    container1.nbins      = nbins;
    container1.nphotons   = nphotons;
    container1.interval   = interval;
    container1.width      = width;
    container1.delay      = delay;
    container1.hyperparam = alpha;

	X[1] = 0.5;
	X[2] = 1.0;

    math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW0W1BinnedData, 2, id, (void *)(&container1), &value, X, (double)precision);

	maxprob  = bayes_MonoExpProbW0W1BinnedData(X[1],X[2],data,nbins,interval,width,delay,alpha,&result);
	
    if (RESULT_INVALID_ARITHMETIC_ERROR == result)
    {
        /* Get the most accurate possible (w0,w1) predictions using the indirect routine... */

        X[1] = X[2];
        math_MinimiseFctDoubleWithGenericContainer(0, bayes_MonoExpMinusLogProbW1BinnedData, 1, id, (void *)(&container1), &value, X, (double)precision);

        *w1  = (float)X[1];
        *w0  = (float)container1.w0;

        return (-1);
    }

    points = MAX((int)(0.5+0.001/precision),1);

    Prob     = Bayes_dmatrix(0,points,0,points);
    w0values = Bayes_dvector(0,points);
	w1values = Bayes_dvector(0,points);

    neww0    = X[1];
	neww1    = X[2];
    *w0      = (float)X[1];
	*w1      = (float)X[2];
    
	if (maxprob == 0.0)
	{
        puts("too many points, cannot calculate error bars");
		return -1;
    }

    /* Estimate the w0-width around the most likely point (w0,w1) */
    range = 0.1;

	while (range > 0.00005)
	{
        neww0 += range;
		
		if (neww0 > 1.0)
			neww0 = 1.0;
		
        newprob  = bayes_MonoExpProbW0W1BinnedData(neww0,neww1,data,nbins,interval,width,delay,alpha,&result);
      
		if (newprob < (0.4*maxprob))
		{
			neww0 -= range;
			range *= 0.1;
		}
    }
   
	neww0 += (5.0*range);
	range  = neww0-X[1];
    w0min  = X[1]-6.0*range;
	w0max  = X[1]+6.0*range;
   
	if (w0min < 0.0)
		w0min = 0.0;
	
	if (w0max > 1.0)
		w0max = 1.0;

    *dw0 = (float)range;

    /* Estimate the w1-width around the most likely point (w0,w1) */
    neww0 = X[1];
    range = 0.1;
	
	while (range > 0.00005)
	{
        neww1   += range;
		newprob  = bayes_MonoExpProbW0W1BinnedData(neww0,neww1,data,nbins,interval,width,delay,alpha,&result);
      
		if (newprob < (0.4*maxprob))
		{
			neww1 -= range;
			range *= 0.1;
		}
    }
   
	neww1 += (5.0*range);
	range  = neww1-X[2];
    w1min  = X[2]-6.0*range;
	w1max  = X[2]+6.0*range;
   
	if (w1min < 0.0)
		w1min = SMALLTIME;
   
	*dw1 = (float)range;

    if(quick)
		return 1;
   
   /*
   		Fill 2D array with p(w|D) values
   
   */
   
    dx   = (w0max-w0min)/(double)points;
	dy   = (w1max-w1min)/(double)points;
	dxdy = dx*dy;
	norm = 0.0;
   
	for (i=0; i<=points; i++)
		w0values[i] = w0min+(double)i*dx;
   
	for (j=0; j<=points; j++)
		w1values[j] = w1min+(double)j*dy;
   
	for (i=0; i<=points; i++)
		for (j=0; j<=points; j++)
		{
            value      = bayes_MonoExpProbW0W1BinnedData(w0values[i],w1values[j],data,nbins,interval,width,delay,alpha,&result);
			Prob[i][j] = value;
			norm      += value;
        }

    norm *= dxdy;

	for (i=0; i<=points; i++)
		for (j=0; j<=points; j++)
			Prob[i][j] = Prob[i][j]/norm;

   // Save data to a file if wanted   
    if (distr)
	{
        sprintf(filename,"prob_density_p%d.txt", nphotons);
		fp = fopen(filename,"w");
      
		//for (i=0; i<=points; i++)
		for (i=0; i<points; i++)
		{
            //for (j=0; j<=points; j++)
            for (j=0; j<points; j++)
			{
                fprintf(fp,"%lf\t%lf\t%g\n",w0values[i],w1values[j],Prob[i][j]);
				//fprintf(fp,"%lf %lf %lf\n",w0values[i],w1values[j],Prob[i][j]);
            }
         
			fprintf(fp,"\n");
        }
      
		fclose(fp);
    }

   // display graphs
#if 0
   if(distr) {
	   for(i=0;i<=points;i++) {
		   sprintf(filename, "w1 at w0=%f, peak w0 is %f", w0values[i], X[1]);
		   XYGraphPopup (filename, w1values, Prob[i], points, VAL_DOUBLE, VAL_DOUBLE);
	   }
   }
#endif
   // Get the averages and errors
   
    average = 0.0;
   
	for (i=0; i<=points; i++)
	{
        marginal = 0.0;
		
		for (j=0; j<=points; j++)
			marginal += Prob[i][j];
      
		marginal *= dy;
		average  +=(w0values[i]*marginal);
    }
    
	average *= dx;
	error    = 0.0;

    for (i=0; i<=points; i++)
	{
        marginal = 0.0;
		
		for (j=0; j<=points; j++)
			marginal += Prob[i][j];
      
		marginal *= dy;
        error    += ((w0values[i]-average)*(w0values[i]-average)*marginal);
    }
    
	error  *= dx;
	error   = sqrt(error);
	*w0     = (float)average;
	*dw0    = (float)error;
    average = 0.0;

    for (j=0; j<=points; j++)
	{
        marginal = 0.0;
		
		for (i=0; i<=points; i++)
			marginal += Prob[i][j];
      
		marginal *= dx;
		average  += (w1values[j]*marginal);
    }
    
	average *= dy;
	error    = 0.0;

    for (j=0; j<=points; j++)
	{
        marginal = 0.0;
		
		for (i=0; i<=points; i++)
			marginal += Prob[i][j];
		
		marginal *= dx;
        error    += ((w1values[j]-average)*(w1values[j]-average)*marginal);
    }
   
	error *= dy;
	error  = sqrt(error);
	*w1    = (float)average;
	*dw1   = (float)error;
   
	free_Bayes_dmatrix(Prob,0,points,0,points);
    free_Bayes_dvector(w0values,0,points);
	free_Bayes_dvector(w1values,0,points);

	return 0;
}
#endif





#if 0
int bayes_AveragesAndErrorBarsBinnedDataOverflowCorrected(  int   *data,
                                                            int    nbins,
                                                            int    nphotons,
						                                    float *w0_mp,
							                                float *w1_mp,
						                                    float *w0_ave,
							                                float *w1_ave,
                                                            float *dw0,
							                                float *dw1,
							                                float  delay,
							                                float  width,
						                                    float  interval,
                                                            float  alpha,
							                                float  precision,
							                                int    p,/* hmm, this is actually nphotons */
							                                int    quick,
							                                int    distr)
{
	double average, error, value, temp, min, dxdy;
	double norm, w0min, w0max, w1min, w1max, X[3], **Prob;
	double *w0values, *w1values, marginal, *marginals, dx, dy, range, maxprob, newprob, neww0, neww1;
	int    i, j, points, maxCount=1000;
    int    id=0, update;
    int    result;
	char   filename[100];
	FILE   *fp;

    float val;

    MonoExpMinusLogProbW0W1BinnedData_t container1;

    /* Determine most accurate most-probable values as the starting point for computation of all distribution parameter statistics... */
    bayes_IndirectMostProbableBinnedDataFast(data,nbins,nphotons,w0_mp,w1_mp,&val,nphotons,delay,width,interval,alpha,precision);

	X[1] = *w0_mp;
	X[2] = *w1_mp;

    /* ...now determine the 'normalization' constant for subsequent probability computations... */
    container1.data       = data;
    container1.nbins      = nbins;
    container1.nphotons   = nphotons;
    container1.interval   = interval;
    container1.width      = width;
    container1.delay      = delay;
    container1.hyperparam = alpha;

    min      = bayes_MonoExpMinusLogProbW0W1BinnedDataOverflowCorrected(X,id,(void *)(&container1));
	maxprob  = bayes_MonoExpProbW0W1BinnedDataOverflowCorrected(X[1],X[2],data,nbins,interval,width,delay,alpha,min,&result);

    /* Now start analysis of the distribution and generate its statistics... */
    points = MAX((int)(0.5+0.001/precision),1);

    Prob     = Bayes_dmatrix(0,points,0,points);
    w0values = Bayes_dvector(0,points);
	w1values = Bayes_dvector(0,points);

    neww0    = X[1];
	neww1    = X[2];
    
    /* Estimate the w0-width around the most likely point (w0,w1) */
    range = 0.1;

    DBG_TIMING_START_TIMING()


    update = 1;

	while (range > 0.00005) /* Determine the width of the background parameter marginal distribution... */
	{
        neww0 += range;
		
		if (neww0 > 1.0)
			neww0 = 1.0;

        /* Only ever one fixed value of w1 used in this loop, use optimised probability function... */
		//newprob  = bayes_MonoExpProbW0W1BinnedDataOverflowCorrected(neww0,neww1,data,nbins,interval,width,delay,alpha,min,&result);

        newprob  = bayes_MonoExpProbW0W1BinnedDataOverflowCorrectedFixedInstrumentFixedLifetime(neww0,neww1,data,nbins,interval,width,delay,update,alpha,min,&result);
        update   = 0;

		if (newprob < (0.4*maxprob))
		{
			neww0 -= range;
			range *= 0.1;
		}
    }

    DBG_TIMING_STOP_TIMING("Check 1")

	neww0 += (5.0*range);
	range  = neww0-X[1];
    w0min  = X[1]-6.0*range;
	w0max  = X[1]+6.0*range;
   
	if (w0min < 0.0)
		w0min = 0.0;
	
	if (w0max > 1.0)
		w0max = 1.0;

    *dw0 = (float)range;

    /* Estimate the w1-width around the most likely point (w0,w1) */
    neww0 = X[1];
    range = 0.1;

    DBG_TIMING_START_TIMING()

	while (range > 0.00005)
	{
        neww1   += range;

        //newprob  = bayes_MonoExpProbW0W1BinnedDataOverflowCorrected(neww0,neww1,data,nbins,interval,width,delay,alpha,min,&result);
		newprob  = bayes_MonoExpProbW0W1BinnedDataOverflowCorrectedFixedInstrument(neww0,neww1,data,nbins,interval,width,delay,1,alpha,min,&result);
      
		if (newprob < (0.4*maxprob))
		{
			neww1 -= range;
			range *= 0.1;
		}
    }

    DBG_TIMING_STOP_TIMING("Check 2")

	neww1 += (5.0*range);
	range  = neww1-X[2];
    w1min  = X[2]-6.0*range;
	w1max  = X[2]+6.0*range;
   
	if (w1min < 0.0)
		w1min = SMALLTIME;
   
	*dw1 = (float)range;

    if(quick)
		return 1;
   
   /*
   		Fill 2D array with p(w|D) values
   
   */
   
    dx   = (w0max-w0min)/(double)points;
	dy   = (w1max-w1min)/(double)points;
	dxdy = dx*dy;
	norm = 0.0;

	for (i=0; i<=points; i++)
    {
		w0values[i] = w0min+(double)i*dx;
		w1values[i] = w1min+(double)i*dy;
    }

    DBG_TIMING_START_TIMING()

	for (j=0; j<=points; j++) /* Run through the lifetime values in outer loop...   */
    {                         /* ...so that inner loop is over a constant lifetime. */
        for (i=0, update=1; i<=points; i++)
		{
            //value = bayes_MonoExpProbW0W1BinnedDataOverflowCorrected(w0values[i],w1values[j],data,nbins,interval,width,delay,alpha,min,&result);
            value      = bayes_MonoExpProbW0W1BinnedDataOverflowCorrectedFixedInstrumentFixedLifetime(w0values[i],w1values[j],data,nbins,interval,width,delay,update,alpha,min,&result);
			update     = 0;
            Prob[i][j] = value;
			norm      += value;
        }
    }

    DBG_TIMING_STOP_TIMING("Check 3")

    norm *= dxdy;

	for (i=0; i<=points; i++)
		for (j=0; j<=points; j++)
			Prob[i][j] = Prob[i][j]/norm;

   // Save data to a file if wanted   
    if (distr)
	{
        sprintf(filename,"prob_density_p%d.txt", nphotons);
		fp = fopen(filename,"w");
      
		//for (i=0; i<=points; i++)
		for (i=0; i<points; i++)
		{
            //for (j=0; j<=points; j++)
            for (j=0; j<points; j++)
			{
                fprintf(fp,"%lf\t%lf\t%g\n",w0values[i],w1values[j],Prob[i][j]);
				//fprintf(fp,"%lf %lf %lf\n",w0values[i],w1values[j],Prob[i][j]);
            }
         
			fprintf(fp,"\n");
        }
      
		fclose(fp);
    }

   // display graphs
#if 0
   if(distr) {
	   for(i=0;i<=points;i++) {
		   sprintf(filename, "w1 at w0=%f, peak w0 is %f", w0values[i], X[1]);
		   XYGraphPopup (filename, w1values, Prob[i], points, VAL_DOUBLE, VAL_DOUBLE);
	   }
   }
#endif
   // Get the averages and errors

    marginals = dvector(0,points);

    average = 0.0;

    DBG_TIMING_START_TIMING()

	for (i=0; i<=points; i++)
	{
        marginal = 0.0;
		
		for (j=0; j<=points; j++)
			marginal += Prob[i][j];
      
        marginal     *= dy;
        average      += w0values[i]*marginal;
        marginals[i]  = marginal;
    }
    
	average *= dx;
	error    = 0.0;

    for (i=0; i<=points; i++)
	{
        temp   = w0values[i]-average;
        error += temp*temp*marginals[i];
    }
    
	error  *= dx;
	error   = sqrt(error);
	*w0_ave = (float)average;
	*dw0    = (float)error;

    DBG_TIMING_STOP_TIMING("Check 4")

    average = 0.0;

    DBG_TIMING_START_TIMING()

    for (j=0; j<=points; j++)
	{
        marginal = 0.0;
		
		for (i=0; i<=points; i++)
			marginal += Prob[i][j];
      
		marginal     *= dx;
		average      += w1values[j]*marginal;
        marginals[j]  = marginal;
    }
    
	average *= dy;
	error    = 0.0;

    for (j=0; j<=points; j++)
    {
        temp   = w1values[j]-average;
        error += temp*temp*marginals[j];
    }

	error  *= dy;
	error   = sqrt(error);
	*w1_ave = (float)average;
	*dw1    = (float)error;

    DBG_TIMING_STOP_TIMING("Check 5")

    free_dvector(marginals,0,points);

	free_Bayes_dmatrix(Prob,0,points,0,points);
    free_Bayes_dvector(w0values,0,points);
	free_Bayes_dvector(w1values,0,points);

	return 0;
}
#endif



















#if 0
int direct_most_probable_wrapped_binned (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, int distr, float *val, float *params, float *fitted, float *residuals)
{
    float timeunit;
    int   i;

    if (w0_ave) *w0_ave = 0.0;
    if (w1_ave) *w1_ave = 0.0;
	if (dw0)    *dw0    = 0.0;
	if (dw1)    *dw1    = 0.0;
	
    DBG_TIMING_START_TIMING()
    bayes_MonoExpDirectMostProbW0W1BinnedData(data, nbins, nphotons, w0_mp, w1_mp, val, p, *delay, *width, interval, alpha, precision);
    DBG_TIMING_STOP_TIMING("bayes_MonoExpDirectMostProbW0W1BinnedData (running time)")

    if (params)
    {
		params[0] = bayes_ComputeRawBackgroundZ(*w0_mp, interval/(float)(nbins), interval, nphotons);                 /* Z   */
		params[1] = bayes_ComputeRawAmplitudeA(*w0_mp, *w1_mp, interval/(float)nbins, interval, *delay, nphotons);    /* A   */
        params[2] = *w1_mp;                                                                                           /* tau */
    }

    if (fitted)
    {
        timeunit = interval/(float)nbins;

        for (i=0; i<nbins; i++)
        {
            fitted[i] = bayes_ModelTransformComputePhotonCount
                            (*delay, *width, interval, timeunit*(float)i, *w0_mp, *w1_mp, timeunit, nphotons);
        }
    }

    if (residuals)
    {
        for (i=0; i<nbins; i++)
            residuals[i] = (float)(data[i])-fitted[i];

        scaleDataAccordingToSignalNoise(residuals, nbins, fitted);    
    }

    return 0;
}
#endif



#if 0
int indirect_most_probable_wrapped_binned (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, int distr, float *val, float *params, float *fitted, float *residuals)
{
    float timeunit;
    int   i;

    if (w0_ave) *w0_ave = 0.0;
    if (w1_ave) *w1_ave = 0.0;
	if (dw0)    *dw0    = 0.0;
	if (dw1)    *dw1    = 0.0;

    DBG_TIMING_START_TIMING()
    bayes_IndirectMostProbableBinnedDataFast(data, nbins, nphotons, w0_mp, w1_mp, val, p, *delay, *width, interval, alpha, precision);
    DBG_TIMING_STOP_TIMING("bayes_IndirectMostProbableBinnedDataFast (running time)")

    if (params)
    {
		params[0] = bayes_ComputeRawBackgroundZ(*w0_mp, interval/(float)(nbins), interval, nphotons);                 /* Z   */
		params[1] = bayes_ComputeRawAmplitudeA(*w0_mp, *w1_mp, interval/(float)nbins, interval, *delay, nphotons);    /* A   */
        params[2] = *w1_mp;                                                                                           /* tau */
    }

    if (fitted)
    {
        timeunit = interval/(float)nbins;

        for (i=0; i<nbins; i++)
            fitted[i] = bayes_ModelTransformComputePhotonCount(*delay, *width, interval, timeunit*(float)i, *w0_mp, *w1_mp, timeunit, nphotons);    
    }

    if (residuals)
    {
        for (i=0; i<nbins; i++)
            residuals[i] = (float)(data[i])-fitted[i];

        scaleDataAccordingToSignalNoise(residuals, nbins, fitted);    
    }
	
    return 0;
}
#endif

#if 0
int averages_and_errorbars_wrapped_binned (int *data, int nbins, int nphotons, float *w0_mp, float *w1_mp, float *w0_ave, float *w1_ave, float *dw0, float *dw1, float *delay, float *width,
						 		  int p, float interval, float alpha, float precision, int quick, int distr, float *val, float *params, float *fitted, float *residuals)
{
    float timeunit;
    int   i, ret;

	if (val) *val=0;
/*
    ret = bayes_AveragesAndErrorBarsBinnedData(data, nbins, nphotons,
                                               w0, w1, dw0, dw1,
                                               *delay, *width, interval, alpha, precision,
                                               p, quick, distr);
*/
    ret = bayes_AveragesAndErrorBarsBinnedDataOverflowCorrected(   data, nbins, nphotons,
                                                                   w0_mp, w1_mp, w0_ave, w1_ave, dw0, dw1,
                                                                   *delay, *width, interval, alpha, precision,
                                                                   p, quick, distr);

    if (-1 == ret) /* Too many points to determine error bars... */
    {
        *w0_ave = *w0_mp;
        *w1_ave = *w1_mp;
        *dw0    = 0.0;
        *dw1    = 0.0;
    }

    if (params)
    {
		params[0] = bayes_ComputeRawBackgroundZ(*w0_ave, interval/(float)(nbins), interval, nphotons);                  /* Z   */
		params[1] = bayes_ComputeRawAmplitudeA(*w0_ave, *w1_ave, interval/(float)nbins, interval, *delay, nphotons);    /* A   */
        params[2] = *w1_ave;                                                                                            /* tau */

		params[3] = bayes_ComputeRawBackgroundZ(*w0_mp, interval/(float)(nbins), interval, nphotons);                   /* Z (MP)   */
		params[4] = bayes_ComputeRawAmplitudeA(*w0_mp, *w1_mp, interval/(float)nbins, interval, *delay, nphotons);      /* A (MP)   */
        params[5] = *w1_mp;                                                                                             /* tau (MP) */
    }

    if (fitted)
    {
        timeunit = interval/(float)nbins;

        for (i=0; i<nbins; i++)
            fitted[i] = bayes_ModelTransformComputePhotonCount(*delay, *width, interval, timeunit*(float)i, *w0_ave, *w1_ave, timeunit, nphotons);    
    }

    if (residuals)
    {
        for (i=0; i<nbins; i++)
            residuals[i] = (float)(data[i])-fitted[i];

        scaleDataAccordingToSignalNoise(residuals, nbins, fitted);    
    }

    return (ret);
}
#endif




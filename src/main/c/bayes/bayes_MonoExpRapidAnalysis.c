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

#include "stdlib.h"
#include "math.h"

#include "extmath.h"
#include "matrices.h"
#include "DTYPE.h"
#include "bayes_Sizes.h"
#include "bayes_Types.h"
#include "bayes_DataManagement.h"
#include "bayes_DistributionFctsBinLikelihoods.h"
#include "bayes_MonoExpAnalysisBinLikelihoods.h"
#include "bayes_ModelTransformTools.h"
#include "bayes_ModelSelection.h"
#include "bayes_Interface.h"


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
       

double bayes_RapidMonoExpMinusLogProbW0W1(double *x, int id, void *container)
{
    int     *data, bin, nbins, fitstart, nphotonsbin;
    double  alpha, value, w0, w1, *likelihoods;

    MonoExpMinusLogProbW0W1_t *params1;

    w0 = x[0];
    w1 = x[1];

    if ((w0 < 0.0) || (w0 > 1.0))
        return (BIG);
    
    if (w1 <= 0.0)
        return (BIG);

    params1     = (MonoExpMinusLogProbW0W1_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    fitstart    = params1->fitstart;
    alpha       = params1->hyperparam;
    likelihoods = params1->likelihoods->logphotonlikelihoodgiventauandw0;

    value       = alpha*w1;

    for (bin=fitstart; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
        {
            value -= (double)(nphotonsbin) * likelihoods[bin];
        }
    }

    return (value);
}


int bayes_RapidMonoExpDirectMostProbW0W1( int                               *data,
                                          int                                nbins,
                                          int                                fitstart,
                                          int                                nphotons,
                                          double                            *w0,
                                          double                            *w1,
                                          float                             *val,
                                          BayesInstrRsp_t                   *instr,
                                          float                              interval,
                                          float                              modulationperiod,
                                          float                              alpha,
                                          BayesRapidMonoExpValueStore_t     *grid,
                                          BayesAveErrDistn_t                *distr)
{
	double value, x[2], valuebest, xbest[2]= { 0, 0 };
    int    i, j, id=0;
    int    type;

    MonoExpMinusLogProbW0W1_t container;

    container.data       = data;
    container.nbins      = nbins;
    container.fitstart   = fitstart;
    container.hyperparam = alpha;

    if (distr)
	{
        distr->x = Bayes_vector(0,grid->settings->nweights-1);
        distr->y = Bayes_vector(0,grid->settings->ntaus-1);
        distr->z = Bayes_matrix(0,grid->settings->nweights-1,0,grid->settings->ntaus-1);
        
        distr->m = grid->settings->nweights;
        distr->n = grid->settings->ntaus;

        for (i=0; i<grid->settings->nweights; i++)
            distr->x[i] = (float)grid->likelihoods[i][1].w0;

        for (j=0; j<grid->settings->ntaus; j++)
            distr->y[j] = (float)grid->likelihoods[1][j].tau;
    }

    for (i=0, valuebest=BIG; i<grid->settings->nweights; i++)
    {
        for (j=0; j<grid->settings->ntaus; j++)
        {
            x[0] = grid->likelihoods[i][j].w0;
            x[1] = grid->likelihoods[i][j].tau;
            container.likelihoods = &(grid->likelihoods[i][j]);

            value = bayes_RapidMonoExpMinusLogProbW0W1(x,id,(void*)(&container));

            if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(value,&type))
            {
                value = BIG;
            }

            if (distr)
                distr->z[i][j] = (float)value;

            if (value<valuebest)
            {
                valuebest = value;

                xbest[0] = x[0];
                xbest[1] = x[1];
            }
        }
    }

    *w0  = (float)xbest[0]; 
    *w1  = (float)xbest[1]; 
    *val = (float)valuebest;

    if (valuebest>=BIG)
        return (-1);

    return (0);
}








#if 1
/*http://www.smipple.net/snippet/Martin%20Ankerl/A%20Fast,%20Compact%20Approximation%20of%20the%20Exponential%20Function*/
double fast_exp(double y)
{
    double d;

    *((int*)(&d) + 0) = 0;
    *((int*)(&d) + 1) = (int)(1512775 * y + 1072632447);

    return (d);
}



int bayes_RapidMonoExpAvgAndErrors( int                           *data,
                                    int                            nbins,
                                    int                            fitstart,
                                    int                            nphotons,
                                    double                        *w0_mp,
	                                double                        *w1_mp,
                                    double                        *w0_ave,
	                                double                        *w1_ave,
                                    double                        *dw0,
	                                double                        *dw1,
                                    BayesInstrRsp_t               *instr,
                                    float                          interval,
                                    float                          modulationperiod,
                                    float                          alpha,
	                                float                          precision,
	                                int                            p,/* hmm, this is actually nphotons */
	                                int                            quick,
	                                BayesAveErrDistn_t            *distr,
                                    BayesRapidMonoExpValueStore_t *grid,
                                    float                         *value)
{
    int    ret = BAYES_AVE_ERRS_ROUTINE_NO_ERRORS;
	float  average, error, temp, min, norm, dx, dy, dxdy, marginal, *marginals, max;
	int    i, j;
    int    id=0;

    BayesAveErrDistn_t logdistr, *probdistr;

    /* Determine most-probable values and compute the -log[p(w0,w1|data)] grid... */
    bayes_RapidMonoExpDirectMostProbW0W1(data,nbins,fitstart,nphotons,w0_mp,w1_mp,&min,
                                         instr,interval,modulationperiod,alpha,grid,&logdistr);

    *value = min;

    if ((*w0_mp<0.0) || (*w0_mp>1.0) || (*w1_mp<0.0))
    {
        *w0_mp  = -1.0; *w1_mp  = -1.0;
        *w0_ave = -1.0; *w1_ave = -1.0;
        *dw0    = -1.0; *dw1    = -1.0;

        return (BAYES_AVE_ERRS_ROUTINE_ERROR);
    }

    if (distr)
        probdistr = distr;
    else
        probdistr = (BayesAveErrDistn_t*)malloc(sizeof(BayesAveErrDistn_t));

    probdistr->x = Bayes_vector(0,grid->settings->nweights);
    probdistr->y = Bayes_vector(0,grid->settings->ntaus);
    probdistr->z = Bayes_matrix(0,grid->settings->nweights,0,grid->settings->ntaus);

    probdistr->manualrange = 2;//for bpda
    probdistr->m = grid->settings->nweights;
    probdistr->n = grid->settings->ntaus;

    for (i=0; i<probdistr->m; i++)
        probdistr->x[i] = (float)grid->likelihoods[i][1].w0;

    for (j=0; j<probdistr->n; j++)
        probdistr->y[j] = (float)grid->likelihoods[1][j].tau;

    for (i=0; i<probdistr->m; i++)
        for (j=0; j<probdistr->n; j++)
            probdistr->z[i][j] = (float)exp((double)min - logdistr.z[i][j]);

    dy   = (float)((grid->settings->tau[grid->settings->ntaus-1]-grid->settings->tau[0])/(double)grid->settings->ntaus);
    dx   = (float)((grid->settings->weight[grid->settings->nweights-1]-grid->settings->weight[0])/(double)grid->settings->nweights);
	dxdy = dx*dy;

	for (i=0,norm=0.0; i<probdistr->m; i++)
		for (j=0; j<probdistr->n; j++)
			norm += probdistr->z[i][j];

    norm *= dxdy;

    for (i=0; i<probdistr->m; i++)
		for (j=0; j<probdistr->n; j++)
            probdistr->z[i][j] /= norm;

    for (i=0,max=0.0; i<probdistr->m; i++)
		for (j=0; j<probdistr->n; j++)
            if (probdistr->z[i][j]>max)
                max = probdistr->z[i][j];

	//for (i = 0; i < probdistr->m; i++) {
	//	printf("%d:\n", i);
	//	for (j = 0; j < probdistr->n; j++)
	//		printf("%11f ", probdistr->z[i][j]);
	//	printf("\n\n");
	//}

    /* Check for insufficient grid... */
//    for (i=0,j=0; j<probdistr->n; j++)
//        if (probdistr->z[i][j]>0.01*max)
//            ret = BAYES_AVE_ERR_RAPID_INSUFFICIENT_GRID;

    if (grid->settings->weight[0]>0.0)
        for (i=0,j=0; j<probdistr->n; j++)
            if (probdistr->z[i][j]>0.02*max)
                ret = BAYES_AVE_ERR_RAPID_INSUFFICIENT_GRID;

    if (grid->settings->weight[grid->settings->nweights-1]<1.0)
        for (i=probdistr->m-1,j=0; j<probdistr->n; j++)
            if (probdistr->z[i][j]>0.02*max)
                ret = BAYES_AVE_ERR_RAPID_INSUFFICIENT_GRID;

    if (grid->settings->tau[0]>0.0)
        for (i=0,j=0; i<probdistr->m; i++)
            if ((probdistr->z[i][j]>0.02*max) || (probdistr->z[i][j]==max))
                ret = BAYES_AVE_ERR_RAPID_INSUFFICIENT_GRID;

	for (i = 0, j = probdistr->n - 1; i < probdistr->m; i++) {
		if (probdistr->z[i][j] > 0.02 * max)
			ret = BAYES_AVE_ERR_RAPID_INSUFFICIENT_GRID;
	}

    /* Determine the average and error for 'w0'... */
    marginals = Bayes_vector(0,probdistr->m);

	for (i=0,average=0.0; i<probdistr->m; i++)
	{
		for (j=0,marginal=0.0; j<probdistr->n; j++)
			marginal += probdistr->z[i][j];
      
        marginal     *= dy;
        average      += probdistr->x[i]*marginal;
        marginals[i]  = marginal;
    }
    
    average *= dx;

    if ((average<0.0) || (average>1.0))
    {
        return (BAYES_AVE_ERRS_ROUTINE_ERROR);
    }

    if (w0_ave)
        *w0_ave  = average;
    
    *w0_mp   = /*average*/*w0_mp;

    for (i=0,error=0.0; i<probdistr->m; i++)
	{
        temp   = probdistr->x[i]-average;
        error += temp*temp*marginals[i];
    }

    if (dw0)
        *dw0 = (float)sqrt((double)error * dx);
	
	free_Bayes_vector(marginals,0,probdistr->m);

    /* Determine the average and error for 'w1'... */
	marginals = Bayes_vector(0,probdistr->n);

    for (j=0,average=0.0; j<probdistr->n; j++)
	{
		for (i=0,marginal=0.0; i<probdistr->m; i++)
			marginal += probdistr->z[i][j];
      
		marginal     *= dx;
        average      += probdistr->y[j]*marginal;
        marginals[j]  = marginal;
    }

    average *= dy;

    if (w1_ave)
    	*w1_ave  = average;
    
    *w1_mp   = /*average*/*w1_mp;

    for (j=0,error=0.0; j<probdistr->n; j++)
    {
        temp   = probdistr->y[j]-average;
        error += temp*temp*marginals[j];
    }

    if (dw1)
        *dw1 = (float)sqrt((double)error * dy);

    free_Bayes_vector(marginals,0,probdistr->n);

    if (!distr)
    {
        free_Bayes_vector(probdistr->x,0,grid->settings->nweights);
        free_Bayes_vector(probdistr->y,0,grid->settings->ntaus);
        free_Bayes_matrix(probdistr->z,0,grid->settings->nweights,0,grid->settings->ntaus);
        free(probdistr);
    }

    free_Bayes_vector(logdistr.x,0,grid->settings->nweights);
    free_Bayes_vector(logdistr.y,0,grid->settings->ntaus);
    free_Bayes_matrix(logdistr.z,0,grid->settings->nweights,0,grid->settings->ntaus);

    return (ret);
}

#endif


#if 1 //hyperparam optimization


double bayes_RapidMonoExpMinusLogProbDataLikelihood(double *x, int id, void *container)
{
    int     *data, bin, nbins, fitstart, nphotonsbin;
    double   value, w0, w1, *likelihoods;

    MonoExpMinusLogProbW0W1_t *params1;

    w0 = x[0];
    w1 = x[1];

    if ((w0 < 0.0) || (w0 > 1.0))
        return (BIG);
    
    if (w1 <= 0.0)
        return (BIG);

    params1     = (MonoExpMinusLogProbW0W1_t *)(container);
    data        = params1->data;
    nbins       = params1->nbins;
    fitstart    = params1->fitstart;
    likelihoods = params1->likelihoods->logphotonlikelihoodgiventauandw0;

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


int bayes_RapidMonoExpPopulateDataLikelihoodGrid(double                       **datalikelihoods,
                                                 int                            nw0s,
                                                 int                            ntaus,
                                                 int                           *data,
                                                 int                            nbins,
                                                 int                            fitstart,
                                                 int                            nphotons,
                                                 double                        *binwalls,
                                                 BayesInstrRsp_t               *instr,
                                                 float                          interval,
                                                 float                          modulationperiod,
                                                 BayesRapidMonoExpValueStore_t *grid)
{
    double x[2], value;
    int    i, j, id=0;
    MonoExpMinusLogProbW0W1_t container;

    container.data             = data;
    container.nbins            = nbins;
    container.fitstart         = fitstart;
    container.nphotons         = nphotons;
    container.instr            = instr;
    container.interval         = interval;
    container.modulationperiod = modulationperiod;
    container.binwalls         = binwalls;

    for (i=0; i<nw0s; i++)
    {
        for (j=0; j<ntaus; j++)
        {
            x[0]                  = grid->likelihoods[i][j].w0;
            x[1]                  = grid->likelihoods[i][j].tau;
            container.likelihoods = &(grid->likelihoods[i][j]);

            value                 = bayes_RapidMonoExpMinusLogProbDataLikelihood(x,id,(void*)(&container));
            datalikelihoods[i][j] = value;
        }
    }

    return (0);
}


int bayes_RapidMonoExpDirectMostProbW0W1PreComputedDataLikelihood(double                         *w0,
                                                                  double                         *w1,
                                                                  float                          *val,
                                                                  float                           alpha,
                                                                  BayesRapidMonoExpValueStore_t  *grid,
                                                                  double                        **datalikelihoods)
{
    double value, valuebest, w0best, w1best;
    int    i, j, id=0;

    for (i=0,valuebest=BIG,w0best=-1.0,w1best=-1.0; i<grid->settings->nweights; i++)
    {
        for (j=0; j<grid->settings->ntaus; j++)
        {
            value  = alpha*grid->likelihoods[i][j].tau;
            value += datalikelihoods[i][j];

            if (value<valuebest)
            {
                valuebest = value;

                w0best = grid->likelihoods[i][j].w0;
                w1best = grid->likelihoods[i][j].tau;
            }
        }
    }

    *w0  = (float)w0best; 
    *w1  = (float)w1best; 
    *val = (float)valuebest;

    if (valuebest>=BIG)
        return (-1);

    return (0);
}


double bayes_RapidMonoExpMinusLogProbAlphaTimesModelEvidence(double *x, int id, void *container)
{
    int     ret, *data, nbins, fitstart, nphotons;
    float   minuslogprob;
    double  alpha, alphamin, logevidence, weights[2], taus[2], interval, modperiod, *binwalls, **datalikelihoods;

    MonoExpMinusLogProbW0W1_t     *paramscontainer;
    BayesRapidMonoExpValueStore_t *grid;
    BayesInstrRsp_t               *instr;

    paramscontainer = (MonoExpMinusLogProbW0W1_t *)(container);
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
    grid            = paramscontainer->grid;
    datalikelihoods = paramscontainer->datalikelihoods;

    /* Determine most-probable parameter (w0,w1) values given the data and the hyperparameter alpha... */
    ret = bayes_RapidMonoExpDirectMostProbW0W1PreComputedDataLikelihood(&weights[0],&taus[1],&minuslogprob,
                                                                        (float)alpha,grid,datalikelihoods);
//    ret = bayes_RapidMonoExpDirectMostProbW0W1(data,nbins,fitstart,nphotons,&weights[0],&taus[1],&minuslogprob,
//                                               instr,interval,modperiod,alpha,grid,NULL);

    if (ret<0)
        return (BIG);

    weights[1] = 1.0-weights[0];

    /* Use the Gaussian approximation to determine the (log of) integral (i.e. log of model evidence)... */
    ret = bayes_DetemineDecayModelEvidence(/*ndecays*/1,weights,taus,NULL,minuslogprob,
                                           nbins,binwalls,data,
                                           interval,modperiod,instr,&logevidence);

    if (ret<0)
        return (BIG);
    else
        return (-(log(alpha)+logevidence));
}



int bayes_RapidMonoExpHyperParamOptimization(int                           *data,
                                             int                            nbins,
                                             int                            fitstart,
                                             int                            nphotons,
                                             double                        *binwalls,
                                             BayesInstrRsp_t               *instr,
                                             float                          interval,
                                             float                          modulationperiod,
                                             float                         *alphastar,
                                             float                          alphamin,
	                                         float                          precision,
                                             float                         *value,
                                             BayesRapidMonoExpValueStore_t *grid)
{
	double x[2], val, valbest, xbest[2] = { 0, 0 }, ** datalikelihoods;
    int    i, id=0, ret;
    double alpha, dalpha=0.05;

    double tempval[100], tempalpha[100];

    MonoExpMinusLogProbW0W1_t container;

    container.data             = data;
    container.nbins            = nbins;
    container.fitstart         = fitstart;
    container.nphotons         = nphotons;
    container.instr            = instr;
    container.interval         = interval;
    container.modulationperiod = modulationperiod;
    container.binwalls         = binwalls;
    container.grid             = grid;
    container.alphamin         = alphamin;
    

    datalikelihoods = Bayes_dmatrix(0,grid->settings->nweights,0,grid->settings->ntaus);

    ret = bayes_RapidMonoExpPopulateDataLikelihoodGrid(datalikelihoods,grid->settings->nweights,grid->settings->ntaus,
                                                       data,nbins,fitstart,nphotons,
                                                       binwalls,instr,interval,modulationperiod,
                                                       grid);

    container.datalikelihoods  = datalikelihoods;

    for (i=0,valbest=BIG; i<50; i++)
    {
        alpha = alphamin+dalpha*(float)i;
        x[1] = (float)alpha;

        val = bayes_RapidMonoExpMinusLogProbAlphaTimesModelEvidence(x,id,(void*)(&container));

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

    free_Bayes_dmatrix(datalikelihoods,0,grid->settings->nweights,0,grid->settings->ntaus);

    return (0);
}




#endif

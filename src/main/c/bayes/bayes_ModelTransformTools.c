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
#include "math.h"
#include "extmath.h"
#include "stdlib.h"
#include "random.h"
#include "matrices.h"
#include "bayes_ModelTransformTools.h"
#include "bayes_Interface.h"
#include "bayes_DistributionFctsBinLikelihoods.h"



/*=================================================================================*/
/*                                                                                 */
/*  ROUTINES TO TRANSFORM BETWEEN:                                                 */
/*   1. BAYES MODEL OF FORM                                                        */
/*        p(w,tau...|data) = w0 + w1.F(tau1|data) + w2.F(tau2|data) +...           */
/*   2.CONVENTIONAL MODEL OF FORM                                                  */
/*        p(Z,A,tau...|data) = Z + A1.F(tau1|data) + A2.F(tau2|data) +...          */
/*                                                                                 */
/*=================================================================================*/


float bayes_ToBayesModelTransformFromParamZ(float Z, int nbins, int nphotons)
{
    return (Z*(float)nbins/(float)nphotons);
}


float bayes_FromBayesModelTransformToParamZ(float w0, int nbins, int nphotons)
{
    return (w0*(float)nphotons/(float)nbins);
}


float bayes_FromBayesModelWeightAndTauToParamA( float w_i,
                                                float tau_i,
                                                float binsize,
                                                float interval,
                                                float delay,
                                                int   nphotons)
{
    double num, den1, R_i, upsilon3, value;

    num  = 1.0/tau_i; 
    den1 = 1.0-exp((delay-interval)/tau_i);

    if (bayes_UseRepetitionEffectsInAnalysis())
    {
        upsilon3  = bayes_ComputeArrBinLikelihoodConstantUpsilon3(bayes_GetConfigParameterValueModulationPeriod(),tau_i);
        R_i       = 1.0+1.0/(exp(bayes_GetConfigParameterValueModulationPeriod()/tau_i)-1.0);
        den1     += upsilon3*(exp(delay/tau_i)-1.0);
        den1     *= R_i;
    }

    value = (w_i)*(num/den1);

    return ((float)((float)nphotons*binsize*value));
}


float bayes_ToBayesModelWeightFromParamAAndTau( float A_i,
                                                float tau_i,
                                                float binsize,
                                                float interval,
                                                float delay,
                                                int   nphotons)
{
    double value;

    value  = 1.0-exp((delay-interval)/tau_i);
    value *= A_i*tau_i;

    return ((float)value/((float)nphotons*binsize));
}


int bayes_CheckForValidBayesModelWeightParamValues(/* Model and params... */
                                                   int     nweights,
                                                   double *weights,
                                                   int    *isweightfixed)
{
    int    i;
    double sum, sumfixed;

    for (i=0, sum=0.0, sumfixed=0.0; i<nweights; i++)
    {
        if (BAYES_PARAM_VALUE_USER_FIXED == isweightfixed[i])
        {
            if ((weights[i]<0.0) || (weights[i]>1.0))
                return (-1); /* User problem, illegal weight value has been fixed... */

            sumfixed += weights[i];
        }

        if ((sumfixed<0.0) || (sumfixed>1.0))
            return (-1); /* User problem, illegal weight value has been fixed... */

        sum += weights[i];

        if ((sum<0.0) || (sum>1.0))
            return (-2);
    }

    return (0);
}


int bayes_CheckForValidBayesModelLifetimeParamValues(/* Model and params... */
                                                     int     nlifetimes,
                                                     double *lifetimes,
                                                     int    *islifetimefixed)
{
    int i;

    for (i=1; i<=nlifetimes; i++)
    {
        if (lifetimes[i]<0.0)
        {
            if (islifetimefixed[i])
                return (-1); /* User problem, illegal lifetime value has been fixed... */
            else
                return (-2);
        }
    }

    return (0);
}


// this function merely converts the given values, doesn't check for errors or weights summing to 1...
int bayes_ToBayesModelParamValuesFromConventionalModelParamValues(/* Model and params... */
                                                                  int    nparams,
                                                                  float *params_bayes,
                                                                  float *params_conventional,
                                                                  /* Instrument response... */
                                                                  float  delay,
                                                                  /* Transient... */
                                                                  int    nbins,
                                                                  int    nphotons,
                                                                  float  interval)
{
    int     ndecays, i;
    double *weights_bayes, *taus_bayes, *weights_conventional, *taus_conventional;

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_conventional,&taus_conventional);
    bayes_PopulateWeightsAndTausVectorsFromParamVector(nparams,params_conventional,weights_conventional,taus_conventional);

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_bayes,&taus_bayes);

    weights_bayes[0] = bayes_ToBayesModelTransformFromParamZ(params_conventional[0],nbins,nphotons);

    for (i=1; i<=ndecays; i++)
    {
        weights_bayes[i] = (double)bayes_ToBayesModelWeightFromParamAAndTau((float)weights_conventional[i],(float)taus_conventional[i],
                                                                            interval/(float)nbins,interval,delay,nphotons);

        taus_bayes[i]    = taus_conventional[i];
    }

    bayes_PopulateParamVectorFromWeightsAndTausVectors(params_bayes,nparams,weights_bayes,taus_bayes);

    bayes_FreeWeightsAndTausVectors(ndecays,weights_conventional,taus_conventional);
    bayes_FreeWeightsAndTausVectors(ndecays,weights_bayes,taus_bayes);

    return (0);
}


// this function merely converts the given values, doesn't check for errors or weights summing to 1...
int bayes_ToBayesModelWeightValuesFromConventionalWeightParamValues(/* Model and params... */
                                                                    int     ndecays,
                                                                    double *weights_conventional,
                                                                    double *weights_bayes,
                                                                    double *taus,
                                                                    /* Instrument response... */
                                                                    double  delay,
                                                                    /* Transient... */
                                                                    int     nbins,
                                                                    int     nphotons,
                                                                    double  interval)
{
    int k;

    weights_bayes[0] = bayes_ToBayesModelTransformFromParamZ((float)weights_conventional[0],nbins,nphotons);

    for (k=1; k<=ndecays; k++)
    {
        weights_bayes[k] = (double)bayes_ToBayesModelWeightFromParamAAndTau((float)weights_conventional[k],(float)taus[k],
                                                                            (float)interval/(float)nbins,(float)interval,(float)delay,nphotons);
    }

    return (0);
}


int bayes_ToConventionalModelParamValuesFromBayesModelParamValues(/* Model and params... */
                                                                  int    nparams,
                                                                  float *params_bayes,
                                                                  float *params_conventional,
                                                                  /* Instrument response... */
                                                                  float  delay,
                                                                  /* Transient... */
                                                                  int    nbins,
                                                                  int    nphotons,
                                                                  float  interval)
{
    int     ndecays, i;
    double *weights_bayes, *taus_bayes, *weights_conventional, *taus_conventional;

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_bayes,&taus_bayes);
    bayes_PopulateWeightsAndTausVectorsFromParamVector(nparams,params_bayes,weights_bayes,taus_bayes);

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_conventional,&taus_conventional);

    weights_conventional[0] = bayes_FromBayesModelTransformToParamZ((float)weights_bayes[0],nbins,nphotons);

    for (i=1; i<=ndecays; i++)
    {
        weights_conventional[i] = (double)bayes_FromBayesModelWeightAndTauToParamA((float)weights_bayes[i],(float)taus_bayes[i],
                                                                            interval/(float)nbins,interval,delay,nphotons);

        taus_conventional[i]    = taus_bayes[i];
    }

    bayes_PopulateParamVectorFromWeightsAndTausVectors(params_conventional,nparams,weights_conventional,taus_conventional);

    bayes_FreeWeightsAndTausVectors(ndecays,weights_conventional,taus_conventional);
    bayes_FreeWeightsAndTausVectors(ndecays,weights_bayes,taus_bayes);

    return (0);
}


int bayes_ToConventionalModelParamValuesFromBayesModelParamValuesFreeParamsOnly(/* Model and params... */
                                                                                int    nparams,
                                                                                int   *isweightfixed,
                                                                                int   *istaufixed,
                                                                                float *params_bayes,
                                                                                float *params_conventional,
                                                                                /* Instrument response... */
                                                                                float  delay,
                                                                                /* Transient... */
                                                                                int    nbins,
                                                                                int    nphotons,
                                                                                float  interval)
{
    int     ndecays, i;
    double *weights_bayes, *taus_bayes, *weights_conventional, *taus_conventional;

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_bayes,&taus_bayes);
    bayes_PopulateWeightsAndTausVectorsFromParamVector(nparams,params_bayes,weights_bayes,taus_bayes);

    bayes_AllocateWeightsAndTausVectors(nparams,&ndecays,&weights_conventional,&taus_conventional);

    if (BAYES_PARAM_VALUE_USER_FIXED != isweightfixed[0])
        weights_conventional[0] = bayes_FromBayesModelTransformToParamZ((float)weights_bayes[0],nbins,nphotons);

    for (i=1; i<=ndecays; i++)
    {
        if (BAYES_PARAM_VALUE_USER_FIXED != isweightfixed[i])
            weights_conventional[i] = (double)bayes_FromBayesModelWeightAndTauToParamA((float)weights_bayes[i],(float)taus_bayes[i],
                                                                                       interval/(float)nbins,interval,delay,nphotons);

        if (BAYES_PARAM_VALUE_USER_FIXED != istaufixed[i])
            taus_conventional[i] = taus_bayes[i];
    }

    bayes_PopulateParamVectorFromWeightsAndTausVectors(params_conventional,nparams,weights_conventional,taus_conventional);

    bayes_FreeWeightsAndTausVectors(ndecays,weights_conventional,taus_conventional);
    bayes_FreeWeightsAndTausVectors(ndecays,weights_bayes,taus_bayes);

    return (0);
}



/*=================================================================================*/
/*                                                                                 */
/*              ROUTINES TO ISOLATE WEIGHT AND LIFETIMES PARAMETERS                */
/*                                                                                 */
/*=================================================================================*/

//this fct takes the vector that indicates the parameters the user has fixed
//it populates the vector that indicates the parameters that are either fixed by the user or the model
//also populates vectors indicating whether the weights and lifetimes are fixed, and why (i.e. user or model)

int bayes_AllocateWeightsAndTausVectors(int       nparams,
                                        int      *ndecays,
                                        double  **weights,
                                        double  **taus)
{
    int ell;

    if (nparams<=0)
        return (-1);

    ell = nparams/2; /* The number of decay components, ell... */

    *ndecays = ell;   

    *weights = Bayes_dvector(0,ell); /* 0-index accounts for constant background, w0... */
    *taus    = Bayes_dvector(1,ell);

    return (0);
}


int bayes_AllocateWeightsAndTausIsFixedVectors(int   ndecays,
                                               int **isweightfixed,
                                               int **istaufixed)
{
    if (ndecays<=0)
        return (-1);

    *isweightfixed = Bayes_ivector(0,ndecays);
    *istaufixed    = Bayes_ivector(1,ndecays);

    return (0);
}


int bayes_PopulateWeightsAndTausVectorsFromParamVector(int     nparams,
                                                       float  *params,
                                                       double *weights,
                                                       double *taus)
{
    int i, j, k;

    if ((nparams<=0) || (!params) || (!weights) || (!taus))
        return (-1);

    i = 0; j = 1; k = 0;

    weights[i] = params[k];
    i++; k++;

    do
    {
        weights[i] = params[k];
        i++; k++;
        taus[j]    = params[k];
        j++; k++;
    }
    while (k<nparams);

    return (0);
}


int bayes_PopulateParamVectorFromWeightsAndTausVectors(float  *params,
                                                       int     nparams,
                                                       double *weights,
                                                       double *taus)
{
    int i, j, k;

    i = 0; j = 1; k = 0;

    params[k] = (float)weights[i];
    i++; k++;

    do
    {
        params[k] = (float)weights[i];
        i++; k++;
        params[k] = (float)taus[j];
        j++; k++;
    }
    while (k<nparams);

    return (0);
}



void bayes_FreeWeightsAndTausVectors(int     ndecays,
                                     double *weights,
                                     double *taus)
{
    free_Bayes_dvector(weights,0,ndecays);
    free_Bayes_dvector(taus,1,ndecays);
}


void bayes_FreeWeightsAndTausFixedVectors(int     ndecays,
                                          int    *isweightfree,
                                          int    *istaufree)
{
    free_Bayes_ivector(isweightfree,0,ndecays);
    free_Bayes_ivector(istaufree,1,ndecays);
}


/*=================================================================================*/
/*                                                                                 */
/*        ROUTINES TO DETERMINE AND ISOLATE FREE AND FIXED MODEL PARAMETERS        */
/*                                                                                 */
/*=================================================================================*/


int bayes_DetermineModelParamsFixed(int   nparams,
                                    int  *isparamuserfixed,
                                    int  *nparamsfixed,
                                    int  *isparamfixed,
                                    int   ndecays,                                    
                                    int  *isweightfixed,
                                    int  *nweightsfixed,
                                    int  *istaufixed,
                                    int  *ntausfixed)
{
    int i, j, k;

    if ((nparams<=0) ||
        (!isparamuserfixed) || (!isparamfixed) ||
        (!isweightfixed) || (!istaufixed))
        return (-1);

    /* First strip out the parameters that are fixed by user choice... */
    i = 0; //weights is 0-based
    j = 1; //taus is 1-based
    k = 0; //params is 0-based

    isweightfixed[i] = isparamuserfixed[k]; //the background fraction, w0...
    i++; k++;

    do
    {
        isweightfixed[i] = isparamuserfixed[k]; //decay component weight...
        i++; k++;
        istaufixed[j]    = isparamuserfixed[k]; //...and associated lifetime
        j++; k++;
    }
    while (k<nparams);

    /* Now update the weights fixed map to incorporate a weight 'fixed' by model definition... */
    for (i=0; i<=ndecays; i++)
    {
        //select the lowest index weight that's not user fixed to be model fixed...
        if (BAYES_PARAM_VALUE_FREE == isweightfixed[i])
        {
            isweightfixed[i] = BAYES_PARAM_VALUE_MODEL_FIXED;
            break;
        }
    }

    /* Now combine the fixed weights and taus map to a fixed parameter map... */
    i = 0; //weights is 0-based
    j = 1; //taus is 1-based
    k = 0; //params is 0-based

    isparamfixed[k] = isweightfixed[i];
    i++; k++;

    do
    {
        isparamfixed[k] = isweightfixed[i]; //decay component weight...
        i++; k++;
        isparamfixed[k] = istaufixed[j]; //...and associated lifetime
        j++; k++;
    }
    while (k<nparams);

    (*nweightsfixed) = 0;
    (*ntausfixed)    = 0;

    for (i=0; i<=ndecays; i++)
        if (BAYES_PARAM_VALUE_FREE != isweightfixed[i])
            (*nweightsfixed)++;

    for (j=1; j<=ndecays; j++)
        if (BAYES_PARAM_VALUE_FREE != istaufixed[j])
            (*ntausfixed)++;

    *nparamsfixed = *nweightsfixed + *ntausfixed;

    return (0);
}


int bayes_AllocateFreeAndFixedParamVectors(int      nxfree,
                                           double **xfree,
                                           int      nxfixed,
                                           double **xfixed)
{
    if (nxfree>0)
        *xfree  = Bayes_dvector(1,nxfree);
    else
        *xfree  = NULL;

    if (nxfixed>0)
        *xfixed = Bayes_dvector(1,nxfixed);
    else
        *xfixed = NULL;

    return (0);
}


void bayes_FreeFreeAndFixedParamVectors(int      nxfree,
                                        double *xfree,
                                        int      nxfixed,
                                        double *xfixed)
{
    free_Bayes_dvector(xfree,1,nxfree);
    free_Bayes_dvector(xfixed,1,nxfixed);
}


int bayes_PopulateParamVectorFromFreeAndFixedVectors(float  *x,
                                                     int     nx,
                                                     double *xfree,
                                                     int     nxfree,
                                                     double *xfixed,
                                                     int     nxfixed,
                                                     int    *isfixed)
{
    int i, j, k;

    if ((!isfixed) || (!x) || (nx<=0) || (nx!=nxfree+nxfixed))
        return (-1);

    if ((nxfree>0 && !xfree) || (nxfixed>0 && !xfixed))
        return (-2);

    for (i=0, j=1, k=1; i<nx; i++)
    {
        if (isfixed[i])
        {
            x[i] = (float)xfixed[j];
            j++;
        }
        else
        {
            x[i] = (float)xfree[k];
            k++;
        } 
    }

    return (0);
}


int bayes_PopulateFreeAndFixedVectorsFromParamVector(float  *x,
                                                     int     nx,
                                                     double *xfree,
                                                     int     nxfree,
                                                     double *xfixed,
                                                     int     nxfixed,
                                                     int    *isfixed)
{
    int i, j, k;

    if ((!isfixed) || (!x) || (nx<=0) || (nx!=nxfree+nxfixed))
        return (-1);

    if ((nxfree>0 && !xfree) || (nxfixed>0 && !xfixed))
        return (-2);

    for (i=0, j=1, k=1; i<nx; i++)
    {
        if (isfixed[i])
        {
            xfixed[j] = x[i];
            j++;
        }
        else
        {
            xfree[k] = x[i];
            k++;
        } 
    }

    return (0);
}


//all but the model defined weight value have been set...
//...update vector to include model defined value
int bayes_UpdateWeightsVectorModelDefinedValue(double *weights,
                                               int     nweights,
                                               int    *isweightfixed)
{
    double remaining;
    int    i, j;

    for (i=0, remaining=1.0, j=0; i<nweights; i++)
    {
        if (BAYES_PARAM_VALUE_MODEL_FIXED != isweightfixed[i])
            remaining -= weights[i];
        else
            j = i;
    }

    weights[j] = remaining;

    return (0);
}



/*=================================================================================*/
/*                                                                                 */
/*             ROUTINES FOR RANDOM INITIALIZATION OF MODEL PARAMETERS              */
/*                                                                                 */
/*=================================================================================*/


int bayes_RandomlyInitWeightsVector(double *weights,
                                    int     nweights,
                                    int     nweightsfixed,
                                    int    *isweightfixed)
{
    double remaining, temp;
    int    i;

    for (i=0, remaining=1.0; i<nweights; i++)
        if (BAYES_PARAM_VALUE_USER_FIXED == isweightfixed[i])
        {
            remaining  -= weights[i];
            //weights[i]  = -1.0;
        }

    rand_InitializeRandomSeed();

    for (i=0; i<nweights; i++)
    {
        if (!isweightfixed[i])
        {
            if (BAYES_PARAM_VALUE_MODEL_FIXED != isweightfixed[i])
            {
                temp = remaining*rand_RandomDouble();
                    
                weights[i]  = temp;
                remaining  -= temp;
            }
            /* else...                                                       */
            /* Invalid value carried through and set appropriately according */
            /* to other weights when parameter vector is reconstructed...    */
        }
    }

    return (0);
}


int bayes_RandomlyInitTausVector(double *taus,
                                 int     ntaus,
                                 int    *istaufixed,
                                 double  interval)
{
    int i;

    rand_InitializeRandomSeed();

    for (i=1; i<=ntaus; i++)
    {
        if (!istaufixed[i])
        {                
            taus[i] = interval*rand_RandomDouble();
        }
    }

    return (0);
}



/*=================================================================================*/
/*                                                                                 */
/*                 MISCELLANEOUS MODEL PARAMETER RELATED TOOLS                     */
/*                                                                                 */
/*=================================================================================*/


int bayes_OrderDecaysByDecreasingLifetimes(int ndecays, double *weights, double *taus)
{
    int    i, j;
    double w, t;

    if ((!weights) || (!taus))
        return (-1);

    /* No need to worry about order of sorting routine here... */
    for (j=2; j<=ndecays; j++)
    {
        t = taus[j];
        w = weights[j];
        i = j-1;

        while ((i>0) && (taus[i]<t))
        {
            taus[i+1]    = taus[i];
            weights[i+1] = weights[i];
            i--;
        }

        taus[i+1]    = t;
        weights[i+1] = w;
    }

    return (0);
}






/////////////////////////////// older code - may look at revising some of this using the generic fcts. above...


float bayes_ComputeRawBackgroundZ(float w0,
                                  float binsize,
                                  float interval,
                                  int nphotons)
{
    double value;

    value = w0/interval;

    return ((float)(value*binsize*(float)nphotons));
}


float bayes_ComputeRawAmplitudeA(float w0,
                                 float w1,
                                 float binsize,
                                 float interval,
                                 float delay,
                                 int nphotons)
{
	double num, den1, value;

	num  = 1.0/w1; 
	den1 = 1.0-exp((delay-interval)/w1);
    
//	value=w0/TIMEWINDOW+(1.0-w0)*(num/den1);  // TC ? - surely this includes the background
	value=(1.0-w0)*(num/den1);

	return ((float)(nphotons*binsize*value));
}


float bayes_ModelTransformComputePhotonCount(float interval,
				                             float w0,
                                             float likelihood,
				                             float binsize,
				                             int   nphotons)
{
	double value, dw0, dp;

    dw0        = (double)w0;
	dp         = (double)nphotons;

    value = dw0/interval + (1.0-dw0)*likelihood;

    return ((float)(dp*binsize*value));
}

#if 0 // not used
//this is a hack for initial integration
double F(double s, double w, double UC, double SIGMA, double TIMEWINDOW)
{
    double num,den1,den2;
	double ROOTTWO;
   
	if(s<0.0) return(0.0);
	if(s>TIMEWINDOW) return(0.0);
   
	if(SIGMA<=0.0)
	{
        if(s<UC) return(0.0);
      
		num  = exp((UC-s)/w)/w;
		den1 = 1.0-exp((UC-TIMEWINDOW)/w);
      
		return(num/den1);
    }
	
	ROOTTWO = sqrt(2.0);
	
    num  = erf((UC/SIGMA+SIGMA/w)/ROOTTWO)-erf(((UC-s)/SIGMA+SIGMA/w)/ROOTTWO);
    num  = num*exp((UC-s)/w)/w;
    den1 = erf((TIMEWINDOW-UC)/(SIGMA*ROOTTWO))+erf(UC/(SIGMA*ROOTTWO));
    den1 = den1*exp(-0.5*(SIGMA/w)*(SIGMA/w));
    den2 = erf((UC/SIGMA+SIGMA/w)/ROOTTWO)-erf(((UC-TIMEWINDOW)/SIGMA+SIGMA/w)/ROOTTWO);
    den2 = den2*exp((UC-TIMEWINDOW)/w);
   /* { printf("F(%lf,%lf): num=%lf, den=%lf",s,w,num,den1-den2);
     puts(""); getchar(); }  */
   return(num/(den1-den2));
}


float photoncount(float uc, float sigma, float T, float t, float w0, float w1, float binsize, int p)
{
    double value, dw0, dw1, dt, dp;
	
    double UC=(double)uc;
	double SIGMA=(double)sigma;
	double TIMEWINDOW=(double)T;
    
	dw0=(double)w0;
	dw1=(double)w1;
	dt=(double)t;
	dp=(double)p;
	
    value=dw0/TIMEWINDOW+(1.0-dw0)*F(dt,dw1,uc,sigma,T);
    
	return ((float)(dp*binsize*value));
}



float photoncount_multiple(float uc, float sigma, float T, float t, float *times, float *weights, int L, float binsize, int p)
{
    double value,dw0=1.0,dt,dp;
    int i;
    double UC=(double)uc;
    double SIGMA=(double)sigma;
    double TIMEWINDOW=(double)T;
   
	for (i=0;i<L;i++)
		dw0 -= (double)weights[i];
	
    if (dw0<0.0)
		dw0=0.0;
	
	if(dw0>1.0)
		dw0=1.0;
	
    dt=(double)t;
	dp=(double)p;
   
	value = dw0/TIMEWINDOW;
    
	for (i=0;i<L;i++)
		value += ((double)weights[i]*F(dt,(double)times[i],uc,sigma,T));
	
    return ((float)(dp*binsize*value));
}
#endif


#if 0
/* Copied from TRFitting.c */
void scaleDataAccordingToSignalNoise (float *data, int n, float *signal)
{
	int i;
	for (i=0; i<n; i++) data[i] /= sqrt(fabs(signal[i]));
}
#endif

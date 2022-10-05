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
#include "stdio.h"
#include "matrices.h"
#include "bayes_Sizes.h"
#include "bayes_Types.h"
#include "bayes_DataManagement.h"
#include "bayes_ModelSelection.h"
#include "bayes_ModelTransformTools.h"
#include "bayes_Interface.h"
#include "bayes_RapidBayesDecayAnalysis.h"
#include "bayes_DistributionFctsBinLikelihoods.h"

#include "bayes_MultiExpAnalysisBinLikelihoods.h"


//#include "nrutil.h"
//#define TINY 1.0e-20;

// for determinant computation
int ludcmp(double **a, int n, int *indx, double *d)
{
	int     i,imax,j,k;
	double  big,dum,sum,temp;
	double *vv;

	vv = Bayes_dvector(1,n);
	
    *d=1.0;
	
    for (i=1;i<=n;i++)
    {
        big=0.0;

        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        
        if (big == 0.0)
            return (-1);//matrices_error("Singular matrix in routine ludcmp");
        
        vv[i]=1.0/big;
	}

    for (j=1;j<=n;j++)
    {
        for (i=1;i<j;i++)
        {
            sum = a[i][j];
        
            for (k=1;k<i;k++)
                sum -= a[i][k]*a[k][j];
        
            a[i][j]=sum;
        }

        big=0.0;
		
        for (i=j;i<=n;i++)
        {
            sum = a[i][j];

            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            
            a[i][j]=sum;
            
            if ( (dum=vv[i]*fabs(sum)) >= big)
            {
                big=dum;
                imax=i;
            }
        }

        if (j != imax)
        {
            for (k=1;k<=n;k++)
            {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }

            *d = -(*d);
            vv[imax]=vv[j];
        }

        indx[j]=imax;

        if (a[j][j] == 0.0)
            a[j][j] = 1.0e-20/*TINY*/;
		
        if (j != n)
        {
			dum=1.0/(a[j][j]);

			for (i=j+1;i<=n;i++)
                a[i][j] *= dum;
		}
	}

	free_Bayes_dvector(vv,1,n);

    return (0);
}
//#undef TINY

double bayes_ComputeDeterminantValue(double **A, int ndim)
{
    int    *index, i;
    double  d;

    index = Bayes_ivector(1,ndim);

    ludcmp(A,ndim,index,&d);

    for (i=1; i<=ndim; i++)
        d *= A[i][i];

	free_Bayes_ivector(index,1,ndim);

    return (d);
}


//direct hessian computation

//uses gaussian approximation
//hessian element variants
// diagonal elements... 
//   d^2/dz^2, 
//   d^2/ds^2 
// off-diagonal elements...
//   d^2/dz_i dz_j
//   d^2/dz_i ds_i
//   d^2/dz_i ds_j
//   d^2/ds_i ds_j





/*===================================================================*/
/*                                                                   */
/*                    First derivative calulators                    */
/*                                                                   */
/*===================================================================*/

int bayes_DataLikelihood1stDerivativesWrtWeight(double          *derivatives,
                                                int              derivativeindex,
                                                int              ndecays,
                                                double          *weights,
                                                double          *taus,
                                                double         **fluorescencephotonlikelihoods,
                                                int              nbins,
                                                double          *binwalls,
                                                int             *data,
                                                double           interval,
                                                double           modperiod,
                                                BayesInstrRsp_t *instr)
{
    int    bin;
    double bjoverT, bL, bH, Fx, Lambda;

    //if (interval==modperiod)
        Lambda = 1.0;
    //else
      //  Lambda = ;

    for (bin=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            bL               = binwalls[bin];
            bH               = binwalls[bin+1];
            bjoverT          = (bH-bL)/interval;
            Fx               = fluorescencephotonlikelihoods[derivativeindex][bin];

            derivatives[bin] = Fx/Lambda-bjoverT;
        }
        else
        {
            derivatives[bin] = 0.0;
        }
    }
    
    return (0);
}


int bayes_DataLikelihood1stDerivativesWrtLifetime(double          *derivatives,
                                              //    int              derivativeindex,
                                                  int              ndecays,
                                                  double           weight,
                                                  double           tau,
                                              //    double         **photonlikelihoods,
                                                  int              nbins,
                                                  double          *binwalls,
                                                  int             *data,
                                                  double           interval,
                                                  double           modperiod,
                                                  BayesInstrRsp_t *instr)
{
    int      bin, valid, i, ell, ellmax;
    double   bL, bH, *chiprimeL, *chiprimeH, valbin, Lambda;
    double   psi, theta, phi, psiprime, thetaprime;
    double   gammaprime[BAYES_INSTR_RSP_MAX_COMPONENTS], delay, width, cutoff;

    ellmax = 2;

    //if (interval==modperiod)
        Lambda = 1.0;
    //else
      //  Lambda = ;

    for (i=0; i<instr->ninstr; i++)
    {
        weight = instr->params[i].weight;
        delay  = instr->params[i].delay;
        width  = instr->params[i].width;
        cutoff = instr->params[i].cutoff;

        gammaprime[i] = weight/(1.0+erf((delay-cutoff)/(width*ROOTTWO)));
    }

    chiprimeL = Bayes_dvector(0,instr->ninstr-1);
    chiprimeH = Bayes_dvector(0,instr->ninstr-1);

    //tau  = taus[derivativeindex];

    for (bin=0,valid=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            bL = binwalls[bin];
            bH = binwalls[bin+1];

            if (valid)
            {
                for (i=0; i<instr->ninstr; i++)
                    chiprimeL[i] = chiprimeH[i];

                /* Summation over 'ell' for each instrument response component... */
                for (i=0; i<instr->ninstr; i++)
                {
                    delay  = instr->params[i].delay;
                    width  = instr->params[i].width;
                    cutoff = instr->params[i].cutoff;

                    for (ell=0,chiprimeH[i]=0.0; ell<ellmax; ell++)
                    {
                        psi        = (delay-modperiod*(double)ell-bH)/tau + width*width/(2.0*tau*tau);
                        theta      = (delay-modperiod*(double)ell-bH+width*width/tau)/(width*ROOTTWO);
                        phi        = (delay-cutoff+width*width/tau)/(width*ROOTTWO);
                        psiprime   = -(delay-modperiod*(double)ell-bH+width*width/tau)/tau/tau;
                        thetaprime = -width/ROOTTWO/tau/tau;

                        chiprimeH[i] += psiprime*exp(psi)*(erf(theta)-erf(phi)) +
                                        (2.0/ROOTPI)*thetaprime*(exp(psi-theta*theta)-exp(psi-phi*phi));
                    }
                }
            }
            else
            {
                for (i=0,valbin=0.0; i<instr->ninstr; i++)
                {
                    delay  = instr->params[i].delay;
                    width  = instr->params[i].width;
                    cutoff = instr->params[i].cutoff;

                    for (ell=0,chiprimeL[i]=0.0,chiprimeH[i]=0.0; ell<ellmax; ell++)
                    {
                        psi        = (delay-modperiod*(double)ell-bL)/tau + width*width/(2.0*tau*tau);
                        theta      = (delay-modperiod*(double)ell-bL+width*width/tau)/(width*ROOTTWO);
                        phi        = (delay-cutoff+width*width/tau)/(width*ROOTTWO);
                        psiprime   = -(delay-modperiod*(double)ell-bL+width*width/tau)/tau/tau;
                        thetaprime = -width/ROOTTWO/tau/tau;
                        
                        chiprimeL[i] += psiprime*exp(psi)*(erf(theta)-erf(phi)) +
                                        (2.0/ROOTPI)*thetaprime*(exp(psi-theta*theta)-exp(psi-phi*phi));

                        psi        = (delay-modperiod*(double)ell-bH)/tau + width*width/(2.0*tau*tau);
                        theta      = (delay-modperiod*(double)ell-bH+width*width/tau)/(width*ROOTTWO);
                        psiprime   = -(delay-modperiod*(double)ell-bH+width*width/tau)/tau/tau;
                        
                        chiprimeH[i] += psiprime*exp(psi)*(erf(theta)-erf(phi)) +
                                        (2.0/ROOTPI)*thetaprime*(exp(psi-theta*theta)-exp(psi-phi*phi));
                    }
                }

                valid = 1;
            }

            /* Now do the sum over all instrument response components... */
            for (i=0,valbin=0.0; i<instr->ninstr; i++)
                valbin += gammaprime[i]*(chiprimeH[i]-chiprimeL[i]);

            derivatives[bin] = weight*valbin/Lambda;
        }
        else
        {
            derivatives[bin] = 0.0;
            valid = 0;
        }
    }

    free_Bayes_dvector(chiprimeL,0,instr->ninstr-1);
    free_Bayes_dvector(chiprimeH,0,instr->ninstr-1);
    
    return (0);
}



/*===================================================================*/
/*                                                                   */
/*                   Second derivative calulators                    */
/*                                                                   */
/*===================================================================*/

int bayes_DataLikelihood2ndDerivativesWrtSingleWeight(double          *derivatives,
                                                      int              derivativeindex,
                                                      int              ndecays,
                                                      double          *weights,
                                                      double          *taus,
                                                      double         **fluorescencephotonlikelihoods,
                                                      int              nbins,
                                                      double          *binwalls,
                                                      int             *data,
                                                      double           interval,
                                                      double           modperiod,
                                                      BayesInstrRsp_t *instr)
{
    int    bin;
    double wx, oneminuswx, Lambda;

//    if (interval==modperiod)
        Lambda = 1.0;
  //  else
    //    Lambda = ;
    
    wx         = weights[derivativeindex];
    oneminuswx = 1.0-wx;

    for (bin=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]) || (1.0==Lambda))
        {
/*            bL      = binwalls[bin];
            bH      = binwalls[bin+1];
            bjoverT = (bH-bL)/interval;

            Fx      = fluorescencephotonlikelihoods[derivativeindex][bin];
            temp1   = wx*(1.0-2.0*wx)*(Fx-bjoverT);

            for (k=1,temp2=0.0; k<=ndecays; k++)
            {
                wk     = weights[k];
                Fk     = fluorescencephotonlikelihoods[k][bin];
                temp2 += wk*(Fk-bjoverT);
            }

            temp2 *= wx*(1.0-2.0*wx);

            derivatives[bin] = temp1-temp2;
*/
            derivatives[bin] = 0.0;//normalisation correction needs to be incorporated
        }
        else
        {
            derivatives[bin] = 0.0;
        }
    }
    
    return (0);
}


int bayes_DataLikelihood2ndDerivativesWrtDiffWeights(double          *derivatives,
                                                     int              weightderivativeindex_x,
                                                     int              weightderivativeindex_y,
                                                     int              ndecays,
                                                     double          *weights,
                                                     double          *taus,
                                                     double         **fluorescencephotonlikelihoods,
                                                     int              nbins,
                                                     double          *binwalls,
                                                     int             *data,
                                                     double           interval,
                                                     double           modperiod,
                                                     BayesInstrRsp_t *instr)
{
    int    bin;
    double wxwy, Lambda;

//    if (interval==modperiod)
        Lambda = 1.0;
  //  else
    //    Lambda = ;

    wxwy = weights[weightderivativeindex_x]*weights[weightderivativeindex_y];

    for (bin=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]) || (1.0==Lambda))
        {
/*
            bL      = binwalls[bin];
            bH      = binwalls[bin+1];
            bjoverT = (bH-bL)/interval;

            Fx    = fluorescencephotonlikelihoods[weightderivativeindex_x][bin];
            Fy    = fluorescencephotonlikelihoods[weightderivativeindex_y][bin];
            temp1 = -wxwy*(Fx+Fy);
            temp2 = 2.0*wxwy*bjoverT;

            derivatives[bin] = temp1+temp2;
*/
            derivatives[bin] = 0.0;//normalisation correction needs to be incorporated
        }
        else
        {
            derivatives[bin] = 0.0;
        }
    }

    return (0);
}


int bayes_DataLikelihood2ndDerivativesWrtSingleLifetime(double          *derivatives,
                                                     //   int              derivativeindex,
                                                        int              ndecays,
                                                        double           weight,
                                                        double           tau,
                                                      //  double         **photonlikelihoods,
                                                        int              nbins,
                                                        double          *binwalls,
                                                        int             *data,
                                                        double           interval,
                                                        double           modperiod,
                                                        BayesInstrRsp_t *instr)
{
    int      bin, valid, i, ell, ellmax;
    double   bL, bH, *chiprimeprimeL, *chiprimeprimeH, valbin, Lambda;
    double   psi, theta, phi, psiprime, psiprimeprime, thetaprime, thetaprimeprime;
    double   temp1, temp2, temp3, twooverrootpi;
    double   gammaprime[BAYES_INSTR_RSP_MAX_COMPONENTS], delay, width, cutoff;

//    if (interval==modperiod)
        Lambda = 1.0;
  //  else
    //    Lambda = ;

    ellmax        = 2;
    twooverrootpi = 2.0/ROOTPI;

    for (i=0; i<instr->ninstr; i++)
    {
        weight = instr->params[i].weight;
        delay  = instr->params[i].delay;
        width  = instr->params[i].width;
        cutoff = instr->params[i].cutoff;

        gammaprime[i] = weight/(1.0+erf((delay-cutoff)/(width*ROOTTWO)));
    }

    chiprimeprimeL = Bayes_dvector(0,instr->ninstr-1);
    chiprimeprimeH = Bayes_dvector(0,instr->ninstr-1);

  //  tau  = taus[derivativeindex];

    for (bin=0,valid=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            bL = binwalls[bin];
            bH = binwalls[bin+1];

            if (valid)
            {
                for (i=0; i<instr->ninstr; i++)
                    chiprimeprimeL[i] = chiprimeprimeH[i];

                /* Summation over 'ell' for each instrument response component... */
                for (i=0; i<instr->ninstr; i++)
                {
                    delay  = instr->params[i].delay;
                    width  = instr->params[i].width;
                    cutoff = instr->params[i].cutoff;

                    for (ell=0,chiprimeprimeH[i]=0.0; ell<ellmax; ell++)
                    {
                        psi             = (delay-modperiod*(double)ell-bH)/tau + width*width/(2.0*tau*tau);
                        theta           = (delay-modperiod*(double)ell-bH+width*width/tau)/(width*ROOTTWO);
                        psiprime        = -(delay-modperiod*(double)ell-bH+width*width/tau)/tau/tau;
                        psiprimeprime   = (2.0*(delay-modperiod*(double)ell-bH)+3.0*width*width/tau)/tau/tau/tau;

                        temp1 = (psiprimeprime+psiprime*psiprime)*(exp(psi)*(erf(theta)-erf(phi)));
                        temp2 = twooverrootpi*(thetaprimeprime+2.0*(thetaprime*(psiprime-theta*thetaprime)))*exp(psi-theta*theta);
                        temp3 = twooverrootpi*(thetaprimeprime+2.0*(thetaprime*(psiprime-phi*thetaprime)))*exp(psi-phi*phi);
                        
                        chiprimeprimeH[i] += temp1+temp2-temp3;
                    }
                }
            }
            else
            {
                for (i=0,valbin=0.0; i<instr->ninstr; i++)
                {
                    delay  = instr->params[i].delay;
                    width  = instr->params[i].width;
                    cutoff = instr->params[i].cutoff;

                    for (ell=0,chiprimeprimeL[i]=0.0,chiprimeprimeH[i]=0.0; ell<ellmax; ell++)
                    {
                        psi             = (delay-modperiod*(double)ell-bL)/tau + width*width/(2.0*tau*tau);
                        theta           = (delay-modperiod*(double)ell-bL+width*width/tau)/(width*ROOTTWO);
                        phi             = (delay-cutoff+width*width/tau)/(width*ROOTTWO);
                        psiprime        = -(delay-modperiod*(double)ell-bL+width*width/tau)/tau/tau;
                        psiprimeprime   = (2.0*(delay-modperiod*(double)ell-bL)+3.0*width*width/tau)/tau/tau/tau;
                        thetaprime      = -width/ROOTTWO/tau/tau;
                        thetaprimeprime = 2.0*width/ROOTTWO/tau/tau/tau;

                        temp1 = (psiprimeprime+psiprime*psiprime)*(exp(psi)*(erf(theta)-erf(phi)));
                        temp2 = twooverrootpi*(thetaprimeprime+2.0*(thetaprime*(psiprime-theta*thetaprime)))*exp(psi-theta*theta);
                        temp3 = twooverrootpi*(thetaprimeprime+2.0*(thetaprime*(psiprime-phi*thetaprime)))*exp(psi-phi*phi);
                        
                        chiprimeprimeL[i] += temp1+temp2-temp3;

                        psi             = (delay-modperiod*(double)ell-bH)/tau + width*width/(2.0*tau*tau);
                        theta           = (delay-modperiod*(double)ell-bH+width*width/tau)/(width*ROOTTWO);
                        psiprime        = -(delay-modperiod*(double)ell-bH+width*width/tau)/tau/tau;
                        psiprimeprime   = (2.0*(delay-modperiod*(double)ell-bH)+3.0*width*width/tau)/tau/tau/tau;

                        temp1 = (psiprimeprime+psiprime*psiprime)*(exp(psi)*(erf(theta)-erf(phi)));
                        temp2 = twooverrootpi*(thetaprimeprime+2.0*(thetaprime*(psiprime-theta*thetaprime)))*exp(psi-theta*theta);
                        temp3 = twooverrootpi*(thetaprimeprime+2.0*(thetaprime*(psiprime-phi*thetaprime)))*exp(psi-phi*phi);
                        
                        chiprimeprimeH[i] += temp1+temp2-temp3;
                    }
                }

                valid = 1;
            }

            /* Now do the sum over all instrument response components... */
            for (i=0,valbin=0.0; i<instr->ninstr; i++)
                valbin += gammaprime[i]*(chiprimeprimeH[i]-chiprimeprimeL[i]);

            derivatives[bin] = weight*valbin/Lambda;
        }
        else
        {
            derivatives[bin] = 0.0;
            valid = 0;
        }
    }

    free_Bayes_dvector(chiprimeprimeL,0,instr->ninstr-1);
    free_Bayes_dvector(chiprimeprimeH,0,instr->ninstr-1);
    
    return (0);
}


int bayes_DataLikelihood2ndDerivativesWrtDiffLifetimes(double *derivatives,
                                                       int     nbins)
{
    int    bin;
    double Lambda;

//    if (interval==modperiod)
        Lambda = 1.0;
  //  else
    //    Lambda = ;

    for (bin=0; bin<nbins; bin++)
    {
        derivatives[bin] = 0.0;//norm correction required
    }

    return (0);
}


int bayes_DataLikelihood2ndDerivativesWrtWeightAndLifetimeSameSignalComponents(double          *derivatives,
                                                                              //    int              derivativeindex,
                                                                               int              ndecays,
                                                                               double           weight,
                                                                               double           tau,
                                                                           //    double         **photonlikelihoods,
                                                                               int              nbins,
                                                                               double          *binwalls,
                                                                               int             *data,
                                                                               double           interval,
                                                                               double           modperiod,
                                                                               BayesInstrRsp_t *instr)
{
    int      bin, valid, i, ell, ellmax;
    double   bL, bH, *chiprimeL, *chiprimeH, valbin, Lambda;
    double   psi, theta, phi, psiprime, thetaprime;
    double   gammaprime[BAYES_INSTR_RSP_MAX_COMPONENTS], delay, width, cutoff;

    ellmax = 2;

    //if (interval==modperiod)
        Lambda = 1.0;
    //else
      //  Lambda = ;

    for (i=0; i<instr->ninstr; i++)
    {
        weight = instr->params[i].weight;
        delay  = instr->params[i].delay;
        width  = instr->params[i].width;
        cutoff = instr->params[i].cutoff;

        gammaprime[i] = weight/(1.0+erf((delay-cutoff)/(width*ROOTTWO)));
    }

    chiprimeL = Bayes_dvector(0,instr->ninstr-1);
    chiprimeH = Bayes_dvector(0,instr->ninstr-1);

    //tau  = taus[derivativeindex];

    for (bin=0,valid=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            bL = binwalls[bin];
            bH = binwalls[bin+1];

            if (valid)
            {
                for (i=0; i<instr->ninstr; i++)
                    chiprimeL[i] = chiprimeH[i];

                /* Summation over 'ell' for each instrument response component... */
                for (i=0; i<instr->ninstr; i++)
                {
                    delay  = instr->params[i].delay;
                    width  = instr->params[i].width;
                    cutoff = instr->params[i].cutoff;

                    for (ell=0,chiprimeH[i]=0.0; ell<ellmax; ell++)
                    {
                        psi        = (delay-modperiod*(double)ell-bH)/tau + width*width/(2.0*tau*tau);
                        theta      = (delay-modperiod*(double)ell-bH+width*width/tau)/(width*ROOTTWO);
                        phi        = (delay-cutoff+width*width/tau)/(width*ROOTTWO);
                        psiprime   = -(delay-modperiod*(double)ell-bH+width*width/tau)/tau/tau;
                        thetaprime = -width/ROOTTWO/tau/tau;

                        chiprimeH[i] += psiprime*exp(psi)*(erf(theta)-erf(phi)) +
                                        (2.0/ROOTPI)*thetaprime*(exp(psi-theta*theta)-exp(psi-phi*phi));
                    }
                }
            }
            else
            {
                for (i=0,valbin=0.0; i<instr->ninstr; i++)
                {
                    delay  = instr->params[i].delay;
                    width  = instr->params[i].width;
                    cutoff = instr->params[i].cutoff;

                    for (ell=0,chiprimeL[i]=0.0,chiprimeH[i]=0.0; ell<ellmax; ell++)
                    {
                        psi        = (delay-modperiod*(double)ell-bL)/tau + width*width/(2.0*tau*tau);
                        theta      = (delay-modperiod*(double)ell-bL+width*width/tau)/(width*ROOTTWO);
                        phi        = (delay-cutoff+width*width/tau)/(width*ROOTTWO);
                        psiprime   = -(delay-modperiod*(double)ell-bL+width*width/tau)/tau/tau;
                        thetaprime = -width/ROOTTWO/tau/tau;
                        
                        chiprimeL[i] += psiprime*exp(psi)*(erf(theta)-erf(phi)) +
                                        (2.0/ROOTPI)*thetaprime*(exp(psi-theta*theta)-exp(psi-phi*phi));

                        psi        = (delay-modperiod*(double)ell-bH)/tau + width*width/(2.0*tau*tau);
                        theta      = (delay-modperiod*(double)ell-bH+width*width/tau)/(width*ROOTTWO);
                        psiprime   = -(delay-modperiod*(double)ell-bH+width*width/tau)/tau/tau;
                        
                        chiprimeH[i] += psiprime*exp(psi)*(erf(theta)-erf(phi)) +
                                        (2.0/ROOTPI)*thetaprime*(exp(psi-theta*theta)-exp(psi-phi*phi));
                    }
                }

                valid = 1;
            }

            /* Now do the sum over all instrument response components... */
            for (i=0,valbin=0.0; i<instr->ninstr; i++)
                valbin += gammaprime[i]*(chiprimeH[i]-chiprimeL[i]);

            derivatives[bin] = valbin/Lambda;
        }
        else
        {
            derivatives[bin] = 0.0;
            valid = 0;
        }
    }

    free_Bayes_dvector(chiprimeL,0,instr->ninstr-1);
    free_Bayes_dvector(chiprimeH,0,instr->ninstr-1);
    
    return (0);
}


int bayes_DataLikelihood2ndDerivativesWrtWeightAndLifetimeDiffSignalComponents(double *derivatives,
                                                                               double  weight,
                                                                               double *likelihood1stderivativeswrttau,
                                                                               int     nbins)
{
    int bin;



    for (bin=0; bin<nbins; bin++)
    {
        derivatives[bin] = 0.0; //lambda correction reqd
    }
    
    return (0);
}


















/*===================================================================*/
/*                                                                   */
/*           Data-likelihood second derivate calculators             */
/*                                                                   */
/*===================================================================*/


double bayes_ComputeHessianElementDataLikelihoodWrtSingleWeight(int              weightderivativeindex,
                                                                int              ndecays,
                                                                double          *weights,
                                                                double          *taus,
                                                                int              nbins,
                                                                double          *binwalls,
                                                                int             *data,
                                                                double          *binlikelihoods,
                                                                double         **fluorescencephotonlikelihoods,
                                                                double           interval,
                                                                double           modperiod,
                                                                BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesx, *secondderivativesxx;

    firstderivativesx   = Bayes_dvector(0,nbins-1);
    secondderivativesxx = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtWeight(firstderivativesx,weightderivativeindex,
                                                ndecays,weights,taus,
                                                fluorescencephotonlikelihoods,
                                                nbins,binwalls,data,
                                                interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtSingleWeight(secondderivativesxx,weightderivativeindex,
                                                      ndecays,weights,taus,
                                                      fluorescencephotonlikelihoods,
                                                      nbins,binwalls,data,
                                                      interval,modperiod,instr);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = (binlikelihoods[bin]*secondderivativesxx[bin])-(firstderivativesx[bin]*firstderivativesx[bin]);
            temp /= (binlikelihoods[bin]*binlikelihoods[bin]);
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(secondderivativesxx,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementDataLikelihoodWrtMixedWeights(int              weightderivativeindex_x,
                                                                int              weightderivativeindex_y,
                                                                int              ndecays,
                                                                double          *weights,
                                                                double          *taus,
                                                                int              nbins,
                                                                double          *binwalls,
                                                                int             *data,
                                                                double          *binlikelihoods,
                                                                double         **fluorescencephotonlikelihoods,
                                                                double           interval,
                                                                double           modperiod,
                                                                BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesx, *firstderivativesy, *secondderivativesxy;

    firstderivativesx   = Bayes_dvector(0,nbins-1);
    firstderivativesy   = Bayes_dvector(0,nbins-1);
    secondderivativesxy = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtWeight(firstderivativesx,weightderivativeindex_x,
                                                ndecays,weights,taus,
                                                fluorescencephotonlikelihoods,
                                                nbins,binwalls,data,
                                                interval,modperiod,instr);

    bayes_DataLikelihood1stDerivativesWrtWeight(firstderivativesy,weightderivativeindex_y,
                                                ndecays,weights,taus,
                                                fluorescencephotonlikelihoods,
                                                nbins,binwalls,data,
                                                interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtDiffWeights(secondderivativesxy,weightderivativeindex_x,weightderivativeindex_y,
                                                     ndecays,weights,taus,
                                                     fluorescencephotonlikelihoods,
                                                     nbins,binwalls,data,
                                                     interval,modperiod,instr);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesxy[bin]-firstderivativesx[bin]*firstderivativesy[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(firstderivativesy,0,nbins-1);
    free_Bayes_dvector(secondderivativesxy,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementDataLikelihoodWrtSingleLifetime(int              ndecays,
                                                                  double           weight,
                                                                  double           tau,
                                                                  int              nbins,
                                                                  double          *binwalls,
                                                                  int             *data,
                                                                  double          *binlikelihoods,
                                                                  double           interval,
                                                                  double           modperiod,
                                                                  BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesx, *secondderivativesxx;

    firstderivativesx   = Bayes_dvector(0,nbins-1);
    secondderivativesxx = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativesx,/*derivativeindex_x,*/
                                                  ndecays,weight,tau,
                                                  /*fluorescencephotonlikelihoods,*/
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtSingleLifetime(secondderivativesxx,/*derivativeindex_x,*/
                                                        ndecays,weight,tau,
                                                        //double         **photonlikelihoods,
                                                        nbins,binwalls,data,
                                                        interval,modperiod,instr);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesxx[bin]-firstderivativesx[bin]*firstderivativesx[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(secondderivativesxx,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementDataLikelihoodWrtMixedLifetimes(int              tauderivativeindex_x,
                                                                  int              tauderivativeindex_y,
                                                                  int              ndecays,
                                                                  double          *weights,
                                                                  double          *taus,
                                                                  double          *binlikelihoods,
                                                                  int              nbins,
                                                                  double          *binwalls,
                                                                  int             *data,
                                                                  double           interval,
                                                                  double           modperiod,
                                                                  BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesx, *firstderivativesy, *secondderivativesxy;

    firstderivativesx   = Bayes_dvector(0,nbins-1);
    firstderivativesy   = Bayes_dvector(0,nbins-1);
    secondderivativesxy = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativesx,
                                                  ndecays,
                                                  weights[tauderivativeindex_x],
                                                  taus[tauderivativeindex_x],
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativesy,
                                                  ndecays,
                                                  weights[tauderivativeindex_y],
                                                  taus[tauderivativeindex_y],
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtDiffLifetimes(secondderivativesxy,nbins);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesxy[bin]-firstderivativesx[bin]*firstderivativesy[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(firstderivativesy,0,nbins-1);
    free_Bayes_dvector(secondderivativesxy,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementDataLikelihoodWrtWeightAndLifetimeSameSignalComponent(int              derivativeindex_x,
                                                                                        int              ndecays,
                                                                                        double          *weights,
                                                                                        double          *taus,
                                                                                        int              nbins,
                                                                                        double          *binwalls,
                                                                                        int             *data,
                                                                                        double          *binlikelihoods,
                                                                                        double         **fluorescencephotonlikelihoods,
                                                                                        double           interval,
                                                                                        double           modperiod,
                                                                                        BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesx, *firstderivativesy, *secondderivativesxy;

    firstderivativesx   = Bayes_dvector(0,nbins-1);
    firstderivativesy   = Bayes_dvector(0,nbins-1);
    secondderivativesxy = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtWeight(firstderivativesx,derivativeindex_x,
                                                ndecays,weights,taus,
                                                fluorescencephotonlikelihoods,
                                                nbins,binwalls,data,
                                                interval,modperiod,instr);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativesy,
                                                  ndecays,weights[derivativeindex_x],taus[derivativeindex_x],
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtWeightAndLifetimeSameSignalComponents(secondderivativesxy,
                                                                               ndecays,weights[derivativeindex_x],taus[derivativeindex_x],
                                                                               nbins,binwalls,data,
                                                                               interval,modperiod,instr);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesxy[bin]-firstderivativesx[bin]*firstderivativesy[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(firstderivativesy,0,nbins-1);
    free_Bayes_dvector(secondderivativesxy,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementDataLikelihoodWrtWeightAndLifetimeDiffSignalComponents(int              weightderivativeindex,
                                                                                         int              tauderivativeindex,
                                                                                         int              ndecays,
                                                                                         double          *weights,
                                                                                         double          *taus,
                                                                                         int              nbins,
                                                                                         double          *binwalls,
                                                                                         int             *data,
                                                                                         double          *binlikelihoods,
                                                                                         double         **fluorescencephotonlikelihoods,
                                                                                         double           interval,
                                                                                         double           modperiod,
                                                                                         BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesweight, *firstderivativestau, *secondderivativesweightandtau;

    firstderivativesweight        = Bayes_dvector(0,nbins-1);
    firstderivativestau           = Bayes_dvector(0,nbins-1);
    secondderivativesweightandtau = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtWeight(firstderivativesweight,weightderivativeindex,
                                                ndecays,weights,taus,
                                                fluorescencephotonlikelihoods,
                                                nbins,binwalls,data,
                                                interval,modperiod,instr);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativestau,
                                                  ndecays,weights[tauderivativeindex],taus[tauderivativeindex],
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);



    bayes_DataLikelihood2ndDerivativesWrtWeightAndLifetimeDiffSignalComponents(secondderivativesweightandtau,
                                                                               weights[weightderivativeindex],
                                                                               firstderivativestau,
                                                                               nbins);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesweightandtau[bin]-firstderivativesweight[bin]*firstderivativestau[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesweight,0,nbins-1);
    free_Bayes_dvector(firstderivativestau,0,nbins-1);
    free_Bayes_dvector(secondderivativesweightandtau,0,nbins-1);

    return (val);
}



/*===================================================================*/
/*                                                                   */
/*                   Hessian element calculators                     */
/*                                                                   */
/*===================================================================*/


double bayes_ComputeHessianElementWrtSingleWeight(int              decayindex_x,
                                                  int              decayindex_y,
                                                  int              ndecays,
                                                  double          *weights,
                                                  double          *taus,
                                            //      double          *hyperparams,
                                                  int              nbins,
                                                  double          *binwalls,
                                                  double          *binlikelihoods,
                                                  double         **fluorescencephotonlikelihoods,
                                                  int             *data,
                                                  double           interval,
                                                  double           modperiod,
                                                  BayesInstrRsp_t *instr)
{
    int    bin, cj, weightderivativeindex;
    double temp, val, *firstderivativesx, *secondderivativesxx;

    firstderivativesx     = Bayes_dvector(0,nbins-1);
    secondderivativesxx   = Bayes_dvector(0,nbins-1);
    weightderivativeindex = decayindex_x;

    bayes_DataLikelihood1stDerivativesWrtWeight(      firstderivativesx,weightderivativeindex,
                                                      ndecays,weights,taus,
                                                      fluorescencephotonlikelihoods,
                                                      nbins,binwalls,data,
                                                      interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtSingleWeight(secondderivativesxx,weightderivativeindex,
                                                      ndecays,weights,taus,
                                                      fluorescencephotonlikelihoods,
                                                      nbins,binwalls,data,
                                                      interval,modperiod,instr);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = (binlikelihoods[bin]*secondderivativesxx[bin])-(firstderivativesx[bin]*firstderivativesx[bin]);
            temp /= (binlikelihoods[bin]*binlikelihoods[bin]);
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(secondderivativesxx,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementWrtMixedWeights(int              decayindex_x,
                                                  int              decayindex_y,
                                                  int              ndecays,
                                                  double          *weights,
                                                  double          *taus,
                                              //    double          *hyperparams,
                                                  int              nbins,
                                                  double          *binwalls,
                                                  double          *binlikelihoods,
                                                  double         **fluorescencephotonlikelihoods,
                                                  int             *data,
                                                  double           interval,
                                                  double           modperiod,
                                                  BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesx, *firstderivativesy, *secondderivativesxy;

    firstderivativesx   = Bayes_dvector(0,nbins-1);
    firstderivativesy   = Bayes_dvector(0,nbins-1);
    secondderivativesxy = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtWeight(     firstderivativesx,decayindex_x,
                                                     ndecays,weights,taus,
                                                     fluorescencephotonlikelihoods,
                                                     nbins,binwalls,data,
                                                     interval,modperiod,instr);

    bayes_DataLikelihood1stDerivativesWrtWeight(     firstderivativesy,decayindex_y,
                                                     ndecays,weights,taus,
                                                     fluorescencephotonlikelihoods,
                                                     nbins,binwalls,data,
                                                     interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtDiffWeights(secondderivativesxy,decayindex_x,decayindex_y,
                                                     ndecays,weights,taus,
                                                     fluorescencephotonlikelihoods,
                                                     nbins,binwalls,data,
                                                     interval,modperiod,instr);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesxy[bin]-firstderivativesx[bin]*firstderivativesy[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(firstderivativesy,0,nbins-1);
    free_Bayes_dvector(secondderivativesxy,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementWrtSingleLifetime(int              decayindex_x,
                                                    int              decayindex_y,
                                                    int              ndecays,
                                                    double          *weights,
                                                    double          *taus,
                                                 //   double          *hyperparams,
                                                    int              nbins,
                                                    double          *binwalls,
                                                    double          *binlikelihoods,
                                                    double         **fluorescencephotonlikelihoods,
                                                    int             *data,
                                                    double           interval,
                                                    double           modperiod,
                                                    BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesx, *secondderivativesxx;

    firstderivativesx   = Bayes_dvector(0,nbins-1);
    secondderivativesxx = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativesx,/*derivativeindex_x,*/
                                                  ndecays,weights[decayindex_x],taus[decayindex_x],
                                                  /*fluorescencephotonlikelihoods,*/
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtSingleLifetime(secondderivativesxx,/*derivativeindex_x,*/
                                                        ndecays,weights[decayindex_x],taus[decayindex_x],
                                                        //double         **photonlikelihoods,
                                                        nbins,binwalls,data,
                                                        interval,modperiod,instr);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesxx[bin]-firstderivativesx[bin]*firstderivativesx[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(secondderivativesxx,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementWrtMixedLifetimes(int              decayindex_x,
                                                    int              decayindex_y,
                                                    int              ndecays,
                                                    double          *weights,
                                                    double          *taus,
                                                   // double          *hyperparams,
                                                    int              nbins,
                                                    double          *binwalls,
                                                    double          *binlikelihoods,
                                                    double         **fluorescencephotonlikelihoods,
                                                    int             *data,
                                                    double           interval,
                                                    double           modperiod,
                                                    BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesx, *firstderivativesy, *secondderivativesxy;

    firstderivativesx   = Bayes_dvector(0,nbins-1);
    firstderivativesy   = Bayes_dvector(0,nbins-1);
    secondderivativesxy = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativesx,
                                                  ndecays,
                                                  weights[decayindex_x],
                                                  taus[decayindex_x],
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativesy,
                                                  ndecays,
                                                  weights[decayindex_y],
                                                  taus[decayindex_y],
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtDiffLifetimes(secondderivativesxy,nbins);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesxy[bin]-firstderivativesx[bin]*firstderivativesy[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(firstderivativesy,0,nbins-1);
    free_Bayes_dvector(secondderivativesxy,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementWrtWeightAndLifetimeSameSignalComponent(int              decayindex_x,
                                                                          int              decayindex_y,
                                                                          int              ndecays,
                                                                          double          *weights,
                                                                          double          *taus,
                                                                        //  double          *hyperparams,
                                                                          int              nbins,
                                                                          double          *binwalls,
                                                                          double          *binlikelihoods,
                                                                          double         **fluorescencephotonlikelihoods,
                                                                          int             *data,
                                                                          double           interval,
                                                                          double           modperiod,
                                                                          BayesInstrRsp_t *instr)
{
    int    bin, cj;
    double temp, val, *firstderivativesx, *firstderivativesy, *secondderivativesxy;

    firstderivativesx   = Bayes_dvector(0,nbins-1);
    firstderivativesy   = Bayes_dvector(0,nbins-1);
    secondderivativesxy = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtWeight(firstderivativesx,decayindex_x,
                                                ndecays,weights,taus,
                                                fluorescencephotonlikelihoods,
                                                nbins,binwalls,data,
                                                interval,modperiod,instr);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativesy,
                                                  ndecays,weights[decayindex_x],taus[decayindex_x],
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);

    bayes_DataLikelihood2ndDerivativesWrtWeightAndLifetimeSameSignalComponents(secondderivativesxy,
                                                                               ndecays,weights[decayindex_x],taus[decayindex_x],
                                                                               nbins,binwalls,data,
                                                                               interval,modperiod,instr);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesxy[bin]-firstderivativesx[bin]*firstderivativesy[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesx,0,nbins-1);
    free_Bayes_dvector(firstderivativesy,0,nbins-1);
    free_Bayes_dvector(secondderivativesxy,0,nbins-1);

    return (val);
}


double bayes_ComputeHessianElementWrtWeightAndLifetimeDiffSignalComponents(int              weightderivativeindex,
                                                                           int              tauderivativeindex,
                                                                           int              ndecays,
                                                                           double          *weights,
                                                                           double          *taus,
                                                                         //  double          *hyperparams,
                                                                           int              nbins,
                                                                           double          *binwalls,
                                                                           double          *binlikelihoods,
                                                                           double         **fluorescencephotonlikelihoods,
                                                                           int             *data,
                                                                           double           interval,
                                                                           double           modperiod,
                                                                           BayesInstrRsp_t *instr)
{
    /* NOTE: Indices are switched prior to calling this function if required...!!! */
    int    bin, cj;
    double temp, val, *firstderivativesweight, *firstderivativestau, *secondderivativesweightandtau;

    firstderivativesweight        = Bayes_dvector(0,nbins-1);
    firstderivativestau           = Bayes_dvector(0,nbins-1);
    secondderivativesweightandtau = Bayes_dvector(0,nbins-1);

    bayes_DataLikelihood1stDerivativesWrtWeight(firstderivativesweight,weightderivativeindex,
                                                ndecays,weights,taus,
                                                fluorescencephotonlikelihoods,
                                                nbins,binwalls,data,
                                                interval,modperiod,instr);

    bayes_DataLikelihood1stDerivativesWrtLifetime(firstderivativestau,
                                                  ndecays,weights[tauderivativeindex],taus[tauderivativeindex],
                                                  nbins,binwalls,data,
                                                  interval,modperiod,instr);



    bayes_DataLikelihood2ndDerivativesWrtWeightAndLifetimeDiffSignalComponents(secondderivativesweightandtau,
                                                                               weights[weightderivativeindex],
                                                                               firstderivativestau,
                                                                               nbins);

    for (bin=0,val=0.0; bin<nbins; bin++)
    {
        cj = data[bin];

        if (cj)
        {
            temp  = binlikelihoods[bin]*secondderivativesweightandtau[bin]-firstderivativesweight[bin]*firstderivativestau[bin];
            temp /= binlikelihoods[bin]*binlikelihoods[bin];
            val  -= (double)cj*temp;
        }
    }

    free_Bayes_dvector(firstderivativesweight,0,nbins-1);
    free_Bayes_dvector(firstderivativestau,0,nbins-1);
    free_Bayes_dvector(secondderivativesweightandtau,0,nbins-1);

    return (val);
}




typedef double(*ptrHessianElementCalculator)(int,int,int,double*,/*double*,*/double*,int,double*,double*,
                                             double**,int*,double,double,BayesInstrRsp_t*);


ptrHessianElementCalculator bayes_DetermineHessianElementCalculatorFromIndices(int i, int j, int k, int *x, int *y)
{
    if ((i<=0) || (j<=0) || (i>2*k) || (j>2*k))
    {
        return (NULL);
    }

    if ((i<=k) && (j<=k)) /* Weights only quadrant... */
    {
        *x = i; *y = j;
        
        if (i==j) /* Diagonal element... */
        {
            return (&bayes_ComputeHessianElementWrtSingleWeight);
        }
        else
        {
            return (&bayes_ComputeHessianElementWrtMixedWeights);
        }
    }
    else if ((i>k) && (j>k)) /* Lifetimes only quadrant... */
    {
        *x = i-k; *y = j-k;

        if (i==j) /* Diagonal element... */
        {
            return (&bayes_ComputeHessianElementWrtSingleLifetime);
        }
        else
        {
            return (&bayes_ComputeHessianElementWrtMixedLifetimes);
        }    
    }
    else /* Mixed weight and lifetime quadrant... */
    {
        if (i<=k)
        {
            *x = i; *y = j-k;
        }
        else
        {
            *x = j; *y = i-k; /* Transpose... */
        }

        if (*x==*y) /* Weight and lifetime variables of same signal component... */
        {
            return (&bayes_ComputeHessianElementWrtWeightAndLifetimeSameSignalComponent);
        }
        else
        {
            return (&bayes_ComputeHessianElementWrtWeightAndLifetimeDiffSignalComponents);
        }    
    }
}




int bayes_PopulateHessianMatrix(double         **hessian,
                                int              ndecays,
                                double          *weights,
                                double          *taus,
                                double          *hyperparams,
                                int              nbins,
                                double          *binwalls,
                                int             *data,
                                double           interval,
                                double           modperiod,
                                BayesInstrRsp_t *instr)
{
    int    nparams, i, j, x, y;
    double value, *binlikelihoods, **fluorescencephotonlikelihoods;
    double (*calculator)(int,int,int,double*,/*double*,*/double*,int,double*,double*,double**,int*,double,double,BayesInstrRsp_t*);

    nparams                       = 2*ndecays;
    fluorescencephotonlikelihoods = Bayes_dmatrix(1,ndecays,0,nbins-1);
    binlikelihoods                = Bayes_dvector(0,nbins-1);

    for (i=1; i<=ndecays; i++)
        bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(fluorescencephotonlikelihoods[i],nbins,binwalls,data,interval,modperiod,instr,taus[i],ndecays,weights,taus);

    bayes_ComputeBinLikelihoodsFromWeightsAndFluorescencePhotonLikelihoods(binlikelihoods,nbins,binwalls,
                                                                           ndecays,fluorescencephotonlikelihoods,weights,
                                                                           instr,interval);
    for (i=1; i<=nparams; i++)
    {
        for (j=i; j<=nparams; j++)
        {
            /* Select the appropriate element calculator... */
            calculator    = bayes_DetermineHessianElementCalculatorFromIndices(i,j,ndecays,&x,&y);

            value         = calculator(x,y,
                                       ndecays,weights,taus,
                                       nbins,binwalls,binlikelihoods,fluorescencephotonlikelihoods,data,
                                       interval,modperiod,instr);
            hessian[i][j] = value;
        }
    }

    for (i=1; i<=nparams; i++)
    {
        for (j=1; j<i; j++)
        {
            hessian[i][j] = hessian[j][i];
        }
    }

    free_Bayes_dvector(binlikelihoods,0,nbins-1);
    free_Bayes_dmatrix(fluorescencephotonlikelihoods,1,ndecays,0,nbins-1);

    return (0);
}


// akaike information criteria simplified analysis needs to be implemented for rapid model selection
int bayes_DetemineDecayModelEvidence(int              ndecays,
                                     double          *weights,
                                     double          *taus,
                                     double          *hyperparams,
                                     double           minuslogprob,
                                     int              nbins,
                                     double          *binwalls,
                                     int             *data,
                                     double           interval,
                                     double           modperiod,
                                     BayesInstrRsp_t *instr,
                                     double          *logmodelevidence)
{
    double **hessian, determinant, logevidence;
    int      nparams;
	
	FILE *fp;
    int   Debug=0, i, j;

    nparams = 2*ndecays;
    hessian = Bayes_dmatrix(1,nparams,1,nparams);

    bayes_PopulateHessianMatrix(hessian,ndecays,weights,taus,hyperparams,nbins,binwalls,data,interval,modperiod,instr);

    if (Debug)
    {
        if (nparams==2)
        {
            if ((fp = fopen("hessian_mono.txt","w")))
                Debug=0;
        }
        else if (nparams==4)
        {
            if ((fp = fopen("hessian_bi.txt","w")))
                Debug=0;
        }
        else
        {
            if ((fp = fopen("hessian.txt","w")))
                Debug=0;        
        }

        if (Debug)
        {
            fprintf(fp, "Most probable values for Hessian computation\n");

            for (i=0; i<=ndecays; i++)
            {
                fprintf(fp, "weights[%d]: %g\n", i, weights[i]);
            }
            
            for (i=1; i<=ndecays; i++)
            {
                fprintf(fp, "taus[%d]: %g\n", i, taus[i]);
            }
            
            fprintf(fp, "minuslogprob: %g\n\n", minuslogprob);

            for (i=1; i<=nparams; i++)
            {
                for (j=1; j<=nparams; j++)
                {
                    fprintf(fp, "%g\t", hessian[i][j]);
                }
                fprintf(fp, "\n");
            }
        }
    }

    determinant = bayes_ComputeDeterminantValue(hessian,nparams);

    if (determinant<=0.0)
    {
        if (Debug)
        {
            fprintf(fp, "\ndeterminant: %g\n", determinant);
            fprintf(fp, "ERROR!");
            fclose(fp);
        }

        *logmodelevidence = 0.0;        

        free_Bayes_dmatrix(hessian,1,nparams,1,nparams);

        return (-1);
    }

    logevidence  = -0.5*log(determinant);
    logevidence += nparams*log(2.0*PI);
    logevidence -= minuslogprob;

    *logmodelevidence = logevidence;

    if (Debug)
    {
        fprintf(fp, "\ndeterminant: %g\n", determinant);
        fprintf(fp, "logevidence: %g\n", logevidence);
        fclose(fp);
    }

    free_Bayes_dmatrix(hessian,1,nparams,1,nparams);

    return (0);
}

int const bayes_NdecaysToModelType[] = {-1,FIT_MONOEXP,FIT_BIEXP,FIT_MULTIEXP};

int bayes_DetemineDecayModelRelativeLikelihoods(/* Data in... */
                                                int                       nbins,
												int                       fitstart,
                                                double                   *binwalls,
                                                int                      *data,
                                                int                       nphotons,
                                                /* Instrumentation... */
                                                double                    interval,
                                                double                    modperiod,
                                                BayesInstrRsp_t          *instr,
                                                /* RLD-derived estimates for search initialisation... */
                                                BayesParamValsAndFit_t   *decayestimates,
                                                /* Model selection data out... */
                                                float                    *decaymodellikelihoods,
                                                BayesParamValsAndFit_t   *paramvalsandfits,
                                                /* Configuration... */
                                                int                       rapidanalysis,
                                                BayesRapidValueStore_t   *rapidgrid)
{
    int     ret=BAYES_ERR_NO_ERROR, i, k, K;
    int     nparams;
    double *weights_mp, *weights_ave, *weights_err, *taus_mp, *taus_ave, *taus_err;
    double *hyperparams, val, alpha, evidence;
    float   valfloat;
    int     paramsfree[] = {1,1,1,1,1,1,1}; //no parameter fixing for model selection
    int     modeltype;
    int     akaike=1;
    BayesUserFixedParams_t paramfixing;

    alpha = 1.0/(double)data_ComputeBinnedDataAverageArrTime(data,nbins,fitstart,nphotons,(float)interval);
    
    for (K=1; K<=2; K++)
    {
        modeltype   = bayes_NdecaysToModelType[K];
        nparams     = 2*K;

        if (paramvalsandfits)
        {
            weights_mp  = paramvalsandfits[K].weights;
            taus_mp     = paramvalsandfits[K].taus;
        }
        else
        {
            weights_mp  = Bayes_dvector(0,K);
            taus_mp     = Bayes_dvector(1,K);        
        }

        /* Load the pre-estimates as starting point in decay parameter value search... */
        for (k=0; k<=K; k++)
            weights_mp[k] = decayestimates[K].weights[k];

        for (k=1; k<=K; k++)
            taus_mp[k] = decayestimates[K].taus[k];

        weights_ave = NULL;
        weights_err = NULL;
        taus_ave    = NULL;
        taus_err    = NULL;
        hyperparams = Bayes_dvector(1,K);

        bayes_CheckParameterValueFixingForBayesFitting(&paramfixing,1+nparams,paramsfree,NULL,
                                                           nbins,fitstart,binwalls,nphotons,interval,modperiod,instr);

        ret = bayes_PerformBayesParameterEstimation(data,nbins,fitstart,binwalls,nphotons,
                                                        instr,(float)modperiod,(float)interval,
                                                        modeltype,K,&paramfixing,alpha,
                                                        weights_mp,taus_mp,weights_ave,taus_ave,weights_err,taus_err,
                                                        &valfloat,NULL,rapidanalysis,rapidgrid);

        bayes_FreeParameterValueFixingForBayesFitting(&paramfixing,K,nbins);

        if (ret>=0) /* Successful parameter estimation required before Gaussian approximated determination of model likelihood... */
        {
            if (akaike) /* Akaike information criterion... */
            {
                /* AICc[K-1]... */
                decaymodellikelihoods[K-1] = (float)(
                    (2.0*(valfloat+(double)nparams)) /* AIC */
                    +
                    ((2.0*(double)nparams*((double)nparams+1.0))/((double)nphotons-(double)nparams-1.0))); /* Correction term */
            } 
            else /* Gaussian approximation... */
            {
				val = valfloat;

				ret = bayes_DetemineDecayModelEvidence(K,weights_mp,taus_mp,hyperparams,val,
													   nbins,binwalls,data,interval,modperiod,instr,&evidence);
				if (ret>=0)
				{
					decaymodellikelihoods[K-1] = (float)evidence;
				}
				else
				{
					ret = BAYES_ERR_MODEL_SEL_HESSIAN_FAILURE;

					if (!paramvalsandfits)
					{
						free_Bayes_dvector(weights_mp,0,K);
						free_Bayes_dvector(taus_mp,1,K);
					}
	                
					free_Bayes_dvector(hyperparams,1,K);

					break;
                }        
            }        
        }
        else
        {
            ret = BAYES_ERR_MODEL_SEL_PARAM_EST_FAILURE;

            if (!paramvalsandfits)
            {
                free_Bayes_dvector(weights_mp,0,K);
                free_Bayes_dvector(taus_mp,1,K);            
            }

            free_Bayes_dvector(hyperparams,1,K);

            break;
        }

        if (!paramvalsandfits)
        {
            free_Bayes_dvector(weights_mp,0,K);
            free_Bayes_dvector(taus_mp,1,K);            
        }

        free_Bayes_dvector(hyperparams,1,K); 
    }

    if (ret>=BAYES_ERR_NO_ERROR)
    {
        if (akaike)
        {
            /* http://en.wikipedia.org/wiki/Akaike_information_criterion         */
            /* Given a set of candidate models for the data, the preferred model */
            /* is the one with the minimum AIC value. Hence AIC not only rewards */
            /* goodness of fit, but also includes a penalty that is an           */
            /* increasing function of the number of estimated parameters. This   */
            /* penalty discourages overfitting.                                  */

            if (decaymodellikelihoods[0]<=decaymodellikelihoods[1]) //mono-exp
            {
                decaymodellikelihoods[3] = 1.0;

                /* Report the relative likelihood of bi-exponential (K=2)...                     */
                /* Denote the AIC values of the candidate models by AIC1, AIC2, AIC3, ..., AICR. */
                /* Let AICmin be the minimum of those values.                                    */
                /* The exp((AICmin-AICi)/2) can be interpreted as the relative probability that  */
                /* the ith model minimizes the (estimated) information loss.                     */
				decaymodellikelihoods[2] = 1.0f/(1.0f+(float)exp((decaymodellikelihoods[0]-decaymodellikelihoods[1])/2.0f));
				decaymodellikelihoods[2] = 1.0f-decaymodellikelihoods[2];
            }
            else //bi-exp
            {
                decaymodellikelihoods[3] = 2.0f;
				decaymodellikelihoods[2] = 1.0f/(1.0f+(float)exp((decaymodellikelihoods[1]-decaymodellikelihoods[0])/2.0f));
			}
        }
        else
        {
			// compute the probability
			//  let d0 = decaymodellikelihoods[0] = ln(p0); sim. for 1
			//	double p0 = exp(d0);
			//	double p1 = exp(d1);
			//  want P0 = p0/(p0+p1)
			//  in logs
			//  ln(P0)= ln(p0) - ln(p0+p1)
			//        = d0     - (ln(p0) + ln{1 + exp[ln(p1)-ln(p0)]})      (http://en.wikipedia.org/wiki/List_of_logarithmic_identities)
			//        = d0     - d0     - ln{1 + exp[d1-d0]}
			//        = -ln{1 + exp[d1 - d0]}
			// So  P0 = 1 / (1 + exp[d1 - d0])
			// and P1 = 1 / (1 + exp[d0 - d1])

		//	decaymodellikelihoods[2] = 1.0f / (1.0f + (float)exp(decaymodellikelihoods[1] - decaymodellikelihoods[0]));  // less complex
			decaymodellikelihoods[2] = 1.0f / (1.0f + (float)exp(decaymodellikelihoods[0] - decaymodellikelihoods[1]));  // more complex

			/* Report the order of the most probable model, i.e. 1=MONO, 2=BI, etc... */
			if (decaymodellikelihoods[0]>decaymodellikelihoods[1])
				decaymodellikelihoods[3] = 1.0;
			else
				decaymodellikelihoods[3] = 2.0;
        }
    }
    else
    {
        for (i=0; i<4; i++)
        {
            decaymodellikelihoods[i] = 0.0;
        }
    }

    return (ret);
}


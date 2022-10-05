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
#include "stdio.h"
#include "matrices.h"
#include "extmath.h"
#include "DTYPE.h"

#include "bayes_Sizes.h"
#include "bayes_DistributionFctsBinLikelihoods.h"
#include "bayes_DataManagement.h"
#include "bayes_Interface.h"
#include "bayes_InstrRspAnalysis.h"



/***********************************************************************************/
/*                                                                                 */
/*              ARRIVAL TIME BIN-LIKELIHOOD CONSTANT COMPUTATION ROUTINES          */
/*                                                                                 */
/***********************************************************************************/


int  bayes_PopulateBinWallsVectorUniformIntervals(double *binwalls, /* Vector running from '0' to 'nbins'... */
                                                  int     nbins,
                                                  double  interval)
{
    int    j;
    double binwidth, t;

    if ((!binwalls) || (nbins<=0) || (interval<=0.0))
        return (-1);

    binwidth = interval/(double)nbins;

    for (j=0, t=0.0; j<=nbins; j++)
    {
        binwalls[j]  = t;
        t           += binwidth;
    }

    return (0);
}


double bayes_LogMultiply(double a, double b)
{
    int sign=1;

    if (a<0.0)
    {
        a     = -a;
        sign *= -1;
    }
    
    if (b<0.0)
    {
        b     = -b;
        sign *= -1;
    }
    
    return ((double)sign*exp(log(a)+log(b)));
}


//may need to introduce error detection here, e.g. for very small tau causing overflow etc...
// 130125 - stability - atomised difference computations
int bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(double          *normalisation,
                                                              double           interval,
                                                              double           modperiod,
                                                              double           dithertime,
                                                              BayesInstrRsp_t *instr,
                                                              int              ndecays,
                                                              double          *weights,
                                                              double          *taus)
{
    int    j, type;
    double T, Tm, Ts, w0dagger, weight, tau;

    int    nbins=1;
    double binwalls[2],likelihoods[1];
    double ditherbinwalls[2],ditherlikelihoods[1];
    double norm;

    binwalls[0] = T  = interval;
    binwalls[1] = Tm = modperiod;

    if (dithertime>0.0)
    {
        ditherbinwalls[0] = 0.0;
        ditherbinwalls[1] = Ts = dithertime;    
    }

    for (j=1,w0dagger=0.0; j<=ndecays; j++)
        w0dagger += weights[j];

    w0dagger = 1.0/w0dagger;

    for (j=1,norm=0.0; j<=ndecays; j++)
    {
        weight = weights[j];
        tau    = taus[j];

        bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(
            likelihoods,nbins,binwalls,NULL,T,Tm,instr,tau,0,NULL,NULL);

        norm += (weight*likelihoods[0]);

        if (dithertime>0.0) /* Normalisation must be corrected as data is being discarded at start of transient due to dither... */
        {
            bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(
                ditherlikelihoods,nbins,ditherbinwalls,NULL,T,Tm,instr,tau,0,NULL,NULL);

            norm += (weight*ditherlikelihoods[0]);
        }
    }

    norm           *= w0dagger;
    *normalisation  = 1.0-norm;

    if ((*normalisation<=0.0) ||
        (*normalisation >1.0) ||
        (BAYES_DM_DOUBLE_TYPE_INVALID==bayes_dm_CheckDoubleValueValid(*normalisation,&type)))
    {
        return (-1);
    } //error case

    return (0);
}


// new model
#if 1 // 130125 - stability - atomised differences computation
int bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(double          *fluorescencephotonlikelihoods,
                                                               int              nbins,
                                                               double          *binwalls,
                                                               int             *data,
                                                               double           interval,
                                                               double           modperiod,
                                                               BayesInstrRsp_t *instr,
                                                               double           tau,
                                                               int              ndecays,
                                                               double          *weights,
                                                               double          *taus)
{
    int    bin, valid, i, ell, ellmax, ret;
    double bL, bH, valbin[BAYES_INSTR_RSP_MAX_COMPONENTS];
    double gammaprime[BAYES_INSTR_RSP_MAX_COMPONENTS], delay, width, cutoff, weight;

    double t_ell_bL, t_ell_bH;
    double arg1_L, arg1_H, erf_arg1_L, erf_arg1_H;
    double arg2, erf_arg2;
    double arg3_L, arg3_H, erf_arg3_L, erf_arg3_H;
    double arg4_L, arg4_H, exp_arg4_L, exp_arg4_H;
    double arg5, exp_arg5;
    int    erf_arg3_L_minus_erf_arg2_non_zero, erf_arg3_H_minus_erf_arg2_non_zero;
    double norm, temp;

    int    DebugTrace, DebugDetail;
    FILE   *fp = NULL;

    DebugTrace  = 0;
    DebugDetail = DebugTrace; //this could be added to ui

    ellmax = bayes_ConfigDecayModelMaxNumberOfModulationPeriodsConsidered();

    /* Compute normalisation if all lifetimes and weights are provided... */
    if ((ndecays>0) && (weights) && (taus))
    {
        ret = bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(
                  &norm,interval,modperiod,0.0,instr,ndecays,weights,taus);

        if (ret<0)
        {
            /*if (DebugTrace)
            {
                if (0 == fopen(&fp, "bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTauDebugOutput.txt", "w"))
                {
                    fprintf(fp,"\nbayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau ==>\n");
                    fprintf(fp,"Error (%d) in normalisation\n",ret);
                    fclose(fp);
                }            
            }*/

            return (ret);
        }
    }
    else
    {
        norm = 1.0; /* Default value if all decay signal parameter values are not supplied... */
    }

    for (i=0; i<instr->ninstr; i++)
    {
        weight = instr->params[i].weight;
        delay  = instr->params[i].delay;
        width  = instr->params[i].width;
        cutoff = instr->params[i].cutoff;

        gammaprime[i] = weight/(1.0+erf((delay-cutoff)/(width*ROOTTWO)));
    }

    /*if (DebugTrace)
    {
        if fp = fopen(("bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTauDebugOutput.txt", "w"))))
        {
            DebugTrace = 0;
        }
        else
        {
            fprintf(fp,"\nbayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau ==>\n");

            fprintf(fp,"\n=============  Settings  =============\n");

            fprintf(fp,"nbins ==> %d\n", nbins);
            fprintf(fp,"interval ==> %g\n", interval);
            fprintf(fp,"modperiod ==> %g\n", modperiod);
            fprintf(fp,"tau ==> %g\n", tau);
            fprintf(fp,"ellmax ==> %d\n", ellmax);

            fprintf(fp,"\n=============  Instrument  =============\n");
            fprintf(fp,"ninstr ==> %d\n", instr->ninstr);

            for (i=0; i<instr->ninstr; i++)
            {
                fprintf(fp,"%d) weight: %g\n", i,instr->params[i].weight);
                fprintf(fp,"%d) delay: %g\n", i,instr->params[i].delay);
                fprintf(fp,"%d) width: %g\n", i,instr->params[i].width);
                fprintf(fp,"%d) cutoff: %g\n", i,instr->params[i].cutoff);
                fprintf(fp,"%d) gammaprime: %g\n", i,gammaprime[i]);
            }

            fprintf(fp,"\n=============  Normalisation  =============\n");
            fprintf(fp,"ndecays ==> %d\n", ndecays);

            if (ndecays>0)
            {
                for (i=1; i<=ndecays; i++)
                {
                    fprintf(fp,"%d) weight: %g\n", i, weights[i]);
                }

                for (i=1; i<=ndecays; i++)
                {
                    fprintf(fp,"%d) tau: %g\n", i, taus[i]);
                }

                fprintf(fp,"norm ==> %g\n", norm);
            }
            else
            {
                fprintf(fp,"norm ==> %g\n", norm);
            }

            if (DebugDetail)
            {
                fprintf(fp,"\n=============  Computation  =============\n");
                fprintf(fp,"arg1 = (delay-t_ell)/(width*ROOTTWO)\n");
                fprintf(fp,"arg2 = ((delay-cutoff)*tau+width*width)/(width*tau*ROOTTWO)\n");
                fprintf(fp,"arg3 = (((delay-t_ell)*tau)+width*width)/(width*tau*ROOTTWO)\n");
                fprintf(fp,"arg4 = -(t_ell-delay)/tau\n");
                fprintf(fp,"arg5 = (0.5*width*width)/(tau*tau)\n");            
            }

            fprintf(fp,"\n========================================\n\n\n");

            fprintf(fp,"bin\t");
            fprintf(fp,"bL\t");
            fprintf(fp,"bH\t");

            if (DebugDetail)
            {
                for (i=0; i<instr->ninstr; i++)
                {
                    for (ell=0; ell<ellmax; ell++)
                    {
                        fprintf(fp,"arg1_L(i=%d,ell=%d)\t", i, ell);
                        fprintf(fp,"arg1_H(%d,%d)\t", i, ell);
                        fprintf(fp,"erf_arg1_L(%d,%d)\t", i, ell);
                        fprintf(fp,"erf_arg1_H(%d,%d)\t", i, ell);
                        fprintf(fp,"arg2(%d,%d)\t", i, ell);
                        fprintf(fp,"erf_arg2(%d,%d)\t", i, ell);
                        fprintf(fp,"arg3_L(%d,%d)\t", i, ell);
                        fprintf(fp,"arg3_H(%d,%d)\t", i, ell);
                        fprintf(fp,"erf_arg3_L(%d,%d)\t", i, ell);
                        fprintf(fp,"erf_arg3_H(%d,%d)\t", i, ell);
                        fprintf(fp,"arg4_L(%d,%d)\t", i, ell);
                        fprintf(fp,"arg4_H(%d,%d)\t", i, ell);
                        fprintf(fp,"exp_arg4_L(%d,%d)\t", i, ell);
                        fprintf(fp,"exp_arg4_H(%d,%d)\t", i, ell);
                        fprintf(fp,"arg5(%d,%d)\t", i, ell);
                        fprintf(fp,"exp_arg5(%d,%d)\t", i, ell);
                        fprintf(fp,"likelihood(%d,%d)\t", i, ell);
                    }            
                }            
            }

            fprintf(fp,"likelihood (total)\n\n");
        }
    }*/

    for (bin=0,valid=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            bL = binwalls[bin];
            bH = binwalls[bin+1];

            /*if (DebugTrace)
            {
                valid =0;
                fprintf(fp,"%d\t%g\t%g\t",bin,bL,bH);
            }*/
#if 0
            if (valid)
            {
                // there is scope to makes things faster by using previously computed values,
                // but need to do this for each argument for each instrument component and for each ell...
            }
            else
#endif 
            {
                for (i=0; i<instr->ninstr; i++)
                {
                    valbin[i] = 0.0;
                    delay     = instr->params[i].delay;
                    width     = instr->params[i].width;
                    cutoff    = instr->params[i].cutoff;

                    arg2      = ((delay-cutoff)*tau+width*width)/(width*tau*ROOTTWO);
                    erf_arg2  = erf(arg2);

                    for (ell=0; ell<ellmax; ell++)
                    {
                        t_ell_bL = modperiod*(double)ell + bL;
                        t_ell_bH = modperiod*(double)ell + bH;

                        if (t_ell_bH>cutoff)
                        {
                            if (t_ell_bL<cutoff)
                                t_ell_bL = cutoff;

                            arg1_L     = (delay-t_ell_bL)/(width*ROOTTWO);
                            arg1_H     = (delay-t_ell_bH)/(width*ROOTTWO);
                            erf_arg1_L = erf(arg1_L);
                            erf_arg1_H = erf(arg1_H);

                            arg3_L     = (((delay-t_ell_bL)*tau)+width*width)/(width*tau*ROOTTWO);
                            arg3_H     = (((delay-t_ell_bH)*tau)+width*width)/(width*tau*ROOTTWO);
                            erf_arg3_L = erf(arg3_L);
                            erf_arg3_H = erf(arg3_H);

                            if ((erf_arg3_L-erf_arg2)!=0.0)
                            {
                                erf_arg3_L_minus_erf_arg2_non_zero = 1;

                                arg4_L     = -(t_ell_bL-delay)/tau;
                                exp_arg4_L = exp(arg4_L);
                            }
                            else
                            {
                                erf_arg3_L_minus_erf_arg2_non_zero = 0;
                            }

                            if ((erf_arg3_H-erf_arg2)!=0.0)
                            {
                                erf_arg3_H_minus_erf_arg2_non_zero = 1;

                                arg4_H     = -(t_ell_bH-delay)/tau;
                                exp_arg4_H = exp(arg4_H);

                            }
                            else
                            {
                                erf_arg3_H_minus_erf_arg2_non_zero = 0;
                            }

                            arg5       = (0.5*width*width)/(tau*tau);
                            exp_arg5   = exp(arg5);

                            temp       = 0.0;

                            if (erf_arg3_H_minus_erf_arg2_non_zero)
                                temp  += -erf_arg2*exp_arg4_H;
                            
                            if (erf_arg3_L_minus_erf_arg2_non_zero)
                                temp  -= -erf_arg2*exp_arg4_L;

                            if (erf_arg3_H_minus_erf_arg2_non_zero)
                                temp  += erf_arg3_H*exp_arg4_H;
                            
                            if (erf_arg3_L_minus_erf_arg2_non_zero)
                                temp  -= erf_arg3_L*exp_arg4_L;
                            
                            temp      *= exp_arg5;
                            temp      += -erf_arg1_H;
                            temp      -= -erf_arg1_L;

                            valbin[i] += temp;

                            if ((DebugTrace) && (DebugDetail))
                            {
                                fprintf(fp,"%g\t%g\t%g\t%g\t",arg1_L,
                                                              arg1_H,
                                                              erf_arg1_L,
                                                              erf_arg1_H);

                                fprintf(fp,"%g\t%g\t",        arg2,
                                                              erf_arg2);

                                fprintf(fp,"%g\t%g\t%g\t%g\t",arg3_L,
                                                              arg3_H,
                                                              erf_arg3_L,
                                                              erf_arg3_H);

                                if (erf_arg3_L_minus_erf_arg2_non_zero)
                                    fprintf(fp,"%g\t",arg4_L);
                                else
                                    fprintf(fp,"\t");

                                if (erf_arg3_H_minus_erf_arg2_non_zero)
                                    fprintf(fp,"%g\t",arg4_H);
                                else
                                    fprintf(fp,"\t");

                                if (erf_arg3_L_minus_erf_arg2_non_zero)
                                    fprintf(fp,"%g\t",exp_arg4_L);
                                else
                                    fprintf(fp,"\t");

                                if (erf_arg3_H_minus_erf_arg2_non_zero)
                                    fprintf(fp,"%g\t",exp_arg4_H);
                                else
                                    fprintf(fp,"\t");

                                fprintf(fp,"%g\t%g\t",        arg5,
                                                              exp_arg5);

                                fprintf(fp,"%g\t",            temp);
                            }                                                
                        }
                        else
                        {
                            if ((DebugTrace) && (DebugDetail))
                            {
                                fprintf(fp,"\t\t\t\t");
                                fprintf(fp,"\t\t");
                                fprintf(fp,"\t\t\t\t");
                                fprintf(fp,"\t\t\t\t");
                                fprintf(fp,"\t\t");
                                fprintf(fp,"\t");
                            }                                                                        
                        }
                    }
                }

                valid = 1;
            }

            /* Now do the sum over all instrument response components... */
            for (i=0,fluorescencephotonlikelihoods[bin]=0.0; i<instr->ninstr; i++)
            {
                fluorescencephotonlikelihoods[bin] += (gammaprime[i]*valbin[i])/norm;
            }

//            if ((valbin<0.0) || (BAYES_DM_DOUBLE_TYPE_INVALID==bayes_dm_CheckDoubleValueValid(valbin,&type)))
//                valbin = valbin;

            /*if (DebugTrace)
            {
                fprintf(fp,"%g\n",fluorescencephotonlikelihoods[bin]);
            }*/
        }
        else
        {
            fluorescencephotonlikelihoods[bin] = 0.0;
            valid = 0;
        }
    }

    /*if (DebugTrace)
    {
        fclose(fp);
    }*/
    
    return (0);
}
#else
int bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(double          *fluorescencephotonlikelihoods,
                                                               int              nbins,
                                                               double          *binwalls,
                                                               int             *data,
                                                               double           interval,
                                                               double           modperiod,
                                                               BayesInstrRsp_t *instr,
                                                               double           tau,
                                                               int              ndecays,
                                                               double          *weights,
                                                               double          *taus)
{
    int     bin, valid, i, ell, ellmax, ret, type;
    double  bL, bH, *chiL, *chiH, valbin;
    double  upsilon1, upsilon2, upsilon3, norm;
    double  gammaprime[BAYES_INSTR_RSP_MAX_COMPONENTS], delay, width, cutoff, weight;

    int     DebugTrace;
    FILE    *fp = NULL;

    DebugTrace = 0;

    ellmax = bayes_ConfigDecayModelMaxNumberOfModulationPeriodsConsidered();

    if ((ndecays>0) && (weights) && (taus))
    {
        ret = bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(
                  &norm,interval,modperiod,0.0,instr,ndecays,weights,taus);

        if (ret<0)
        {
            if (DebugTrace)
            {
                if (0 == fopen(&fp, "bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTauDebugOutput.txt", "w"))
                {
                    fprintf(fp,"\nbayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau ==>\n");
                    fprintf(fp,"Error (%d) in normalisation\n",ret);
                    fclose(fp);
                }            
            }

            return (ret);
        }
    }
    else
    {
        norm = 1.0; /* Default value if all decay signal parameter values are not supplied... */
    }

    for (i=0; i<instr->ninstr; i++)
    {
        weight = instr->params[i].weight;
        delay  = instr->params[i].delay;
        width  = instr->params[i].width;
        cutoff = instr->params[i].cutoff;

        gammaprime[i] = weight/(1.0+erf((delay-cutoff)/(width*ROOTTWO)));
    }

    if (DebugTrace)
    {
        if fp = fopen(("bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTauDebugOutput.txt", "w"))))
        {
            DebugTrace = 0;
        }
        else
        {
            fprintf(fp,"\nbayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau ==>\n");

            fprintf(fp,"\n=============  Settings  =============\n");

            fprintf(fp,"nbins ==> %d\n", nbins);
            fprintf(fp,"interval ==> %g\n", interval);
            fprintf(fp,"modperiod ==> %g\n", modperiod);
            fprintf(fp,"tau ==> %g\n", tau);
            fprintf(fp,"ellmax ==> %d\n", ellmax);

            fprintf(fp,"\n=============  Instrument  =============\n");
            fprintf(fp,"ninstr ==> %d\n", instr->ninstr);

            for (i=0; i<instr->ninstr; i++)
            {
                fprintf(fp,"%d) weight: %g\n", i,instr->params[i].weight);
                fprintf(fp,"%d) delay: %g\n", i,instr->params[i].delay);
                fprintf(fp,"%d) width: %g\n", i,instr->params[i].width);
                fprintf(fp,"%d) cutoff: %g\n", i,instr->params[i].cutoff);
                fprintf(fp,"%d) gammaprime: %g\n", i,gammaprime[i]);
            }

            fprintf(fp,"\n=============  Normalisation  =============\n");
            fprintf(fp,"ndecays ==> %d\n", ndecays);

            if (ndecays>0)
            {
                for (i=1; i<=ndecays; i++)
                {
                    fprintf(fp,"%d) weight: %g\n", i, weights[i]);
                }

                for (i=1; i<=ndecays; i++)
                {
                    fprintf(fp,"%d) tau: %g\n", i, taus[i]);
                }

                fprintf(fp,"norm ==> %g\n", norm);
            }
            else
            {
                fprintf(fp,"norm ==> %g (default value)\n", norm);
            }


            fprintf(fp,"\n========================================\n\n\n");

            fprintf(fp,"bin\t");
            fprintf(fp,"bL\t");
            fprintf(fp,"bH\t");
            fprintf(fp,"upsilon1\t");
            fprintf(fp,"upsilon2\t");
            fprintf(fp,"upsilon3\t");
            fprintf(fp,"upsilon4\t");
            fprintf(fp,"chiL\t");
            fprintf(fp,"upsilon1\t");
            fprintf(fp,"upsilon2\t");
            fprintf(fp,"upsilon3\t");
            fprintf(fp,"upsilon4\t");
            fprintf(fp,"chiH\t");
            fprintf(fp,"likelihood\n\n");
        }
    }

    chiL = Bayes_dvector(0,instr->ninstr-1);
    chiH = Bayes_dvector(0,instr->ninstr-1);

    for (bin=0,valid=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            bL = binwalls[bin];
            bH = binwalls[bin+1];

            if (DebugTrace)
            {
                fprintf(fp,"%d\t%g\t%g\t",bin,bL,bH);
            }

            if (valid)
            {
                for (i=0; i<instr->ninstr; i++)
                    chiL[i] = chiH[i];

                /* Summation over 'ell' for each instrument response component... */
                for (i=0; i<instr->ninstr; i++)
                {
                    delay  = instr->params[i].delay;
                    width  = instr->params[i].width;
                    cutoff = instr->params[i].cutoff;

                    for (ell=0,chiH[i]=0.0; ell<ellmax; ell++)
                    {
                        if (DebugTrace)
                        {
                            fprintf(fp,"%g\t%g\t%g\t%g\t%g\t",upsilon1,
                                                              upsilon2,
                                                              upsilon3,
                                                              (-(modperiod*(double)ell+bL-delay)/tau)+((width*width)/(2.0*tau*tau)),
                                                              chiL[i]);
                        }

                        upsilon1 = erfc((delay-modperiod*(double)ell-bH)/(width*ROOTTWO));
                        upsilon2 = erfc(((delay-cutoff)*tau+width*width)/(width*tau*ROOTTWO));
                        upsilon3 = erfc(((delay-modperiod*(double)ell-bH)*tau+width*width)/(width*tau*ROOTTWO));

                        chiH[i] += upsilon1;
                        
                        if ((upsilon2-upsilon3)>0.0)
                            chiH[i] += exp((-(modperiod*(double)ell+bH-delay)/tau)+((width*width)/(2.0*tau*tau))+log(upsilon2-upsilon3));
                        else
                            chiH[i] -= exp((-(modperiod*(double)ell+bH-delay)/tau)+((width*width)/(2.0*tau*tau))+log(upsilon3-upsilon2));

                        if (DebugTrace)
                        {
                            fprintf(fp,"%g\t%g\t%g\t%g\t%g\t",upsilon1,
                                                              upsilon2,
                                                              upsilon3,
                                                              (-(modperiod*(double)ell+bH-delay)/tau)+((width*width)/(2.0*tau*tau)),
                                                              chiH[i]);
                        }
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

                    for (ell=0,chiL[i]=0.0,chiH[i]=0.0; ell<ellmax; ell++)
                    {
                        upsilon1 = erfc((delay-modperiod*(double)ell-bL)/(width*ROOTTWO));
                        upsilon2 = erfc(((delay-cutoff)*tau+width*width)/(width*tau*ROOTTWO));
                        upsilon3 = erfc(((delay-modperiod*(double)ell-bL)*tau+width*width)/(width*tau*ROOTTWO));

                        chiL[i] += upsilon1;

                        if ((upsilon2-upsilon3)>0.0)
                            chiL[i] += exp((-(modperiod*(double)ell+bL-delay)/tau)+((width*width)/(2.0*tau*tau))+log(upsilon2-upsilon3));
                        else
                            chiL[i] -= exp((-(modperiod*(double)ell+bL-delay)/tau)+((width*width)/(2.0*tau*tau))+log(upsilon3-upsilon2));

                        if (DebugTrace)
                        {
                            fprintf(fp,"%g\t%g\t%g\t%g\t%g\t",upsilon1,
                                                              upsilon2,
                                                              upsilon3,
                                                              (-(modperiod*(double)ell+bL-delay)/tau)+((width*width)/(2.0*tau*tau)),
                                                              chiL[i]);
                        }

                        upsilon1 = erfc((delay-modperiod*(double)ell-bH)/(width*ROOTTWO));
                        //upsilon2 = erfc(((delay-cutoff)*tau+width*width)/(width*tau*ROOTTWO));
                        upsilon3 = erfc(((delay-modperiod*(double)ell-bH)*tau+width*width)/(width*tau*ROOTTWO));

                        chiH[i] += upsilon1;
                        
                        if ((upsilon2-upsilon3)>0.0)
                            chiH[i] += exp((-(modperiod*(double)ell+bH-delay)/tau)+((width*width)/(2.0*tau*tau))+log(upsilon2-upsilon3));
                        else
                            chiH[i] -= exp((-(modperiod*(double)ell+bH-delay)/tau)+((width*width)/(2.0*tau*tau))+log(upsilon3-upsilon2));
                        
                        if (DebugTrace)
                        {
                            fprintf(fp,"%g\t%g\t%g\t%g\t%g\t",upsilon1,
                                                              upsilon2,
                                                              upsilon3,
                                                              (-(modperiod*(double)ell+bH-delay)/tau)+((width*width)/(2.0*tau*tau)),
                                                              chiH[i]);
                        }
                    }
                }

                valid = 1;
            }

            /* Now do the sum over all instrument response components... */
            for (i=0,valbin=0.0; i<instr->ninstr; i++)
            {
                weight = instr->params[i].weight;
                delay  = instr->params[i].delay;
                width  = instr->params[i].width;
                cutoff = instr->params[i].cutoff;

                valbin += bayes_LogMultiply(weight/(2.0-erfc((delay-cutoff)/(width*ROOTTWO))),chiH[i]-chiL[i]);
            }

            fluorescencephotonlikelihoods[bin] = valbin/norm;

            if ((fluorescencephotonlikelihoods[bin]<0.0) ||
                (BAYES_DM_DOUBLE_TYPE_INVALID==bayes_dm_CheckDoubleValueValid(fluorescencephotonlikelihoods[bin],&type)))
            {
                free_Bayes_dvector(chiL,0,instr->ninstr-1);
                free_Bayes_dvector(chiH,0,instr->ninstr-1);

                if (DebugTrace)
                {
                    fprintf(fp,"****  ERROR: INVALID LIKELIHOOD DETECTED  ****\n\n");
                    fprintf(fp,"****  ERROR: bin %d, likelihood %g\n  ****",bin,fluorescencephotonlikelihoods[bin]);
                    fclose(fp);
                }

                return (-2);
            }

            if (DebugTrace)
            {
                fprintf(fp,"%g\n",fluorescencephotonlikelihoods[bin]);
            }
        }
        else
        {
            fluorescencephotonlikelihoods[bin] = 0.0;
            valid = 0;
        }
    }

    free_Bayes_dvector(chiL,0,instr->ninstr-1);
    free_Bayes_dvector(chiH,0,instr->ninstr-1);

    if (DebugTrace)
    {
        fclose(fp);
    }
    
    return (0);
}
#endif

#if 1 //remove all old code
#if 1

#if 0
void bayes_ComputeArrBinLikelihoodConstantUpsilon1(double **upsilon1,/* array of length nbins+1 */
                                                   int     *data,
                                                   int      nbins,
                                                   double   interval,
                                                   double   width,
                                                   double   delay)
{
    double oneoverwidthroottwo, binwidth;
    int    bin, valid;

    if (width>0.0)
    {
        binwidth            = interval/(double)nbins;
        oneoverwidthroottwo = 1.0/(width*ROOTTWO);

        /* Lifetime independent bin endpoint values (occupied bins only)... */
        for (bin=0,valid=0; bin<nbins; bin++)
        {
            if ((!data) || (data[bin])) /* If data is not provided compute all endpoint values... */
            {
                if (valid)
                {
                    (*upsilon1)[bin+1] = erf((delay-(bin+1)*binwidth)*oneoverwidthroottwo);
                    valid              = 1;
                }
                else
                {
                    (*upsilon1)[bin]   = erf((delay-bin*binwidth)*oneoverwidthroottwo);
                    (*upsilon1)[bin+1] = erf((delay-(bin+1)*binwidth)*oneoverwidthroottwo);
                    valid              = 1;
                }
            }
            else
            {
                valid = 0;
            }
        }    
    }
}
#endif

void bayes_CreateAndPopulateVectorInstrRspConstantGammaTilde(double **gammatilde, BayesInstrRsp_t *instr)
{
    int    i;
    double weight, delay, width, cutoff;

    *gammatilde = Bayes_dvector(0,instr->ninstr-1);

    for (i=0; i<instr->ninstr; i++)
    {
        weight = instr->params[i].weight;
        delay  = instr->params[i].delay;
        width  = instr->params[i].width;
        cutoff = instr->params[i].cutoff;

        (*gammatilde)[i] = weight/(1.0+erf((delay-cutoff)/(width*ROOTTWO)));
    }
}

double bayes_ComputeArrBinLikelihoodConstantUpsilon2(double delay,
                                                     double width,
                                                     double tau)
{
    if (width>0.0)
    {
        return (erf((delay/(width*ROOTTWO))+(width/(tau*ROOTTWO))));
    }
    else
    {
        return (0.0);
    }
}


double bayes_ComputeArrBinLikelihoodConstantUpsilon3(double repetitonperiod,
                                                     double tau)
{
    return (exp(-repetitonperiod/tau));
}


double bayes_ComputeArrBinLikelihoodConstantUpsilon4(double width,
                                                     double tau)
{
    if (width>0.0)
    {
        return (exp(-0.5*width*width/(tau*tau)));
    }
    else
    {
        return (0.0);
    }
}


double bayes_ArrBinLikelihoodEndpointValue(double x,
                                           double upsilon1i,
                                           double upsilon2,
                                           double upsilon3,
                                           double upsilon4,
                                           double interval,
                                           double width,
                                           double delay,
                                           double tau)
{
    double val, delayminusx, oneminusupsilon3;

    if (width <= 0.0)
    {
        if (delay <= 0.0)
        {
            return (-exp(-x/tau));
        }
        else /* Delay only... */
        {
            if (x < delay)
                return (upsilon3*exp(-(x-delay)/tau));
            else
                return (exp(-(x-delay)/tau));
        }
    }
    else
    {
        delayminusx      = delay-x;
        oneminusupsilon3 = 1.0-upsilon3;

        val  = upsilon2;
        val -= oneminusupsilon3*erf((delayminusx/(width*ROOTTWO))+(width/(tau*ROOTTWO)));
        val += upsilon3;
        val *= exp(delayminusx/tau);
        val += oneminusupsilon3*upsilon4*upsilon1i;
    
        return (val);
    }
}


double bayes_ArrBinLikelihoodGivenTauNormConstant(double upsilon2,
                                                  double upsilon3,
                                                  double upsilon4,
                                                  double interval,
                                                  double width,
                                                  double delay,
                                                  double tau)
{
    double upsilon1, norm;

    if (width <= 0.0)
    {
        if (delay <= 0.0)
        {
            return (1.0-exp((delay-interval)/tau));
        }
        else
        {
            return (-1.0*(1.0-exp((delay-interval)/tau)+(upsilon3*(exp(delay/tau)-1))));
        }
    }
    else
    {
        upsilon1  = erf((delay-interval)/(width*ROOTTWO));
        norm      = bayes_ArrBinLikelihoodEndpointValue(interval,upsilon1,upsilon2,upsilon3,upsilon4,interval,width,delay,tau);
        upsilon1  = erf((delay)/(width*ROOTTWO));
        norm     -= bayes_ArrBinLikelihoodEndpointValue(0.0,upsilon1,upsilon2,upsilon3,upsilon4,interval,width,delay,tau);

        return (norm);
    }
}




/***********************************************************************************/
/*                                                                                 */
/*                         LIKELIHOOD / PROBABILITY ROUTINES                       */
/*                                                                                 */
/***********************************************************************************/


#if 1 //new model


#if 0
//new form, takes a vector which defines m+1 bin wall locations (makes it easier in future to have bins of unequal widths)
void bayes_ComputeArrBinLikelihoodConstantUpsilon1(double  *binwalls,
                                                   double **upsilon1,/* array of length nbins+1 */
                                                   int     *data,
                                                   int      nbins,
                                                   double   interval,
                                                   double   width,
                                                   double   delay)
{
    double oneoverwidthroottwo, binwidth;
    int    bin, valid;

    if (width>0.0)
    {
        binwidth            = interval/(double)nbins;
        oneoverwidthroottwo = 1.0/(width*ROOTTWO);

        /* Lifetime independent bin endpoint values (occupied bins only)... */
        for (bin=0,valid=0; bin<nbins; bin++)
        {
            if ((!data) || (data[bin])) /* If data is not provided compute all endpoint values... */
            {
                if (valid)
                {
                    (*upsilon1)[bin+1] = erf((delay-(bin+1)*binwidth)*oneoverwidthroottwo);
                    valid              = 1;
                }
                else
                {
                    (*upsilon1)[bin]   = erf((delay-bin*binwidth)*oneoverwidthroottwo);
                    (*upsilon1)[bin+1] = erf((delay-(bin+1)*binwidth)*oneoverwidthroottwo);
                    valid              = 1;
                }
            }
            else
            {
                valid = 0;
            }
        }    
    }
}
#endif

//new form, takes a vector which defines m+1 bin wall locations (makes it easier in future to have bins of unequal widths)
void bayes_CreateAndPopulateMatrixArrBinLikelihoodConstantUpsilon1(double        ***upsilon1,
                                                                   double          *binwalls,
                                                                   int              nbins,
                                                                   int             *data,
                                                                   int              ellmax,
                                                                   double           modperiod,
                                                                   BayesInstrRsp_t *instr)
{
    int     bin, ell, i;
    double  val, delay, width;

    *upsilon1 = Bayes_dmatrix(0,instr->ninstr-1,0,nbins);

    for (i=0; i<instr->ninstr; i++)
    {
        delay = instr->params[i].delay;
        width = instr->params[i].width;

        for (bin=0; bin<=nbins; bin++)
        {
            for (ell=0,val=0.0; ell<ellmax; ell++)
                val += erf((modperiod*(double)ell+binwalls[bin]-delay)/(width*ROOTTWO));            

            (*upsilon1)[i][bin] = val;
            
        }    
    }
}


int  bayes_ArrBinLikelihoodsGivenTau(double          *likelihoods,
                                     double          *binwalls_g,
                                     double          *upsilon1_g,
                                    // ArrLikelihoodConstants_t *constants_g,
                                     int             *data,
                                     int              nbins,
                                     double           interval,
                                     double           modperiod,
                                     BayesInstrRsp_t *instr,
                                     double           tau)
{
    int      bin, valid, i, ell, ellmax;
    double  *binwalls, bL, bH, *chiL, *chiH, valbin;
    double   upsilon1, upsilon2, upsilon3, upsilon4;
    double   gammaprime[BAYES_INSTR_RSP_MAX_COMPONENTS], delay, width, cutoff, weight;

    ellmax = bayes_ConfigDecayModelMaxNumberOfModulationPeriodsConsidered();

    if (binwalls_g)
    {
        binwalls = binwalls_g;
    }
    else /* Default to uniform width bins spanning the measurement interval... */
    {
        binwalls = Bayes_dvector(0,nbins);
        bayes_PopulateBinWallsVectorUniformIntervals(binwalls,nbins,interval);
    }

    for (i=0; i<instr->ninstr; i++)
    {
        weight = instr->params[i].weight;
        delay  = instr->params[i].delay;
        width  = instr->params[i].width;
        cutoff = instr->params[i].cutoff;

        gammaprime[i] = weight/(1.0+erf((delay-cutoff)/(width*ROOTTWO)));
    }

    chiL = Bayes_dvector(0,instr->ninstr-1);
    chiH = Bayes_dvector(0,instr->ninstr-1);

    for (bin=0,valid=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            bL = binwalls[bin];
            bH = binwalls[bin+1];

            if (valid)
            {
                for (i=0; i<instr->ninstr; i++)
                    chiL[i] = chiH[i];

                /* Summation over 'ell' for each instrument response component... */
                for (i=0; i<instr->ninstr; i++)
                {
                    delay  = instr->params[i].delay;
                    width  = instr->params[i].width;
                    cutoff = instr->params[i].cutoff;

                    for (ell=0,chiH[i]=0.0; ell<ellmax; ell++)
                    {
                        upsilon1 = erf((modperiod*(double)ell+bH-delay)/(width*ROOTTWO));
                        upsilon2 = erf(((delay-cutoff)*tau+width*width)/(width*tau*ROOTTWO));
                        upsilon3 = erf(((delay-modperiod*(double)ell-bH)*tau+width*width)/(width*tau*ROOTTWO));
                        upsilon4 = exp((-(modperiod*(double)ell+bH-delay)/tau)+((width*width)/(2.0*tau*tau)));

                        chiH[i] += upsilon1+upsilon4*(upsilon3-upsilon2);
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

                    for (ell=0,chiL[i]=0.0,chiH[i]=0.0; ell<ellmax; ell++)
                    {
                        upsilon1 = erf((modperiod*(double)ell+bL-delay)/(width*ROOTTWO));
                        upsilon2 = erf(((delay-cutoff)*tau+width*width)/(width*tau*ROOTTWO));
                        upsilon3 = erf(((delay-modperiod*(double)ell-bL)*tau+width*width)/(width*tau*ROOTTWO));
                        upsilon4 = exp((-(modperiod*(double)ell+bL-delay)/tau)+((width*width)/(2.0*tau*tau)));

                        chiL[i] += upsilon1+upsilon4*(upsilon3-upsilon2);

                        upsilon1 = erf((modperiod*(double)ell+bH-delay)/(width*ROOTTWO));
                        upsilon2 = erf(((delay-cutoff)*tau+width*width)/(width*tau*ROOTTWO));
                        upsilon3 = erf(((delay-modperiod*(double)ell-bH)*tau+width*width)/(width*tau*ROOTTWO));
                        upsilon4 = exp((-(modperiod*(double)ell+bH-delay)/tau)+((width*width)/(2.0*tau*tau)));

                        chiH[i] += upsilon1+upsilon4*(upsilon3-upsilon2);
                    }
                }

                valid = 1;
            }

            /* Now do the sum over all instrument response components... */
            for (i=0,valbin=0.0; i<instr->ninstr; i++)
                valbin += gammaprime[i]*(chiH[i]-chiL[i]);

            likelihoods[bin] = valbin;
        }
        else
        {
            likelihoods[bin] = 0.0;
            valid = 0;
        }
    }

    free_Bayes_dvector(chiL,0,instr->ninstr-1);
    free_Bayes_dvector(chiH,0,instr->ninstr-1);
    if (!binwalls_g) free_Bayes_dvector(binwalls,0,nbins);
    
    return (0);
}
#if 0
int  bayes_ArrBinLikelihoodsGivenTau(double          *likelihoods,
                                     double          *binwalls_g,
                                     double          *upsilon1_g,
                                     int             *data,
                                     int              nbins,
                                     double           interval,
                                     double           modperiod,
                                     BayesInstrRsp_t *instr,
                                     double           tau)
{
    int    bin, valid, ret;
    double *binwalls, xl, xh, numl, numh, den, val;
    double *upsilon1, upsilon2, upsilon3, upsilon4;
    int    type;

    if (binwalls_g)
    {
        binwalls = binwalls_g;
    }
    else /* Default to uniform width bins spanning the measurement interval... */
    {
        binwalls = Bayes_dvector(0,nbins);
        bayes_PopulateBinWallsVectorUniformIntervals(binwalls,nbins,interval);
    }

    if (upsilon1_g) /* Bin-likelihood quantity that does not depend on the lifetime... */
    {
        upsilon1 = upsilon1_g;
    }
    else
    {
        upsilon1 = Bayes_dvector(0,nbins);
        bayes_ComputeArrBinLikelihoodConstantUpsilon1(&upsilon1,data,nbins,interval,width,delay);
    }
    
    if (width > 0.0)
    {
        /* New tau, compute new tau-dependent 'constants' for subsequent routines... */
        upsilon2 = bayes_ComputeArrBinLikelihoodConstantUpsilon2(delay,width,tau);
        
        if (bayes_UseRepetitionEffectsInAnalysis())
            upsilon3 = bayes_ComputeArrBinLikelihoodConstantUpsilon3(modperiod,tau);
        else
            upsilon3 = 0.0;

        upsilon4 = bayes_ComputeArrBinLikelihoodConstantUpsilon4(width,tau);
    }
    else
    {
        upsilon2 = 0.0;

        if (bayes_UseRepetitionEffectsInAnalysis())
            upsilon3 = bayes_ComputeArrBinLikelihoodConstantUpsilon3(modperiod,tau);
        else
            upsilon3 = 0.0;

        upsilon4 = 0.0;    
    }

    den = bayes_ArrBinLikelihoodGivenTauNormConstant(upsilon2,upsilon3,upsilon4,interval,width,delay,tau);

    for (bin=0,valid=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            xl = (double)bin*xincr;
            xh = xl+xincr;

            if (valid) /* If width<=0.0 and valid has been set, already at bin where the delay no longer needs to be considered... */
            {
                numl = numh;
                numh = bayes_ArrBinLikelihoodEndpointValue
                        (xh,upsilon1[bin+1],upsilon2,upsilon3,upsilon4,interval,width,delay,tau);
            }
            else
            {
                numl = bayes_ArrBinLikelihoodEndpointValue
                         (xl,upsilon1[bin],upsilon2,upsilon3,upsilon4,interval,width,delay,tau);
                numh = bayes_ArrBinLikelihoodEndpointValue
                         (xh,upsilon1[bin+1],upsilon2,upsilon3,upsilon4,interval,width,delay,tau);

                valid = 1;
            }

            if ((width<=0.0) && (delay>0.0) && (xl<delay) && (delay<xh))
                numl  += 1.0-upsilon3;

            val = (numh-numl)/(den*xincr);

            /* An error has occured (likely to underflow/overflow)... */
            if ((val<0.0) || (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(val,&type)))
            {
                val = 0.0;
            }

            likelihoods[bin] = val;
        }
        else
        {
            likelihoods[bin] = 0.0;
            valid = 0;
        }
    }

    if (!upsilon1_g)
        free_Bayes_dvector(upsilon1,0,nbins);
    
    return (0);
}
#endif
#else
int  bayes_ArrBinLikelihoodsGivenTau(double *likelihoods,
                                     double *upsilon1_g,
                                     int    *data,
                                     int     nbins,
                                     double  interval,
                                     double  width,
                                     double  delay,
                                     double  tau,
                                     double  modperiod)
{
    int    bin, valid, ret;
    double xincr, xl, xh, numl, numh, den, val;
    double *upsilon1, upsilon2, upsilon3, upsilon4;
    int    type;
    
#ifdef BAYES_DIST_FCTS_DEBUG_OUTPUT
    FILE   *fp;
    int     i, error=0;
#endif

    xincr = interval/(double)nbins;

    if (upsilon1_g) /* This may have been pre-computed... */
    {
        upsilon1 = upsilon1_g;
    }
    else
    {
        upsilon1 = Bayes_dvector(0,nbins);
        bayes_ComputeArrBinLikelihoodConstantUpsilon1(&upsilon1,data,nbins,interval,width,delay);
    }
    
    if (width > 0.0)
    {
        /* New tau, compute new tau-dependent 'constants' for subsequent routines... */
        upsilon2 = bayes_ComputeArrBinLikelihoodConstantUpsilon2(delay,width,tau);
        
        if (bayes_UseRepetitionEffectsInAnalysis())
            upsilon3 = bayes_ComputeArrBinLikelihoodConstantUpsilon3(modperiod,tau);
        else
            upsilon3 = 0.0;

        upsilon4 = bayes_ComputeArrBinLikelihoodConstantUpsilon4(width,tau);
    }
    else
    {
        upsilon2 = 0.0;

        if (bayes_UseRepetitionEffectsInAnalysis())
            upsilon3 = bayes_ComputeArrBinLikelihoodConstantUpsilon3(modperiod,tau);
        else
            upsilon3 = 0.0;

        upsilon4 = 0.0;    
    }

    den = bayes_ArrBinLikelihoodGivenTauNormConstant(upsilon2,upsilon3,upsilon4,interval,width,delay,tau);

    for (bin=0,valid=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            xl = (double)bin*xincr;
            xh = xl+xincr;

            if (valid) /* If width<=0.0 and valid has been set, already at bin where the delay no longer needs to be considered... */
            {
                numl = numh;
                numh = bayes_ArrBinLikelihoodEndpointValue
                        (xh,upsilon1[bin+1],upsilon2,upsilon3,upsilon4,interval,width,delay,tau);
            }
            else
            {
                numl = bayes_ArrBinLikelihoodEndpointValue
                         (xl,upsilon1[bin],upsilon2,upsilon3,upsilon4,interval,width,delay,tau);
                numh = bayes_ArrBinLikelihoodEndpointValue
                         (xh,upsilon1[bin+1],upsilon2,upsilon3,upsilon4,interval,width,delay,tau);

                valid = 1;
            }

            if ((width<=0.0) && (delay>0.0) && (xl<delay) && (delay<xh))
                numl  += 1.0-upsilon3;

            val = (numh-numl)/(den*xincr);

            /* An error has occured (likely to underflow/overflow)... */
            if ((val<0.0) || (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(val,&type)))
            {
                val = 0.0;

#ifdef BAYES_DIST_FCTS_DEBUG_OUTPUT
                error = 1;
#endif
            }

            likelihoods[bin] = val;
        }
        else
        {
            likelihoods[bin] = 0.0;
            valid = 0;
        }
    }

#ifdef BAYES_DIST_FCTS_DEBUG_OUTPUT
    if (error)
    {
        /* Debug output */
        if fp = fopen(("BayesDistFctsDebugOutput.txt", "w"))))
            return (-1);

        fprintf(fp,"Upsilon1\n\n");
        for (i=0; i<nbins; i++)
            fprintf(fp,"%g\t%g\n", upsilon1[i], likelihoods[i]);

        fclose(fp);
        /* Debug output */    

        if (!upsilon1_g)
            free_Bayes_dvector(upsilon1,0,nbins);

        return (-1);
    }
#endif
    if (!upsilon1_g)
        free_Bayes_dvector(upsilon1,0,nbins);
    
    return (0);
}
#endif
#else
void bayes_ComputeArrBinLikelihoodConstantUpsilon1(double **upsilon1,/* array of length nbins+1 */
                                                   int     *data,
                                                   int      nbins,
                                                   double   interval,
                                                   double   width,
                                                   double   delay)
{
    double oneoverwidthroottwo, binwidth;
    int    bin, valid;

    if (width>0.0)
    {
        binwidth            = interval/(double)nbins;
        oneoverwidthroottwo = 1.0/(width*ROOTTWO);

        /* Lifetime independent bin endpoint values (occupied bins only)... */
        for (bin=0,valid=0; bin<nbins; bin++)
        {
            if ((!data) || (data[bin])) /* If data is not provided compute all endpoint values... */
            {
                if (valid)
                {
                    (*upsilon1)[bin+1] = erf((delay-(bin+1)*binwidth)*oneoverwidthroottwo);
                    valid              = 1;
                }
                else
                {
                    (*upsilon1)[bin]   = erf((delay-bin*binwidth)*oneoverwidthroottwo);
                    (*upsilon1)[bin+1] = erf((delay-(bin+1)*binwidth)*oneoverwidthroottwo);
                    valid              = 1;
                }
            }
            else
            {
                valid = 0;
            }
        }    
    }
}

//used to be 3, now 2
double bayes_ComputeArrBinLikelihoodConstantUpsilon2(double delay,
                                                     double width,
                                                     double tau)
{
    if (width>0.0)
    {
        return (erf((delay/(width*ROOTTWO))+(width/(tau*ROOTTWO))));
    }
    else
    {
        return (0.0);
    }
}

//new upsilon 3
double bayes_ComputeArrBinLikelihoodConstantUpsilon3(double repetitonperiod,
                                                     double tau)
{
    return (exp(-repetitonperiod/tau));
}

//this has not changed
double bayes_ComputeArrBinLikelihoodConstantUpsilon4(double width,
                                                     double tau)
{
    if (width>0.0)
    {
        return (exp(-0.5*width*width/(tau*tau)));
    }
    else
    {
        return (0.0);
    }
}

double bayes_ArrBinLikelihoodEndpointValue(double x,
                                           double upsilon1i,
                                           double upsilon2,
                                           double upsilon3,
                                           double upsilon4,
                                           double interval,
                                           double width,
                                           double delay,
                                           double tau)
{
    double val, delayminusx, oneminusupsilon3; int type;

    if (width <= 0.0)
    {
        if (delay <= 0.0)
        {
            return (-exp(-x/tau));
        }
        else /* Delay only... */
        {
            if (x < delay)
                return (upsilon3*exp(-(x-delay)/tau));
            else
                return (exp(-(x-delay)/tau));
        }
    }
    else
    {
        delayminusx      = delay-x;
        oneminusupsilon3 = 1.0-upsilon3;

        val  = upsilon2;
        val -= oneminusupsilon3*erf((delayminusx/(width*ROOTTWO))+(width/(tau*ROOTTWO)));
        val += upsilon3;
        val *= exp(delayminusx/tau);
        val += oneminusupsilon3*upsilon4*upsilon1i;

        return (val);
    }
}

double bayes_ArrBinLikelihoodGivenTauNormConstant(double upsilon2,
                                                  double upsilon3,
                                                  double upsilon4,
                                                  double interval,
                                                  double width,
                                                  double delay,
                                                  double tau)
{
    double upsilon1, norm;

    if (width <= 0.0)
    {
        if (delay <= 0.0)
        {
            return (1.0-exp((delay-interval)/tau));
        }
        else
        {
            return (-1.0*(1.0-exp((delay-interval)/tau)+(upsilon3*(exp(delay/tau)-1))));
        }
    }
    else
    {
        upsilon1  = erf((delay-interval)/(width*ROOTTWO));
        norm      = bayes_ArrBinLikelihoodEndpointValue(interval,upsilon1,upsilon2,upsilon3,upsilon4,interval,width,delay,tau);
        upsilon1  = erf((delay)/(width*ROOTTWO));
        norm     -= bayes_ArrBinLikelihoodEndpointValue(0.0,upsilon1,upsilon2,upsilon3,upsilon4,interval,width,delay,tau);

        return (norm);
    }
}

int bayes_ArrBinLikelihoodsGivenTau(double *likelihoods,
                                     double *upsilon1_g,
                                     int    *data,
                                     int     nbins,
                                     double  interval,
                                     double  width,
                                     double  delay,
                                     double  tau,
                                     double  modperiod)
{
    int    bin, valid;
    double xincr, xl, xh, numl, numh, den, val;
    double *upsilon1, upsilon2, upsilon3, upsilon4;

    xincr = interval/(double)nbins;

    if (upsilon1_g) /* This may have been pre-computed... */
    {
        upsilon1 = upsilon1_g;
    }
    else
    {
        upsilon1 = Bayes_dvector(0,nbins);
        bayes_ComputeArrBinLikelihoodConstantUpsilon1(&upsilon1,data,nbins,interval,width,delay);
    }
    
    if (width > 0.0)
    {
        /* New tau, compute new tau-dependent 'constants' for subsequent routines... */
        upsilon2 = bayes_ComputeArrBinLikelihoodConstantUpsilon2(delay,width,tau);
        
        if (bayes_GetIncludeRepetitionEffectsInAnalysisFlag())
            upsilon3 = bayes_ComputeArrBinLikelihoodConstantUpsilon3(modperiod,tau);
        else
            upsilon3 = 0.0;

        upsilon4 = bayes_ComputeArrBinLikelihoodConstantUpsilon4(width,tau);
    }
    else
    {
        upsilon2 = 0.0;

        if (bayes_GetIncludeRepetitionEffectsInAnalysisFlag())
            upsilon3 = bayes_ComputeArrBinLikelihoodConstantUpsilon3(modperiod,tau);
        else
            upsilon3 = 0.0;

        upsilon4 = 0.0;    
    }

    den = bayes_ArrBinLikelihoodGivenTauNormConstant(upsilon2,upsilon3,upsilon4,interval,width,delay,tau);

    for (bin=0,valid=0; bin<nbins; bin++)
    {
        if ((!data) || (data[bin]))
        {
            xl = (double)bin*xincr;
            xh = xl+xincr;

            if (valid) /* If width<=0.0 and valid has been set, already at bin where the delay no longer needs to be considered... */
            {
                numl = numh;
                numh = bayes_ArrBinLikelihoodEndpointValue
                        (xh,upsilon1[bin+1],upsilon2,upsilon3,upsilon4,interval,width,delay,tau);
            }
            else
            {
                numl = bayes_ArrBinLikelihoodEndpointValue
                         (xl,upsilon1[bin],upsilon2,upsilon3,upsilon4,interval,width,delay,tau);
                numh = bayes_ArrBinLikelihoodEndpointValue
                         (xh,upsilon1[bin+1],upsilon2,upsilon3,upsilon4,interval,width,delay,tau);

                valid = 1;
            }

            if ((width<=0.0) && (delay>0.0) && (xl<delay) && (delay<xh))
                numl  += 1.0-upsilon3;

            val = (numh-numl)/(den*xincr);

            if ((val<0.0) || (numh>numl))
                return (-1);

            likelihoods[bin] = val;
        }
        else
        {
            likelihoods[bin] = 0.0;
            valid = 0;
        }
    }

    return (0);
}
#endif









int bayes_ComputeBinLikelihoodsFromWeightsAndFluorescencePhotonLikelihoods(double           *binlikelihoods,
                                                                           int               nbins,
                                                                           double           *binwalls,
                                                                           int               ndecays,
                                                                           double          **fluorescencephotonlikelihoods,
                                                                           double           *weights,
                                                                           BayesInstrRsp_t  *instr,
                                                                           double            interval)
{
    double bL, bH, bjoverT, value;
    int    bin, k;    

    if ((!binlikelihoods) || (!binwalls) || (nbins<=1))
        return (-1);
    
    if ((!instr) || (interval<=0.0))
        return (-2);

    if ((ndecays<1) || (!fluorescencephotonlikelihoods) || (!weights))
        return (-3);

    for (k=1; k<=ndecays; k++)
        if ((weights[k]<0.0) || (weights[k]>1.0))
            return (-4);

    for (bin=0; bin<nbins; bin++)
    {
        bL      = binwalls[bin];
        bH      = binwalls[bin+1];
        bjoverT = (bH-bL)/interval;

        for (k=1,value=bjoverT; k<=ndecays; k++)
            value += weights[k]*(fluorescencephotonlikelihoods[k][bin]-bjoverT);

        binlikelihoods[bin] = value;
    }

    return (0);
}

#endif //remove all old code








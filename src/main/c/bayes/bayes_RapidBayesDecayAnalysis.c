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
#include "math.h"
#include "stdlib.h"
#include "extmath.h"
#include "matrices.h"
#include "random.h"
#include "bayes_DistributionFctsBinLikelihoods.h"
#include "bayes_InstrRspAnalysis.h"
#include "bayes_Sizes.h"
#include "bayes_Types.h"
#include "bayes_DataManagement.h"

#include "bayes_RapidBayesDecayAnalysis.h"



/*=======================================================*/
/*                     LOCAL DEFINES                     */
/*=======================================================*/

#define BAYES_RAPID_NTHREADS 8

typedef struct
{
    BayesRapidLikelihoodValues_t       *fluorescencelikelihoods;
    int                                *low;
    int                                *high;
    int                                 ntaus;
    double                             *taus;
    int                                 nbins;
    double                             *binwalls;
    BayesInstrRsp_t                    *instr;
    double                              interval;
    double                              modulationperiod;    
    int                                 nvalid[1];
    int                                 ninvalid[1];
}
RapidLikelihoodsDiscreteThreadData;


typedef struct
{
    BayesRapidMonoExpDiscreteValues_t **likelihoods;
    int                                *low;
    int                                *high;
    int                                 ntaus;
    double                             *taus;
    int                                 nweights;
    double                             *weights;
    double                              backgroundmax;
    int                                 nbins;
    double                             *binwalls;
    BayesInstrRsp_t                    *instr;
    double                              interval;
    double                              modulationperiod;
    BayesRapidLikelihoodValues_t       *fluorescencelikelihoods;
    int                                 nvalid[1];
    int                                 ninvalid[1];
}
RapidMonoExpDiscreteThreadData;


typedef struct
{
    BayesRapidBiExpDiscreteValues_t *****likelihoods;
    int                                 *low;
    int                                 *high;
    int                                  ntaus;
    double                              *taus;
    int                                  nweights;
    double                              *weights;
    double                               backgroundmin;
    double                               backgroundmax;
    int                                  nbins;
    double                              *binwalls;
    BayesInstrRsp_t                     *instr;
    double                               interval;
    double                               modulationperiod;
    BayesRapidLikelihoodValues_t        *fluorescencelikelihoods;
    int                                  nvalid[1];
    int                                  ninvalid[1];
}
RapidBiExpDiscreteThreadData;


#define THRPTRLIKELIHOOD ((RapidLikelihoodsDiscreteThreadData *)p)
#define THRPTRMONO       ((RapidMonoExpDiscreteThreadData *)p)
#define THRPTRBI         ((RapidBiExpDiscreteThreadData *)p)


/*=======================================================*/
/*          RAPID VALUE STORE ACCESS FUNCTIONS           */
/*=======================================================*/

BayesRapidValueStore_t *bayes_GetRapidValueStorePtr(void)
{
	return bayes_GetRapidValueStorePtrSafe();
}

/*=======================================================*/
/*              GENERAL UTILITY FUNCTIONS                */
/*=======================================================*/

int bayes_MapWeightValueToClosestRapidGridPoint(double value, int nweights, double *weights)
{
    int    i, j;
    double d;

    for (i=0, j=nweights/2, d=BAYES_SIZE_DOUBLE_HUGE; i<nweights; i++)
    {
        if (mod(value-weights[i]) < d)
        {
            d = mod(value-weights[i]);
            j = i;
        }
    }

    return (j);
}


int bayes_MapLifetimeValueToClosestRapidGridPoint(double value, int ntaus, double *taus)
{
    int    i, j;
    double d;

    for (i=0, j=ntaus/2, d=BAYES_SIZE_DOUBLE_HUGE; i<ntaus; i++)
    {
        if (mod(value-taus[i]) < d)
        {
            d = mod(value-taus[i]);
            j = i;
        }
    }

    return (j);
}



/*=======================================================*/
/*              MEMORY MANAGEMENT FUNCTIONS              */
/*=======================================================*/

BayesRapidLikelihoodValues_t *bayes_AllocateRapidLikelihoodValuesVector(int nentries)
{
    BayesRapidLikelihoodValues_t *m;

    m = (BayesRapidLikelihoodValues_t*)malloc((size_t)(nentries*sizeof(BayesRapidLikelihoodValues_t)));

    return (m);
}


void bayes_FreeRapidLikelihoodValuesVector(BayesRapidLikelihoodValues_t *m)
{
    free((char*)(m));
}


#define NR_END   1
BayesRapidMonoExpDiscreteValues_t **bayes_AllocateRapidMonoExpDiscreteValuesMatrix(int nrl,int nrh,int ncl,int nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long                                i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	BayesRapidMonoExpDiscreteValues_t **m;

	/* allocate pointers to rows */
	m = (BayesRapidMonoExpDiscreteValues_t **)malloc((size_t)((nrow+NR_END)*sizeof(BayesRapidMonoExpDiscreteValues_t*)));
	
    if (!m)
        return (NULL);

	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (BayesRapidMonoExpDiscreteValues_t *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(BayesRapidMonoExpDiscreteValues_t)));
	
    if (!m[nrl])
        return (NULL);

	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1; i<=nrh; i++)
        m[i] = m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return (m);
}


void bayes_FreeRapidMonoExpDiscreteValuesMatrix(BayesRapidMonoExpDiscreteValues_t **m,
                                                int nrl, int nrh,
                                                int ncl, int nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	free((char*)(m[nrl]+ncl-NR_END));
	free((char*)(m+nrl-NR_END));
}


int bayes_PopulateRapidMonoExpDiscreteValueMatrix(BayesRapidMonoExpDiscreteValues_t **likelihoods,
                                                  int                                *low,
                                                  int                                *high,
                                                  int                                 ntaus,
                                                  double                             *taus,
                                                  int                                 nweights,
                                                  double                             *weights,
                                                  double                              backgroundmax,
                                                  int                                 nbins,
                                                  double                             *binwalls,
                                                  BayesInstrRsp_t                    *instr,
                                                  double                              interval,
                                                  double                              modulationperiod,
                                                  BayesRapidLikelihoodValues_t       *fluorescencelikelihoods,
                                                  int                                *nvalid,
                                                  int                                *ninvalid)
{
    int    i, j, bin, ndecays, ret, type;
    double norm, w[2], t[2], tau;
    double w0, oneminusw0, bL, bH, bjoverT;
    double photonlikelihoodgiventauandw0, *binlikelihood, *loglikelihoods;

    int    DebugTrace=0;

    for (j=low[2]; j<=high[2]; j++) /* Loop over lifetime values... */
    {
        if (fluorescencelikelihoods[j].valid)
        {
            tau           = fluorescencelikelihoods[j].tau;
            binlikelihood = fluorescencelikelihoods[j].fluorescencedecayphotonlikelihoodsgiventau;

            for (i=low[1]; i<=high[1]; i++) /* Loop over the 'w0' values... */
            {
                ndecays    = 1;
                w0         = weights[i];
                oneminusw0 = 1.0-w0;
                w[0]       = w0;
                w[1]       = oneminusw0;
                t[1]       = tau;
                
                likelihoods[i][j].w0  = w0;
                likelihoods[i][j].tau = tau;
                
                ret = bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(
                          &norm,interval,modulationperiod,/*binwalls[fitstart]*/0.0,instr,ndecays,w,t);

                if ((ret<0) && (DebugTrace)) /* Create error file... */
                {
                    likelihoods[i][j].w0  = -1.0;
                    likelihoods[i][j].tau = -1.0;
                }

                
                loglikelihoods = likelihoods[i][j].logphotonlikelihoodgiventauandw0;

                for (bin=0; bin<nbins; bin++)
                {
                    bH                            = binwalls[bin+1];
                    bL                            = binwalls[bin];
                    bjoverT                       = (bH-bL)/interval;
                    photonlikelihoodgiventauandw0 = w0*bjoverT + oneminusw0*binlikelihood[bin]/norm;

                    /* NOTE: Dither correction is not applied in the dedicated         */
                    /* mono-exponential grid, the psuedo rapid grid is used instead... */

                    loglikelihoods[bin] = log(photonlikelihoodgiventauandw0);

                    if (BAYES_DM_DOUBLE_TYPE_INVALID == bayes_dm_CheckDoubleValueValid(loglikelihoods[bin],&type))
                    {
                        loglikelihoods[bin] = -BIG;
                    }
                }

                (*nvalid)++;
            }
        }
        else /* Likelihoods have not been computed succesfully...*/
        {
            (*ninvalid)++;                    
        }
    }

    return (0);
}

/*
int CVICALLBACK bayes_PopulateRapidMonoExpDiscreteValueMatrix_ThreadFunction (void *p)
{
    int ret;

    ret = bayes_PopulateRapidMonoExpDiscreteValueMatrix(
              THRPTRMONO->likelihoods,THRPTRMONO->low,THRPTRMONO->high,
              THRPTRMONO->ntaus,THRPTRMONO->taus,THRPTRMONO->nweights,THRPTRMONO->weights,THRPTRMONO->backgroundmax,
              THRPTRMONO->nbins,THRPTRMONO->binwalls,THRPTRMONO->instr,THRPTRMONO->interval,THRPTRMONO->modulationperiod,
              THRPTRMONO->fluorescencelikelihoods,
              THRPTRMONO->nvalid,THRPTRMONO->ninvalid);

	return (ret);
}


void bayes_LoadRapidMonoExpDiscreteThreadData(RapidMonoExpDiscreteThreadData     *data,
                                              BayesRapidMonoExpDiscreteValues_t **likelihoods,
                                              int                                *low,
                                              int                                *high,
                                              int                                 ntaus,
                                              double                             *taus,
                                              int                                 nweights,
                                              double                             *weights,
                                              double                              backgroundmax,
                                              int                                 nbins,
                                              double                             *binwalls,
                                              BayesInstrRsp_t                    *instr,
                                              double                              interval,
                                              double                              modulationperiod,
                                              BayesRapidLikelihoodValues_t       *fluorescencelikelihoods)
{
	data->likelihoods             = likelihoods;
    data->low                     = low;
    data->high                    = high;
	data->ntaus	                  = ntaus;
	data->taus                    = taus;
	data->nweights                = nweights;
	data->weights                 = weights;
	data->backgroundmax	          = backgroundmax;
	data->nbins                   = nbins;
	data->binwalls                = binwalls;
	data->instr                   = instr;
	data->interval                = interval; 
	data->modulationperiod	      = modulationperiod; 
	data->fluorescencelikelihoods = fluorescencelikelihoods;
	*(data->nvalid)               = 0;
    *(data->ninvalid)             = 0;
}


int bayes_PopulateRapidMonoExpDiscreteValueMatrixMultiThreaded(BayesRapidMonoExpDiscreteValues_t **likelihoods,
                                                               int                                *low,
                                                               int                                *high,
                                                               int                                 ntaus,
                                                               double                             *taus,
                                                               int                                 nweights,
                                                               double                             *weights,
                                                               double                              backgroundmax,
                                                               int                                 nbins,
                                                               double                             *binwalls,
                                                               BayesInstrRsp_t                    *instr,
                                                               double                              interval,
                                                               double                              modulationperiod,
                                                               BayesRapidLikelihoodValues_t       *fluorescencelikelihoods,
                                                               int                                 nstates,
                                                               int                                *nvalid,
                                                               int                                *ninvalid)
{
	int complete[BAYES_RAPID_NTHREADS], nComplete, ThreadPoolHandle, ThreadID[BAYES_RAPID_NTHREADS], i, states;
    RapidMonoExpDiscreteThreadData   data[BAYES_RAPID_NTHREADS];
    CmtThreadFunctionExecutionStatus status;

    int ret=BAYES_RAPID_RESULT_NO_ERROR, threadlow[BAYES_RAPID_NTHREADS][3], threadhigh[BAYES_RAPID_NTHREADS][3], low1, high1, diff;

    diff = (high[1]-low[1]+1)/BAYES_RAPID_NTHREADS;

    for (i=0; i<BAYES_RAPID_NTHREADS; i++)
    {
            threadlow[i][2]  = low[2];
            threadhigh[i][2] = high[2];
    }

    CmtNewThreadPool(BAYES_RAPID_NTHREADS, &ThreadPoolHandle);

    // start nThreads-1 threads
    for (i=0,low1=low[1],high1=low[1]+diff-1; i<BAYES_RAPID_NTHREADS-1; i++)
    {
        threadlow[i][1]  = low1;
        threadhigh[i][1] = high1;

        bayes_LoadRapidMonoExpDiscreteThreadData(&data[i],
                                                 likelihoods,threadlow[i],threadhigh[i],
                                                 ntaus,taus,nweights,weights,backgroundmax,
                                                 nbins,binwalls,instr,interval,modulationperiod,
                                                 fluorescencelikelihoods);

        low1  = high1+1;
        high1 = low1+diff-1;

        CmtScheduleThreadPoolFunction(ThreadPoolHandle,bayes_PopulateRapidMonoExpDiscreteValueMatrix_ThreadFunction,(void*)&data[i],&ThreadID[i]);
        	
        complete[i]=0;
    }

    threadlow[i][1]   = low1;
    threadhigh[i][1]  = high[1];

    // last thread gets the remainder of the grid
    bayes_LoadRapidMonoExpDiscreteThreadData(&data[i],
                                             likelihoods,threadlow[i],threadhigh[i],
                                             ntaus,taus,nweights,weights,backgroundmax,
                                             nbins,binwalls,instr,interval,modulationperiod,
                                             fluorescencelikelihoods);
    
    CmtScheduleThreadPoolFunction(ThreadPoolHandle,bayes_PopulateRapidMonoExpDiscreteValueMatrix_ThreadFunction,(void*)&data[i],&ThreadID[i]);
    complete[i]=0;

    nComplete=0;

    while (nComplete < BAYES_RAPID_NTHREADS)
    {
        // this loop queries each thread to see if it has finished
        // Add up all the progress values from the threads to report an overall progress against the number of rows
        // Checks for an interupt, each thread does this also and will exit if it is detected.

        for (i=0; i<BAYES_RAPID_NTHREADS; i++)
        {
            if (complete[i]==0)
            {
                CmtGetThreadPoolFunctionAttribute(ThreadPoolHandle,ThreadID[i],ATTR_TP_FUNCTION_EXECUTION_STATUS,&status);
                
                if (status == kCmtThreadFunctionComplete)
                {
                    complete[i]=1;
                    nComplete++;
                    CmtReleaseThreadPoolFunctionID(ThreadPoolHandle,ThreadID[i]);
                }
            }
        }

        Delay(0.5);

        for (i=0,states=0; i<BAYES_RAPID_NTHREADS; i++)
        {
            states += (*(data[i].nvalid))+(*(data[i].ninvalid));
        }

        if (UpdateProgressDialog(progress,100*states/nstates,1))
		{   // user cancelled
            ret = BAYES_RAPID_RESULT_USER_CANCEL;
			break;
		}
    }

    for (i=0,*nvalid=0,*ninvalid=0; i<BAYES_RAPID_NTHREADS; i++)
    {
        *nvalid   += *data[i].nvalid;
        *ninvalid += *data[i].ninvalid;
    }

    CmtDiscardThreadPool(ThreadPoolHandle);

    return (ret);
}
*/

BayesRapidBiExpDiscreteValues_t *****bayes_AllocateRapidBiExpDiscreteValuesMatrix(int                            *low,
                                                                                  int                            *high,
                                                                                  int                            *nstates,
                                                                                  int                             ntaus,
                                                                                  double                         *taus,
                                                                                  int                             nweights,
                                                                                  double                         *weights,
                                                                                  int                             nbins,
                                                                                  BayesRapidBiExpMemoryControl_t *memory)
{
	long                                 i, j, k, ell, ndim[5], nvalid;
	BayesRapidBiExpDiscreteValues_t *****m;
    double                               sum, w0, w0low, w0high;

    for (i=1; i<=4; i++)
    {
        ndim[i] = high[i]-low[i]+1;
    }

    for (i=1,*nstates=1; i<=4; i++)
        *nstates *= ndim[i];

	/* Allocate the grid... */
	m = (BayesRapidBiExpDiscreteValues_t*****)malloc((size_t)(ndim[1]*sizeof(BayesRapidBiExpDiscreteValues_t****)));
	
    if (!m)
        return (NULL);

    m -= low[1];

    for (i=low[1]; i<=high[1]; i++)
    {
        m[i] = (BayesRapidBiExpDiscreteValues_t****)malloc((size_t)(ndim[2]*sizeof(BayesRapidBiExpDiscreteValues_t***)));

        if (!m[i])
            return (NULL);

	    m[i] -= low[2];

        for (j=low[2]; j<=high[2]; j++)
        {
            m[i][j] = (BayesRapidBiExpDiscreteValues_t***)malloc((size_t)(ndim[3]*sizeof(BayesRapidBiExpDiscreteValues_t**)));

            if (!m[i][j])
                return (NULL);

	        m[i][j] -= low[3];

            for (k=low[3]; k<=high[3]; k++)
            {
                m[i][j][k] = (BayesRapidBiExpDiscreteValues_t**)malloc((size_t)(ndim[4]*sizeof(BayesRapidBiExpDiscreteValues_t*)));

                if (!m[i][j][k])
                    return (NULL);

	            m[i][j][k] -= low[4];
            }
        }
    }

    /* Determine the size of the memory pools and allocate for valid grid elements... */
    memory->npools = high[1]-low[1]+1;
    memory->pools  = (BayesRapidBiExpMemoryPool_t*)malloc(memory->npools*sizeof(BayesRapidBiExpMemoryPool_t));

    w0low  = weights[low[0]];
    w0high = weights[high[0]];

    for (i=low[1]; i<=high[1]; i++)
    {
        for (j=low[2],nvalid=0; j<=high[2]; j++)
        {
            sum = weights[i]+weights[j];
            w0  = 1.0-sum;

            if ((sum>=0.0) && (sum<=1.0) && (w0>=w0low) && (w0<=w0high))
            {
                for (k=low[3]; k<=high[3]; k++)
                {
                    for (ell=low[4]; ell<=high[4]; ell++)
                    {
                        if (taus[k]>taus[ell])
                            nvalid++;
                    }
                }            
            }
        }

        memory->pools[i-low[1]].nvalid         = nvalid;
        memory->pools[i-low[1]].likelihoods    = (BayesRapidBiExpDiscreteValues_t*)malloc(nvalid*sizeof(BayesRapidBiExpDiscreteValues_t));
        memory->pools[i-low[1]].loglikelihoods = Bayes_dmatrix(0,nvalid,0,nbins-1);
    }

    /* Use pooled memory for the valid grid entries... */
    for (i=low[1]; i<=high[1]; i++)
    {
        for (j=low[2],nvalid=0; j<=high[2]; j++)
        {
            sum = weights[i]+weights[j];
            w0  = 1.0-sum;

            if ((sum>=0.0) && (sum<=1.0) && (w0>=w0low) && (w0<=w0high))
            {
                for (k=low[3]; k<=high[3]; k++)
                {
                    for (ell=low[4]; ell<=high[4]; ell++)
                    {
                        if (taus[k]>taus[ell])
                        {
                            m[i][j][k][ell] = &memory->pools[i-low[1]].likelihoods[nvalid];
                            m[i][j][k][ell]->logphotonlikelihoodgiventausandweights = memory->pools[i-low[1]].loglikelihoods[nvalid];
                            nvalid++;
                        }
                    }
                }            
            }
        }
    }

    return (m);
}


int bayes_FreeRapidBiExpDiscreteValuesMatrix(BayesRapidBiExpDiscreteValues_t *****m,
                                             int                                 *low,
                                             int                                 *high,
                                             BayesRapidBiExpMemoryControl_t      *memory,
                                             int                                  nbins)
{
	long i, j, k;

    for (i=low[1]; i<=high[1]; i++)
    {
        for (j=low[2]; j<=high[2]; j++)
        {
            for (k=low[3]; k<=high[3]; k++)
            {
                if (m[i][j][k]+low[4])
                    free(m[i][j][k]+low[4]);
            }

            if (m[i][j]+low[3])
                free(m[i][j]+low[3]);
        }

        if (m[i]+low[2])
            free(m[i]+low[2]);
    }

    if (m+low[1])
        free(m+low[1]);

    /* Now free the memory pools... */
    for (i=0; i<memory->npools; i++)
    {
        free(memory->pools[i].likelihoods);
        free_Bayes_dmatrix(memory->pools[i].loglikelihoods,0,memory->pools[i].nvalid,0,nbins);
    }

    free(memory->pools);
    free(memory);

    return (0);
}


int bayes_ComputeRapidBiExpDiscreteValuesNumOfValidEntries(int    *low,
                                                           int    *high,
                                                           int     ntaus,
                                                           double *taus,
                                                           int     nweights,
                                                           double *weights,
                                                           double  backgroundmin,
                                                           double  backgroundmax)
{
    int    nvalid, i, j, k, ell;
    double sum, w0;

    for (i=low[1],nvalid=0; i<=high[1]; i++)
    {
        for (j=low[2]; j<=high[2]; j++)
        {
            sum = weights[i]+weights[j];
            w0  = 1.0-sum;

            if ((sum>=0.0) && (sum<=1.0) && (w0>=/*backgroundmin*/weights[low[0]]) && (w0<=/*backgroundmax*/weights[high[0]]))
            {
                for (k=low[3]; k<=high[3]; k++)
                {
                    for (ell=low[4]; ell<=high[4]; ell++)
                    {
                        if (taus[k]>taus[ell])
                            nvalid++;
                    }
                }            
            }
        }
    }

    return (nvalid);
}

#define BAYES_RAPID_MEGABYTE 1048576 /* Bytes, http://en.wikipedia.org/wiki/Megabyte */

double bayes_ComputeRapidBiExpDiscreteValuesNumOfMegaBytesReqd(int npts,
                                                               int nptsvalid,
                                                               int nbins)
{
    int bytes;

    bytes  = npts*sizeof(BayesRapidBiExpDiscreteValues_t*);
    bytes += nptsvalid*sizeof(BayesRapidBiExpDiscreteValues_t);
    bytes += nptsvalid*nbins*sizeof(double);

    return ((double)(bytes)/(double)(BAYES_RAPID_MEGABYTE));
}


int bayes_PopulateRapidBiExpDiscreteValueMatrixElement( BayesRapidBiExpDiscreteValues_t *likelihood,
                                                        int                             *x,
                                                        int                              ntaus,
                                                        double                          *taus,
                                                        int                              nweights,
                                                        double                          *weights,
                                                        double                           backgroundmax,
                                                        int                              nbins,
                                                        double                          *binwalls,
                                                        BayesInstrRsp_t                 *instr,
                                                        double                           interval,
                                                        double                           modulationperiod,
                                                        BayesRapidLikelihoodValues_t    *fluorescencelikelihoods)
{
    int    ret;
    int    i, bin, ndim=4, ndecays=2;
    double bL, bH, bjoverT, norm, photonlikelihood, sum;

    sum = weights[x[1]]+weights[x[2]];

    if ((sum>=0.0) && (sum<=1.0) && ((1.0-sum)<=backgroundmax) && (taus[x[3]]>taus[x[4]]))
    {
        likelihood = (BayesRapidBiExpDiscreteValues_t*)malloc(sizeof(BayesRapidBiExpDiscreteValues_t));
        likelihood->logphotonlikelihoodgiventausandweights = Bayes_dvector(0,nbins-1);

        ret = bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(&norm,
                                                                        interval,modulationperiod,/*binwalls[fitstart]*/0.0,instr,
                                                                        ndecays,weights,taus);

        if (ret<0)
        {
            free_Bayes_dvector(likelihood->logphotonlikelihoodgiventausandweights,0,nbins-1);
            free(likelihood);
            likelihood = NULL;
            return (-1);
        }

        for (bin=0; bin<nbins; bin++)
        {
            bL      = binwalls[bin];
            bH      = binwalls[bin+1];
            bjoverT = (bH-bL)/interval;

            for (i=1, photonlikelihood=(1.0-sum)*bjoverT; i<=ndecays; i++)
                photonlikelihood += (weights[x[i]]*fluorescencelikelihoods[x[2+i]].fluorescencedecayphotonlikelihoodsgiventau[bin])/norm;

            likelihood->logphotonlikelihoodgiventausandweights[bin] = (float)(log(photonlikelihood));
        }

        likelihood->weights[1] = weights[x[1]];
        likelihood->weights[2] = weights[x[2]];
        likelihood->taus[1]    = taus[x[3]];
        likelihood->taus[2]    = taus[x[4]];
        likelihood->valid      = 1;
    }
    else
    {
        likelihood = NULL;
    }

    return (0);
}


int bayes_PopulateRapidBiExpDiscreteValueMatrix(BayesRapidBiExpDiscreteValues_t *****likelihoods,
                                                int                                 *low,
                                                int                                 *high,
                                                int                                  ntaus,
                                                double                              *taus,
                                                int                                  nweights,
                                                double                              *weights,
                                                double                               backgroundmin,
                                                double                               backgroundmax,
                                                int                                  nbins,
                                                double                              *binwalls,
                                                BayesInstrRsp_t                     *instr,
                                                double                               interval,
                                                double                               modulationperiod,
                                                BayesRapidLikelihoodValues_t        *fluorescencelikelihoods,
                                                int                                 *nvalid,
                                                int                                 *ninvalid)
{
    int    ret;
    int    x[5], *xmin, *xmax, i, j, k, ell;
    int    ndim=4, *weightindexes, *tauindexes, decay, ndecays=2, nstates;
    BayesRapidBiExpDiscreteValues_t *likelihood;
    
    double bL, bH, bjoverT, norm, photonlikelihood;
    int    bin;


    *nvalid   = 0;
    *ninvalid = 0;

    for (i=1,nstates=1; i<=4; i++)
        nstates *= (high[i]-low[i]+1);

    xmin = low;
    xmax = high;

    weightindexes = x;
    tauindexes    = &x[2];

    for (i=xmin[1]; i<=xmax[1]; i++)
    {
        x[1] = i;

        for (j=xmin[2]; j<=xmax[2]; j++)
        {
            x[2] = j;

            for (k=xmin[3]; k<=xmax[3]; k++)
            {
                x[3] = k;

                for (ell=xmin[4]; ell<=xmax[4]; ell++)
                {
                    x[4] = ell;

                    if (((weights[x[1]]+weights[x[2]]) >= 0.0) &&
                        ((weights[x[1]]+weights[x[2]]) <= 1.0) &&
                        (1.0-(weights[x[1]]+weights[x[2]]) >= /*backgroundmin*/weights[low[0]]) &&
                        (1.0-(weights[x[1]]+weights[x[2]]) <= /*backgroundmax*/weights[high[0]]) &&
                        (taus[x[3]] > taus[x[4]]))
                    {
                        (*nvalid)++;

                        likelihood = likelihoods[x[1]][x[2]][x[3]][x[4]];

                        likelihood->weights[0] = 1.0-weights[weightindexes[1]]-weights[weightindexes[2]];
                        likelihood->weights[1] = weights[/*x*/weightindexes[1]];
                        likelihood->weights[2] = weights[/*x*/weightindexes[2]];
                        likelihood->taus[1]    = taus[/*x*/tauindexes[1/*3*/]];
                        likelihood->taus[2]    = taus[/*x*/tauindexes[2/*4*/]];

                        ret = bayes_ComputeFluorescenceDecayPhotonNormalisationConstant(&norm,
                                                                                        interval,modulationperiod,/*binwalls[fitstart]*/0.0,instr,
                                                                                        ndecays,likelihood->weights,likelihood->taus);

                        if (ret<0)
                        {
                        }

                        for (bin=0; bin<nbins; bin++)
                        {
                            bL      = binwalls[bin];
                            bH      = binwalls[bin+1];
                            bjoverT = (bH-bL)/interval;

                            for (decay=1, photonlikelihood=likelihood->weights[0]*bjoverT; decay<=ndecays; decay++)
                                photonlikelihood += (likelihood->weights[decay]*fluorescencelikelihoods[tauindexes[decay]].fluorescencedecayphotonlikelihoodsgiventau[bin])/norm;

                            likelihood->logphotonlikelihoodgiventausandweights[bin] = (float)(log(photonlikelihood));
                        }

                        likelihood->valid = 1;
                    }
                    else
                    {
                        (*ninvalid)++;
                        likelihoods[x[1]][x[2]][x[3]][x[4]] = NULL;
                    }
                }
            }
        }
    }

    return (0);
}

/*
int CVICALLBACK bayes_PopulateRapidBiExpDiscreteValueMatrix_ThreadFunction (void *p)
{
    int ret;

    ret = bayes_PopulateRapidBiExpDiscreteValueMatrix(
              THRPTRBI->likelihoods,THRPTRBI->low,THRPTRBI->high,
              THRPTRBI->ntaus,THRPTRBI->taus,THRPTRBI->nweights,THRPTRBI->weights,THRPTRBI->backgroundmin,THRPTRBI->backgroundmax,
              THRPTRBI->nbins,THRPTRBI->binwalls,THRPTRBI->instr,THRPTRBI->interval,THRPTRBI->modulationperiod,
              THRPTRBI->fluorescencelikelihoods,
              THRPTRBI->nvalid,THRPTRBI->ninvalid);

	return (ret);
}



void bayes_LoadRapidBiExpDiscreteThreadData(RapidBiExpDiscreteThreadData        *data,
                                            BayesRapidBiExpDiscreteValues_t *****likelihoods,
                                            int                                 *low,
                                            int                                 *high,
                                            int                                  ntaus,
                                            double                              *taus,
                                            int                                  nweights,
                                            double                              *weights,
                                            double                               backgroundmin,
                                            double                               backgroundmax,
                                            int                                  nbins,
                                            double                              *binwalls,
                                            BayesInstrRsp_t                     *instr,
                                            double                               interval,
                                            double                               modulationperiod,
                                            BayesRapidLikelihoodValues_t        *fluorescencelikelihoods)
{
	data->likelihoods             = likelihoods;
	data->low                     = low;
	data->high                    = high;
	data->ntaus	                  = ntaus;
	data->taus                    = taus;
	data->nweights                = nweights;
	data->weights                 = weights;
    data->backgroundmin           = backgroundmin;
	data->backgroundmax	          = backgroundmax;
	data->nbins                   = nbins;
	data->binwalls                = binwalls;
	data->instr                   = instr;
	data->interval                = interval; 
	data->modulationperiod	      = modulationperiod; 
	data->fluorescencelikelihoods = fluorescencelikelihoods;
	*(data->nvalid)               = 0;
    *(data->ninvalid)             = 0;
}

int bayes_PopulateRapidBiExpDiscreteValueMatrixMultiThreaded(BayesRapidBiExpDiscreteValues_t *****likelihoods,
                                                             int                                 *low,
                                                             int                                 *high,
                                                             int                                  ntaus,
                                                             double                              *taus,
                                                             int                                  nweights,
                                                             double                              *weights,
                                                             double                               backgroundmin,
                                                             double                               backgroundmax,
                                                             int                                  nbins,
                                                             double                              *binwalls,
                                                             BayesInstrRsp_t                     *instr,
                                                             double                               interval,
                                                             double                               modulationperiod,
                                                             BayesRapidLikelihoodValues_t        *fluorescencelikelihoods,
                                                             int                                  nstates,
                                                             int                                 *nvalid,
                                                             int                                 *ninvalid)
{
	int complete[BAYES_RAPID_NTHREADS], nComplete, ThreadPoolHandle, ThreadID[BAYES_RAPID_NTHREADS], progress, i, j, states;
    RapidBiExpDiscreteThreadData     data[BAYES_RAPID_NTHREADS];
    CmtThreadFunctionExecutionStatus status;

    int ret=BAYES_RAPID_RESULT_NO_ERROR, threadlow[BAYES_RAPID_NTHREADS][5], threadhigh[BAYES_RAPID_NTHREADS][5], low1, high1, diff;

    diff = (high[1]-low[1]+1)/BAYES_RAPID_NTHREADS;

    for (i=0; i<BAYES_RAPID_NTHREADS; i++)
    {
        for (j=0; j<=4; j++)
        {
            threadlow[i][j]  = low[j];
            threadhigh[i][j] = high[j];
        }
    }

    CmtNewThreadPool(BAYES_RAPID_NTHREADS, &ThreadPoolHandle);

    // start nThreads-1 threads
    for (i=0,low1=low[1],high1=low[1]+diff-1; i<BAYES_RAPID_NTHREADS-1; i++)
    {
        threadlow[i][1]  = low1;
        threadhigh[i][1] = high1;

        bayes_LoadRapidBiExpDiscreteThreadData(&data[i],
                                               likelihoods,threadlow[i],threadhigh[i],
                                               ntaus,taus,nweights,weights,backgroundmin,backgroundmax,
                                               nbins,binwalls,instr,interval,modulationperiod,
                                               fluorescencelikelihoods);

        low1  = high1+1;
        high1 = low1+diff-1;

        CmtScheduleThreadPoolFunction(ThreadPoolHandle,bayes_PopulateRapidBiExpDiscreteValueMatrix_ThreadFunction,(void*)&data[i],&ThreadID[i]);
        	
        complete[i]=0;
    }

    threadlow[i][1]   = low1;
    threadhigh[i][1]  = high[1];

    // last thread gets the remainder of the grid
    bayes_LoadRapidBiExpDiscreteThreadData(&data[i],
                                           likelihoods,threadlow[i],threadhigh[i],
                                           ntaus,taus,nweights,weights,backgroundmin,backgroundmax,
                                           nbins,binwalls,instr,interval,modulationperiod,
                                           fluorescencelikelihoods);
    
    CmtScheduleThreadPoolFunction(ThreadPoolHandle,bayes_PopulateRapidBiExpDiscreteValueMatrix_ThreadFunction,(void*)&data[i],&ThreadID[i]);
    complete[i]=0;

    nComplete=0;

    while (nComplete < BAYES_RAPID_NTHREADS)
    {
        // this loop queries each thread to see if it has finished
        // Add up all the progress values from the threads to report an overall progress against the number of rows
        // Checks for an interupt, each thread does this also and will exit if it is detected.

        for (i=0; i<BAYES_RAPID_NTHREADS; i++)
        {
            if (complete[i]==0)
            {
                CmtGetThreadPoolFunctionAttribute(ThreadPoolHandle,ThreadID[i],ATTR_TP_FUNCTION_EXECUTION_STATUS,&status);
                
                if (status == kCmtThreadFunctionComplete)
                {
                    complete[i]=1;
                    nComplete++;
                    CmtReleaseThreadPoolFunctionID(ThreadPoolHandle,ThreadID[i]);
                }
            }
        }

        Delay(0.5);

        for (i=0,states=0; i<BAYES_RAPID_NTHREADS; i++)
        {
            states += (*(data[i].nvalid))+(*(data[i].ninvalid));
        }

    }

    for (i=0,*nvalid=0,*ninvalid=0; i<BAYES_RAPID_NTHREADS; i++)
    {
        *nvalid   += *data[i].nvalid;
        *ninvalid += *data[i].ninvalid;
    }

    CmtDiscardThreadPool(ThreadPoolHandle);

    return (ret);
}
*/

int bayes_AllocateForRapidValueStore(BayesRapidValueStore_t   *store,
                                     int                       updatetype,
                                     /* Mono-exponential grid settings... */
                                     int                       ntaus,
                                     int                       nweights,
                                     /* Bi-exponential grid settings... */
                                     int                       ntaus_bi,
                                     double                   *taus_bi,
                                     int                       nweights_bi,
                                     double                   *weights_bi,
                                     int                      *low,
                                     int                      *high,
                                     /* General settings... */
                                     int                       nbins)
{
    int i, j;

    /* Psuedo-rapid fluorescence decay likelihood given tau vector  */

    /*==================================================================================*/
    /*                             Mono-exponential value store                         */
    /*==================================================================================*/
    if ((updatetype==BAYES_RAPID_GRID_MONO) ||
        (updatetype==BAYES_RAPID_GRID_MONO_AND_BI))
    {
        store->monoexpvaluestore = (BayesRapidMonoExpValueStore_t *)malloc(sizeof(BayesRapidMonoExpValueStore_t));

        /* Settings */
        store->monoexpvaluestore->settings         = (BayesRapidGridSettings_t *)malloc(sizeof(BayesRapidGridSettings_t));
        store->monoexpvaluestore->settings->tau    = Bayes_dvector(0,ntaus-1);
        store->monoexpvaluestore->settings->weight = Bayes_dvector(0,nweights-1);

        /* Fluorescence likelihoods... */
        store->monoexpvaluestore->fluorescencelikelihoods = bayes_AllocateRapidLikelihoodValuesVector(ntaus);

        for (j=0; j<ntaus; j++)
        {
            store->monoexpvaluestore->fluorescencelikelihoods[j].fluorescencedecayphotonlikelihoodsgiventau
                = (double*)malloc(sizeof(double)*(nbins));
        }

        /* Likelihoods... */
        store->monoexpvaluestore->likelihoods = bayes_AllocateRapidMonoExpDiscreteValuesMatrix(0,nweights-1,0,ntaus-1);

        for (i=0; i<nweights; i++)
        {
            for (j=0; j<ntaus; j++)
            {
                store->monoexpvaluestore->likelihoods[i][j].logphotonlikelihoodgiventauandw0
                    = (double *)malloc(sizeof(double)*(nbins));
            }
        }    
    }

    /*==================================================================================*/
    /*                              Bi-exponential value store                          */
    /*==================================================================================*/
    if ((updatetype==BAYES_RAPID_GRID_BI) ||
        (updatetype==BAYES_RAPID_GRID_MONO_AND_BI))
    {
        /* Dedicated bi-exponential log bin likelihood given taus and weights matrix... */
        if ((ntaus_bi>0) && (taus_bi) && (nweights_bi>0) && (weights_bi) && (low) && (high))
        {
            store->biexpvaluestore                          = (BayesRapidBiExpValueStore_t *)malloc(sizeof(BayesRapidBiExpValueStore_t));
            store->biexpvaluestore->settings                = (BayesRapidGridSettings_t *)malloc(sizeof(BayesRapidGridSettings_t));
            store->biexpvaluestore->settings->tau           = Bayes_dvector(0,ntaus_bi-1);
            store->biexpvaluestore->settings->weight        = Bayes_dvector(0,nweights_bi-1);
            store->biexpvaluestore->low                     = Bayes_ivector(0,4);
            store->biexpvaluestore->high                    = Bayes_ivector(0,4);
            store->biexpvaluestore->nstates                 = 0;
            store->biexpvaluestore->fluorescencelikelihoods = bayes_AllocateRapidLikelihoodValuesVector(ntaus_bi);

            //using the number of 'w1' entries for memory management...
            store->biexpvaluestore->memory                  = (BayesRapidBiExpMemoryControl_t *)malloc(sizeof(BayesRapidBiExpMemoryControl_t));
            
            for (j=0; j<ntaus_bi; j++)
            {
                store->biexpvaluestore->fluorescencelikelihoods[j].fluorescencedecayphotonlikelihoodsgiventau
                                                            = (double*)malloc(sizeof(double)*(nbins));
            }
#if 0
            store->biexpvaluestore->likelihoods             = bayes_AllocateRapidBiExpDiscreteValuesMatrix(low,high,&store->biexpvaluestore->nstates);
#else
            store->biexpvaluestore->likelihoods = bayes_AllocateRapidBiExpDiscreteValuesMatrix(low,high,&store->biexpvaluestore->nstates,
                                                                                               ntaus_bi,taus_bi,nweights_bi,weights_bi,nbins,
                                                                                               store->biexpvaluestore->memory);

#endif    
        }
    }

    return (BAYES_RAPID_RESULT_NO_ERROR);
}


int bayes_FreeForRapidValueStore(BayesRapidValueStore_t *store, int gridtype)
{
    int nweights, ntaus, i, j;

    if (store)
    {
        if ((gridtype==BAYES_RAPID_GRID_MONO) ||
            (gridtype==BAYES_RAPID_GRID_MONO_AND_BI))
        {
            /* Free the dedicated mono-exponential value store... */
            if (store->monoexpvaluestore)
            {
                ntaus    = store->monoexpvaluestore->settings->ntaus;
                nweights = store->monoexpvaluestore->settings->nweights;
                //nweights = store->monoexpvaluestore->nw0s;
                
                /* Free the likelihoods vector... */
                if (store->monoexpvaluestore->fluorescencelikelihoods)
                {
                    for (i=0; i<ntaus; i++)
                    {
                        if (store->monoexpvaluestore->fluorescencelikelihoods[i].fluorescencedecayphotonlikelihoodsgiventau)
                        {
                            free(store->monoexpvaluestore->fluorescencelikelihoods[i].fluorescencedecayphotonlikelihoodsgiventau);
                            store->monoexpvaluestore->fluorescencelikelihoods[i].fluorescencedecayphotonlikelihoodsgiventau = NULL;
                        }
                    }

                    bayes_FreeRapidLikelihoodValuesVector(store->monoexpvaluestore->fluorescencelikelihoods);
                    store->monoexpvaluestore->fluorescencelikelihoods = NULL;
                }

                /* Free the likelihoods matrix... */
                for (j=0; j<ntaus; j++)
                {
                    for (i=0; i<nweights; i++)
                    {
                        if (store->monoexpvaluestore->likelihoods[i][j].logphotonlikelihoodgiventauandw0)
                        {
                            free(store->monoexpvaluestore->likelihoods[i][j].logphotonlikelihoodgiventauandw0);
                            store->monoexpvaluestore->likelihoods[i][j].logphotonlikelihoodgiventauandw0 = NULL;
                        }
                    }
                }

                bayes_FreeRapidMonoExpDiscreteValuesMatrix(store->monoexpvaluestore->likelihoods,0,nweights-1,0,ntaus-1);
                free(store->monoexpvaluestore);
                store->monoexpvaluestore = NULL;
            }        
        }

        if ((gridtype==BAYES_RAPID_GRID_BI) ||
            (gridtype==BAYES_RAPID_GRID_MONO_AND_BI))
        {
            /* Free dedicated bi-exponential store... */
            if (store->biexpvaluestore)
            {
                bayes_FreeRapidBiExpDiscreteValuesMatrix(store->biexpvaluestore->likelihoods,
                                                         store->biexpvaluestore->low,
                                                         store->biexpvaluestore->high,
                                                         store->biexpvaluestore->memory,
                                                         store->biexpvaluestore->nbins);

                free_Bayes_ivector(store->biexpvaluestore->low,0,4);
                free_Bayes_ivector(store->biexpvaluestore->high,0,4);

                ntaus = store->biexpvaluestore->settings->ntaus;

                if (store->biexpvaluestore->fluorescencelikelihoods)
                {
                    for (j=0; j<ntaus; j++)
                    {
                        if (store->biexpvaluestore->fluorescencelikelihoods[j].fluorescencedecayphotonlikelihoodsgiventau)
                            free(store->biexpvaluestore->fluorescencelikelihoods[j].fluorescencedecayphotonlikelihoodsgiventau);
                    }
                }

                free(store->biexpvaluestore->fluorescencelikelihoods);

                free_Bayes_dvector(store->biexpvaluestore->settings->tau,0,store->biexpvaluestore->settings->ntaus);
                free_Bayes_dvector(store->biexpvaluestore->settings->weight,0,store->biexpvaluestore->settings->nweights);
                free(store->biexpvaluestore->settings);
                store->biexpvaluestore->settings = NULL;

                free(store->biexpvaluestore);
                store->biexpvaluestore = NULL;
            }        
        }
    }

    return (0);
}


/*=======================================================*/
/*     FLUORESCENCE LIKELIHOOD POPULATION FUNCTIONS      */
/*=======================================================*/

void bayes_LoadRapidLikelihoodsDiscreteThreadData(RapidLikelihoodsDiscreteThreadData *data,
                                                  BayesRapidLikelihoodValues_t       *fluorescencelikelihoods,
                                                  int                                *low,
                                                  int                                *high,
                                                  int                                 ntaus,
                                                  double                             *taus,
                                                  int                                 nbins,
                                                  double                             *binwalls,
                                                  BayesInstrRsp_t                    *instr,
                                                  double                              interval,
                                                  double                              modulationperiod)
{
	data->fluorescencelikelihoods = fluorescencelikelihoods;
    data->low                     = low;
    data->high                    = high;
	data->ntaus	                  = ntaus;
	data->taus                    = taus;
	data->nbins                   = nbins;
	data->binwalls                = binwalls;
	data->instr                   = instr;
	data->interval                = interval; 
	data->modulationperiod	      = modulationperiod; 
	*(data->nvalid)               = 0;
    *(data->ninvalid)             = 0;
}


int bayes_PopulateRapidLikelihoodsDiscreteValueVector(BayesRapidLikelihoodValues_t *fluorescencelikelihoods,
                                                      int                          *low,
                                                      int                          *high,
                                                      int                           ntaus,
                                                      double                       *taus,
                                                      int                           nbins,
                                                      double                       *binwalls,
                                                      BayesInstrRsp_t              *instr,
                                                      double                        interval,
                                                      double                        modulationperiod,
                                                      int                          *nvalid,
                                                      int                          *ninvalid)
{
    int j, ret;

    /* Fluorescence likelihoods given tau... */
    for (j=low[1]; j<=high[1]; j++)
    {
        ret = bayes_ComputeFluorescenceDecayPhotonBinLikelihoodsGivenTau(fluorescencelikelihoods[j].fluorescencedecayphotonlikelihoodsgiventau,
                                                                         nbins,binwalls,NULL,
                                                                         interval,modulationperiod,instr,
                                                                         taus[j],
                                                                         0,NULL,NULL);//not normalised

        if (ret>=0)
        {
            /* Populate likelihood container... */
            fluorescencelikelihoods[j].tau   = taus[j];
            fluorescencelikelihoods[j].valid = 1;

            (*nvalid)++;
        }
        else /* Likelihoods have not been computed succesfully...*/
        {
            fluorescencelikelihoods[j].tau   = -1.0;
            fluorescencelikelihoods[j].valid = 0; 

            (*ninvalid)++;
        }
    }

    return (0);
}


/*

	THIS MULTITHREADED VERSION IS BUILT ON THE LabWindows/CVI THREADING FUNCTIONS.
	Could be re-written to be more portable.

int CVICALLBACK bayes_PopulateRapidLikelihoodsDiscreteValueVector_ThreadFunction (void *p)
{
    int ret;

    ret = bayes_PopulateRapidLikelihoodsDiscreteValueVector(
              THRPTRLIKELIHOOD->fluorescencelikelihoods,THRPTRLIKELIHOOD->low,THRPTRLIKELIHOOD->high,
              THRPTRLIKELIHOOD->ntaus,THRPTRLIKELIHOOD->taus,
              THRPTRLIKELIHOOD->nbins,THRPTRLIKELIHOOD->binwalls,THRPTRLIKELIHOOD->instr,THRPTRLIKELIHOOD->interval,THRPTRLIKELIHOOD->modulationperiod,
              THRPTRLIKELIHOOD->nvalid,THRPTRLIKELIHOOD->ninvalid);

	return (ret);
}

int bayes_PopulateRapidLikelihoodsDiscreteValueVectorMultiThreaded(BayesRapidLikelihoodValues_t *fluorescencelikelihoods,
                                                                   int                           ntaus,
                                                                   double                       *taus,
                                                                   int                           nbins,
                                                                   double                       *binwalls,
                                                                   BayesInstrRsp_t              *instr,
                                                                   double                        interval,
                                                                   double                        modulationperiod,
                                                                   int                           nstates,
                                                                   int                          *nvalid,
                                                                   int                          *ninvalid)
{
	int complete[BAYES_RAPID_NTHREADS], nComplete, ThreadPoolHandle, ThreadID[BAYES_RAPID_NTHREADS], i, states;
    RapidLikelihoodsDiscreteThreadData   data[BAYES_RAPID_NTHREADS];
    CmtThreadFunctionExecutionStatus status;

    int low[2], high[2];

    int ret=BAYES_RAPID_RESULT_NO_ERROR, threadlow[BAYES_RAPID_NTHREADS][2], threadhigh[BAYES_RAPID_NTHREADS][2], low1, high1, diff;

    low[1]  = 0;
    high[1] = ntaus-1;
    diff = (high[1]-low[1]+1)/BAYES_RAPID_NTHREADS;

    for (i=0; i<BAYES_RAPID_NTHREADS; i++)
    {
            threadlow[i][1]  = low[1];
            threadhigh[i][1] = high[1];
    }

    CmtNewThreadPool(BAYES_RAPID_NTHREADS, &ThreadPoolHandle);

    // start nThreads-1 threads
    for (i=0,low1=low[1],high1=low[1]+diff-1; i<BAYES_RAPID_NTHREADS-1; i++)
    {
        threadlow[i][1]  = low1;
        threadhigh[i][1] = high1;

        bayes_LoadRapidLikelihoodsDiscreteThreadData(&data[i],
                                                     fluorescencelikelihoods,threadlow[i],threadhigh[i],
                                                     ntaus,taus,
                                                     nbins,binwalls,instr,interval,modulationperiod);

        low1  = high1+1;
        high1 = low1+diff-1;

        CmtScheduleThreadPoolFunction(ThreadPoolHandle,bayes_PopulateRapidLikelihoodsDiscreteValueVector_ThreadFunction,(void*)&data[i],&ThreadID[i]);
        	
        complete[i]=0;
    }

    threadlow[i][1]   = low1;
    threadhigh[i][1]  = high[1];

    // last thread gets the remainder of the grid
    bayes_LoadRapidLikelihoodsDiscreteThreadData(&data[i],
                                                 fluorescencelikelihoods,threadlow[i],threadhigh[i],
                                                 ntaus,taus,
                                                 nbins,binwalls,instr,interval,modulationperiod);
    
    CmtScheduleThreadPoolFunction(ThreadPoolHandle,bayes_PopulateRapidLikelihoodsDiscreteValueVector_ThreadFunction,(void*)&data[i],&ThreadID[i]);
    complete[i]=0;

    nComplete=0;

    while (nComplete < BAYES_RAPID_NTHREADS)
    {
        // this loop queries each thread to see if it has finished
        // Add up all the progress values from the threads to report an overall progress against the number of rows
        // Checks for an interupt, each thread does this also and will exit if it is detected.

        for (i=0; i<BAYES_RAPID_NTHREADS; i++)
        {
            if (complete[i]==0)
            {
                CmtGetThreadPoolFunctionAttribute(ThreadPoolHandle,ThreadID[i],ATTR_TP_FUNCTION_EXECUTION_STATUS,&status);
                
                if (status == kCmtThreadFunctionComplete)
                {
                    complete[i]=1;
                    nComplete++;
                    CmtReleaseThreadPoolFunctionID(ThreadPoolHandle,ThreadID[i]);
                }
            }
        }

        Delay(0.05);

        for (i=0,states=0; i<BAYES_RAPID_NTHREADS; i++)
        {
            states += (*(data[i].nvalid))+(*(data[i].ninvalid));
        }

    }

    for (i=0,*nvalid=0,*ninvalid=0; i<BAYES_RAPID_NTHREADS; i++)
    {
        *nvalid   += *data[i].nvalid;
        *ninvalid += *data[i].ninvalid;
    }

    CmtDiscardThreadPool(ThreadPoolHandle);

    return (ret);
}
*/

/*=======================================================*/
/*         GRID MANAGEMENT/POPULATION FUNCTIONS          */
/*=======================================================*/

void bayes_InvalidateRapidValueStore(BayesRapidValueStore_t *store)
{
    store->validlikelihoodsgrid = 0;
    store->validmonoexpgrid     = 0;
    store->validbiexpgrid       = 0;
}


void bayes_InvalidateRapidBiExpValueStore(BayesRapidValueStore_t *store)
{
    if (store->biexpvaluestore)
        store->biexpvaluestore->validgrid = 0;
}


int bayes_CheckForValidRapidValueStore(BayesRapidValueStore_t *store, int gridtype)
{
    if ((gridtype==BAYES_RAPID_GRID_MONO) ||
        (gridtype==BAYES_RAPID_GRID_MONO_AND_BI))
    {
        if ((0==store->validmonoexpgrid) ||
            (!store->monoexpvaluestore))
            return (0);
    }
    
    if ((gridtype==BAYES_RAPID_GRID_BI) ||
        (gridtype==BAYES_RAPID_GRID_MONO_AND_BI))
    {
        if ((0==store->validbiexpgrid) ||
            (!store->biexpvaluestore))
            return (0);
    }

    return (1);
}


int bayes_InitializeRapidValueStore(BayesRapidValueStore_t *store)
{
    store->validlikelihoodsgrid  = 0;
    store->likelihoodsvaluestore = NULL;
    store->validmonoexpgrid      = 0;
    store->monoexpvaluestore     = NULL;
    store->validbiexpgrid        = 0;
    store->biexpvaluestore       = NULL;

    return (0);
}

#if 1
int bayes_CheckForDifferentTauVectorValues(int ntaus, double *stored, double *userinterface)
{
    int i;

    for (i=0; i<ntaus; i++)
    {
        if (stored[i] != userinterface[i])
        {
            return (1);
        }
    }

    return (0);
}
int bayes_CheckForDifferentWeightVectorValues(int nweights, double *stored, double *userinterface)
{
    int i;

    for (i=0; i<nweights; i++)
    {
        if (stored[i] != userinterface[i])
        {
            return (1);
        }
    }

    return (0);
}
int bayes_CheckForDifferentLowIndexVectorValues(int nlow, int *stored, int *userinterface)
{
    int i;

    for (i=0; i<nlow; i++)
    {
        if (stored[i] != userinterface[i])
        {
            return (1);
        }
    }

    return (0);
}
int bayes_CheckForDifferentHighIndexVectorValues(int nhigh, int *stored, int *userinterface)
{
    int i;

    for (i=0; i<nhigh; i++)
    {
        if (stored[i] != userinterface[i])
        {
            return (1);
        }
    }

    return (0);
}

int bayes_DetermineIfBayesGridUpdateReqd(BayesRapidValueStore_t *store,
                                         int                     gridtype,
                                         /* Mono-exp... */
                                         int                     ntaus,
                                         double                 *taus,
                                         int                     nweights,
                                         double                 *weights,
                                         double                  backgroundmin,
                                         double                  backgroundmax,
                                         /* Bi-exp... */
                                         int                     ntausbi,
                                         double                 *tausbi,
                                         int                     nweightsbi,
                                         double                 *weightsbi,
                                         double                  backgroundminbi,
                                         double                  backgroundmaxbi,
                                         int                    *lowbi,
                                         int                    *highbi,
                                         /* Instrumentation and data... */
                                         int                     nbins,
                                         BayesInstrRsp_t        *instr,
                                         double                  interval,
                                         double                  modulationperiod)
{
    int updatemono=0, updatebi=0;

    /* NOTE: The 'store->validgrid' flag is set invalid when an event    */
    /*       occurs (e.g. configuration panel opened, parameter value    */
    /*       changed, new image opened, etc) that MAY invalidate grid... */

    /* Check for a valid mono-exponential grid... */
    if ((gridtype==BAYES_RAPID_GRID_MONO) ||
        (gridtype==BAYES_RAPID_GRID_MONO_AND_BI))
    {
        if (!store->monoexpvaluestore) /* Check whether a mono-exponential grid even exists... */
        {
            updatemono = 1;
        }
        else /* General settings (instrumentation and rebinning) */
        if ((nbins            != store->monoexpvaluestore->nbins) ||
            (interval         != store->monoexpvaluestore->interval) ||
            (modulationperiod != store->monoexpvaluestore->modulationperiod) ||
            (bayes_CheckForDifferentInstrRspConfigParams(&store->monoexpvaluestore->instr,instr)))
        {
            updatemono = 1;
        }
        else /* Mono-exp grid... */
        if ((ntaus         != store->monoexpvaluestore->settings->ntaus) ||
            (nweights      != store->monoexpvaluestore->settings->nweights) ||
            (backgroundmin != store->monoexpvaluestore->settings->backgroundmin) ||
            (backgroundmax != store->monoexpvaluestore->settings->backgroundmax) ||
            (bayes_CheckForDifferentTauVectorValues(ntaus,store->monoexpvaluestore->settings->tau,taus)) ||
            (bayes_CheckForDifferentWeightVectorValues(nweights,store->monoexpvaluestore->settings->weight,weights)))
        {
            updatemono = 1;
        }
        else
        {
            /* No changes detected... */
            store->validmonoexpgrid = 1;
        }
    }

    /* Check for a valid bi-exponential grid... */
    if ((gridtype==BAYES_RAPID_GRID_BI) ||
        (gridtype==BAYES_RAPID_GRID_MONO_AND_BI))
    {
        /* Check whether a bi-exponential grid even exists... */
        if (!store->biexpvaluestore)
        {
            updatebi = 1; /* Bi-exp grid doesn't exist... */
        }
        else /* General settings (instrumentation and rebinning) */
        if ((nbins            != store->biexpvaluestore->nbins) ||
            (interval         != store->biexpvaluestore->interval) ||
            (modulationperiod != store->biexpvaluestore->modulationperiod) ||
            (bayes_CheckForDifferentInstrRspConfigParams(&store->biexpvaluestore->instr,instr)))
        {
            updatebi = 1;
        }
        else /* Bi-exp grid... */
        if ((ntausbi         != store->biexpvaluestore->settings->ntaus) ||
            (nweightsbi      != store->biexpvaluestore->settings->nweights) ||
            (backgroundminbi != store->biexpvaluestore->settings->backgroundmin) ||
            (backgroundmaxbi != store->biexpvaluestore->settings->backgroundmax) ||
            (bayes_CheckForDifferentTauVectorValues(ntausbi,store->biexpvaluestore->settings->tau,tausbi)) ||
            (bayes_CheckForDifferentWeightVectorValues(nweightsbi,store->biexpvaluestore->settings->weight,weightsbi)) ||
            (bayes_CheckForDifferentHighIndexVectorValues(4,store->biexpvaluestore->high,highbi)) ||
            (bayes_CheckForDifferentHighIndexVectorValues(4,store->biexpvaluestore->low,lowbi)))
        {
            updatebi = 1;
        }
        else
        {
            /* No changes detected... */
            store->validbiexpgrid = 1;                    
        }
    }

    if (gridtype==BAYES_RAPID_GRID_MONO)
    {
        if (updatemono)
            return (BAYES_RAPID_GRID_MONO);
        else
            return (0);
    }
    else if (gridtype==BAYES_RAPID_GRID_BI)
    {
        if (updatebi)
            return (BAYES_RAPID_GRID_BI);
        else
            return (0);
    }
    else if (gridtype==BAYES_RAPID_GRID_MONO_AND_BI)
    {
        if ((updatemono) && (updatebi))
            return (BAYES_RAPID_GRID_MONO_AND_BI);
        else if (updatemono)
            return (BAYES_RAPID_GRID_MONO);
        else if (updatebi)
            return (BAYES_RAPID_GRID_BI);
        else
            return (0);
    }
    else
    {
        return (BAYES_RAPID_GRID_MONO_AND_BI); //default case, update everything...
    }
}
#else
int bayes_DetermineIfBayesGridUpdateReqd(BayesRapidValueStore_t *store,
                                         /* Mono-exp... */
                                         int                     ntaus,
                                         int                     nweights,
                                         double                  backgroundmax,
                                         /* Bi-exp... */
                                         int                     ntausbi,
                                         double                  taulowbi,
                                         double                  tauhighbi,
                                         int                     nweightsbi,
                                         double                  weightlowbi,
                                         double                  weighthighbi,
                                         double                  backgroundminbi,
                                         double                  backgroundmaxbi,
                                         double                  w0lowbi,
                                         double                  w0highbi,
                                         double                  w1lowbi,
                                         double                  w1highbi,
                                         double                  w2lowbi,
                                         double                  w2highbi,
                                         double                  tau1lowbi,
                                         double                  tau1highbi,
                                         double                  tau2lowbi,
                                         double                  tau2highbi,
                                         /* Instrumentation and data... */
                                         int                     nbins,
                                         BayesInstrRsp_t        *instr,
                                         double                  interval,
                                         double                  modulationperiod)
{
    /* Flags first...                                                                 */
    /* NOTE: The flags are set invalid when an event occurs (e.g. configuration panel */
    /* opened, parameter value changed) that MAY invalidate grid...                   */

	if (!store->validgrid)
    {
        return (1);
    }

        /* General settings (instrumentation and rebinning) */
        if ((nbins            != store->nbins) ||
            (interval         != store->interval) ||
            (modulationperiod != store->modulationperiod) ||
            (bayes_CheckForDifferentInstrRspConfigParams(&store->instr,instr)))
        {
            return (1);
        }

        /* Mono-exp grid... */
        if ((ntaus            != store->settings->ntaus) ||
            (nweights         != store->settings->nweights) ||
            (backgroundmax    != store->settings->backgroundmax))
        {
            return (1);
        }

        /* Bi-exp grid... */
        if (!store->biexpvaluestore)
        {
            return (1);
        }

        if ((ntausbi         != store->biexpvaluestore->settings->ntaus) ||
            (taulowbi        != store->biexpvaluestore->settings->tau[0]) ||
            (tauhighbi       != store->biexpvaluestore->settings->tau[ntausbi-1]) ||
            (nweightsbi      != store->biexpvaluestore->settings->nweights) ||
            (weightlowbi     != store->biexpvaluestore->settings->weight[0]) ||
            (weighthighbi    != store->biexpvaluestore->settings->weight[nweightsbi-1]) ||
            (w0lowbi         != store->biexpvaluestore->settings->weight[store->biexpvaluestore->low[0]]) ||
            (w0highbi        != store->biexpvaluestore->settings->weight[store->biexpvaluestore->high[0]]) ||
            (w1lowbi         != store->biexpvaluestore->settings->weight[store->biexpvaluestore->low[1]]) ||
            (w1highbi        != store->biexpvaluestore->settings->weight[store->biexpvaluestore->high[1]]) ||
            (w2lowbi         != store->biexpvaluestore->settings->weight[store->biexpvaluestore->low[2]]) ||
            (w2highbi        != store->biexpvaluestore->settings->weight[store->biexpvaluestore->high[2]]) ||
            (tau1lowbi       != store->biexpvaluestore->settings->weight[store->biexpvaluestore->low[3]]) ||
            (tau1highbi      != store->biexpvaluestore->settings->weight[store->biexpvaluestore->high[3]]) ||
            (tau2highbi      != store->biexpvaluestore->settings->weight[store->biexpvaluestore->low[4]]) ||
            (tau2highbi      != store->biexpvaluestore->settings->weight[store->biexpvaluestore->high[4]]) ||
            (backgroundminbi != store->biexpvaluestore->settings->backgroundmin) ||
            (backgroundmaxbi != store->biexpvaluestore->settings->backgroundmax))
        {
            return (1);
        }


    /* No changes detected... */
    //store->validgrid = 1;

    return (0);
}
#endif


int bayes_DestroyRapidValueStore(BayesRapidValueStore_t *store, int gridtype)
{
    bayes_FreeForRapidValueStore(store,gridtype);
    //bayes_InitializeRapidValueStore(gridtype);

    return (BAYES_RAPID_RESULT_NO_ERROR);
}


int bayes_PopulateRapidValueStore(BayesRapidValueStore_t   *store,
                                  int                       updatetype,
                                  /* Mono-exponential grid settings... */
                                  int                       ntaus,
                                  double                   *taus,
                                  int                       nweights,
                                  double                   *weights,
                                  double                    backgroundmin,
                                  double                    backgroundmax,
                                  /* Bi-exponential grid settings... */
                                  int                       ntaus_bi,
                                  double                   *taus_bi,
                                  int                       nweights_bi,
                                  double                   *weights_bi,
                                  double                    backgroundmin_bi,
                                  double                    backgroundmax_bi,
                                  int                      *low,
                                  int                      *high,
                                  /* General (common) settings... */
                                  int                       nbins,
                                  double                   *binwalls,
                                  BayesInstrRsp_t          *instr,
                                  double                    interval,
                                  double                    modulationperiod)
{
    int     i, ret;

    int     nw0s, npts, nvalid, ninvalid;
    int     lowmono[3], highmono[3];
    double  w0low, w0high;
    FILE   *fp=NULL;
	int single_low[2] = { 0, 0 };
	int single_high[2] = { 0, ntaus - 1 };


    /* Mono-exponential grid fluorescence likelihoods population... */
    if ((updatetype==BAYES_RAPID_GRID_MONO) ||
        (updatetype==BAYES_RAPID_GRID_MONO_AND_BI))
    {
#if 0 
		// Multi-threaded option:
        bayes_PopulateRapidLikelihoodsDiscreteValueVectorMultiThreaded(store->monoexpvaluestore->fluorescencelikelihoods,
                                                                       ntaus,taus,
                                                                       nbins,binwalls,instr,interval,modulationperiod,
                                                                       ntaus,&nvalid,&ninvalid);
#else
		// Single Thread:
		bayes_PopulateRapidLikelihoodsDiscreteValueVector(store->monoexpvaluestore->fluorescencelikelihoods,
															single_low, single_high,
															ntaus, taus,
															nbins, binwalls, instr, interval, modulationperiod,
															&nvalid, &ninvalid);
#endif


        store->monoexpvaluestore->nlikelihoods        = ntaus;
        store->monoexpvaluestore->nvalidlikelihoods   = nvalid;
        store->monoexpvaluestore->ninvalidlikelihoods = ninvalid;

        /* Dedicated mono-exponential settings... */
        w0low  = weights[0];
        w0high = backgroundmax;
        nw0s   = nweights;

        store->monoexpvaluestore->nstates = nw0s*ntaus;

        lowmono[1]  = 0;
        lowmono[2]  = 0;
        highmono[1] = nw0s-1;
        highmono[2] = ntaus-1;

#if 0
		// Multi-threaded option:
		ret = bayes_PopulateRapidMonoExpDiscreteValueMatrixMultiThreaded(store->monoexpvaluestore->likelihoods,lowmono,highmono,
                                                                         ntaus,taus,nweights,weights,backgroundmax,
                                                                         nbins,binwalls,instr,interval,modulationperiod,
                                                                         store->monoexpvaluestore->fluorescencelikelihoods,
                                                                         store->monoexpvaluestore->nstates,&nvalid,&ninvalid);
#else
		// Single Threaded:
		ret = bayes_PopulateRapidMonoExpDiscreteValueMatrix(store->monoexpvaluestore->likelihoods, lowmono, highmono,
															ntaus, taus, nweights, weights, backgroundmax,
															nbins, binwalls, instr, interval, modulationperiod,
															store->monoexpvaluestore->fluorescencelikelihoods,
															&nvalid, &ninvalid);
#endif

        store->monoexpvaluestore->settings->ntaus         = ntaus;
        
        for (i=0; i<ntaus; i++)
            store->monoexpvaluestore->settings->tau[i]    = taus[i];

        store->monoexpvaluestore->settings->nweights         = nweights;

        for (i=0; i<nweights; i++)
            store->monoexpvaluestore->settings->weight[i] = weights[i];

        store->monoexpvaluestore->settings->backgroundmin = backgroundmin;
        store->monoexpvaluestore->settings->backgroundmax = backgroundmax;
        

        store->monoexpvaluestore->nbins            = nbins;
        store->monoexpvaluestore->interval         = interval;
        store->monoexpvaluestore->modulationperiod = modulationperiod;
        bayes_CopyInstrRspConfigParams(instr,&store->monoexpvaluestore->instr);
        store->monoexpvaluestore->nvalidstates     = nvalid;
        store->monoexpvaluestore->ninvalidstates   = ninvalid;
        store->monoexpvaluestore->validgrid        = 1;

        store->validmonoexpgrid                    = 1;
    }

    /* Bi-exponential grid... */
    if ((updatetype==BAYES_RAPID_GRID_BI) ||
        (updatetype==BAYES_RAPID_GRID_MONO_AND_BI))
    {
        if ((ntaus_bi>0) && (taus_bi) && (nweights_bi>0) && (weights_bi) && (low) && (high))
        {
            for (i=0; i<=4; i++)
            {
                store->biexpvaluestore->low[i]  = low[i];
                store->biexpvaluestore->high[i] = high[i];
            }

            /* Bi-exp fluorescence likelihoods population */
#if 0
			// Multi-thread option:
            bayes_PopulateRapidLikelihoodsDiscreteValueVectorMultiThreaded(store->biexpvaluestore->fluorescencelikelihoods,
                                                                           ntaus_bi,taus_bi,
                                                                           nbins,binwalls,instr,interval,modulationperiod,
                                                                           ntaus_bi,&nvalid,&ninvalid);
#else
			// Single Thread:
			bayes_PopulateRapidLikelihoodsDiscreteValueVector(store->biexpvaluestore->fluorescencelikelihoods,
				single_low, single_high,
				ntaus_bi, taus_bi,
				nbins, binwalls, instr, interval, modulationperiod,
				&nvalid, &ninvalid);
#endif

            store->biexpvaluestore->nlikelihoods        = ntaus_bi;
            store->biexpvaluestore->nvalidlikelihoods   = nvalid;
            store->biexpvaluestore->ninvalidlikelihoods = ninvalid;

#if 0
            ret = bayes_PopulateRapidBiExpDiscreteValueMatrixMultiThreaded(store->biexpvaluestore->likelihoods,low,high,
                                                                           ntaus_bi,taus_bi,nweights_bi,weights_bi,
                                                                           backgroundmin_bi,backgroundmax_bi,
                                                                           nbins,binwalls,instr,interval,modulationperiod,
                                                                           store->biexpvaluestore->fluorescencelikelihoods,
                                                                           store->biexpvaluestore->nstates,&nvalid,&ninvalid);
#else

			ret = bayes_PopulateRapidBiExpDiscreteValueMatrix(store->biexpvaluestore->likelihoods, low, high,
				ntaus_bi, taus_bi, nweights_bi, weights_bi,
				backgroundmin_bi, backgroundmax_bi,
				nbins, binwalls, instr, interval, modulationperiod,
				store->biexpvaluestore->fluorescencelikelihoods,
				&nvalid, &ninvalid);
#endif



            if (ret>=BAYES_RAPID_RESULT_NO_ERROR)
            {
                /* Bi-exponential grid settings... */
                store->biexpvaluestore->settings->ntaus         = ntaus_bi;
                store->biexpvaluestore->settings->nweights      = nweights_bi;

                for (i=0; i<ntaus_bi; i++)
                    store->biexpvaluestore->settings->tau[i]    = taus_bi[i];

                for (i=0; i<nweights_bi; i++)
                    store->biexpvaluestore->settings->weight[i] = weights_bi[i];

                store->biexpvaluestore->settings->backgroundmax = backgroundmax_bi;
                store->biexpvaluestore->settings->backgroundmin = backgroundmin_bi;
                store->biexpvaluestore->nbins                   = nbins;
                bayes_CopyInstrRspConfigParams(instr,&store->biexpvaluestore->instr);
                store->biexpvaluestore->interval                = interval;
                store->biexpvaluestore->modulationperiod        = modulationperiod;
                store->biexpvaluestore->nvalidstates            = nvalid;
                store->biexpvaluestore->ninvalidstates          = ninvalid;
                store->biexpvaluestore->validgrid               = 1;

                bayes_RapidBiExpDetermineGridSizeAdv(ntaus_bi,taus_bi,nweights_bi,weights_bi,
                                                  backgroundmin_bi,backgroundmax_bi,low,high,
                                                  nbins,&npts,&nvalid,&store->biexpvaluestore->megabytes);
            }

            store->validbiexpgrid = 1;
        }    
    }

    if (ret==BAYES_RAPID_RESULT_USER_CANCEL)
        return (BAYES_RAPID_RESULT_USER_CANCEL);
    else
        return (BAYES_RAPID_RESULT_NO_ERROR);
}


int bayes_CreateRapidValueStore(BayesRapidValueStore_t *store,
                                int                     updatetype,
                                /* Mono-exponential settings... */
                                int                     ntaus,
                                double                 *taus,
                                int                     nweights,
                                double                 *weights,
                                double                  backgroundmax,
                                double                  backgroundmin,
                                /* Bi-exponential settings... */
                                int                     ntaus_bi,
                                double                 *taus_bi,
                                int                     nweights_bi,
                                double                 *weights_bi,
                                double                  backgroundmin_bi,
                                double                  backgroundmax_bi,
                                int                    *low,
                                int                    *high,
                                /* General (common) settings... */
                                int                     nbins,
                                double                 *binwalls,
                                BayesInstrRsp_t        *instr,
                                double                  interval,
                                double                  modulationperiod)
{
    int ret, nw0s=0, nw0s_bi=0;

//    if ((nweights>0) && (weights))
//        nw0s    = (int)((double)nweights*backgroundmax/(weights[nweights-1]-weights[0]));
    
//    if ((nweights_bi>0) && (weights_bi))
//        nw0s_bi = (int)((double)nweights_bi*backgroundmax_bi/(weights_bi[nweights_bi-1]-weights_bi[0]));

    ret = bayes_AllocateForRapidValueStore(store,updatetype,
                                           ntaus,nweights,
                                           ntaus_bi,taus_bi,nweights_bi,weights_bi,low,high,
                                           nbins);
    
    if (ret<0)
        return (ret);

    ret = bayes_PopulateRapidValueStore(store,updatetype,ntaus,taus,nweights,weights,
                                                         backgroundmax,backgroundmin,
                                                         ntaus_bi,taus_bi,nweights_bi,weights_bi,
                                                         backgroundmin_bi,backgroundmax_bi,low,high,
                                                         nbins,binwalls,instr,interval,modulationperiod);

    if (ret==BAYES_RAPID_RESULT_USER_CANCEL)
    {
        bayes_DestroyRapidValueStore(store,updatetype);
        return (BAYES_RAPID_RESULT_USER_CANCEL);
    }

//    bayes_OutputGridToFile("grid.txt","",store);

    return (ret);
}


int bayes_OutputGridToFile(char                   *filename,
                           char                   *comment,
                           BayesRapidValueStore_t *grid)
{
	FILE *fp;
	int i, j, k, ell, bin;
        
    if ((fp = fopen(filename,"w")))
        return (-1);
    
    fprintf(fp,"#==========================================================#\n");
    fprintf(fp,"#               BAYESIAN RAPID GRID DUMP                   #\n");
    fprintf(fp,"#==========================================================#\n\n");

    fprintf(fp,"Comments:\n%s\n\n",comment);

    fprintf(fp,"#==========================================================#\n");
    fprintf(fp,"#                 MONO-EXPONENTIAL GRID                    #\n");
    fprintf(fp,"#==========================================================#\n\n");

    if (grid->monoexpvaluestore)
    {
        /* Grid settings... */
        fprintf(fp,"Settings ==>\n\n");
        fprintf(fp,"ntaus: %d\n",grid->monoexpvaluestore->settings->ntaus);
        fprintf(fp,"nweights: %d\n",grid->monoexpvaluestore->settings->nweights);
        fprintf(fp,"backgroundmin: %g\n",grid->monoexpvaluestore->settings->backgroundmin);
        fprintf(fp,"backgroundmax: %g\n\n",grid->monoexpvaluestore->settings->backgroundmax);
        fprintf(fp,"taus ==>\n");
        for (i=0; i<grid->monoexpvaluestore->settings->ntaus; i++) fprintf(fp,"%g\t",grid->monoexpvaluestore->settings->tau[i]);
        fprintf(fp,"\n\nweights ==>\n");
        for (i=0; i<grid->monoexpvaluestore->settings->nweights; i++) fprintf(fp,"%g\t",grid->monoexpvaluestore->settings->weight[i]);
        fprintf(fp,"\n\n\n");

        /* Instrument and data... */
        fprintf(fp,"Instrument and data ==>\n\n");
        fprintf(fp,"nbins: %d\n",grid->monoexpvaluestore->nbins);
        fprintf(fp,"interval: %g\n",grid->monoexpvaluestore->interval);
        fprintf(fp,"modulationperiod: %g\n\n\n",grid->monoexpvaluestore->modulationperiod);

        /* Fluorescence likelihoods... */
        fprintf(fp,"Fluorescence likelihoods (total: %d, valid: %d, invalid: %d) ==>\n",grid->monoexpvaluestore->nlikelihoods,
                                                                                        grid->monoexpvaluestore->nvalidlikelihoods,
                                                                                        grid->monoexpvaluestore->ninvalidlikelihoods);

        fprintf(fp,"valid\t");
        fprintf(fp,"tau\t");
        for (bin=0; bin<grid->monoexpvaluestore->nbins; bin++) fprintf(fp,"%d\t",bin);
        fprintf(fp,"\n");

        for (i=0; i<grid->monoexpvaluestore->nlikelihoods; i++)
        {
            fprintf(fp,"%d\t",grid->monoexpvaluestore->fluorescencelikelihoods[i].valid);
            fprintf(fp,"%g\t",grid->monoexpvaluestore->fluorescencelikelihoods[i].tau);

            if (grid->monoexpvaluestore->fluorescencelikelihoods[i].valid==1)
                for (bin=0; bin<grid->monoexpvaluestore->nbins; bin++)
                    fprintf(fp,"%g\t",grid->monoexpvaluestore->fluorescencelikelihoods[i].fluorescencedecayphotonlikelihoodsgiventau[bin]);
            
            fprintf(fp,"\n");
        }

        fprintf(fp,"\n\n");
    }
    else
    {
        fprintf(fp,"No mono-exponential grid available.\n\n\n");
    }


    fprintf(fp,"#==========================================================#\n");
    fprintf(fp,"#                  BI-EXPONENTIAL GRID                     #\n");
    fprintf(fp,"#==========================================================#\n\n");

    if (grid->biexpvaluestore)
    {
        /* Grid settings... */
        fprintf(fp,"Settings ==>\n\n");
        fprintf(fp,"ntaus: %d\n",grid->biexpvaluestore->settings->ntaus);
        fprintf(fp,"nweights: %d\n",grid->biexpvaluestore->settings->nweights);
        fprintf(fp,"backgroundmin: %g\n",grid->biexpvaluestore->settings->backgroundmin);
        fprintf(fp,"backgroundmax: %g\n\n",grid->biexpvaluestore->settings->backgroundmax);
        fprintf(fp,"taus ==>\n");
        for (i=0; i<grid->biexpvaluestore->settings->ntaus; i++) fprintf(fp,"%g\t",grid->biexpvaluestore->settings->tau[i]);
        fprintf(fp,"\n\nweights ==>\n");
        for (i=0; i<grid->biexpvaluestore->settings->nweights; i++) fprintf(fp,"%g\t",grid->biexpvaluestore->settings->weight[i]);
        fprintf(fp,"\n\n\n");

        /* Instrument and data... */
        fprintf(fp,"Instrument and data ==>\n\n");
        fprintf(fp,"nbins: %d\n",grid->biexpvaluestore->nbins);
        fprintf(fp,"interval: %g\n",grid->biexpvaluestore->interval);
        fprintf(fp,"modulationperiod: %g\n\n\n",grid->biexpvaluestore->modulationperiod);

        /* Fluorescence likelihoods... */
        fprintf(fp,"Fluorescence likelihoods (total: %d, valid: %d, invalid: %d) ==>\n",grid->biexpvaluestore->nlikelihoods,
                                                                                        grid->biexpvaluestore->nvalidlikelihoods,
                                                                                        grid->biexpvaluestore->ninvalidlikelihoods);

        fprintf(fp,"valid\t");
        fprintf(fp,"tau\t");
        for (bin=0; bin<grid->biexpvaluestore->nbins; bin++) fprintf(fp,"%d\t",bin);
        fprintf(fp,"\n");

        for (i=0; i<grid->biexpvaluestore->nlikelihoods; i++)
        {
            fprintf(fp,"%d\t",grid->biexpvaluestore->fluorescencelikelihoods[i].valid);
            fprintf(fp,"%g\t",grid->biexpvaluestore->fluorescencelikelihoods[i].tau);

            if (grid->biexpvaluestore->fluorescencelikelihoods[i].valid==1)
                for (bin=0; bin<grid->biexpvaluestore->nbins; bin++)
                    fprintf(fp,"%g\t",grid->biexpvaluestore->fluorescencelikelihoods[i].fluorescencedecayphotonlikelihoodsgiventau[bin]);
            
            fprintf(fp,"\n");
        }

        /* Bi-exponential pre-computed likelihoods... */
        fprintf(fp,"\n\nBi-exponential likelihoods ==>\n\n");
        fprintf(fp,"memory (MB): %g\n\n",grid->biexpvaluestore->megabytes);
        fprintf(fp,"low (indices):\t");
        for (i=1; i<=4; i++) fprintf(fp,"%d\t",grid->biexpvaluestore->low[i]);
        fprintf(fp,"\n");
        fprintf(fp,"high (indices):\t");
        for (i=1; i<=4; i++) fprintf(fp,"%d\t",grid->biexpvaluestore->high[i]);
        fprintf(fp,"\n\n");

        fprintf(fp,"low (values):\t");
        for (i=1; i<=2; i++) fprintf(fp,"%g\t",grid->biexpvaluestore->settings->weight[grid->biexpvaluestore->low[i]]);
        for (; i<=4; i++)    fprintf(fp,"%g\t",grid->biexpvaluestore->settings->tau[grid->biexpvaluestore->low[i]]);
        fprintf(fp,"\n");
        fprintf(fp,"high (values):\t");
        for (i=1; i<=2; i++) fprintf(fp,"%g\t",grid->biexpvaluestore->settings->weight[grid->biexpvaluestore->high[i]]);
        for (; i<=4; i++)    fprintf(fp,"%g\t",grid->biexpvaluestore->settings->tau[grid->biexpvaluestore->high[i]]);
        fprintf(fp,"\n\n");

        fprintf(fp,"nstates: %d\n",grid->biexpvaluestore->nstates);
        fprintf(fp,"nvalidstates: %d\n",grid->biexpvaluestore->nvalidstates);
        fprintf(fp,"ninvalidstates: %d\n\n",grid->biexpvaluestore->ninvalidstates);

        fprintf(fp,"valid\t");
        fprintf(fp,"w0\t");
        fprintf(fp,"w1\t");
        fprintf(fp,"w2\t");
        fprintf(fp,"tau1\t");
        fprintf(fp,"tau2\n");

        for (i=grid->biexpvaluestore->low[1]; i<=grid->biexpvaluestore->high[1]; i++)
        {
            for (j=grid->biexpvaluestore->low[2]; j<=grid->biexpvaluestore->high[2]; j++)
            {
                for (k=grid->biexpvaluestore->low[3]; k<=grid->biexpvaluestore->high[3]; k++)
                {
                    for (ell=grid->biexpvaluestore->low[4]; ell<=grid->biexpvaluestore->high[4]; ell++)
                    {
                        if (grid->biexpvaluestore->likelihoods[i][j][k][ell])
                        {
                            fprintf(fp,"%d\t",1);
                            fprintf(fp,"%g\t",grid->biexpvaluestore->likelihoods[i][j][k][ell]->weights[0]);
                            fprintf(fp,"%g\t",grid->biexpvaluestore->likelihoods[i][j][k][ell]->weights[1]);
                            fprintf(fp,"%g\t",grid->biexpvaluestore->likelihoods[i][j][k][ell]->weights[2]);
                            fprintf(fp,"%g\t",grid->biexpvaluestore->likelihoods[i][j][k][ell]->taus[1]);
                            fprintf(fp,"%g\t",grid->biexpvaluestore->likelihoods[i][j][k][ell]->taus[2]);
                            fprintf(fp,"\n");
                        }
                        /*else
                        {
                            fprintf(fp,"%d\n",0);
                        }*/
                    }
                }
            }
        }
    }
    else
    {
        fprintf(fp,"No bi-exponential grid available.\n\n");
    }

    fclose(fp);

    return (0);
}


int bayes_RapidBiExpDetermineGridSizeAdv(/* Input... */
                                      int     ntaus,
                                      double *taus,
                                      int     nweights,
                                      double *weights,
                                      double  backgroundmin,
                                      double  backgroundmax,
                                      int    *low,
                                      int    *high,
                                      int     nbins,
                                      /* Output... */
                                      int    *npts,
                                      int    *nptsvalid,
                                      double *megabytesreqd)
{
    int i;

    for (i=1,*npts=1; i<=4; i++)
        *npts *= (high[i]-low[i]+1);

    *nptsvalid = bayes_ComputeRapidBiExpDiscreteValuesNumOfValidEntries(low,high,
                                                                        ntaus,taus,nweights,weights,
                                                                        backgroundmin,backgroundmax);

    *megabytesreqd = bayes_ComputeRapidBiExpDiscreteValuesNumOfMegaBytesReqd(*npts,*nptsvalid,nbins);

    return (0);
}

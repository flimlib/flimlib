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
/* mrowley - creation 070904 */


#include "stdio.h"
#include "math.h"

#include "extmath.h"
#include "random.h"
#include "matrices.h"
#include "bayes_Sizes.h"

#include "bayes_MultiExpRapidAnalysis.h" // for extended debug traces only
#include "bayes_ModelTransformTools.h" // for extended debug traces only


/*--------------- mathematical manipulations ------------------------*/

// https://www.johndcook.com/blog/cpp_erf/
// Public Domain erf code
// (Previously code was tested with the GSL error function gsl_sf_erf_e)

/*
double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
	double t, y;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    t = 1.0/(1.0 + p*x);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}
*/

double mod(double x)
{
    if (x >= 0.0)
		return (x);

	return (0.0-x);
}


/* Minimisation Routines */

#define MAX_N 1000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

double math_AmotryFctDoubleWithGenericContainer(double **p,
                                                double  *y,
                                                double  *psum,
                                                int      ndim,
                                                double  (*func)(double *, int, void *),
                                                int      id,
                                                void    *container,
                                                int      high_index,
                                                int     *nfunc,
                                                double   f)
{
    int j, count=*nfunc;
    double f1, f2, new_y, *ptry=NULL;

    ptry = Bayes_dvector(1,ndim);
	if(!ptry) return 0.0;

    f1 = (1.0-f)/ndim;
    f2 = f1-f;

    for (j=1; j<=ndim; j++)
        ptry[j] = (f1*psum[j]) - (f2*p[high_index][j]);

    new_y = (*func)(ptry, id, container);  // Get new y
	count++;

    if (new_y < y[high_index]) 
    {
        y[high_index] = new_y;

        for (j=1; j<=ndim; j++)
        {
            psum[j]  += ptry[j] - p[high_index][j];
            p[high_index][j] = ptry[j];
        }
    }

    free_Bayes_dvector(ptry, 1, ndim);

	*nfunc = count;

    return (new_y);
}



int math_MinimiseFctDoubleWithGenericContainer(double (*func)(double *, int, void *),
                                               int    id,
                                               void   *container,
                                               int    ndim,
                                               double *w,
                                               double *value,
                                               void   *config)
{
//* Based on "amoeba" for Downhill Simplex Method
//*    MONITOR = 0,1     whether or not to show intermediate results   
//*    func()            the function to be minimised                  
//*    ndim = 1,2,...    the number of arguments of func()             
//*    *value            the actual minimum value of func() found      
//*    w[]           the vector w the minimum is located       
//*                      w[1],...,w[ndim]                      
//*                      must be initialised to serve as a starting    
//*                      position for amoeba                           
//*    tolerance         the relative accuracy in the minimum search          

    int    ret=MATH_MINIMISATION_RESULT_SUCCESS;
    int    nfunc, i, j, lo_index, high_index, high_n, mpts=ndim+1;
    double **p=NULL, *deltas;
    double *x=NULL, *y=NULL, new_y, ysave, sum, rtol, *psum;

    double tolerance;
    int    MONITOR;

    MONITOR   = ((AmoebaConfigParams_t *)config)->monitor;
    tolerance = ((AmoebaConfigParams_t *)config)->tolerance;
    deltas    = ((AmoebaConfigParams_t *)config)->deltas;

    x = Bayes_dvector(1, ndim);
    y = Bayes_dvector(1, mpts);
    p = Bayes_dmatrix(1, mpts, 1, ndim);

	if(!x || !y || !p)
		return -1;
        
    for (i=1; i<=mpts; i++)
    {
        for (j=1; j<=ndim; j++)
        {
            x[j]  = p[i][j] = w[j];
        }

        if (i < mpts)
        {
            x[i]    += deltas[i];
            p[i][i] += deltas[i];
        }

        y[i] = func(x, id, container);
    }

    psum = Bayes_dvector(1, ndim);
	if(!psum)
		return -1;

    nfunc = 0;

	for (j=1; j<=ndim; j++)
    {
		for (i=1,sum=0.0; i<=mpts; i++)
		{
		    sum += p[i][j];
		}
        psum[j]=sum;
    }

	while(TRUE)
    {
	    
		lo_index = 1;

		if (y[1] > y[2])
		{
			high_n=2;
			high_index=1;
		}
		else
		{
			high_n=1;
			high_index=2;
		}

        for (i=1; i<=mpts; i++)
        {
            if (y[i] < y[lo_index])
                lo_index=i;

            if (y[i] > y[high_index])
            {
                high_n = high_index;
                high_index  = i;
            }
            else if (y[i] > y[high_n])
                if (i != high_index)
                    high_n = i;
        }

        rtol = 2.0 * fabs(y[high_index]-y[lo_index])/(fabs(y[high_index])+fabs(y[lo_index]));

        if (MONITOR)
        { 
            printf(" >>amoeba: min=%lf max=%lf",y[lo_index],y[high_index]);
            puts("");
        } 

        if (rtol < tolerance)
            break;

        if (nfunc >= MAX_N)
        { 
            ret = MATH_MINIMISATION_RESULT_MAX_FCT_CALLS_REACHED;
            break;
        }

        new_y = math_AmotryFctDoubleWithGenericContainer(p,y,psum,ndim,func,id,container,high_index,&nfunc,-ALPHA);

        if (new_y <= y[lo_index])
        {
            new_y = math_AmotryFctDoubleWithGenericContainer(p,y,psum,ndim,func,id,container,high_index,&nfunc,GAMMA);
        }
        else if (new_y >= y[high_n])
        {
            ysave = y[high_index];
            new_y  = math_AmotryFctDoubleWithGenericContainer(p,y,psum,ndim,func,id,container,high_index,&nfunc,BETA);

            if (new_y >= ysave)
            {
                for (i=1; i<=mpts; i++)
                {
                    if (i != lo_index)
                    {
                        for (j=1; j<=ndim; j++)
                        {
                            psum[j] = 0.5*(p[i][j]+p[lo_index][j]);
                            p[i][j] = psum[j];
                        }

                        y[i] = (*func)(psum,id,container);
                    }
                }

                nfunc += ndim;

				for (j=1; j<=ndim; j++)
				{
					for (i=1,sum=0.0; i<=mpts; i++)
					{
						sum += p[i][j];
					}
					psum[j]=sum;
				}
            }
        }
    }

    *value = y[lo_index]; 

    for (j=1; j<=ndim; j++)
        w[j] = p[lo_index][j];

    free_Bayes_dvector(psum,1,ndim);
    free_Bayes_dvector(x,1,ndim);
    free_Bayes_dvector(y,1,mpts);
    free_Bayes_dmatrix(p,1,mpts,1,ndim);

    return (ret);
}

/*
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
 
/** 
 * FLIMLib - Phasor analysis.
 * Classic Phasor or Polar approach to FLIM.
 * See Clayton 2004 or Leray 2008.
 *
 * \file GCI_Phasor.c
 */

#include "GCI_Phasor.h"
#include "GCI_PhasorInternal.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


#ifndef NULL
#define NULL 0
#endif

#define PHASOR_ERR_NO_ERROR                         0
#define PHASOR_ERR_INVALID_DATA                    -1
#define PHASOR_ERR_INVALID_WINDOW                  -2
#define PHASOR_ERR_INVALID_MODEL                   -3
#define PHASOR_ERR_FUNCTIONALITY_NOT_SUPPORTED     -4

/* Internal function prototypes */
void setPhasorPeriod(double period);


static double phasorPeriod;

void setPhasorPeriod(double period)
{
	phasorPeriod=period;
}

double GCI_Phasor_getPeriod()
{	// calculate and return the used phasor period
	// value is only valid after GCI_Phasor has been called
	// returns the period in current time units
	return phasorPeriod;
}

/** GCI_Phasor.

 See Clayton 2004 or Leray 2008
 Classic Phasor or Polar approach to FLIM

 u = integral(data(t) * cos(wt)) dt / integral(data(t)) dt
 v = integral(data(t) * sin(wt)) dt / integral(data(t)) dt
 
 tau phi = taup = 1/w * (v/u)
 tau mod = taum = 1/w * sqrt(1/(u^2 + v^2) - 1)
 tau average = tau = (taup + taum) / 2

 Z must have been estimated previously so that it can be subtracted from the data here
 A is estimated by making the photon count in the fit the same as the data

 chisq is calculated for comparison with other methods but is not used as there is no optimisation with this method
*/

int GCI_Phasor(float xincr, float y[], int fit_start, int fit_end,
							  const float *Z, float *U, float *V, float *taup, float *taum, float *tau, float *fitted, float *residuals,
							  float *chisq)
{
	int nBins = (fit_end - fit_start);
	if (nBins < 0)
		return (PHASOR_ERR_INVALID_WINDOW);
	float* cosine = malloc((size_t)nBins * sizeof(float));
	float* sine = malloc((size_t)nBins * sizeof(float));
	createSinusoids(nBins, cosine, sine);
	int ret = GCI_Phasor_compute(xincr, y, fit_start, fit_end, Z, cosine, sine, U, V, taup, taum, tau, fitted, residuals, chisq);
	free(cosine);
	free(sine);
	return ret;
}

void createSinusoids(int nBins, float* cosine, float* sine) {
	float w = 2.0f * 3.1415926535897932384626433832795028841971f / (float)nBins; //2.0*PI/(float)nBins;
	// Take care that values correspond to the centre of the bin, hence i+0.5
	for (int i = 0; i < nBins; i++) {
		cosine[i] = cosf(w * ((float)i + 0.5f));
		sine[i] = sinf(w * ((float)i + 0.5f));
	}
}

int GCI_Phasor_compute(float xincr, float y[], int fit_start, int fit_end,
	const float* Z, float* cosine, float* sine, float* U, float* V, float* taup, float* taum, float* tau, float* fitted, float* residuals,
	float* chisq)
{
    // Z must contain a bg estimate
	// fitted and residuals must be arrays big enough to hold possibly fit_end floats.

	int   i, ret = PHASOR_ERR_NO_ERROR, nBins;
	float *data, u, v, A, w, I, Ifit, bg, chisq_local, res, sigma2, *validFittedArray;

	// we require residuals or chisq but have not supplied a "fitted" array so must malloc one
	if (fitted==NULL && (residuals!=NULL || chisq!=NULL)){
		if ((validFittedArray = (float *)malloc((long unsigned int)fit_end * sizeof(float)))== NULL) return (-1);
	}
	else validFittedArray = fitted;

	data = &(y[fit_start]);	
	nBins = (fit_end - fit_start);
	bg = *Z;
    if (!data)
        return (PHASOR_ERR_INVALID_DATA);
    if (nBins<0)
        return (PHASOR_ERR_INVALID_WINDOW);

	// rep frequency, lets use the period of the measurement, but we can stay in the units of bins
	w = 2.0f*3.1415926535897932384626433832795028841971f/(float)nBins; //2.0*PI/(float)nBins;
	setPhasorPeriod((float)nBins*xincr); // store the real phasor period used for future external use.

	// integral over data
	for (i=0, I=0.0f; i<nBins; i++) 
		I += (data[i]-bg);

	// Phasor coords
	for (i = 0, u = 0.0f, v = 0.0f; i < nBins; i++) {
		u += (data[i] - bg) * cosine[i];
		v += (data[i] - bg) * sine[i];
	}
	u /= I;
	v /= I;

	// taus, convert now to real time with xincr
	*taup = (xincr/w) * (v/u);
	*taum = (xincr/w) * sqrtf(1.0f/(u*u + v*v) - 1.0f);

	*tau = ((*taup) + (*taum))/2.0f;

	*U = u;
	*V = v;

	/* Now calculate the fitted curve and chi-squared if wanted. */
	/* if validFittedArray is NULL, malloc was not used */
	if (validFittedArray == NULL)
		return ret;
	memset(validFittedArray, 0, (size_t)fit_end * sizeof(float));
	if (residuals != NULL)
		memset(residuals, 0, (size_t)fit_end * sizeof(float));
	// Madison report some "memory issue", and replaced the 2 line above with new arrays.
	// Not sure what that was but I breaks the filling of the fitted array, probably not allocating the arrays before calling

	// integral over nominal fit data
	for (Ifit=0.0f, i=fit_start; i<fit_end; i++) 
		Ifit += expf((float)(-(i-fit_start))*xincr/(*tau));
	// Estimate A
	A = I / Ifit;

	// Calculate fit
	for (i=fit_start; i<fit_end; i++){
		validFittedArray[i] = bg + A * expf((float)(-(i-fit_start))*xincr/(*tau));
	}
	// OK, so now fitted contains our data for the timeslice of interest.
	// We can calculate a chisq value and plot the graph, along with
	// the residuals.

	/* if both residuals and chisq are NULL, malloc was not used for validFittedArray */
	if (residuals == NULL && chisq == NULL)
		return ret;

	chisq_local = 0.0f;
	for (i=0; i<fit_start; i++) {
		res = y[i]-validFittedArray[i];
		if (residuals != NULL)
			residuals[i] = res;
	}


//	case NOISE_POISSON_FIT:
		/* Summation loop over all data */
		for (i=fit_start ; i<fit_end; i++) {
			res = y[i] - validFittedArray[i];
			if (residuals != NULL)
				residuals[i] = res;
			/* don't let variance drop below 1 */
			sigma2 = (validFittedArray[i] > 1 ? 1.0f/validFittedArray[i] : 1.0f);
			chisq_local += res * res * sigma2;
		}

	if (chisq != NULL)
		*chisq = chisq_local;

	if (fitted==NULL){
		free (validFittedArray);
	}

	return (ret);
}

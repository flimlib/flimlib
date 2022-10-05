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

 /* This is GCI_PhasorInternal.h, the header file for internal functions in
	the 2022 version of the GCI_Phasor.c file. */

#ifndef _PHASOR_FITTING_INTERNAL_H
#define _PHASOR_FITTING_INTERNAL_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Generates one period of both a cosine and a sine. Values are computed from
 * the center of the bins.
 * 
 * \param[in] nBins The number of data points in the resulting sinusoids.
 * \param[out] cosine
 * \param[out] sine
 */
void createSinusoids(int nBins, float* cosine, float* sine);

/**
 * Handles the computation for phasor analysis utilizing a provided sine and cosine
 */
int GCI_Phasor_compute(float xincr, float y[], int fit_start, int fit_end,
	const float* Z, float* cosine, float* sine, float* U, float* V, float* taup, 
	float* taum, float* tau, float* fitted, float* residuals, float* chisq);

#ifdef __cplusplus
}
#endif

#endif /* _PHASOR_FITTING_INTERNAL_H */

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

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PHASOR_FITTING_H
#define PHASOR_FITTING_H

/** 
 * FLIMLib - Phasor analysis Header.
 * Classic Phasor or Polar approach to FLIM.
 * See Clayton 2004 or Leray 2008.
 *
 * \file GCI_Phasor.h
 */

/**
 * Take transient data and perform phasor analysis, returning u, v, and estimated lifetime.
 *
 * u = integral(data(t) * cos(wt)) dt / integral(data(t)) dt
 * v = integral(data(t) * sin(wt)) dt / integral(data(t)) dt
 *
 * tau phi = taup = 1/w * (v/u)
 * tau mod = taum = 1/w * sqrt(1/(u^2 + v^2) - 1)
 * tau average = tau = (taup + taum) / 2
 *
 * \param[in] xincr The time increment inbetween the values in the y array.
 * \param[in] y The transient (time resolved) signal to be analysed, the 'data'.
 * \param[in] fit_start The index into the y array marking the start to the data to be used in the fit.
 * \param[in] fit_end The index into the y array marking the end of the data to be used in the fit.
 * \param[in] Z must have been estimated previously so that it can be subtracted from the data here.
 * \param[out] u The 'horizontal' phasor coordinate.
 * \param[out] v The 'vertical' phasor coordinate.
 * \param[out] taup The lifetime calculated from the phase change.
 * \param[out] taum The lifetime calculated from the amplitude change (the demodulation).
 * \param[out] tau The average of the other taus.
 * \param[out] fitted An array containing values fitted to the data, the 'fit'. Fit points are coincident in time with the data points.
 * \param[out] residuals An array containing the difference between the data and the fit.
 * \param[out] chisq The resulting reduced chi squared value of the fit.
 * \return An error code, 0 = success.
 */
int    GCI_Phasor(float xincr, float y[], int fit_start, int fit_end,
							  const float *Z, float *u, float *v, float *taup, float *taum, float *tau, float *fitted, float *residuals,
							  float *chisq);

/**
 * Get the phasor period that was calculated and used in the last call to GCI_Phasor.
 *
 * \return The period in the time units of xinc.
 *
 */
double GCI_Phasor_getPeriod();

#ifdef __cplusplus
}
#endif

#endif /* PHASOR_FITTING_H */

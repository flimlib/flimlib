/* 
This file is part of the SLIM-curve package for exponential curve fitting of spectral lifetime data.

Copyright (c) 2010-2014, Gray Institute University of Oxford & UW-Madison LOCI.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PHASOR_FITTING_H
#define PHASOR_FITTING_H

int    GCI_Phasor(float xincr, float y[], int fit_start, int fit_end,
							  float *Z, float *u, float *v, float *taup, float *taum, float *tau, float *fitted, float *residuals,
							  float *chisq);

double GCI_Phasor_getPeriod();

#endif /* PHASOR_FITTING_H */

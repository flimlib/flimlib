/*-
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2025 University of Oxford and Board of Regents of the
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
#pragma once

typedef void (*fitfunc)(float, float [], float *, float [], int);

class FitFunc {
public:
	virtual float fit(float x, float param[], float dy_dparam[], int nparam) {
		float y = 0;
		if (this->func_ptr)
			this->func_ptr(x, param, &y, dy_dparam, nparam);
		else
			fprintf(stderr, "Warning: FitFunc->func_ptr invalid.");
		return y;
	}

	FitFunc() : func_ptr(NULL) {}

	FitFunc(fitfunc func_ptr) : func_ptr(func_ptr) {}

	virtual ~FitFunc() {}

	// needed for creating array to feed into java callback
	// int nparam;

private:
	const fitfunc func_ptr;
};

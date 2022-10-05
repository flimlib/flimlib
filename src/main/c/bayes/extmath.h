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
#ifndef BAYES_EXT_MATH_H
#define BAYES_EXT_MATH_H

#include "bayes_Interface.h"
#include "math.h"

#define MATH_MINIMISATION_RESULT_MAX_FCT_CALLS_REACHED  +1  //warning rather than an error
#define MATH_MINIMISATION_RESULT_SUCCESS                  0
#define MATH_MINIMISATION_RESULT_ERROR_INVALID_INPUT    -1
#define MATH_MINIMISATION_RESULT_ERROR_INVALID_SETTINGS -2
#define MATH_MINIMISATION_RESULT_USERCANCEL             -99

struct MathContainerAmoebaConfigParams
{
    double  tolerance;
    int     monitor;
    double *deltas;
};
typedef struct MathContainerAmoebaConfigParams AmoebaConfigParams_t;

int math_MinimiseFctDoubleWithGenericContainer(double (*funk)(double *, int, void *),
                                               int    id,
                                               void   *container,
                                               int    ndim,
                                               double *where,
                                               double *value,
                                               void   *config);


//double erf(double Y);

double mod(double x);

#endif /* BAYES_EXT_MATH_H */

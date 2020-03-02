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

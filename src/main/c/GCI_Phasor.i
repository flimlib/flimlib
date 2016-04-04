/* File: GCI_Phasor.i */
%module GCI_Phasor 
/* For things like float* */
/*%include "typemaps.i"*/
%include "cpointer.i"
%include "arrays_java.i"

%{
#include "GCI_Phasor.h"
%}

/*%apply float[] {float*} */
%pointer_functions(float, floatp);

#define PHASOR_ERR_NO_ERROR                         0
#define PHASOR_ERR_INVALID_DATA                    -1
#define PHASOR_ERR_INVALID_WINDOW                  -2
#define PHASOR_ERR_INVALID_MODEL                   -3
#define PHASOR_ERR_FUNCTIONALITY_NOT_SUPPORTED     -4

/*extern double phasorPeriod;*/
extern double GCI_Phasor_getPeriod();
extern int GCI_Phasor(float xincr, float y[], int fit_start, int fit_end, float *Z, float *U, float *V, float *taup, float *taum, float *tau, float *fitted, float *residuals, float *chisq);

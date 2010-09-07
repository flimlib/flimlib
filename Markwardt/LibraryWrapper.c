#include "LibraryWrapper.h"
#include "mpfit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* This is the private data structure which contains the data points
   and their uncertainties */
struct vars_struct {
    int fit_start;
    int fit_end;
    double x_incr;
    double *y;
};

 double fun(double x, double *param, int n) {
    if (n != 3) {
        printf("Can't handle %d parameters\n", n);
        return 0.0;
    }
    return param[0] * exp (-param[1] * x) + param[2];
}

 /*
 * "The user must define a function which computes the appropriate values as specified above.
 * The user function may also optionally compute explicit derivatives (see below).
 * The function must compute the weighted deviations between the model and the data.
 * The user function should be of the following form:"
 *
 * "The data points should be passed in the private parameter.
 * This parameter is passed but not modified by mpfit(), and can have any form desired by the user.
 *
 * int m     - number of data points
 * int n     - number of parameters
 * double *p - array of n parameters
 * double *deviates - array of m deviates to be returned by myfunct()
 * double **derivs - used for user-computed derivatives (see below)
 *                   (= 0  when automatic finite differences are computed)
 *
 * User functions may also indicate a fatal error condition by returning a negative error code.
 *  Error codes from -15 to -1 are reserved for the user function."
 */
int myfunct(int m, int n, double *params, double *deviates,
        double **derivs, void *private)
 {
   // printf("m is %d n is %d\n", m, n);
   // return 0;
   printf("params %g %g %g\n", params[0], params[1], params[2]);
   // printf("estimated lambda is %g\n", params[1]);

    // retrieve private data
    struct vars_struct *v = (struct vars_struct *) private;

  //  return 0;
    
    int i;
    double x = v->fit_start * v->x_incr;
    /* Compute function deviates */
    for (i = 0; i < m; i++) {
        deviates[i] = v->y[i] - fun(x, params, n);
        x += v->x_incr;
        printf(" y %g dev %g", v->y[i], deviates[i]);
    }
    printf("\n");

    /* If derivs is non-zero then user-computed derivatives are //TODO compute deviates and derivs in a single loop, to reuse exp(-params[1]*x) term.
       requested */
    if (derivs) {
        x = v->fit_start * v->x_incr;
        for (i = 0; i < m; i++) {
            if (derivs[0]) {
                derivs[0][i] = exp(-params[1] * x);
            }
            if (derivs[1]) {
                derivs[1][i] = -params[0] * x * exp(-params[1] * x);
            }
            if (derivs[2]) {
                derivs[2][i] = 1;
            }
            printf("derivs for %d %g %g %g\n", i, derivs[0][i], derivs[1][i], derivs[2][i]); //TODO would bomb if not all derivatives requested
            x += v->x_incr;
        }
        printf("\n");
    }
    else printf("no derivs\n");

    /*
          if (parameterIndex == 2 * numExp) return 1; // c term
      int e = parameterIndex / 2;
      int off = parameterIndex % 2;
      double aTerm = a[2 * e];
      double bTerm = a[2 * e + 1];
      switch (off) {
        case 0:
          return Math.exp(-bTerm * x); // a term
        case 1:
          return -aTerm * x * Math.exp(-bTerm * x); // b term
      }
 */

    return 0;
 }

int markwardt_fit(double x_incr, double y[], int fit_start, int fit_end,
        double param[], int param_free[], int n_param) {

    // save incoming info in our private data structure
    struct vars_struct v;
    v.fit_start = fit_start;
    v.fit_end   = fit_end;
    v.x_incr    = x_incr;
    v.y         = y;
    //printf("okay %g\n", v.x_incr);

    printf("INCOMING Y");
    int ii;
    for (ii = 0; ii < fit_end; ++ii) {
        printf(" %g ", y[ii]);
    }
    printf("\n");
    //int j;
    //for (j = 0; j < n_param; ++j) {
    //    printf(" %g ", param[j]);
    //}
    //printf("\n");

    // allocate result structure
    struct mp_result_struct result;
    result.covar = 0;
    result.xerror = 0;
    result.resid = 0;

    struct mp_par_struct par[n_param];
    int i;
    for (i = 0; i < n_param; ++i) {
        par[i].fixed = 0;
        par[i].limited[0] = 0;
        par[i].limited[1] = 0;
        par[i].parname = 0;
        par[i].step = 0;
        par[i].relstep = 0;
        par[i].side = 3;
        par[i].deriv_debug = 1;
    }

    // call fitting function
    //TODO need to acount for param_free!
    //TODO first 0 is mp_par second 0 is mp_config
    int m = fit_end - fit_start + 1;
    int status = mpfit(myfunct, m, n_param, param, par, 0, (void *) &v, &result);

    printf("status %d\n", status);
    int k;
    for (k = 0; k < n_param; ++k) {
        printf(" %g", param[k]);
    }
    printf("\n");
    return 0;
}




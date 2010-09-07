#include "EcfWrapper.h"
#include <stdio.h>

int RLD_fit(
        double x_inc,
        double y[],
        int fit_start,
        int fit_end,
        double instr[],
        int n_instr,
        double sig[],
        double *z,
        double *a,
        double *tau,
        double fitted[],
        double *chi_square,
        double chi_square_target
        ) {

    printf("RLD_fit\n");

    int noise = 0;
    int n_data = fit_end + 1;
    float y_float[n_data];
    float sig_float[n_data];
    int i;
    for (i = 0; i < n_data; ++i) {
       // y_float[i] = (float) y[i];
       // sig_float[i] = (float) sig[i];
    }
    float residuals[n_data];
    /////float fitted_float[]; //TODO how to handle these outgoing values; c/b null
    /////float instr_float[n_instr];

    int returnValue =  GCI_triple_integral_fitting_engine(
            (float) x_inc,
            y_float,
            fit_start,
            fit_end,
            instr,
            n_instr,
            noise,
            sig_float,
            z,
            a,
            tau,
            fitted,
            residuals,
            chi_square,
            chi_square_target
            );

    printf("returns %d z %g a %g tau %g\n", returnValue, *z, *a, *tau);
    return returnValue;
}

int LMA_fit(
        double x_inc,
        double y[],
        int fit_start,
        int fit_end,
        double instr[],
        int n_instr,
        double sig[],
        double param[],
        int param_free[],
        int n_param,
        double fitted[],
        double *chi_square,
        double chi_square_target
        ) {
    return 0;
}


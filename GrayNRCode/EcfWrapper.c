#include "EcfWrapper.h"
#include "ecf.h"
#include <stdio.h>
#include <stdlib.h>

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

    int noise = 0;
    int n_data = fit_end + 1;

    float *y_float = (float *)malloc(n_data * sizeof(float));
    float *sig_float = (float *)malloc(n_data * sizeof(float));
    int i;
    for (i = 0; i < n_data; ++i) {
       y_float[i] = (float) y[i];
       sig_float[i] = ( sig ? (float) sig[i] : 1.0f );
    }

    float *instr_float = 0;
    if (instr) {
        instr_float = (float *)malloc(n_instr * sizeof(float));
        int i;
        for (i = 0; i < n_instr; ++i) {
            instr_float[i] = (float) instr[i];
        }
    }
    float *fitted_float = 0;
    if (fitted) {
        fitted_float = (float*)malloc(n_data * sizeof(float));
    }
    float *residuals_float = (float *)malloc(n_data * sizeof(float));

    float x_inc_float      = (float) x_inc;
    float z_float          = (float) *z;
    float a_float          = (float) *a;
    float tau_float        = (float) *tau;
    float chi_square_float = (float) *chi_square;

    int returnValue =  GCI_triple_integral_fitting_engine(
            x_inc_float,
            y_float,
            fit_start,
            fit_end + 1, // want exclusive, not inclusive end
            instr_float,
            n_instr,
            noise,
            sig_float,
            &z_float,
            &a_float,
            &tau_float,
            fitted_float,
            residuals_float,
            &chi_square_float,
            chi_square_target
            );

    *z          = (double) z_float;
    *a          = (double) a_float;
    *tau        = (double) tau_float;
    *chi_square = (double) chi_square_float;
    
    free(y_float);
    free(sig_float);

    if (instr_float) {
        free(instr_float);
    }
    if (fitted_float) {
        int i;
        for (i = 0; i < n_data; ++i) {
            fitted[i] = (double) fitted_float[i];
        }
        free(fitted_float);
    }
    free(residuals_float);

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

    int noise = 0;
    int restrain = 0;
    float chi_square_delta = 0;
    float chi_square_percent = 500;
    int n_data = fit_end + 1;

    float *y_float = (float *)malloc(n_data * sizeof(float));
    float *sig_float = (float *)malloc(n_data * sizeof(float));
    int i;
    for (i = 0; i < n_data; ++i) {
       y_float[i] = (float) y[i];
       sig_float[i] = ( sig ? (float) sig[i] : 1.0f );
    }

    float *instr_float = 0;
    if (instr) {
        instr_float = (float *)malloc(n_instr * sizeof(float));
        int i;
        for (i = 0; i < n_instr; ++i) {
            instr_float[i] = (float) instr[i];
        }
    }
    float *fitted_float = 0;
    if (fitted) {
        fitted_float = (float*)malloc(n_data * sizeof(float));
    }
    float *param_float = (float *)malloc(n_param * sizeof(float));
    if (3 == n_param) {
        param_float[0] = (float) param[2]; // z
        param_float[1] = (float) param[0]; // a1
        param_float[2] = (float) param[1]; // tau1
    }
    if (5 == n_param) {
        param_float[0] = (float) param[4]; // z
        param_float[1] = (float) param[0]; // a1
        param_float[2] = (float) param[1]; // tau1
        param_float[3] = (float) param[2]; // a2
        param_float[4] = (float) param[3]; // tau2
    }
    if (7 == n_param) {
        param_float[0] = (float) param[6]; // z
        param_float[1] = (float) param[0]; // a1
        param_float[2] = (float) param[1]; // tau1
        param_float[3] = (float) param[2]; // a2
        param_float[4] = (float) param[3]; // tau2
        param_float[5] = (float) param[4]; // a3
        param_float[6] = (float) param[5]; // tau3
    }
    //for (i = 0; i < n_param; ++i) {
    //    param_float[i] = (float) param[i];
    //}
    float *residuals_float = (float *)malloc(n_data * sizeof(float));
    float **covar_float    = GCI_ecf_matrix(n_data, n_data);
    float **alpha_float    = GCI_ecf_matrix(n_data, n_data);
    float **err_axes_float = GCI_ecf_matrix(n_data, n_data);

    float x_inc_float      = (float) x_inc;
    float chi_square_float = (float) *chi_square;
    
    fit_type fit = FIT_GLOBAL_MULTIEXP;
    void (*fitfunc)(float, float [], float *, float[], int) = NULL;
    switch (fit) {
        case FIT_GLOBAL_MULTIEXP:
            fitfunc = GCI_multiexp_tau;
            break;
        case FIT_GLOBAL_STRETCHEDEXP:
            fitfunc = GCI_stretchedexp;
            break;
    }

 /*   printf("x inc %g\n", x_inc_float);
    printf("param ");
    for (i = 0; i < n_param; ++i) {
        printf(" %g", param_float[i]);
        printf(" %d", param_free[i]);
    }
    printf("\n");
    printf("y ");
    for (i = 0; i < n_data; ++i) {
        printf(" %g", y_float[i]);
    }
    printf("\n");
    printf("instr %d\n", instr_float);
    printf("n instr %d\n", n_instr);
    printf("noise %d\n", noise);
    printf("sig ");
    for (i = 0; i < n_data; ++i) {
        printf(" %g", sig_float[i]);
    }
    printf("\n");*/

    int returnValue = GCI_marquardt_fitting_engine(
            x_inc_float, //TODO not necessary, just cast it, I think
            y_float,
            n_data,
            fit_start,
            fit_end + 1, // want exclusive, not inclusive end
	    instr_float,
            n_instr,
            noise,
            sig_float,
	    param_float,
            param_free,
	    n_param,
            restrain,
	    fitfunc,
	    fitted_float,
            residuals_float,
            &chi_square_float,
	    covar_float,
            alpha_float,
            err_axes_float,
	    chi_square_target,
            chi_square_delta,
            chi_square_percent);

    *chi_square = (double) chi_square_float;

    free(y_float);
    free(sig_float);

    if (instr_float) {
        free(instr_float);
    }
    if (fitted_float) {
        int i;
        for (i = 0; i < n_data; ++i) {
            fitted[i] = (double) fitted_float[i];
        }
        free(fitted_float);
    }
    //for (i = 0; i < n_param; ++i) {
    //    param[i] = (double) param_float[i];
    //}
    if (3 == n_param) {
        param[2] = (double) param_float[0]; // z
        param[0] = (double) param_float[1]; // a1
        param[1] = (double) param_float[2]; // tau1
    }
    if (5 == n_param) {
        param[4] = (double) param_float[0]; // z
        param[0] = (double) param_float[1]; // a1
        param[1] = (double) param_float[2]; // tau1
        param[2] = (double) param_float[3]; // a2
        param[3] = (double) param_float[4]; // tau2
    }
    if (7 == n_param) {
        param[6] = (double) param_float[0]; // z
        param[0] = (double) param_float[1]; // a1
        param[1] = (double) param_float[2]; // tau1
        param[2] = (double) param_float[3]; // a2
        param[3] = (double) param_float[4]; // tau2
        param[4] = (double) param_float[5]; // a3
        param[5] = (double) param_float[6]; // tau3
    }
    free(param_float);
    
    free(residuals_float);
    GCI_ecf_free_matrix(covar_float);
    GCI_ecf_free_matrix(alpha_float);
    GCI_ecf_free_matrix(err_axes_float);

    //printf("return value is %d\n", returnValue);
    return returnValue;
}


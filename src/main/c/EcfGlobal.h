


#ifdef __cplusplus
extern "C" {
#endif

int GCI_marquardt_global_exps_instr(float xincr, float **trans,
                    int ndata, int ntrans, int fit_start, int fit_end,
                    float instr[], int ninstr,
                    noise_type noise, float sig[], int ftype,
                    float **param, int paramfree[], int nparam,
                    restrain_type restrain, float chisq_delta,
                    float **fitted, float **residuals,
                    float chisq_trans[], float *chisq_global, int *df,
                    int drop_bad_transients);
int GCI_marquardt_global_generic_instr(float xincr, float **trans,
                     int ndata, int ntrans, int fit_start, int fit_end,
                     float instr[], int ninstr,
                     noise_type noise, float sig[],
                     float **param, int paramfree[], int nparam, int gparam[],
                     restrain_type restrain, float chisq_delta,
                     void (*fitfunc)(float, float [], float *, float [], int),
                     float **fitted, float **residuals,
                     float chisq_trans[], float *chisq_global, int *df);


#ifdef __cplusplus
}
#endif

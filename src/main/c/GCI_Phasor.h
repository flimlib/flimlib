#ifndef PHASOR_FITTING_H
#define PHASOR_FITTING_H



#define PHASOR_MONO_EXP_MODEL    1


int    GCI_Phasor(float xincr, float y[], int fit_start, int fit_end,
							  float *Z, float *u, float *v, float *taup, float *taum, float *tau, float *fitted, float *residuals,
							  float *chisq);


#endif /* PHASOR_FITTING_H */

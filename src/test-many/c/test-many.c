#include <stdio.h>
#include "EcfMultiple.h"
#include <math.h>
#define NDATA 256
#define NPARAM 3
#define NROWS 2

void createExponential(float* time, float* photonCount, float a, float tau, float period, size_t count);

int main() {
	float a = 10.0;
	float tauIn = 2.5;
	float period = 0.039;
	float time[NDATA];
	float photonCount[NDATA];
	float scaledPhotonCount[NDATA];
	

	float covar_data[NPARAM * NPARAM * NROWS] = { 0 };
	float alpha_data[NPARAM * NPARAM * NROWS] = { 0 };
	float erraxes_data[NPARAM * NPARAM * NROWS] = { 0 };
	struct array3d covar = { covar_data, {NROWS, NPARAM, NPARAM}, {NPARAM * NPARAM * sizeof(float), NPARAM * sizeof(float), sizeof(float)} };
	struct array3d alpha = { NULL }; // passing null should be fine. skips copying step
	struct array3d erraxes = { NULL };

	float exponential2d[NDATA * NROWS] = { 0 };
	// fill with two experimental exponential curves
	createExponential(time, exponential2d, a, tauIn, period, NDATA);
	createExponential(time, exponential2d + NDATA, a * 2, tauIn * 2, period, NDATA); // double a and tau to see 2 different "pixels"
	struct array2d photonCount2d = { exponential2d, {NROWS, NDATA}, {NDATA * sizeof(float) , sizeof(float)} };

	float param2d_data[NPARAM * NROWS] = { 0.1, a + 0.1, tauIn + 0.1, 0.1, a + 0.1, tauIn + 0.1 }; // add 0.1 to make sure fitting works
	struct array2d param2d = { param2d_data, {NROWS, NPARAM}, {NPARAM * sizeof(float), sizeof(float)} };

	float fitted2d_data[NDATA * NROWS];
	struct array2d fitted2d = { fitted2d_data, {NROWS, NDATA}, {NDATA * sizeof(float) , sizeof(float)} };

	float residuals2d_data[NDATA * NROWS];
	struct array2d residuals2d = { residuals2d_data, {NROWS, NDATA}, {NDATA * sizeof(float) , sizeof(float)} };

	float one_arr[] = { 1.0 };
	struct array1d paramfree = {one_arr, 1, 0};

	float *chisq_data[NROWS];
	struct array1d chisq = {chisq_data, NROWS, sizeof(float)};

	struct array1d fit_mask = { one_arr, 1, 0 };

	struct common_params common_in = {period, &photonCount2d, 0, NDATA, &fitted2d, &residuals2d, &chisq, NULL};

	struct marquardt_params marquardt_in = {.instr=NULL, .noise=NOISE_POISSON_FIT, .sig=NULL, .param=&param2d, .paramfree=&paramfree, 
		.restrain=ECF_RESTRAIN_DEFAULT, .fitfunc=GCI_multiexp_tau, .covar=&covar, .alpha=&alpha, .erraxes=&erraxes, .chisq_target=1.1, .chisq_delta=1E-5, .chisq_percent=95 };

	struct flim_params flim_in = { &common_in, &marquardt_in };

	GCI_marquardt_fitting_engine_many(&flim_in);

	printf("params\n");
	for (int i = 0; i < NPARAM * NROWS; i++)
		printf("%f\n", param2d.data[i]);

	printf("\n3d covariance matrix\n");
	//for (int i = 0; i < NPARAM * NPARAM * NROWS; i++)
	//	printf("%f\n", covar.data[i]);
	print_array3d(&covar, 5);

	// what happens if we only look at every other element in photonCount?
	photonCount2d.sizes[1] = NDATA/2;
	photonCount2d.strides[1] = 2 * sizeof(float);
	common_in.xincr = period * 2;

	// also lets reverse the order of param
	float param2d_data_rev[NPARAM * NROWS] = { tauIn + 0.1, a + 0.1, 0.1, tauIn + 0.1, a + 0.1, 0.1 };
	param2d.strides[1] = -sizeof(float);
	param2d.data = param2d_data_rev + 2;

	common_in.fit_end /= 2;

	GCI_marquardt_fitting_engine_many(&flim_in);

	printf("params\n");
	for (int i = 0; i < NPARAM * NROWS; i++)
		printf("%f\n", param2d_data_rev[i]);

	print_common(&common_in);
}

void createExponential(float* time, float* photonCount, float a, float tau, float period, size_t count) {
	for (int i = 0; i < count; i++) {
		float t = period * i;
		time[i] = t;
		photonCount[i] = a * exp(-t / tau);
	}
}




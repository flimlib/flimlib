#include <stdio.h>
#include "Ecf.h"
#include <math.h>
#define NUM 256

void createExponential(float* time, float* photonCount, float a, float tau, float period, size_t count);

int main() {
	float a = 10.0;
	float tauIn = 2.5;
	float period = 0.039;
	float time[NUM];
	float photonCount[NUM];
	float scaledPhotonCount[NUM];

	struct array2d photonCount2d;
	struct array2d param2d;
	struct array2d fitted2d;
	struct array2d residuals2d;

	float covar_data[18] = { 0 };
	float alpha_data[18] = { 0 };
	float erraxes_data[18] = { 0 };
	struct array3d covar = { covar_data, 3, 3, 2 };
	struct array3d alpha = { NULL }; // passing null should be fine. skips copying step
	struct array3d erraxes = { NULL };

	float exponential2d[NUM * 2] = { 0 };
	// fill with two experimental exponential curves
	createExponential(time, exponential2d, a, tauIn, period, NUM);
	createExponential(time, exponential2d + NUM, a * 2, tauIn * 2, period, NUM);
	photonCount2d.data = exponential2d;
	photonCount2d.x_size = NUM;
	photonCount2d.y_size = 2;
	photonCount2d.x_stride_bytes = NUM * sizeof(float);

	float param2d_data[6] = { 0.1, a + 0.1, tauIn + 0.1, 0.1, a + 0.1, tauIn + 0.1 }; // add 0.1 to make sure fitting works
	param2d.data = param2d_data;
	param2d.x_size = 3;
	param2d.y_size = 2;
	param2d.x_stride_bytes = 3 * sizeof(float);

	float fitted2d_data[NUM * 2];
	fitted2d.data = fitted2d_data;
	fitted2d.x_size = NUM;
	fitted2d.y_size = 2;
	fitted2d.x_stride_bytes = NUM * sizeof(float);

	float residuals2d_data[NUM * 2];
	residuals2d.data = residuals2d_data;
	residuals2d.x_size = NUM;
	residuals2d.y_size = 2;
	residuals2d.x_stride_bytes = NUM * sizeof(float);

	float paramfree[3] = { 1, 1, 1 };

	float chisq[2];

	GCI_marquardt_fitting_engine_many(period, photonCount2d, 0, NUM, NULL, NULL, NOISE_POISSON_FIT, NULL, param2d, paramfree, 3, ECF_RESTRAIN_DEFAULT, GCI_multiexp_tau, fitted2d, residuals2d, chisq, covar, alpha, erraxes, 1.1, 1E-5, 95);

	printf("params\n");
	for (int i = 0; i < 6; i++)
		printf("%f\n", param2d_data[i]);

	printf("\n3d covariance matrix\n");
	for (int i = 0; i < 18; i++)
		printf("%f\n", covar_data[i]);
}

void createExponential(float* time, float* photonCount, float a, float tau, float period, size_t count) {
	for (int i = 0; i < count; i++) {
		float t = period * i;
		time[i] = t;
		photonCount[i] = a * exp(-t / tau);
	}
}
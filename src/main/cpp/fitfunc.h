#pragma once

typedef void (*fitfunc)(float, float [], float *, float [], int);

class FitFunc {
public:
	virtual float fit(float x, float param[], float dy_dparam[], int nparam) {
		float y = 0;
		if (this->func_ptr)
			this->func_ptr(x, param, &y, dy_dparam, nparam);
		else
			fprintf(stderr, "Warning: FitFunc->func_ptr invalid.");
		return y;
	}

	FitFunc() : func_ptr(NULL) {}

	FitFunc(fitfunc func_ptr) : func_ptr(func_ptr) {}

	virtual ~FitFunc() {}

	// needed for creating array to feed into java callback
	// int nparam;

private:
	const fitfunc func_ptr;
};

#pragma once

typedef void (*fitfunc)(float, float [], float *, float [], int);

class FitFunc {
public:
	virtual void fit(float x, float param[], float y[], float dy_dparam[]) {
		if (this->func_ptr)
			this->func_ptr(x, param, y, dy_dparam, nparam);
	}

	FitFunc() : func_ptr(NULL) {}

	FitFunc(fitfunc func_ptr) : func_ptr(func_ptr) {}

	virtual ~FitFunc() {}

	int nparam; // hidden from java side

private:
	const fitfunc func_ptr;
};

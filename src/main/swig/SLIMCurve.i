%module SLIMCurve

%{
#include <algorithm>
#include <iostream>
#include "Ecf.h"
#include "EcfGlobal.h"
#include "EcfWrapper.h"
#include "GCI_Phasor.h"
#define PKG_NAME "slim"
%}

%include "enums.swg"
// for custom typemaps
%include "SLIMCurve_1DArray.i" // arrays (with length parameter)
%include "SLIMCurve_2DMatrix.i" // 2D arrays (with/out length parameter)
%include "SLIMCurve_FittingFunc.i" // fitting function pointer
%include "SLIMCurve_DMSPVAF.i" // struct used by mode selection engine

%javaconst(1);
// rename enums to meet java naming conventions
%rename(NoiseType) noise_type;
%rename(RestrainType) restrain_type;
%rename(FitType) fit_type;
%rename(FitFunc) fitfunc;

// input 1d array maps
ARRMAP(INTARR_in, int, Int, JNI_ABORT)
ARRMAP(FLTARR_in, float, Float, JNI_ABORT)
ARRMAP(DBLARR_in, double, Double, JNI_ABORT)

// input 2d array maps
MATMAP(F2D_in, float, Float, F, Float2DMatrix)
MATMAP(I2D_in, int, Int, I, Int2DMatrix)

// Tell swig to use corresponding typemaps (OUTPUT defined in typemaps.i)
%apply int *INPUT { int* };
%apply float *OUTPUT { float * };
%apply double *OUTPUT { double * };
%apply F2D_in {
	(float **trans, int ndata, int ntrans)
}
%apply INTARR_in {
	(int paramfree[], int nparam),
	(int param_free[], int n_param)
}
%apply FLTARR_in {
	(float instr[], int ninstr),
	(float y[], int ndata)
}
%apply DBLARR_in {
	(double instr[], int n_instr)
}

// Grab functions from header files as class methods
%include "../c/Ecf.h"
%include "../c/EcfGlobal.h"
%include "../c/EcfWrapper.h"
%include "../c/GCI_Phasor.h"

%pragma(java) jniclasscode=%{
	static {
		
		try {
			System.loadLibrary("slim-curve");
			System.loadLibrary("slim-curve-java");
			//NativeLoader.loadLibrary("slim-curve");
			//NativeLoader.loadLibrary("slim-curve-java");
		} catch (Exception e) {
			//System.err.println("Native library failed to load. Exiting.\n" + e);
			System.err.println("Cannot extract native library. Exiting.\n" + e);
			System.exit(1);
		}
	}
%}

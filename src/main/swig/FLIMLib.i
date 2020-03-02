%module(directors="1") FLIMLib

%{
#include <algorithm>
#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

#include "Ecf.h"
#include "EcfGlobal.h"
#include "EcfWrapper.h"
#include "GCI_Phasor.h"
#include "BayesAnalysis.h"

#ifdef __cplusplus
}
#endif

#define PKG_NAME "flimlib"
%}

// custom typemaps
%include "FLIMLib_1DArray.i" // arrays (with length parameter)
%include "FLIMLib_2DMatrix.i" // 2D arrays (with/out length parameter)
%include "FLIMLib_FittingFunc.i" // fitting function pointer
%include "FLIMLib_DMSPVAF.i" // struct used by mode selection engine
%include "FLIMLib_Enums.i" // all of the enums
%include "FLIMLib_ProgressFunc.i" // progress function for SPA

%javaconst(1);
// rename enums to meet java naming conventions
ENUMMAP(FITTYPEENUM, fit_type, FitType);
%rename(NoiseType) noise_type;
%rename(RestrainType) restrain_type;

// input 1d array (with length) maps
ARRMAP(BOLARRIN_LEN, 1, 1, boolean, Boolean, JNI_ABORT, false)
ARRMAP(FLTARRIN_LEN, 1, 1, float, Float, JNI_ABORT, false)
ARRMAP(FLTPTRIN_LEN, 0, 1, float, Float, JNI_ABORT, false)
ARRMAP(FLTARRIN_LUL, 1, 1, float, Float, JNI_ABORT, true)
ARRMAP(FLTARRIN_NUL, 1, 0, float, Float, JNI_ABORT, true)
ARRMAP(DBLARRIN_LUL, 1, 1, double, Double, JNI_ABORT, true)
ARRMAP(DBLARRIN_NUL, 1, 0, double, Double, JNI_ABORT, true)

// input 2d array maps
MATMAP(F2D_in, float, Float, F, Float2DMatrix)
MATMAP(I2D_in, int, Int, I, Int2DMatrix)

// Tell swig to use corresponding typemaps (OUTPUT defined in typemaps.i)
%apply int *INOUT { int* };
%apply float *INOUT { float *Z,  float *A, float *tau, float *residuals, float *fitted, float *error };
%apply int *OUTPUT { int *nphotons };
%apply float *OUTPUT { float * };
%apply double *INOUT { double * };
%apply FITTYPEENUM { int ftype };
%apply F2D_in {
	(float **trans, int ndata, int ntrans)
}
%apply BOLARRIN_LEN {
	(int paramfree[], int nparam),
	(int paramfree[], int nparamfree),
	(int param_free[], int n_param),
	(int restrain[], int nparam)
}
%apply FLTARRIN_LEN {
	(float params[], int nparam),
	(float y[], int ndata)
}
%apply FLTPTRIN_LEN {
	(float *trans, int ndata)
}
%apply FLTARRIN_LUL {
	(float instr[], int ninstr)
}
%apply FLTARRIN_NUL {
	float sig[]
}
%apply DBLARRIN_LUL {
	(double instr[], int n_instr)
}
%apply DBLARRIN_NUL {
	double sig[]
}

// Grab functions from header files as class methods
%include "../c/Ecf.h"
%include "../c/EcfGlobal.h"
%include "../c/EcfWrapper.h"
%include "../c/GCI_Phasor.h"
%include "../c/BayesAnalysis.h"

%pragma(java) jniclassimports=%{
  import org.scijava.nativelib.NativeLoader;
  import java.io.IOException;
%}

%pragma(java) jniclasscode=%{
	static {
		try {
			NativeLoader.loadLibrary("flimlib");
			NativeLoader.loadLibrary("flimlib-jni");
		} catch (IOException e) {
			throw new ExceptionInInitializerError(e);
		}
	}
%}

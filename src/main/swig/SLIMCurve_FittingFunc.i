%module(directors="1") SLIMCurve

%rename(FitFuncNative) FitFunc;
%typemap(javaclassmodifiers) FitFunc "class"
%typemap(javainterfaces) FitFunc "FitFunc"

// ignore global variables
%rename("$ignore", regextarget=1) "^t_";
%rename("$ignore", regextarget=1) "^do_";
%rename("$ignore", regextarget=1) "^NTV_";
%rename("$ignore", regextarget=1) "^GCI_(multi|stretched)exp";
%ignore FitFunc::FitFunc(fitfunc func_ptr);
%ignore FitFunc::func_ptr;
%ignore FitFunc::nparam;

%javaconst(0);

// Builtin fitting functions as constants on java side
%constant FitFunc *GCI_MULTIEXP_LAMBDA = &NTV_GCI_MULTIEXP_LAMBDA;
%constant FitFunc *GCI_MULTIEXP_TAU = &NTV_GCI_MULTIEXP_TAU;
%constant FitFunc *GCI_STRETCHEDEXP = &NTV_GCI_STRETCHEDEXP;

%inline %{
#include "fitfunc.h"

// The thread-local global variable is defined so that multiple fitting
// functions can coexist in different threads safely.
#ifndef THREAD_FITFUNC
#define THREAD_FITFUNC
// At any given time, no more than two c instance of FitFunc should exist
// (only one of them should be stored at [0] here just before its execution).
// The only exception is .
static thread_local FitFunc *t_fitfunc[2] = { 0 };
static thread_local jfloatArray t_param = NULL;
static thread_local jfloatArray t_dy_dparam = NULL;
#endif

// Static fitting worker for executing the cached fitfunc
static void do_fit0(float x, float param[], float *y, float dy_dparam[], int nparam) {
	*y = t_fitfunc[0]->fit(x, param, dy_dparam, nparam);
}
// Used in GCI_EcfModelSelectionEngine. There is probably no way of knowing
// the FitFunc of which model is being called during the execution, so it is
// hard-coded.
static void do_fit1(float x, float param[], float *y, float dy_dparam[], int nparam) {
	*y = t_fitfunc[1]->fit(x, param, dy_dparam, nparam);
}

const FitFunc& NTV_GCI_MULTIEXP_LAMBDA = FitFunc(GCI_multiexp_lambda);
const FitFunc& NTV_GCI_MULTIEXP_TAU = FitFunc(GCI_multiexp_tau);
const FitFunc& NTV_GCI_STRETCHEDEXP = FitFunc(GCI_stretchedexp);
%}

///////////////////////////////////////////////////////////////////////////////

%define fitfunc_ptr
void (*fitfunc)(float, float [], float *, float [], int)
%enddef
%typemap(jstype) fitfunc_ptr "FitFunc"
%typemap(javain,pgcppname="n",
		pre="    FitFuncNative n = FitFuncNative.makeNative($javainput);")
		fitfunc_ptr  "FitFuncNative.getCPtr(n)"
%typemap(jtype) fitfunc_ptr "long"
%typemap(jni) fitfunc_ptr "jlong"
%typemap(in) fitfunc_ptr %{
	t_fitfunc[0] = (FitFunc *) $input;
	// feed the static fitting function instead
	$1 = &do_fit0;
%}
%typemap(freearg) fitfunc_ptr %{
	// in case of reuse
	t_fitfunc[0] = t_fitfunc[1] = NULL;
%}

%define fitfunc_ref FitFunc & %enddef
%clear fitfunc_ref;
%typemap(jstype) fitfunc_ref "FitFunc"
%typemap(javain,pgcppname="n",
		pre="    FitFuncNative n = FitFuncNative.makeNative($javainput);")
		fitfunc_ref  "FitFuncNative.getCPtr(n)"
%typemap(jtype) fitfunc_ref "long"
%typemap(jni) fitfunc_ref "jlong"
%typemap(in) fitfunc_ref %{
	$1 = (FitFunc *)$input;
%}
%typemap(freearg) fitfunc_ref %{
%}
// TODO: out

// directorin copies c array into java array before the jni call
// directorout copies back to c array after that
%typemap(directorin,descriptor="[F") float param[] {
	// don't bother allocating a new array if length is the same
	// $1_name = param or dy_dparam depending on which is being mapped
	if (!t_$1_name || jenv->GetArrayLength(t_$1_name) != nparam) {
		if (t_$1_name)
			jenv->DeleteGlobalRef(t_$1_name);
		t_$1_name = (jfloatArray)jenv->NewGlobalRef(jenv->NewFloatArray(nparam));
	}
	$input = t_$1_name;
	jenv->SetFloatArrayRegion($input, 0, nparam, $1);
}
%typemap(directorin,descriptor="[F") (float dy_dparam[], int nparam) {
	// don't bother allocating a new array if length is the same
	// $1_name = param or dy_dparam depending on which is being mapped
	if (!t_$1_name || jenv->GetArrayLength(t_$1_name) != nparam) {
		if (t_$1_name)
			jenv->DeleteGlobalRef(t_$1_name);
		t_$1_name = (jfloatArray)jenv->NewGlobalRef(jenv->NewFloatArray(nparam));
	}
	$input = t_$1_name;
	jenv->SetFloatArrayRegion($input, 0, nparam, $1);
}
%typemap(directorargout) (float dy_dparam[], int nparam) {
	jsize sz = jenv->GetArrayLength($input);
	jenv->GetFloatArrayRegion($input, 0, nparam, $1);
}

%feature("director") FitFunc;

// make sure swig knows float[] type and don't mess up in SLIMCurveJNI
ARRMAP(FLTARRIN, 1, 0, float, Float, 0, false)
%apply FLTARRIN { float param[], float dy_dparam[] };
// override
%typemap(in) float dy_dparam[] (jfloat *jarr, jfloatArray, bool do_clean) %{
#define MAP_1D_ARR_dy_dparam
	// local reference to the java array and the length of it
	jsize len_dy_dparam;
	do_clean = true;

	len_dy_dparam = (jsize) jenv->GetArrayLength($input);
	if (!SWIG_JavaArrayInFloat(jenv, &jarr, (float **)&$1, $input))
		exit(1);
	t_fitfunc[0] = (FitFunc*) jarg1;
	if (len_param != len_dy_dparam)
		SWIG_JavaThrowException(jenv, SWIG_JavaIllegalArgumentException,
			"param number does not match dy_dparam number");
%}

%typemap(javacode) FitFunc %{

	private static class FitFuncNativeProxy extends FitFuncNative {
		private FitFunc delegate;
		public FitFuncNativeProxy(FitFunc i) {
			delegate = i;
		}

		public float fit(final float x, final float[] param, final float[] dy_dparam) {
			return delegate.fit(x, param, dy_dparam);
		}
	}

	static FitFuncNative makeNative(FitFunc i) {
		if (i instanceof FitFuncNative) {
			// If it already *is* a FitFuncNative don't bother wrapping it again
			return (FitFuncNative)i;
		}
		return new FitFuncNativeProxy(i);
	}

	public float fit(float x, float[] param, float[] dy_dparam) {
		return fit(x, param, dy_dparam, param.length);
	}
%}

%include "../cpp/fitfunc.h"

%module SLIMCurve

%include "enums.swg"
%javaconst(1);

%rename(FitFunc) fit_funcs;

%{
typedef enum { GCI_MULTIEXP_LAMBDA, GCI_MULTIEXP_TAU, GCI_STRETCHEDEXP } fit_funcs;
%}

// Conversion: FitFunc(J) -> void (*fitfunc)(float, float [], float *, float [], int)(C) in arguments
%define fFunction
void (*fitfunc)(float, float [], float *, float [], int)
%enddef
%typemap(jstype) fFunction "FitFunc"
%typemap(javain) fFunction "$javainput.swigValue()"
%typemap(jtype) fFunction "int"
%typemap(jni) fFunction "jint"
%typemap(in) fFunction {
	switch($input) {
	case GCI_MULTIEXP_LAMBDA:
		$1 = GCI_multiexp_lambda;
		break;
	case GCI_MULTIEXP_TAU:
		$1 = GCI_multiexp_tau;
		break;
	case GCI_STRETCHEDEXP:
		$1 = GCI_stretchedexp;
		break;
	}
}

%ignore GCI_multiexp_lambda;
%ignore GCI_multiexp_tau;
%ignore GCI_stretchedexp;

typedef enum { GCI_MULTIEXP_LAMBDA, GCI_MULTIEXP_TAU, GCI_STRETCHEDEXP } fit_funcs;

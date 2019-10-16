%module FLIMLib

%ignore _DMSPVAF;
%ignore DecayModelSelParamValuesAndFit;

// let swig know about the struct before defining DecayModel
typedef struct _DMSPVAF {
    void          (*fitfunc)(float, float [], float *, float [], int);
    int             nparam;
    float           params[MAXFIT];
    int             nparamfree;
    int             paramfree[MAXFIT];
    restrain_type   restrain;
    float          *fitted;
    float          *residuals;
    float           chisq_target;
    float           chisq_delta;
    int             chisq_percent;
    float           chisq;
    float         **covar;
    float         **alpha;
    float         **erraxes;
} DecayModelSelParamValuesAndFit;

%inline %{
#include "decmod.h"

// The thread-local global variable is defined so that multiple fitting
// functions can coexist in different threads safely.
#ifndef THREAD_FITFUNC
#define THREAD_FITFUNC
// At any given time, no more than two c instance of FitFunc should exist
// (only one of them should be stored at [0] here just before its execution).
// The only exception is .
static thread_local FitFunc *t_fitfunc[2] = { 0 };
#endif
%}

// numinputs=0 ignores the jni argument.
// this map feeds the constructor the JNI environment
%typemap(in, numinputs=0) JNIEnv * {
	$1 = jenv;
}

// Conversion: DecayModelSelParamValuesAndFit[2](J) -> (C) */
// %define DMSPVAF _DMSPVAF *paramsandfits %enddef
%define DMSPVAF DecayModelSelParamValuesAndFit* paramsandfits %enddef
%typemap(jstype) DMSPVAF "DecayModel[]"
%typemap(javain, pre="
	if($javainput.length != 2)
		throw new IllegalArgumentException(\"Requires 2 models.\");
	/*if($javainput[0].params.length > 20 || $javainput[1].params.length > 20)
		throw new IllegalArgumentException(\"More than 20 parameters.\");
	if($javainput[0].paramfree.length != $javainput[0].paramfree.length ||
		$javainput[1].paramfree.length != $javainput[1].paramfree.length)
		throw new IllegalArgumentException(\"Lengths of param and paramfree disagree.\");*/
	if($javainput[0] == null || $javainput[1] == null)
		throw new NullPointerException(\"Array contains null\");
") DMSPVAF "$javainput"
%typemap(jtype) DMSPVAF "DecayModel[]"
%typemap(jni) DMSPVAF "jobjectArray"
%typemap(in) DMSPVAF (DecayModel *decmods[2]) {
	// $input = DecayModel[]
	// $1 = DecayModelSelParamValuesAndFit*
	// This conversion recovers the objects from $input.cPtr (jlong)
	jclass DMSClass = JCALL1(FindClass, jenv, PKG_NAME"/DecayModel");
	jfieldID ptrID = JCALL3(GetFieldID, jenv, DMSClass, "swigCPtr", "J");
	jobject dms_0 = JCALL2(GetObjectArrayElement, jenv, $input, 0);
	decmods[0] = (DecayModel*)JCALL2(GetLongField, jenv, dms_0, ptrID);
	jobject dms_1 = JCALL2(GetObjectArrayElement, jenv, $input, 1);
	decmods[1] = (DecayModel*)JCALL2(GetLongField, jenv, dms_1, ptrID);

	// fill in the fitfunc cache and DecayModelSelParamValuesAndFit's fitfunc pointer
	t_fitfunc[0] = decmods[0]->getFitfunc();
	decmods[0]->fitfunc = &do_fit0;
	t_fitfunc[1] = decmods[1]->getFitfunc();
	decmods[1]->fitfunc = &do_fit1;

	// index starts at 1 in GCI_EcfModelSelectionEngine
	$1 = new DecayModelSelParamValuesAndFit[3];
	$1[1] = *dynamic_cast<DecayModelSelParamValuesAndFit*>(decmods[0]);
	$1[2] = *dynamic_cast<DecayModelSelParamValuesAndFit*>(decmods[1]);
}

%typemap(freearg) DMSPVAF {
	// copy back results
	std::copy($1[1].params, $1[1].params + $1[1].nparam, decmods$argnum[0]->params);
	std::copy($1[2].params, $1[2].params + $1[2].nparam, decmods$argnum[1]->params);
	decmods$argnum[0]->chisq = $1[1].chisq;
	decmods$argnum[1]->chisq = $1[2].chisq;
	// in case of reuse
	t_fitfunc[0] = t_fitfunc[1] = 0;

	delete[] $1;
}

%include "../cpp/decmod.h"

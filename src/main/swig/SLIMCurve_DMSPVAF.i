%module SLIMCurve

// struct name of DecayModelSelParamValuesAndFit
%extend _DMSPVAF {
	// Add constructor
	_DMSPVAF(
		void (*fitfunc)(float, float [], float *, float [], int),
		int nparam,
		float params[],
		int nparamfree,
		int paramfree[],
		restrain_type restrain,
		float *fitted,
		float *residuals,
		float chisqtarget,
		float chisqdelta,
		int chisqpercent,
		float chisq,
		float **covar,
		float **alpha,
		float **erraxes) {
		DecayModelSelParamValuesAndFit* ret = new DecayModelSelParamValuesAndFit();
		ret->fitfunc = fitfunc;
		ret->nparam = nparam;
		std::copy(params, params + nparam, ret->params);
		ret->nparamfree = nparamfree;
		std::copy(paramfree, paramfree + nparam, ret->paramfree);
		ret->restrain = restrain;
		ret->fitted = fitted;
		ret->residuals = residuals;
		ret->chisq_target = chisqtarget;
		ret->chisq_delta = chisqdelta;
		ret->chisq_percent = chisqpercent;
		ret->chisq = chisq;
		ret->covar = covar;
		ret->alpha = alpha;
		ret->erraxes = erraxes;
		return ret;
	}
}

// Conversion: DecayModelSelParamValuesAndFit[2](J) -> (C) */
%define DMSPVAF DecayModelSelParamValuesAndFit* paramsandfits %enddef
%typemap(jstype) DMSPVAF "DecayModelSelParamValuesAndFit[]"
%typemap(javain, pre="
	if($javainput.length != 2)
		throw new IllegalArgumentException(\"Requires 2 models.\");
	if($javainput[0] == null || $javainput[1] == null)
		throw new NullPointerException(\"Array contains null\");
") DMSPVAF "$javainput"
%typemap(jtype) DMSPVAF "DecayModelSelParamValuesAndFit[]"
%typemap(jni) DMSPVAF "jobjectArray"
%typemap(in) DMSPVAF {
	// $input = DecayModelSelParamValuesAndFit[]
	// $1 = DecayModelSelParamValuesAndFit*
	// This conversion recovers the objects from $input.cPtr (jlong)
	$1 = new DecayModelSelParamValuesAndFit[2];
	jclass DMSClass = JCALL1(FindClass, jenv, PKG_NAME"/DecayModelSelParamValuesAndFit");
	jfieldID ptrID = JCALL3(GetFieldID, jenv, DMSClass, "cPtr", "J");
	jobject dms_0 = JCALL2(GetObjectArrayElement, jenv, $input, 0);
	$1[0] = *(DecayModelSelParamValuesAndFit*)JCALL2(GetLongField, jenv, dms_0, ptrID);
	jobject dms_1 = JCALL2(GetObjectArrayElement, jenv, $input, 1);
	$1[1] = *(DecayModelSelParamValuesAndFit*)JCALL2(GetLongField, jenv, dms_1, ptrID);
}

// do not generate getter/setter for members in DecayModelSelParamValuesAndFit (one-time use)
%rename("$ignore", fullname=1, regextarget=1, %$isvariable)"_DMSPVAF::.*";
%ignore _DMSPVAF::fitfunc; // not sure why I have to ignore it again

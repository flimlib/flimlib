/*-
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2022 University of Oxford and Board of Regents of the
 * University of Wisconsin-Madison.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
#pragma once
#include <jni.h>
#include "Ecf.h"

// the same version as in swig-generated Director class
#define GETJNIENV JNIEnv *jenv; jvm->GetEnv((void**)&jenv, JNI_VERSION_1_2);

/* 
 * This is a wrapper for DecayModelSelParamValuesAndFit struct.
 * All the paramenter members can only be set in constructor, and all the "output"
 * members can only be retrieved.
 */
struct DecayModel : public _DMSPVAF {
public:
	FitFunc *getFitfunc() { return this->fitfuncobj; }

	jfloatArray getParams();

	jfloatArray getFitted();

	jfloatArray getResiduals();

	float getChisq() { return this->chisq; }

	ParamMatrix<float> *getCovar() { return this->covarobj; }

	ParamMatrix<float> *getAlpha() { return this->alphaobj; }

	ParamMatrix<float> *getErraxes() { return this->erraxesobj; }

	DecayModel() : _DMSPVAF() {}

	DecayModel(JNIEnv *jenv, FitFunc &fitfunc, jfloatArray params, jbooleanArray paramfree, restrain_type restrain,
		int ndata, float chisq_target, float chisq_delta, int chisq_percent);

	~DecayModel();

private:
	// cached for getMethods and desctructor
	JavaVM *jvm;

	// corresponding objects/references
	FitFunc *fitfuncobj;
	int ndata;
	jfloatArray paramsref;
	jfloatArray fittedref;
	jfloatArray residualsref;
	ParamMatrix<float> *covarobj;
	ParamMatrix<float> *alphaobj;
	ParamMatrix<float> *erraxesobj;
};

jfloatArray DecayModel::getParams() {
	GETJNIENV;
	// update the java array object
	// this method is different from the other two below
	// because param lives within a Decmod whereas fitted and residuals
	// are allocated
	jenv->SetFloatArrayRegion(paramsref, 0, nparam, params);
	return paramsref;
}

jfloatArray DecayModel::getFitted() {
	GETJNIENV;
	// update the java array object
	jenv->ReleaseFloatArrayElements(fittedref, fitted, JNI_COMMIT);
	return fittedref;
}

jfloatArray DecayModel::getResiduals() {
	GETJNIENV;
	// update the java array object
	jenv->ReleaseFloatArrayElements(residualsref, residuals, JNI_COMMIT);
	return residualsref;
}

DecayModel::DecayModel(JNIEnv *jenv, FitFunc &fitfunc, jfloatArray params, jbooleanArray paramfree, restrain_type restrain,
	int ndata, float chisq_target, float chisq_delta, int chisq_percent) : _DMSPVAF() {
	// save for destructor and getFitted/getResiduals
	jenv->GetJavaVM(&this->jvm);
	// fitfunc (do_fit<0/1>) will be determined at fit-time
	this->fitfuncobj = &fitfunc; this->fitfunc = NULL;
	this->nparam = std::min<jsize>(jenv->GetArrayLength(params), MAXFIT);
	// create a GlobalRef so that the array cannot be GC-ed
	this->paramsref = (jfloatArray)jenv->NewGlobalRef(jenv->NewFloatArray(nparam));
	jenv->GetFloatArrayRegion(params, 0, nparam, (jfloat*)this->params);
	bool paramfree_tmp[20];
	jenv->GetBooleanArrayRegion(paramfree, 0, nparam, (jboolean *)paramfree_tmp);
	this->nparamfree = 0;
	for (int i = 0; i < nparam; i++) {
		this->paramfree[i] = (paramfree_tmp[i] ? 1 : 0);
		this->nparamfree += this->paramfree[i];
	}
	this->restrain = restrain;
	this->ndata = ndata;
	this->fittedref = (jfloatArray)jenv->NewGlobalRef(jenv->NewFloatArray(ndata));
	this->fitted = jenv->GetFloatArrayElements(fittedref, NULL);
	this->residualsref = (jfloatArray)jenv->NewGlobalRef(jenv->NewFloatArray(ndata));
	this->residuals = jenv->GetFloatArrayElements(residualsref, NULL);
	this->chisq_target = chisq_target;
	this->chisq_delta = chisq_delta;
	this->chisq_percent = chisq_percent;
	this->covarobj = new Float2DMatrix(nparam, nparam); covar = covarobj->arr;
	this->alphaobj = new Float2DMatrix(nparam, nparam); alpha = alphaobj->arr;
	this->erraxesobj = new Float2DMatrix(nparam, nparam); erraxes = erraxesobj->arr;
}

DecayModel::~DecayModel() {
	GETJNIENV;
	// after GlobalRef is stripped, the refs will be GC-ed if no java ref to them exists
	if (!jenv->IsSameObject(paramsref, NULL))
		jenv->DeleteGlobalRef(paramsref);
	if (!jenv->IsSameObject(fittedref, NULL)) {
		jenv->ReleaseFloatArrayElements(fittedref, fitted, JNI_ABORT);
		jenv->DeleteGlobalRef(fittedref);
	}
	if (!jenv->IsSameObject(residualsref, NULL)) {
		jenv->ReleaseFloatArrayElements(residualsref, residuals, JNI_ABORT);
		jenv->DeleteGlobalRef(residualsref);
	}
	delete covarobj;
	delete alphaobj;
	delete erraxesobj;
}

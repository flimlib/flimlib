/*
 * #%L
 * SLIM Curve package for exponential curve fitting of spectral lifetime data.
 * %%
 * Copyright (C) 2010 - 2014 Gray Institute University of Oxford and Board of
 * Regents of the University of Wisconsin-Madison.
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

#include "loci_slim_SLIMCurve.h"
#include "EcfWrapper.h"

/*
 * Class:     loci_slim_SLIMCurve
 * Method:    RLD_fit
 * Signature: (D[DII[DII[D[D[D[D[D[DD)I
 */
JNIEXPORT jint JNICALL Java_loci_slim_SLIMCurve_RLD_1fit
  (JNIEnv *env, jobject obj, jdouble x_inc, jdoubleArray y,
        jint fit_start, jint fit_end, jdoubleArray instr, jint n_instr,
        jint noise, jdoubleArray sig,
        jdoubleArray z, jdoubleArray a, jdoubleArray tau,
        jdoubleArray fitted, jdoubleArray chi_square, jdouble chi_square_target) {

    jdouble *y_array;
    jdouble *instr_array;
    jdouble *sig_array;
    jdouble *fitted_array;
    jdouble *z_ref;
    jdouble *a_ref;
    jdouble *tau_ref;
    jdouble *chi_square_ref;
    int return_value;

    // convert the arrays
    y_array = (*env)->GetDoubleArrayElements(env, y, 0);
    instr_array = NULL;
    if (NULL != instr) {
        instr_array  = (*env)->GetDoubleArrayElements(env, instr, 0);
    }
    sig_array = NULL;
    if (NULL != sig) {
        sig_array    = (*env)->GetDoubleArrayElements(env, sig, 0);
    }
    fitted_array = (*env)->GetDoubleArrayElements(env, fitted, 0);

    // also handle double by reference as arrays
    z_ref = (*env)->GetDoubleArrayElements(env, z, 0);
    a_ref = (*env)->GetDoubleArrayElements(env, a, 0);
    tau_ref = (*env)->GetDoubleArrayElements(env, tau, 0);
    chi_square_ref = (*env)->GetDoubleArrayElements(env, chi_square, 0);

    return_value = RLD_fit(x_inc, y_array, fit_start, fit_end,
            instr_array, n_instr, noise, sig_array, z_ref, a_ref, tau_ref,
            fitted_array, chi_square_ref, chi_square_target);

    // pass back the arrays
    (*env)->ReleaseDoubleArrayElements(env, y, y_array, 0);
    if (NULL != instr) {
        (*env)->ReleaseDoubleArrayElements(env, instr, instr_array, 0);
    }
    if (NULL != sig) {
        (*env)->ReleaseDoubleArrayElements(env, sig, sig_array, 0);
    }
    (*env)->ReleaseDoubleArrayElements(env, fitted, fitted_array, 0);

    // pass back the double by reference values
    (*env)->ReleaseDoubleArrayElements(env, z, z_ref, 0);
    (*env)->ReleaseDoubleArrayElements(env, a, a_ref, 0);
    (*env)->ReleaseDoubleArrayElements(env, tau, tau_ref, 0);
    (*env)->ReleaseDoubleArrayElements(env, chi_square, chi_square_ref, 0);

    return return_value;
}

/*
 * Class:     loci_slim_SLIMCurve
 * Method:    LMA_fit
 * Signature: (D[DII[DII[D[D[II[D[DDD)I
 */
JNIEXPORT jint JNICALL Java_loci_slim_SLIMCurve_LMA_1fit
  (JNIEnv *env, jobject obj, jdouble x_inc, jdoubleArray y,
        jint fit_start, jint fit_end, jdoubleArray instr, jint n_instr,
        jint noise, jdoubleArray sig,
        jdoubleArray param, jintArray param_free, jint n_param,
        jdoubleArray fitted, jdoubleArray chi_square,
        jdouble chi_square_target, jdouble chi_square_delta) {

    jdouble *y_array;
    jdouble *instr_array;
    jdouble *sig_array;
    jdouble *param_array;
    jint *param_free_array;
    jdouble *fitted_array;
    jdouble *chi_square_ref;
    int return_value;

    // convert the arrays
    y_array = (*env)->GetDoubleArrayElements(env, y, 0);
    instr_array = NULL;
    if (NULL != instr) {
        instr_array  = (*env)->GetDoubleArrayElements(env, instr, 0);
    }
    sig_array = NULL;
    if (NULL != sig) {
        sig_array    = (*env)->GetDoubleArrayElements(env, sig, 0);
    }
    param_array = (*env)->GetDoubleArrayElements(env, param, 0);
    param_free_array = (*env)->GetIntArrayElements(env, param_free, 0);
    fitted_array = (*env)->GetDoubleArrayElements(env, fitted, 0);

    // also handle double by reference as array
    chi_square_ref = (*env)->GetDoubleArrayElements(env, chi_square, 0);

    return_value = LMA_fit(x_inc, y_array, fit_start, fit_end,
            instr_array, n_instr, noise, sig_array,
            param_array, param_free_array, n_param,
            fitted_array, chi_square_ref, chi_square_target, chi_square_delta);

    // pass back the arrays
    (*env)->ReleaseDoubleArrayElements(env, y, y_array, 0);
    if (NULL != instr) {
        (*env)->ReleaseDoubleArrayElements(env, instr, instr_array, 0);
    }
    if (NULL != sig) {
        (*env)->ReleaseDoubleArrayElements(env, sig, sig_array, 0);
    }
    (*env)->ReleaseDoubleArrayElements(env, param, param_array, 0);
    (*env)->ReleaseIntArrayElements(env, param_free, param_free_array, 0);
    (*env)->ReleaseDoubleArrayElements(env, fitted, fitted_array, 0);

    // pass back the double by reference value
    (*env)->ReleaseDoubleArrayElements(env, chi_square, chi_square_ref, 0);

    return return_value;
  }

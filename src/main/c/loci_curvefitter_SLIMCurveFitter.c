#include "loci_curvefitter_SLIMCurveFitter.h"
#include "EcfWrapper.h"

/*
 * Class:     loci_curvefitter_SLIMCurveFitter
 * Method:    RLD_fit
 * Signature: (D[DII[DI[D[D[D[D[D[DD)I
 */
JNIEXPORT jint JNICALL Java_loci_curvefitter_SLIMCurveFitter_RLD_1fit
  (JNIEnv *env, jobject obj, jdouble x_inc, jdoubleArray y,
        jint fit_start, jint fit_end, jdoubleArray instr, jint n_instr,
        jdoubleArray sig, jdoubleArray z, jdoubleArray a, jdoubleArray tau,
        jdoubleArray fitted, jdoubleArray chi_square, jdouble chi_square_target) {

    // convert the arrays
    jdouble *y_array = (*env)->GetDoubleArrayElements(env, y, 0);
    jdouble *instr_array = NULL;
    if (NULL != instr) {
        instr_array  = (*env)->GetDoubleArrayElements(env, instr, 0);
    }
    jdouble *sig_array = NULL;
    if (NULL != sig) {
        sig_array    = (*env)->GetDoubleArrayElements(env, sig, 0);
    }
    jdouble *fitted_array = (*env)->GetDoubleArrayElements(env, fitted, 0);

    // also handle double by reference as arrays
    jdouble *z_ref = (*env)->GetDoubleArrayElements(env, z, 0);
    jdouble *a_ref = (*env)->GetDoubleArrayElements(env, z, 0);
    jdouble *tau_ref = (*env)->GetDoubleArrayElements(env, tau, 0);
    jdouble *chi_square_ref = (*env)->GetDoubleArrayElements(env, chi_square, 0);

    int return_value = RLD_fit(x_inc, y_array, fit_start, fit_end,
            instr_array, n_instr, sig_array, z_ref, a_ref, tau_ref,
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
 * Class:     loci_curvefitter_SLIMCurveFitter
 * Method:    LMA_fit
 * Signature: (D[DII[DI[D[D[II[D[DD)I
 */
JNIEXPORT jint JNICALL Java_loci_curvefitter_SLIMCurveFitter_LMA_1fit
  (JNIEnv *env, jobject obj, jdouble x_inc, jdoubleArray y,
        jint fit_start, jint fit_end, jdoubleArray instr, jint n_instr,
        jdoubleArray sig, jdoubleArray param, jintArray param_free, jint n_param,
        jdoubleArray fitted, jdoubleArray chi_square, jdouble chi_square_target) {

    // convert the arrays
    jdouble *y_array = (*env)->GetDoubleArrayElements(env, y, 0);
    jdouble *instr_array = NULL;
    if (NULL != instr) {
        instr_array  = (*env)->GetDoubleArrayElements(env, instr, 0);
    }
    jdouble *sig_array = NULL;
    if (NULL != sig) {
        sig_array    = (*env)->GetDoubleArrayElements(env, sig, 0);
    }
    jdouble *param_array = (*env)->GetDoubleArrayElements(env, param, 0);
    jint *param_free_array = (*env)->GetIntArrayElements(env, param_free, 0);
    jdouble *fitted_array = (*env)->GetDoubleArrayElements(env, fitted, 0);

    // also handle double by reference as array
    jdouble *chi_square_ref = (*env)->GetDoubleArrayElements(env, chi_square, 0);

    int return_value = LMA_fit(x_inc, y_array, fit_start, fit_end,
            instr_array, n_instr, sig_array, param_array, param_free_array, n_param,
            fitted_array, chi_square_ref, chi_square_target);

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

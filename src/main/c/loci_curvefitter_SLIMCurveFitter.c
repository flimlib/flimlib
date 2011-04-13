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

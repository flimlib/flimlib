import ctypes
import numpy as np
import math
import warnings

_flimlib = ctypes.cdll.flimlib

_GCI_triple_integral_fitting_engine = _flimlib.GCI_triple_integral_fitting_engine

_noise_types = {'NOISE_CONST': 0, 'NOISE_GIVEN': 1, 'NOISE_POISSON_DATA': 2,
                'NOISE_POISSON_FIT': 3, 'NOISE_GAUSSIAN_FIT': 4, 'NOISE_MLE': 5}  # 4 and 5 should be an error

def _prep_sig(noise_type, sig):
    """
    Helper function to make sure noise_type is valid and sig is the correct type 
    for the noise_type given
    """
    if not noise_type in _noise_types.keys():
        raise ValueError(
            "invalid noise type. The valid types are: ", _noise_types.keys)
    elif noise_type == 'NOISE_GIVEN':
        sig = np.asarray(sig, dtype=np.float32)
        if sig.ndim != 1:
            raise ValueError("sig must be 1 dimensional")
        sig = np.ctypeslib.as_ctypes(sig)
    elif noise_type == 'NOISE_CONST':
        sig = float(np.asarray(sig))  # convert to float
        sig = ctypes.c_float(sig)
    elif sig is not None:
        warnings.warn("Expected sig=None for noise type " +
                      str(noise_type) + ". The given value of sig will be ignored")
        return None
    return sig

def _prep_instr_ninstr(instr):
    ninstr = 0
    if instr is not None:
        instr = np.asarray(instr, dtype=np.float32)
        if instr.ndim != 1:
            raise ValueError("instr must be 1 dimensional")
        ninstr = instr.shape[0] # presumably shorter than photon_count
        instr = np.ctypeslib.as_ctypes(instr)
    return instr, ninstr

def _prep_common_fit_params(photon_count):
    photon_count = np.asarray(photon_count, dtype=np.float32)
    if photon_count.ndim != 1:
        raise ValueError("photon_count must be a 1 dimensional")
    
    fit_start = 0
    fit_end = photon_count.shape[0]-1

    photon_count = np.ctypeslib.as_ctypes(photon_count)
    fitted = np.ctypeslib.as_ctypes(np.empty(fit_end, dtype=np.float32))
    residuals = np.ctypeslib.as_ctypes(np.empty(fit_end, dtype=np.float32))

    return photon_count, fit_start, fit_end, fitted, residuals

# C header

# int GCI_triple_integral_fitting_engine(float xincr, float y[], int fit_start, int fit_end,
#                                        float instr[], int ninstr, noise_type noise, float sig[],
#                                        float *Z, float *A, float *tau, float *fitted, float *residuals,
#                                        float *chisq, float chisq_target);

_GCI_triple_integral_fitting_engine.argtypes = [
    ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_int,
    ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_int,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float), ctypes.c_float]


def GCI_triple_integral_fitting_engine(period, photon_count,
                                       instr=None, noise_type='NOISE_CONST', sig=1.0,
                                       chisq_target=1.1, output_fitted_and_residuals=False):
    """
    put documentation here
    """
    period = ctypes.c_float(period)

    photon_count, fit_start, fit_end, fitted, residuals = _prep_common_fit_params(photon_count)

    instr, ninstr = _prep_instr_ninstr(instr)

    # this lack of implementation is unique to triple integral
    if noise_type == 'NOISE_GAUSSIAN_FIT' or noise_type == 'NOISE_MLE':
        raise ValueError(
            "Noise types 'NOISE_GAUSSIAN_FIT' and 'NOISE_MLE' are currently unimplemented")

    sig = _prep_sig(noise_type, sig)

    chisq_target = ctypes.c_float(chisq_target)

    Z, A, tau, chisq = ctypes.c_float(), ctypes.c_float(
    ), ctypes.c_float(), ctypes.c_float()  # output values

    tries = _GCI_triple_integral_fitting_engine(
        period, photon_count, fit_start, fit_end, instr, ninstr, 
        _noise_types[noise_type], sig, Z, A, tau, fitted, 
        residuals, chisq, chisq_target)
    
    if(output_fitted_and_residuals):
        return (tries, Z.value, A.value, tau.value, chisq.value, 
                np.asarray(fitted), np.asarray(residuals))
    return tries, Z.value, A.value, tau.value, chisq.value


_GCI_Phasor = _flimlib.GCI_Phasor

# C header
#int    GCI_Phasor(float xincr, float y[], int fit_start, int fit_end,
#							  float *Z, float *u, float *v, float *taup, float *taum, 
#                             float *tau, float *fitted, float *residuals, float *chisq);

_GCI_Phasor.argtypes = [
    ctypes.c_float,ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_int,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float), 
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.POINTER(ctypes.c_float)]

def GCI_Phasor(period, photon_count, output_fitted_and_residuals=False):
    """
    put documentation here
    """
    period = ctypes.c_float(period)

    photon_count, fit_start, fit_end, fitted, residuals = _prep_common_fit_params(photon_count)
    
    Z, u, v, taup, taum, tau, chisq = (ctypes.c_float(), 
    ctypes.c_float(), ctypes.c_float(), ctypes.c_float(), ctypes.c_float(), 
    ctypes.c_float(), ctypes.c_float())  # output values

    error_code = _GCI_Phasor(period, photon_count, fit_start, fit_end, Z, u, v, taup, taum, tau, fitted, residuals, chisq)

    if(output_fitted_and_residuals):
        return (error_code, Z.value, u.value, v.value, taup.value, taum.value, 
                tau.value, chisq.value, np.asarray(fitted), np.asarray(residuals))
    return error_code, Z.value, u.value, v.value, taup.value, taum.value, tau.value, chisq.value


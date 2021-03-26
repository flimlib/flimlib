import ctypes
import numpy as np
import math
import warnings

_flimlib = ctypes.cdll.flimlib

_GCI_triple_integral_fitting_engine = _flimlib.GCI_triple_integral_fitting_engine

_noise_types = {'NOISE_CONST':0, 'NOISE_GIVEN':1, 'NOISE_POISSON_DATA':2,
                    'NOISE_POISSON_FIT':3, 'NOISE_GAUSSIAN_FIT':4, 'NOISE_MLE':5} #4 and 5 should be an error

#C header

#int GCI_triple_integral_fitting_engine(float xincr, float y[], int fit_start, int fit_end,
#                                        float instr[], int ninstr, noise_type noise, float sig[],
#                                        float *Z, float *A, float *tau, float *fitted, float *residuals,
#                                        float *chisq, float chisq_target);

_GCI_triple_integral_fitting_engine.argtypes = [ctypes.c_float,ctypes.POINTER(ctypes.c_float),ctypes.c_int,ctypes.c_int,
                                                ctypes.POINTER(ctypes.c_float),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_float),
                                                ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),
                                                ctypes.POINTER(ctypes.c_float),ctypes.c_float]


def GCI_triple_integral_fitting_engine(period, photonCount, fit_start=0, fit_end=None, 
                                       instr=None, noise_type='NOISE_CONST', sig=1.0, 
                                       chisq_target=1.1, output_fitted_and_residuals=False):
    """
    put documentation here
    """
    
    try:
        photonCount = np.asarray(photonCount)
    except ValueError:
        print("photonCount must be numpy array or array-like")
        raise
    else:
        if photonCount.ndim != 1:
            raise ValueError("photonCount must be a 1D array")
    
    if fit_end is None:
        fit_end = photonCount.shape[0]-1
    
    ninstr = 0
    if instr is not None:
        try:
            instr = np.asarray(instr)
        except ValueError:
            print("photonCount must be numpy array or array-like")
            raise
        instr = np.ctypeslib.as_ctypes(instr) #presumably shorter than photonCount
        ninstr = len(instr)
    
    if noise_type in _noise_types.keys():
        if noise_type == 'NOISE_GAUSSIAN_FIT' or noise_type == 'NOISE_MLE':
            raise ValueError("Noise types 'NOISE_GAUSSIAN_FIT' and 'NOISE_MLE' are currently unimplemented")
    else:
        raise ValueError("invalid noise type. The valid types are:", _noise_types.keys)
    
    if noise_type == 'NOISE_GIVEN':
        try:
            sig = np.asarray(sig)
        except ValueError:
            print("sig must be array like for type NOISE_GIVEN")
            raise
        sig = np.ctypeslib.as_ctypes(sig)
    elif noise_type == 'NOISE_CONST':
        try:
            sig = np.float(np.asarray(sig)) # convert to float
        except(TypeError, ValueError):
            print("sig must be float or length 1 float array for noise type NOISE_CONST")
            raise
        sig = ctypes.c_float(sig)
    elif sig is not None:
        warnings.warn("Expected sig=None for noise type", noise_type, ". Value passed will be ignored")
        sig = None
    
    samples = fit_end-fit_start #exclusive? the fitted and residuals had strange final values if included last index
    
    period = ctypes.c_float(period)
    photonCount = np.ctypeslib.as_ctypes(photonCount)
    fit_start = ctypes.c_int(fit_start)
    fit_end = ctypes.c_int(fit_end)
    chisq_target = ctypes.c_float(chisq_target)
    
    Z, A, tau, chisq = ctypes.c_float(), ctypes.c_float(), ctypes.c_float(), ctypes.c_float() #output values
    
    
    fitted = np.empty(samples,dtype=np.float32)
    residuals = np.empty(samples,dtype=np.float32)
    
    tries = _GCI_triple_integral_fitting_engine(period,photonCount,fit_start,fit_end,
                                                 instr,ninstr,_noise_types[noise_type],sig,
                                                 Z,A,tau,np.ctypeslib.as_ctypes(fitted),np.ctypeslib.as_ctypes(residuals),
                                                 chisq,chisq_target)
    if(output_fitted_and_residuals):
        return (tries,Z.value,A.value,tau.value,chisq.value,fitted,residuals)

    return (tries,Z.value,A.value,tau.value,chisq.value)
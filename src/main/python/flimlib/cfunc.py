import ctypes
import glob
import os
import warnings
from collections import namedtuple

import numpy as np

folder = os.path.dirname(os.path.abspath(__file__))
dll_candidates = (glob.glob(os.path.join(folder, '_flimlib.*.pyd')) +
                  glob.glob(os.path.join(folder, '_flimlib.*.so')))
if len(dll_candidates) > 1:
    raise RuntimeError("More than one flimlib extension found???")
if not dll_candidates:
    raise RuntimeError("flimlib extension missing")
dll_path = dll_candidates[0]

# 0x8 = LOAD_WITH_ALTERED_SEARCH_PATH, allowing absolute path loading
_flimlib = ctypes.CDLL(dll_path, winmode=0x8)

def _prep_sig(noise_type, sig, len):
    """
    Helper function to make sure noise_type is valid and sig is the correct type 
    for the noise_type given
    """
    if not noise_type in _noise_types.keys():
        raise ValueError(
            "invalid noise type. The valid types are: ", _noise_types.keys)
    elif noise_type == 'NOISE_GIVEN':
        sig = np.asarray(sig, dtype=np.float32)
        if sig.shape != (len,):
            raise ValueError("incorrect shape of sig")
        sig = np.ctypeslib.as_ctypes(sig)
    elif noise_type == 'NOISE_CONST':
        sig = float(np.asarray(sig))  # convert to float
        sig = ctypes.c_float(sig)
    elif sig is not None:
        message = "Expected sig=None for noise type " + str(noise_type) + ". The given value of sig will be ignored"
        warnings.warn(message)
        return None
    return sig

def _prep_instr_ninstr(instr):
    """
    Helper function to handle instr input and generate a value for ninstr
    """
    ninstr = 0
    if instr is not None:
        instr = np.asarray(instr, dtype=np.float32)
        if instr.ndim != 1:
            raise ValueError("instr must be 1 dimensional")
        ninstr = instr.shape[0] # presumably shorter than photon_count
        instr = np.ctypeslib.as_ctypes(instr)
    return instr, ninstr

def _prep_common_fit_params(photon_count):
    """
    Helper function to generate common flimlib inputs given the shape of photon_count
    """
    photon_count = np.asarray(photon_count, dtype=np.float32)
    if photon_count.ndim != 1:
        raise ValueError("photon_count must be a 1 dimensional")
    
    fit_start = 0
    fit_end = photon_count.shape[0] # should it be len or len-1?

    photon_count = np.ctypeslib.as_ctypes(photon_count)
    fitted = np.ctypeslib.as_ctypes(np.empty(fit_end, dtype=np.float32))
    residuals = np.ctypeslib.as_ctypes(np.empty(fit_end, dtype=np.float32))

    return photon_count, fit_start, fit_end, fitted, residuals

_GCI_ecf_matrix = _flimlib.GCI_ecf_matrix
_GCI_ecf_matrix.argtypes = [ctypes.c_int, ctypes.c_int]
_GCI_ecf_matrix.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_float))

_GCI_ecf_free_matrix = _flimlib.GCI_ecf_free_matrix
_GCI_ecf_free_matrix.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_float))]


class _EcfMatrix:
    """
    A class used to represent an allocated matrix used by flimlib

    ...

    Attributes
    ----------
    matrix : LP_LP_c_float
        a ctypes object representing a float** matrix
    nrows : int
        the number of rows
    ncols : int
        the number of columns

    Methods
    -------
    asarray()
        Returns a numpy.ndarray copy of the matrix
    """
    def __init__(self, nrows, ncols):
        self.matrix = _GCI_ecf_matrix(nrows, ncols)
        self.nrows = nrows
        self.ncols = ncols

    def __del__(self):
        _GCI_ecf_free_matrix(self.matrix)

    def asarray(self):
        """Returns a numpy.ndarray copy of the matrix"""
        # find the address of the start of the pointers at the beginning of the matrix
        # the first pointer can be found by indexing
        ptr_addr = ctypes.addressof(self.matrix[0])
        # perform pointer arithmatic to skip the pointers at the beginning
        # this gives us the address of the beginning of the data we are interested in
        data_addr = ptr_addr + self.nrows*ctypes.sizeof(ctypes.c_void_p)
        # apparently we must first cast to c_void_p to properly access the address that we calculated
        data_void_ptr = ctypes.c_void_p(data_addr)
        # cast the void pointer into the desired type (a normal 2d float array)
        data_ptr = ctypes.cast(data_void_ptr, ctypes.POINTER(ctypes.c_float))

        return np.ctypeslib.as_array(data_ptr, shape=(self.nrows, self.ncols)).copy()
    
    def fill(self):
        """for debug purposes only"""
        for r in range(self.nrows):
            for c in range(self.ncols):
                self.matrix[r][c] = r*self.ncols+c # enumerate

# Noise types used by flimlib
# 4 and 5 should raise an error for GCI_triple_integral_fitting_engine
_noise_types = {'NOISE_CONST': 0, 'NOISE_GIVEN': 1, 'NOISE_POISSON_DATA': 2,
                'NOISE_POISSON_FIT': 3, 'NOISE_GAUSSIAN_FIT': 4, 'NOISE_MLE': 5}

# restrain types used by flimlib
_restrain_types = {'ECF_RESTRAIN_DEFAULT': 0, 'ECF_RESTRAIN_USER': 1}

_GCI_triple_integral_fitting_engine = _flimlib.GCI_triple_integral_fitting_engine
TripleIntegralResult = namedtuple('TripleIntegralResult', 'tries Z A tau fitted residuals chisq')
_GCI_triple_integral_fitting_engine.argtypes = [
    ctypes.c_float,                 # float xincr
    ctypes.POINTER(ctypes.c_float), # float y[]
    ctypes.c_int,                   # int fit_start
    ctypes.c_int,                   # int fit_end
    ctypes.POINTER(ctypes.c_float), # float instr[]
    ctypes.c_int,                   # int ninstr
    ctypes.c_int,                   # noise_type noise
    ctypes.POINTER(ctypes.c_float), # float sig[]
    ctypes.POINTER(ctypes.c_float), # float *Z
    ctypes.POINTER(ctypes.c_float), # float *A
    ctypes.POINTER(ctypes.c_float), # float *tau
    ctypes.POINTER(ctypes.c_float), # float *fitted
    ctypes.POINTER(ctypes.c_float), # float *residuals
    ctypes.POINTER(ctypes.c_float), # float *chisq
    ctypes.c_float                  # float chisq_target
    ]

def GCI_triple_integral_fitting_engine(period, photon_count,
                                       instr=None, noise_type='NOISE_POISSON_FIT', sig=None,
                                       chisq_target=1.1):
    """
    Performs an exponential fit on the data using Rapid Lifetime Determination

    Parameters
    ----------
    period : float
        The time between samples in `photon_count`
    photon_count : array_like
        The data to be fit. the length of this array determines the length of the fit
    instr : {None, array_like}, optional
        The instrument response. If `instr` is None, no instrument response will be used
        (default is None)
    noise_type : str, optional
        The noise type to use. Valid values are: 'NOISE_CONST', 'NOISE_GIVEN', 
        'NOISE_POISSON_DATA', 'NOISE_POISSON_FIT' (default is 'NOISE_POISSON_FIT')
    sig : {None, float, array_like}, optional
        The noise. A float array the same length as `photon_count` if `noise_type` 
        is 'NOISE_GIVEN', a float if `noise_type` is 'NOISE_CONST' and None otherwise 
        (default is None)
    chisq_target : float, optional
        The target chi squared value for the fit (default is 1.1)

    Returns
    -------
    TripleIntegralResult
        A namedtuple containing values in order: tries, Z, A, tau, fitted, residuals, chisq

    Raises
    ------
    ValueError
        If noise type passed is not implemented yet
    """
    period = ctypes.c_float(period)

    photon_count, fit_start, fit_end, fitted, residuals = _prep_common_fit_params(photon_count)

    instr, ninstr = _prep_instr_ninstr(instr)

    # this lack of implementation is unique to triple integral
    if noise_type == 'NOISE_GAUSSIAN_FIT' or noise_type == 'NOISE_MLE':
        raise ValueError(
            "Noise types 'NOISE_GAUSSIAN_FIT' and 'NOISE_MLE' are currently unimplemented")

    sig = _prep_sig(noise_type, sig, np.asarray(photon_count).shape[0])

    chisq_target = ctypes.c_float(chisq_target)

    Z, A, tau, chisq = ctypes.c_float(), ctypes.c_float(
    ), ctypes.c_float(), ctypes.c_float()  # output values

    tries = _GCI_triple_integral_fitting_engine(
        period, photon_count, fit_start, fit_end, instr, ninstr, 
        _noise_types[noise_type], sig, Z, A, tau, fitted, 
        residuals, chisq, chisq_target)
    
    return TripleIntegralResult(tries, Z.value, A.value, tau.value, 
        np.asarray(fitted), np.asarray(residuals), chisq.value)


_GCI_Phasor = _flimlib.GCI_Phasor # C function
PhasorResult = namedtuple('PhasorResult', 'error_code Z u v taup taum tau fitted residuals chisq')
_GCI_Phasor.argtypes = [
    ctypes.c_float,                 # float xincr
    ctypes.POINTER(ctypes.c_float), # float y[]
    ctypes.c_int,                   # int fit_start
    ctypes.c_int,                   # int fit_end
    ctypes.POINTER(ctypes.c_float), # float *Z
    ctypes.POINTER(ctypes.c_float), # float *u
    ctypes.POINTER(ctypes.c_float), # float *v
    ctypes.POINTER(ctypes.c_float), # float *taup
    ctypes.POINTER(ctypes.c_float), # float *taum
    ctypes.POINTER(ctypes.c_float), # float *tau
    ctypes.POINTER(ctypes.c_float), # float *fitted
    ctypes.POINTER(ctypes.c_float), # float *residuals
    ctypes.POINTER(ctypes.c_float)  # float *chisq
    ]

def GCI_Phasor(period, photon_count):
    """
    Performs an exponential fit on the data using Phasors

    Parameters
    ----------
    period : float
        The time between samples in `photon_count`
    photon_count : array_like
        The data to be fit. the length of this array determines the length of the fit

    Returns
    -------
    TripleIntegralResult
        A namedtuple containing values in order: error_code, Z, u, v, taup, taum, tau, fitted, residuals, chisq
    """
    period = ctypes.c_float(period)

    photon_count, fit_start, fit_end, fitted, residuals = _prep_common_fit_params(photon_count)
    
    Z, u, v, taup, taum, tau, chisq = (ctypes.c_float(), 
    ctypes.c_float(), ctypes.c_float(), ctypes.c_float(), ctypes.c_float(), 
    ctypes.c_float(), ctypes.c_float())  # output values

    error_code = _GCI_Phasor(period, photon_count, fit_start, fit_end, Z, u, v, taup, taum, tau, fitted, residuals, chisq)
    return PhasorResult(error_code, Z.value, u.value, v.value, taup.value, taum.value, 
        tau.value, np.asarray(fitted), np.asarray(residuals), chisq.value)


class FitFunc:
    """
    A class used to represent a fitfunc used by flimlib

    ...

    Attributes
    ----------
    python_func : function
        A python function passed by the user to wrap as a ctypes function pointer
    c_func : ctypes._CFuncPtr
        A ctypes function pointer for the fitfunc for flimlib to use
    nparam_predicate : function
        A python function that returns True if the correct number of parameters
        has been passed for the given fitfunc

    Methods
    -------
    get_c_func(nparam_in)
        Returns ctypes function pointer if `nparam_in` is valid
    """
    def __init__(self, python_func, nparam_predicate=None):
        """
        Parameters
        ----------
        python_func : {function, ctypes._CFuncPtr)
            A python function passed by the user to wrap as a ctypes function pointer. 
            If a ctypes function pointer is passed, it will be used directly instead
        nparam_predicate : {None, function}, optional
            A python function that returns True if the correct number of parameters
            has been passed for the given fitfunc (default is None)
        """
        if isinstance(python_func, ctypes._CFuncPtr):
            # a ctypes c function was passed instead of a python function. We just use it directly
            self.c_func = python_func
        else:
            self.python_func = python_func
            CFUNC = ctypes.CFUNCTYPE(None, ctypes.c_float, ctypes.POINTER(ctypes.c_float), 
                ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.c_int)
            self.c_func = CFUNC(self._wrapped_python_func)
        self.nparam_predicate=nparam_predicate

    def _wrapped_python_func(self, x, param, y, dy_dparam, nparam):
        param_in = np.ctypeslib.as_array(param, shape=(nparam,))
        y_out, dy_dparam_out = self.python_func(x, param_in)
        y.contents.value = ctypes.c_float(y_out).value
        dy_dparam = np.ctypeslib.as_array(dy_dparam,shape=(nparam,))
        dy_dparam[:] = np.asarray(dy_dparam_out, dtype=np.float32)

    def get_c_func(self, nparam_in):
        if self.nparam_predicate is None or self.nparam_predicate(nparam_in):
            return self.c_func
        else:
            raise TypeError("Incorrect size of param")


def _multiexp_predicate(n_param):
    # 3 or greater and odd
    return n_param >= 3 and n_param % 2 == 1

def _stretchedexp_predicate(n_param):
    return n_param == 4

GCI_multiexp_tau = FitFunc(_flimlib.GCI_multiexp_tau, nparam_predicate=_multiexp_predicate)
GCI_multiexp_lambda = FitFunc(_flimlib.GCI_multiexp_lambda, nparam_predicate=_multiexp_predicate)
GCI_stretchedexp = FitFunc(_flimlib.GCI_stretchedexp, nparam_predicate=_stretchedexp_predicate)

_GCI_set_restrain_limits = _flimlib.GCI_set_restrain_limits
_GCI_set_restrain_limits.argtypes= [
    ctypes.POINTER(ctypes.c_int),   # int restrain[]
    ctypes.c_int,                   # int nparam
    ctypes.POINTER(ctypes.c_float), # float minval[]
    ctypes.POINTER(ctypes.c_float)  # float maxval[]
]

def GCI_set_restrain_limits(restrain, minval, maxval):
    """
    Sets global constraints on fit parameters

    Parameters
    ----------
    restrain : array_like
        Array of 1s and 0s. 1 indicates that the parameter will be constrained, 
        0 indicates it will not be constrained
    minval : array_like
        The minimum acceptable value of each parameter
    maxval : array_like
        The maximum acceptable value of each parameter

    Returns
    -------
    None
    """
    minval = np.asarray(minval,dtype=np.float32)
    maxval = np.asarray(maxval,dtype=np.float32)
    restrain = np.asarray(restrain,dtype=np.intc)

    if not(restrain.ndim == minval.ndim == maxval.ndim == 1):
        raise ValueError("restrain, minval and maxval must be 1 dimentional!")
    if not(restrain.shape == minval.shape == maxval.shape):
        raise ValueError("restrain, minval and maxval must have the same shape!")

    nparam = restrain.shape[0]
    minval = np.ctypeslib.as_ctypes(minval)
    maxval = np.ctypeslib.as_ctypes(maxval)
    restrain = np.ctypeslib.as_ctypes(restrain)

    error_code = _GCI_set_restrain_limits(restrain, nparam, minval, maxval)

    if(error_code == -1):
        raise ValueError("invalid size of restrain")
    if(error_code == -2):
        raise ValueError('maxval must be element-wise greater than minval')

_GCI_marquardt_fitting_engine = _flimlib.GCI_marquardt_fitting_engine # C function
MarquardtResult = namedtuple('MarquardtResult', 'error_code param fitted residuals chisq covar alpha erraxes')
_GCI_marquardt_fitting_engine.argtypes= [
    ctypes.c_float,                                 # float xincr
    ctypes.POINTER(ctypes.c_float),                 # float *trans
    ctypes.c_int,                                   # int ndata
    ctypes.c_int,                                   # int fit_start
    ctypes.c_int,                                   # int fit_end
    ctypes.POINTER(ctypes.c_float),                 # float instr[]
    ctypes.c_int,                                   # int ninstr
    ctypes.c_int,                                   # noise_type noise
    ctypes.POINTER(ctypes.c_float),                 # float sig[]
    ctypes.POINTER(ctypes.c_float),                 # float param[]
    ctypes.POINTER(ctypes.c_int),                   # int paramfree[]
    ctypes.c_int,                                   # int nparam 
    ctypes.c_int,                                   # restrain_type restrain
    ctypes.c_void_p,                                # void (*fitfunc)(float, float [], float *, float [], int)
    ctypes.POINTER(ctypes.c_float),                 # float *fitted
    ctypes.POINTER(ctypes.c_float),                 # float *residuals
    ctypes.POINTER(ctypes.c_float),                 # float *chisq
    ctypes.POINTER(ctypes.POINTER(ctypes.c_float)), # float **covar
    ctypes.POINTER(ctypes.POINTER(ctypes.c_float)), # float **alpha
    ctypes.POINTER(ctypes.POINTER(ctypes.c_float)), # float **erraxes
    ctypes.c_float,                                 # float chisq_target
    ctypes.c_float,                                 # float chisq_delta
    ctypes.c_int,                                   # int chisq_percent
    ]

def GCI_marquardt_fitting_engine(period, photon_count, param, paramfree=None, restrain_type='ECF_RESTRAIN_DEFAULT',
                                       fitfunc=GCI_multiexp_tau, instr=None, noise_type='NOISE_POISSON_FIT', sig=None,
                                       chisq_target=1.1, chisq_delta=1E-5, chisq_percent=95):
    """
    Performs an exponential fit on the data using the Levenberg–Marquardt Algorithm

    Parameters
    ----------
    period : float
        The time between samples in `photon_count`
    photon_count : array_like
        The data to be fit. the length of this array determines the length of the fit
    param : array_like
        The initial parameters for the fit. The "guess"
    paramfree : array_like, optional
        Array of 1s and 0s. 1 indicates that the parameter is free, 0 indicates it is fixed.
        If None is passed, all parameters will be free. (default is None)
    restrain_type : str, optional
        restrain type to use. Valid values are 'ECF_RESTRAIN_DEFAULT' and 'ECF_RESTRAIN_USER' (default is 'ECF_RESTRAIN_DEFAULT')
    fitfunc : FitFunc, optional
        A FitFunc object that contains the fit function to be used in the fit. (default is GCI_multiexp_tau)
    instr : {None, array_like}, optional
        The instrument response. If `instr` is None, no instrument response will be used
        (default is None)
    noise_type : str, optional
        The noise type to use. Valid values are: 'NOISE_CONST', 'NOISE_GIVEN', 
        'NOISE_POISSON_DATA', 'NOISE_POISSON_FIT' (default is 'NOISE_POISSON_FIT')
    sig : {None, float, array_like}, optional
        The noise. A float array the same length as `photon_count` if `noise_type` 
        is 'NOISE_GIVEN', a float if `noise_type` is 'NOISE_CONST' and None otherwise 
        (default is None)
    chisq_target : float, optional
        The target chi squared value for the fit (default is 1.1)
    chisq_delta : float, optional
        An individual fit will continue if the chisq value changes by more then this amount (default is 1E-5)
    chisq_percent : int, optional
        Defines the confidence interval when calculating the error axes (default is 95)

    Returns
    -------
    MarquardtResult
        A namedtuple containing values in order: error_code, param, fitted, residuals, chisq, covar, alpha, erraxes
    """
    period = ctypes.c_float(period)

    photon_count, fit_start, fit_end, fitted, residuals = _prep_common_fit_params(photon_count)

    ndata = np.asarray(photon_count).shape[0]

    instr, ninstr = _prep_instr_ninstr(instr)

    sig = _prep_sig(noise_type, sig, np.asarray(photon_count).shape[0])

    param = np.asarray(param,dtype=np.float32)
    nparam = param.shape[0]
    param = np.ctypeslib.as_ctypes(param)

    if paramfree is None:
        paramfree = np.ones(nparam,dtype=np.intc) # default all parameters are free
    else:
        paramfree = np.asarray(paramfree,dtype=np.intc) # TODO should I check to make sure it's ones and zeros?
    paramfree = np.ctypeslib.as_ctypes(paramfree)

    chisq = ctypes.c_float()

    chisq_target = ctypes.c_float(chisq_target)
    chisq_delta = ctypes.c_float(chisq_delta)
    chisq_percent = ctypes.c_int(chisq_percent)

    # allocate the specific style of 2d array used by flimlib using EcfMatrix class
    covar = _EcfMatrix(nparam, nparam)
    alpha = _EcfMatrix(nparam, nparam)
    erraxes = _EcfMatrix(nparam, nparam)

    # cast ctypes function pointer to c_void_p
    fitfunc = ctypes.cast(fitfunc.get_c_func(nparam), ctypes.c_void_p)

    error_code = _GCI_marquardt_fitting_engine(period, photon_count, ndata, fit_start, fit_end, 
        instr, ninstr, _noise_types[noise_type], sig, param, paramfree, nparam, _restrain_types[restrain_type],
        fitfunc, fitted, residuals, chisq,
        covar.matrix, alpha.matrix, erraxes.matrix, chisq_target, chisq_delta, chisq_percent)
    return MarquardtResult(error_code, np.asarray(param), np.asarray(fitted), 
        np.asarray(residuals), chisq.value, covar.asarray(), alpha.asarray(), erraxes.asarray())


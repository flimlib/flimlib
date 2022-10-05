###
# #%L
# FLIMLib package for exponential curve fitting of fluorescence lifetime data.
# %%
# Copyright (C) 2010 - 2022 University of Oxford and Board of Regents of the
# University of Wisconsin-Madison.
# %%
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public
# License along with this program.  If not, see
# <http://www.gnu.org/licenses/gpl-3.0.html>.
# #L%
###
import ctypes
import glob
import os
import warnings
from dataclasses import dataclass

import numpy as np

folder = os.path.dirname(os.path.abspath(__file__))
dll_candidates = glob.glob(os.path.join(folder, "_flimlib.*.pyd")) + glob.glob(
    os.path.join(folder, "_flimlib.*.so")
)
if len(dll_candidates) > 1:
    raise RuntimeError("More than one flimlib extension found???")
if not dll_candidates:
    raise RuntimeError("flimlib extension missing")
dll_path = dll_candidates[0]

# 0x8 = LOAD_WITH_ALTERED_SEARCH_PATH, allowing absolute path loading
_flimlib = ctypes.CDLL(dll_path, winmode=0x8)

# Noise types used by flimlib
# 4 and 5 should raise an error for GCI_triple_integral_fitting_engine
_noise_types = {
    "NOISE_CONST": 0,
    "NOISE_GIVEN": 1,
    "NOISE_POISSON_DATA": 2,
    "NOISE_POISSON_FIT": 3,
    "NOISE_GAUSSIAN_FIT": 4,
    "NOISE_MLE": 5,
}

# restrain types used by flimlib
_restrain_types = {"ECF_RESTRAIN_DEFAULT": 0, "ECF_RESTRAIN_USER": 1}


class FitFunc:
    """
    A class used to represent a fitfunc used by flimlib

    Attributes
    ----------
    python_func : function
        A python function passed by the user to wrap as a ctypes function
        pointer
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
        python_func : {function, ctypes._CFuncPtr}
            A python function passed by the user to wrap as a ctypes function
            pointer. If a ctypes function pointer is passed, it will be used
            directly instead
        nparam_predicate : {None, function}, optional
            A python function that returns True if the correct number of
            parameters has been passed for the given fitfunc (default is None)
        """
        if isinstance(python_func, ctypes._CFuncPtr):
            # a ctypes c function was passed instead of a python function. We
            # just use it directly
            self.c_func = python_func
        else:
            self.python_func = python_func
            CFUNC = ctypes.CFUNCTYPE(
                None,
                ctypes.c_float,
                ctypes.POINTER(ctypes.c_float),
                ctypes.POINTER(ctypes.c_float),
                ctypes.POINTER(ctypes.c_float),
                ctypes.c_int,
            )
            self.c_func = CFUNC(self._wrapped_python_func)
        self.nparam_predicate = nparam_predicate

    def _wrapped_python_func(self, x, param, y, dy_dparam, nparam):
        param_in = np.ctypeslib.as_array(param, shape=(nparam,))
        y_out, dy_dparam_out = self.python_func(x, param_in)
        y.contents.value = ctypes.c_float(y_out).value
        dy_dparam = np.ctypeslib.as_array(dy_dparam, shape=(nparam,))
        dy_dparam[:] = np.asarray(dy_dparam_out, dtype=np.float32)

    def get_c_func(self, nparam_in):
        if self.nparam_predicate is None or self.nparam_predicate(nparam_in):
            return self.c_func
        else:
            raise TypeError("Incorrect size of param: " + str(nparam_in))


def _multiexp_predicate(n_param):
    # 3 or greater and odd
    return n_param >= 3 and n_param % 2 == 1


def _stretchedexp_predicate(n_param):
    return n_param == 4


GCI_multiexp_tau = FitFunc(
    _flimlib.GCI_multiexp_tau, nparam_predicate=_multiexp_predicate
)
GCI_multiexp_lambda = FitFunc(
    _flimlib.GCI_multiexp_lambda, nparam_predicate=_multiexp_predicate
)
GCI_stretchedexp = FitFunc(
    _flimlib.GCI_stretchedexp, nparam_predicate=_stretchedexp_predicate
)

_GCI_set_restrain_limits = _flimlib.GCI_set_restrain_limits
_GCI_set_restrain_limits.argtypes = [
    ctypes.POINTER(ctypes.c_int),  # int restrain[]
    ctypes.c_int,  # int nparam
    ctypes.POINTER(ctypes.c_float),  # float minval[]
    ctypes.POINTER(ctypes.c_float),  # float maxval[]
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
    minval = np.asarray(minval, dtype=np.float32)
    maxval = np.asarray(maxval, dtype=np.float32)
    restrain = np.asarray(restrain, dtype=np.intc)

    if not (restrain.ndim == minval.ndim == maxval.ndim == 1):
        raise ValueError("restrain, minval and maxval must be 1 dimentional!")
    if not (restrain.shape == minval.shape == maxval.shape):
        raise ValueError(
            "restrain, minval and maxval must have the same shape!"
        )

    nparam = restrain.shape[0]
    minval = np.ctypeslib.as_ctypes(minval)
    maxval = np.ctypeslib.as_ctypes(maxval)
    restrain = np.ctypeslib.as_ctypes(restrain)

    error_code = _GCI_set_restrain_limits(restrain, nparam, minval, maxval)

    if error_code == -1:
        raise ValueError("invalid size of restrain")
    if error_code == -2:
        raise ValueError("maxval must be element-wise greater than minval")


# wrappers for C structs
class _Array1D(ctypes.Structure):
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.c_float)),
        ("sizes", ctypes.c_size_t * 1),
        ("strides", ctypes.c_ssize_t * 1),
    ]


class _Array2D(ctypes.Structure):
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.c_float)),
        ("sizes", ctypes.c_size_t * 2),
        ("strides", ctypes.c_ssize_t * 2),
    ]


class _Array3D(ctypes.Structure):
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.c_float)),
        ("sizes", ctypes.c_size_t * 3),
        ("strides", ctypes.c_ssize_t * 3),
    ]


class _Array1DInt8(ctypes.Structure):
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.c_int8)),
        ("sizes", ctypes.c_size_t * 1),
        ("strides", ctypes.c_ssize_t * 1),
    ]


class _CommonParams(ctypes.Structure):
    _fields_ = [
        ("xincr", ctypes.c_float),
        ("trans", ctypes.POINTER(_Array2D)),
        ("fit_start", ctypes.c_int),
        ("fit_end", ctypes.c_int),
        ("fitted", ctypes.POINTER(_Array2D)),
        ("residuals", ctypes.POINTER(_Array2D)),
        ("chisq", ctypes.POINTER(_Array1D)),
        ("fit_mask", ctypes.POINTER(_Array1DInt8)),
    ]


class _MarquardtParams(ctypes.Structure):
    _fields_ = [
        ("instr", ctypes.POINTER(_Array1D)),
        ("noise", ctypes.c_int),
        ("sig", ctypes.POINTER(_Array1D)),
        ("param", ctypes.POINTER(_Array2D)),
        ("paramfree", ctypes.POINTER(_Array1DInt8)),
        ("restrain", ctypes.c_int),
        ("fitfunc", ctypes.c_void_p),
        ("covar", ctypes.POINTER(_Array3D)),
        ("alpha", ctypes.POINTER(_Array3D)),
        ("erraxes", ctypes.POINTER(_Array3D)),
        ("chisq_target", ctypes.c_float),
        ("chisq_delta", ctypes.c_float),
        ("chisq_percent", ctypes.c_int),
    ]


class _TripleIntegralParams(ctypes.Structure):
    _fields_ = [
        ("instr", ctypes.POINTER(_Array1D)),
        ("noise", ctypes.c_int),
        ("sig", ctypes.POINTER(_Array1D)),
        ("Z", ctypes.POINTER(_Array1D)),
        ("A", ctypes.POINTER(_Array1D)),
        ("tau", ctypes.POINTER(_Array1D)),
        ("chisq_target", ctypes.c_float),
    ]


class _PhasorParams(ctypes.Structure):
    _fields_ = [
        ("Z", ctypes.POINTER(_Array1D)),
        ("u", ctypes.POINTER(_Array1D)),
        ("v", ctypes.POINTER(_Array1D)),
        ("taup", ctypes.POINTER(_Array1D)),
        ("taum", ctypes.POINTER(_Array1D)),
        ("tau", ctypes.POINTER(_Array1D)),
    ]


class _SpecificParams(ctypes.Union):
    _fields_ = [
        ("marquardt", ctypes.POINTER(_MarquardtParams)),
        ("triple_integral", ctypes.POINTER(_TripleIntegralParams)),
        ("phasor", ctypes.POINTER(_PhasorParams)),
    ]


class _FlimParams(ctypes.Structure):
    _anonymous_ = [
        "specific",
    ]
    _fields_ = [
        ("common", ctypes.POINTER(_CommonParams)),
        ("specific", _SpecificParams),
    ]


def _equal_shapes(shape1: tuple, shape2: tuple):
    """shape comparison with wildcard is Ellipses"""
    if shape1 is Ellipsis or shape2 is Ellipsis:  # any shape
        return True
    if len(shape1) != len(shape2):  # different ndims
        return False
    else:
        for i in range(len(shape1)):
            if shape1[i] != shape2[i] and not (
                shape1[i] is Ellipsis or shape2[i] is Ellipsis
            ):
                return False
    return True


def _as_strided_array(
    array_in,
    check_shape,
    ctypes_shape,
    ctypes_type=ctypes.c_float,
    numpy_type=np.float32,
    shape_override=None,
    strides_override=None,
):
    arr_in = np.asarray(array_in, dtype=numpy_type)
    if not _equal_shapes(arr_in.shape, check_shape):
        raise ValueError(
            "expected array with shape="
            + str(check_shape)
            + ", got shape="
            + str(arr_in.shape)
        )
    arr = arr_in.reshape(ctypes_shape)

    if arr.ndim == 1:
        if ctypes_type == ctypes.c_int8:
            result = _Array1DInt8()
        else:
            result = _Array1D()
    elif arr.ndim == 2:
        result = _Array2D()
    elif arr.ndim == 3:
        result = _Array3D()
    else:
        raise ValueError(
            "Invalid ndim. This is an internal bug with the flimlib python wrapper"
        )
    result.data = arr.ctypes.data_as(ctypes.POINTER(ctypes_type))
    # must be unsigned
    result.sizes = (
        (ctypes.c_size_t * arr.ndim)(*arr.shape)
        if shape_override is None
        else (ctypes.c_size_t * arr.ndim)(*shape_override)
    )
    result.strides = (
        arr.ctypes.strides
        if strides_override is None
        else (ctypes.c_ssize_t * arr.ndim)(*strides_override)
    )

    return ctypes.pointer(result), arr


def _prep_optional_output(input, check_shape, ctypes_shape, fit_mask, flag=True):
    if input is None:
        if flag:
            temp_in = np.empty(check_shape, dtype=np.float32)
            if fit_mask is not None:
                temp_in.fill(np.nan)
            return _as_strided_array(
                temp_in,
                check_shape,
                ctypes_shape,
            )
        else:
            return None, None
    elif type(input) is np.ndarray:
        if input.dtype == np.float32:
            return _as_strided_array(input, check_shape, ctypes_shape)
        elif np.issubdtype(input.dtype, np.floating):
            message = "Pre-allocated output is not of type float32. An internal copy will be made."
            warnings.warn(message)
            if fit_mask is None:
                temp_in = np.empty(input.shape, dtype=np.float32)
            else:
                temp_in = input.astype(np.float32) # always makes a copy
            return _as_strided_array(
                temp_in,
                check_shape,
                ctypes_shape,
            )
        else:
            raise TypeError(
                "pre-allocated outputs must be of type float"
            )
    else:
        raise TypeError("pre-allocated outputs must be numpy.ndarray")


def _copy_to_provided_output(user_arr, c_arr, shape):
    """
    user_arr
        The provided output (eg. fitted). It may have more dimensions than was
        expected by internal C code
    c_arr
        The reshaped array that was passed into the C code. If the data had the
        correct dtype and was contiguous no copy is made
    """
    if c_arr is None:
        return None
    reshaped_c_arr = c_arr.reshape(shape)
    if user_arr is not None:
        if user_arr.dtype != np.float32:  # copy was made to change dtype
            np.copyto(user_arr, reshaped_c_arr)
        else:
            try:
                user_shape = user_arr.shape
                # if error, copy was made for reshape (non-contiguous data)
                user_arr.shape = np.prod(user_arr.shape)
                user_arr.shape = user_shape  # undo change
            except AttributeError:
                np.copyto(user_arr, reshaped_c_arr)
        return user_arr
    else:
        return reshaped_c_arr


def _prep_strided_sig(noise_type, sig, check_shape, ctypes_shape):
    """
    Helper function to make sure noise_type is valid and sig is the correct
    type for the noise_type given
    """
    if not noise_type in _noise_types.keys():
        raise ValueError(
            "invalid noise type. The valid types are: ", _noise_types.keys
        )
    elif noise_type == "NOISE_GIVEN":
        return _as_strided_array(sig, check_shape, ctypes_shape)
    elif noise_type == "NOISE_CONST":
        return _as_strided_array(np.asarray(sig).flatten(), (1,), (1,))
    elif sig is not None:
        message = (
            "Expected sig=None for noise type "
            + str(noise_type)
            + ". The given value of sig will be ignored"
        )
        warnings.warn(message)
    return None, None


def _prep_common_params(
    period,
    photon_count,
    fit_start,
    fit_end,
    fit_mask,
    fitted,
    residuals,
    chisq,
    compute_fitted,
    compute_residuals,
    compute_chisq,
    ndata_known,
):
    """
    Helper function to generate common flimlib inputs
    """
    common = _CommonParams()

    # marquardt expects fitted to be size of trans
    # triple integral expects fitted to be size = fit_end
    # phasors expects fitted to be size = fit_end

    dshape = np.asarray(photon_count).shape  # shape of the input data
    if dshape == ():
        raise TypeError("photon_count must be array-like")
    fstart = 0 if fit_start is None else fit_start
    fend = dshape[-1] if fit_end is None else fit_end
    assert 0 <= fstart <= fend <= dshape[-1], "invalid fit range!"
    common.fit_start = fstart
    common.fit_end = fend
    # triple_integral and phasor use fit_end as the size of the data
    data_shape = dshape if ndata_known else (*dshape[0:-1], fend)
    npixels = np.prod(data_shape[0:-1], dtype=int)
    assert period > 0.0, "period must be greater than 0"
    common.xincr = period
    # the shape of trans is wildcard
    common.trans, referenced_trans = _as_strided_array(
        photon_count, ..., (npixels, dshape[-1])
    )
    common.fitted, fitted_out = _prep_optional_output(
        fitted, data_shape, (npixels, data_shape[-1]), fit_mask, flag=compute_fitted
    )
    common.residuals, residuals_out = _prep_optional_output(
        residuals,
        data_shape,
        (npixels, data_shape[-1]),
        fit_mask,
        flag=compute_residuals,
    )
    common.chisq, chisq_out = _prep_optional_output(
        chisq, data_shape[0:-1], (npixels,), fit_mask, flag=compute_chisq
    )
    common.fit_mask, referenced_fit_mask = (
        (None, None)
        if fit_mask is None
        else _as_strided_array(
            fit_mask,
            data_shape[0:-1],
            (npixels,),
            ctypes_type=ctypes.c_int8,
            numpy_type=np.int8,
        )
    )
    # stuff we want to prevent from getting garbage collected
    referenced_objects = (
        referenced_trans,
        referenced_fit_mask,
    )
    return (
        common,
        fitted_out,
        residuals_out,
        chisq_out,
        data_shape,
        npixels,
        referenced_objects,
    )


@dataclass
class MarquardtResult:
    """
    A dataclass containing the outputs of GCI_marquardt_fitting_engine

    Attributes
    ----------
    param : numpy.ndarray
        An array containing arrays of parameters, if the input `param` was
        dtype=float32, the input is modified and returned.
    fitted : {None, numpy.ndarray}.
        An array containing arrays containing values fitted to the data for
        each fit. Fit points are coincident in time with the data points.
    residuals : {None, numpy.ndarray}
        An array containing arrays containing the difference between the data
        and the fit for each fit.
    chisq : {None, numpy.ndarray}
        An array containing the resulting raw chi squared value of each fit. To
        get the reduced chisq, divide by the degrees of freedom (fit_start -
        fit_end - nparam)
    covar : numpy.ndarray
        An array containing the covariance matrix from each fit.
    alpha : numpy.ndarray
        An array containing the alpha matrix from each fit.
    erraxes : numpy.ndarray
        An array containing the dimensions of the confidence ellipsoid of the
        chisq from each fit.
    """

    param: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    chisq: np.ndarray
    covar: np.ndarray
    alpha: np.ndarray
    erraxes: np.ndarray


# C function
_GCI_marquardt_fitting_engine = _flimlib.GCI_marquardt_fitting_engine_many
_GCI_marquardt_fitting_engine.argtypes = [ctypes.POINTER(_FlimParams)]


def GCI_marquardt_fitting_engine(
    period,
    photon_count,
    param,
    fit_start=None,
    fit_end=None,
    instr=None,
    noise_type="NOISE_POISSON_FIT",
    sig=None,
    paramfree=None,
    restrain_type="ECF_RESTRAIN_DEFAULT",
    fitfunc=GCI_multiexp_tau,
    fitted=None,
    residuals=None,
    chisq=None,
    covar=None,
    alpha=None,
    erraxes=None,
    chisq_target=1.1,
    chisq_delta=1e-5,
    chisq_percent=95,
    fit_mask=None,
    compute_fitted=True,
    compute_residuals=True,
    compute_chisq=True,
    compute_covar=True,
    compute_alpha=True,
    compute_erraxes=True,
) -> MarquardtResult:
    """
    Performs a multipixel exponential fits on the data using the
    Levenberg-Marquardt Algorithm

    Parameters
    ----------
    period : float
        The time between samples in `photon_count`
    photon_count : array_like
        An N-dimensional array containing the data to be fit. The last axis is
        the time axis over which the data is fit. Any preceding axes will be
        preserved in the outputs.
    param : array_like
        An N-dimensional array of parameters. Provide parameter estimates,
        these are overridden with the fitted values. The shape must match
        `photon_count` except the last axis, which must match the parameters of
        `fitfunc`
    fit_start : {None, int}, optional
        The index of the start of the fit. Some data before this start index is
        required if convolving with the prompt. If is None, the fit will begin
        at index 0 (default is None)
    fit_end : {None, int}, optional
        The index of the end of the fit. If is None, the fit will cover the
        entire time axis of `photon_count`
    instr : {None, array_like}, optional
        A 1D array containing the instrument response (IRF) or prompt signal.
        If is None, no instrument response will be used (default is None)
    noise_type : str, optional
        The noise type to use. Valid values are: 'NOISE_CONST', 'NOISE_GIVEN',
        'NOISE_POISSON_DATA', 'NOISE_POISSON_FIT' (default is
        'NOISE_POISSON_FIT')
    sig : {None, float, array_like}, optional
        The standard deviation at each data point. A 1D float array the same
        length as the time axis of `photon_count` if `noise_type` is
        'NOISE_GIVEN', a float if `noise_type` is 'NOISE_CONST' and None
        otherwise (default is None)
    paramfree : {None, array_like}, optional
        A 1D array indicating which parameters are free (1), fixed (0)
        If is None, all parameters will be free. (default is None)
    restrain_type : str, optional
        restrain type to use. Valid values are 'ECF_RESTRAIN_DEFAULT' and
        'ECF_RESTRAIN_USER' (default is 'ECF_RESTRAIN_DEFAULT')
    fitfunc : FitFunc, optional
        A FitFunc object that contains the fit function to be used in the fit.
        (default is GCI_multiexp_tau)
    fitted : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed fitted plot for
        each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    residuals : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed residuals plot
        for each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    chisq : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed raw chi squared
        value for each fit. To avoid copying, use dtype=np.float32. If is None,
        a new array will be created (default is None)
    covar : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed covariance matrix
        for each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    alpha : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed alpha matrix for
        each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    erraxes : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed dimensions of the
        confidence ellipsoid of the chisq for each fit. To avoid copying, use
        dtype=np.float32. If is None, a new array will be created (default is
        None)
    chisq_target : float, optional
        chisq_target A raw chi squared value to aim for. If this value is
        reached fitting will stop. If you want to aim for a reduced chisq (say
        1.1 or 1.0) you must multiply by the degree of freedom. (default is
        1.1)
    chisq_delta : float, optional
        An individual fit will continue if the chisq value changes by more then
        this amount (default is 1E-5)
    chisq_percent : int, optional
        Defines the confidence interval when calculating the error axes
        (default is 95)
    fit_mask : {None, array_like}, optional
        An N-dimensional array of bool or 1s and 0s to select which pixels to
        fit. If is None, all pixels will be fit (default is None)
    compute_fitted : bool, optional
        If True, the fitted plot for each fit is kept in memory and returned.
        Ignored if `fitted` is not None (default is True)
    compute_residuals : bool, optional
        If True, the residuals plot for each fit is kept in memory and
        returned. Ignored if `residuals` is not None (default is True)
    compute_chisq : bool, optional
        If True, the raw chi squared value for each fit is kept in memory and
        returned. Ignored if `chisq` is not None (default is True)
    compute_covar : bool, optional
        If True, the covariance matrix for each fit is kept in memory and
        returned. Ignored if `covar` is not None (default is True)
    compute_alpha : bool, optional
        If True, the alpha matrix for each fit is kept in memory and
        returned. Ignored if `alpha` is not None (default is True)
    compute_erraxes : bool, optional
        If True, the dimensions of the confidence ellipsoid of the chisq for
        each fit is kept in memory and returned. Ignored if `erraxes` is not
        None (default is True)

    Returns
    -------
    MarquardtResult
        A dataclass containing fields to be acessed using dot (.) notation:
        param, fitted, residuals, chisq, covar, alpha, erraxes. For each fit
        that fails, its corresponding outputs are filled with `NaN`
    """

    (
        common_in,
        fitted_out,
        residuals_out,
        chisq_out,
        data_shape,
        npixels,
        referenced_objects,
    ) = _prep_common_params(
        period,
        photon_count,
        fit_start,
        fit_end,
        fit_mask,
        fitted,
        residuals,
        chisq,
        compute_fitted,
        compute_residuals,
        compute_chisq,
        True,
    )

    marquardt_in = _MarquardtParams()

    # must pass unit instr because of a bug with flimlib
    marquardt_in.instr, referenced_instr = (
        _as_strided_array([1.0], (1,), (1,))
        if instr is None
        else _as_strided_array(instr, (...,), (-1,))
    )
    # same size as second axis of fitted and residuals
    marquardt_in.sig, referenced_sig = _prep_strided_sig(
        noise_type, sig, (data_shape[-1],), (data_shape[-1],)
    )
    marquardt_in.noise = _noise_types[noise_type]

    # make a copy of param. this is less than ideal performance-wise but in
    # python it is expected that inputs are not modified
    param_in = np.asarray(param).copy()
    nparam = param_in.shape[-1]
    # Does it make sense to pass None and start all guesses at 0 or something?
    marquardt_in.param, param_out = _as_strided_array(
        param_in, (*data_shape[0:-1], ...), (npixels, nparam)
    )
    marquardt_in.paramfree, referenced_paramfree = (
        (None, None)
        if paramfree is None
        else _as_strided_array(
            paramfree,
            (nparam,),
            (nparam,),
            ctypes_type=ctypes.c_int8,
            numpy_type=np.int8,
        )
    )
    marquardt_in.restrain = _restrain_types[restrain_type]
    if all(data_shape):  # is there any data
        # nparam will be checked by the fitfunc predicate
        marquardt_in.fitfunc = ctypes.cast(
            fitfunc.get_c_func(nparam), ctypes.c_void_p
        )
    marquardt_in.covar, covar_out = _prep_optional_output(
        covar,
        (*data_shape[0:-1], nparam, nparam),
        (npixels, nparam, nparam),
        fit_mask,
        flag=compute_covar,
    )
    marquardt_in.alpha, alpha_out = _prep_optional_output(
        alpha,
        (*data_shape[0:-1], nparam, nparam),
        (npixels, nparam, nparam),
        fit_mask,
        flag=compute_alpha,
    )
    marquardt_in.erraxes, erraxes_out = _prep_optional_output(
        erraxes,
        (*data_shape[0:-1], nparam, nparam),
        (npixels, nparam, nparam),
        fit_mask,
        flag=compute_erraxes,
    )
    marquardt_in.chisq_target = chisq_target
    marquardt_in.chisq_delta = chisq_delta
    marquardt_in.chisq_percent = chisq_percent

    flim_in = _FlimParams()

    flim_in.common = ctypes.pointer(common_in)
    flim_in.marquardt = ctypes.pointer(marquardt_in)

    # ctypes automatically gets the pointer
    _GCI_marquardt_fitting_engine(flim_in)

    # Verify that reference was held until after above call
    referenced_instr
    referenced_sig
    referenced_paramfree
    referenced_objects

    # print(data_shape)
    # DO NOT copy to provided param since it is also an input
    param_out = _copy_to_provided_output(None, param_out, param_in.shape)
    fitted_out = _copy_to_provided_output(fitted, fitted_out, data_shape)
    residuals_out = _copy_to_provided_output(
        residuals, residuals_out, data_shape
    )
    chisq_out = _copy_to_provided_output(chisq, chisq_out, data_shape[0:-1])
    covar_out = _copy_to_provided_output(
        covar, covar_out, (*data_shape[0:-1], nparam, nparam)
    )
    alpha_out = _copy_to_provided_output(
        alpha, alpha_out, (*data_shape[0:-1], nparam, nparam)
    )
    erraxes_out = _copy_to_provided_output(
        erraxes, erraxes_out, (*data_shape[0:-1], nparam, nparam)
    )

    return MarquardtResult(
        param_out,
        fitted_out,
        residuals_out,
        chisq_out,
        covar_out,
        alpha_out,
        erraxes_out,
    )


@dataclass
class TripleIntegralResult:
    """
    A dataclass containing the outputs of GCI_triple_integral_fitting_engine

    Attributes
    ----------
    Z : np.ndarray
        An array containing the returned background value from each fit.
    A : np.ndarray
        An array containing the returned amplitude value from each fit.
    tau : np.ndarray
        An array containing the returned lifetime value from each fit.
    fitted : {None, numpy.ndarray}
        An array containing arrays containing values fitted to the data for
        each fit. Fit points are coincident in time with the data points.
    residuals : {None, numpy.ndarray}
        An array containing arrays containing the difference between the data
        and the fit for each fit.
    chisq : np.ndarray
        An array containing the resulting raw chi squared value of each fit. To
        get the reduced chisq, divide by the degrees of freedom (fit_start -
        fit_end - nparam)
    """

    Z: np.ndarray
    A: np.ndarray
    tau: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    chisq: np.ndarray


# C function
_GCI_triple_integral_fitting_engine = (
    _flimlib.GCI_triple_integral_fitting_engine_many
)
_GCI_triple_integral_fitting_engine.argtypes = [ctypes.POINTER(_FlimParams)]


def GCI_triple_integral_fitting_engine(
    period,
    photon_count,
    fit_start=None,
    fit_end=None,
    instr=None,
    noise_type="NOISE_POISSON_FIT",
    sig=None,
    Z=None,
    A=None,
    tau=None,
    fitted=None,
    residuals=None,
    chisq=None,
    chisq_target=-1.0,
    fit_mask=None,
    compute_fitted=True,
    compute_residuals=True,
    compute_chisq=True,
):
    """
    Performs a multipixel exponential fits on the data using Rapid Lifetime
    Determination

    Parameters
    ----------
    period : float
        The time between samples in `photon_count`
    photon_count : array_like
        An N-dimensional array containing the data to be fit. The last axis is
        the time axis over which the data is fit. Any preceding axes will be
        preserved in the outputs.
    fit_start : {None, int}, optional
        The index of the start of the fit. Some data before this start index is
        required if convolving with the prompt.
        If is None, the fit will begin at index 0 (default is None)
    fit_end : {None, int}, optional
        The index of the end of the fit. If is None, the fit will cover the
        entire time axis of `photon_count`
    instr : {None, array_like}, optional
        A 1D array containing the instrument response (IRF) or prompt signal.
        If is None, no instrument response will be used (default is None)
    noise_type : str, optional
        The noise type to use. Valid values are: 'NOISE_CONST', 'NOISE_GIVEN',
        'NOISE_POISSON_DATA', 'NOISE_POISSON_FIT' (default is
        'NOISE_POISSON_FIT')
    sig : {None, float, array_like}, optional
        The standard deviation at each data point. A 1D float array the same
        length as the time axis of `photon_count` if `noise_type` is
        'NOISE_GIVEN', a float if `noise_type` is 'NOISE_CONST' and None
        otherwise (default is None)
    Z : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed background value
        for each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    A : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed amplitude value
        for each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    tau : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed lifetime value
        for each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    fitted : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed fitted plot for
        each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    residuals : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed residuals plot
        for each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    chisq : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed raw chi squared
        value for each fit. To avoid copying, use dtype=np.float32. If is None,
        a new array will be created (default is None)
    chisq_target : float, optional
        chisq_target A raw chi squared value to aim for. If this value is
        reached fitting will stop. Retries will shorten the fit range (without
        changing fit_start). If you want to aim for a reduced chisq (say 1.1 or
        1.0) you must multiply by the degree of freedom. A negative value will
        lead to only a single iteration (default is -1.0)
    fit_mask : {None, array_like}, optional
        An N-dimensional array of bool or 1s and 0s to select which pixels to
        fit. If is None, all pixels will be fit (default is None)
    compute_fitted : bool, optional
        If True, the fitted plot for each fit is computed. Ignored if `fitted` is
        not None (default is True)
    compute_residuals : bool, optional
        If True, the residuals plot for each fit is computed. Ignored if
        `residuals` is not None (default is True)
    compute_chisq : bool, optional
        If True, the raw chi squared value for each fit is computed. Ignored if
        `chisq` is not None (default is True)

    Returns
    -------
    TripleIntegralResult
        A dataclass containing fields to be acessed using dot (.) notation: Z,
        A, tau, fitted, residuals, chisq. For each fit that fails, its
        corresponding outputs are filled with `NaN`
    """
    (
        common_in,
        fitted_out,
        residuals_out,
        chisq_out,
        data_shape,
        npixels,
        referenced_objects,
    ) = _prep_common_params(
        period,
        photon_count,
        fit_start,
        fit_end,
        fit_mask,
        fitted,
        residuals,
        chisq,
        compute_fitted,
        compute_residuals,
        compute_chisq,
        False,
    )
    triple_integral_in = _TripleIntegralParams()
    triple_integral_in.instr, referenced_instr = (
        (None, None)
        if instr is None
        else _as_strided_array(instr, (...,), (-1,))
    )
    # this lack of implementation is unique to triple integral
    if noise_type == "NOISE_GAUSSIAN_FIT" or noise_type == "NOISE_MLE":
        raise ValueError(
            "Noise types 'NOISE_GAUSSIAN_FIT' and 'NOISE_MLE' are currently unimplemented for GCI_triple_integral"
        )
    # same size as second axis of fitted and residuals
    triple_integral_in.sig, referenced_sig = _prep_strided_sig(
        noise_type, sig, (data_shape[-1],), (data_shape[-1],)
    )
    triple_integral_in.noise = _noise_types[noise_type]
    triple_integral_in.Z, Z_out = _prep_optional_output(
        Z, data_shape[0:-1], (npixels,), fit_mask
    )
    triple_integral_in.A, A_out = _prep_optional_output(
        A, data_shape[0:-1], (npixels,), fit_mask
    )
    triple_integral_in.tau, tau_out = _prep_optional_output(
        tau, data_shape[0:-1], (npixels,), fit_mask
    )
    triple_integral_in.chisq_target = chisq_target

    flim_in = _FlimParams()

    flim_in.common = ctypes.pointer(common_in)
    flim_in.triple_integral = ctypes.pointer(triple_integral_in)

    _GCI_triple_integral_fitting_engine(flim_in)

    # Verify that reference was held until after above call
    referenced_instr
    referenced_sig
    referenced_objects

    fitted_out = _copy_to_provided_output(fitted, fitted_out, data_shape)
    residuals_out = _copy_to_provided_output(
        residuals, residuals_out, data_shape
    )
    chisq_out = _copy_to_provided_output(chisq, chisq_out, data_shape[0:-1])
    Z_out = _copy_to_provided_output(Z, Z_out, data_shape[0:-1])
    A_out = _copy_to_provided_output(A, A_out, data_shape[0:-1])
    tau_out = _copy_to_provided_output(tau, tau_out, data_shape[0:-1])

    return TripleIntegralResult(
        Z_out, A_out, tau_out, fitted_out, residuals_out, chisq_out
    )


@dataclass
class PhasorResult:
    """
    A dataclass containing the outputs of GCI_Phasor

    Attributes
    ----------
    u : numpy.ndarray
        An array containing the 'horizontal' phasor coordinate for each fit.
    v : numpy.ndarray
        An array containing the 'vertical' phasor coordinate for each fit.
    taup : float
        An array containg the lifetime calculated from the phase change for
        each fit.
    taum : numpy.ndarray
        An array containing the lifetimes calculated from the amplitude change
        (the demodulation) for each fit.
    tau : numpy.ndarray
        An array containing the averages of the other taus for each fit
    fitted : {None, numpy.ndarray}
        An array containing arrays containing values fitted to the data for
        each fit. Fit points are coincident in time with the data points.
    residuals : {None, numpy.ndarray}
        An array containing arrays containing the difference between the data
        and the fit for each fit.
    chisq : np.ndarray
        An array containing the resulting reduced chi squared value of each
        fit. To get the raw chisq, multiply by the degrees of freedom
        (fit_start - fit_end - nparam)
    """

    u: np.ndarray
    v: np.ndarray
    taup: np.ndarray
    taum: np.ndarray
    tau: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    chisq: np.ndarray


_GCI_Phasor = _flimlib.GCI_Phasor_many  # C function
_GCI_Phasor.argtypes = [ctypes.POINTER(_FlimParams)]


def GCI_Phasor(
    period,
    photon_count,
    fit_start=None,
    fit_end=None,
    Z=0.0,
    u=None,
    v=None,
    taup=None,
    taum=None,
    tau=None,
    fitted=None,
    residuals=None,
    chisq=None,
    fit_mask=None,
    compute_fitted=True,
    compute_residuals=True,
    compute_chisq=True,
):
    """
    Performs a multipixel exponential fits on the data using Phasors

    Parameters
    ----------
    period : float
        The time between samples in `photon_count`
    photon_count : array_like
        An N-dimensional array containing the data to be fit. The last axis is
        the time axis over which the data is fit. Any preceding axes will be
        preserved in the outputs.
    fit_start : {None, int}, optional
        The index of the start of the fit. Some data before this start index is
        required if convolving with the prompt.
        If is None, the fit will begin at index 0 (default is None)
    fit_end : {None, int}, optional
        The index of the end of the fit. If is None, the fit will cover the
        entire time axis of `photon_count`
    Z : {float, array_like}, optional
        The background to be subtracted from the data. If is a float, it will
        be constant for all pixels. (default is 0.0)
    u : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed 'horizontal'
        phasor coordinate. for each fit. To avoid copying, use
        dtype=np.float32. If is None, a new array will be created (default is
        None)
    v : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed 'vertical' phasor
        coordinate. for each fit. To avoid copying, use dtype=np.float32. If is
        None, a new array will be created (default is None)
    taup : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the lifetime calculated from
        the phase change for each fit. To avoid copying, use dtype=np.float32.
        If is None, a new array will be created (default is None)
    taum : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the lifetime calculated from
        the amplitude change for each fit. To avoid copying, use
        dtype=np.float32. If is None, a new array will be created (default is
        None)
    tau : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the average of the other taus
        for each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    fitted : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed fitted plot for
        each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    residuals : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed residuals plot
        for each fit. To avoid copying, use dtype=np.float32. If is None, a new
        array will be created (default is None)
    chisq : {None, numpy.ndarray}, optional
        An N-dimensional array to be filled with the computed reduced chi
        squared value for each fit. To avoid copying, use dtype=np.float32. If
        is None, a new array will be created (default is None)
    fit_mask : {None, array_like}, optional
        An N-dimensional array of bool or 1s and 0s to select which pixels to
        fit. If is None, all pixels will be fit (default is None)
    compute_fitted : bool, optional
        If True, the fitted plot for each fit is computed. Ignored if `fitted` is
        not None (default is True)
    compute_residuals : bool, optional
        If True, the residuals plot for each fit is computed. Ignored if
        `residuals` is not None (default is True)
    compute_chisq : bool, optional
        If True, the reduced chi squared value for each fit is computed.
        Ignored if `chisq` is not None (default is True)

    Returns
    -------
    PhasorResult
        A dataclass containing fields to be acessed using dot (.) notation: u,
        v, taup, taum, tau, fitted, residuals, chisq. For each fit that fails,
        its corresponding outputs are filled with `NaN`
    """

    (
        common_in,
        fitted_out,
        residuals_out,
        chisq_out,
        data_shape,
        npixels,
        referenced_objects,
    ) = _prep_common_params(
        period,
        photon_count,
        fit_start,
        fit_end,
        fit_mask,
        fitted,
        residuals,
        chisq,
        compute_fitted,
        compute_residuals,
        compute_chisq,
        False,
    )
    phasor_in = _PhasorParams()
    try:
        Zf = float(Z)
        # stride 0 array!
        phasor_in.Z, referenced_Z = _as_strided_array(
            [Zf], (1,), (1,), shape_override=(npixels,), strides_override=(0,)
        )
    except TypeError:
        phasor_in.Z, referenced_Z = _as_strided_array(
            Z, data_shape[0:-1], (npixels,)
        )

    phasor_in.u, u_out = _prep_optional_output(u, data_shape[0:-1], (npixels,), fit_mask)
    phasor_in.v, v_out = _prep_optional_output(v, data_shape[0:-1], (npixels,), fit_mask)
    phasor_in.taup, taup_out = _prep_optional_output(
        taup, data_shape[0:-1], (npixels,), fit_mask
    )
    phasor_in.taum, taum_out = _prep_optional_output(
        taum, data_shape[0:-1], (npixels,), fit_mask
    )
    phasor_in.tau, tau_out = _prep_optional_output(
        tau, data_shape[0:-1], (npixels,), fit_mask
    )

    flim_in = _FlimParams()

    flim_in.common = ctypes.pointer(common_in)
    flim_in.phasor = ctypes.pointer(phasor_in)

    _GCI_Phasor(flim_in)

    # Verify that reference was held until after above call
    referenced_Z
    referenced_objects

    fitted_out = _copy_to_provided_output(fitted, fitted_out, data_shape)
    residuals_out = _copy_to_provided_output(
        residuals, residuals_out, data_shape
    )
    chisq_out = _copy_to_provided_output(chisq, chisq_out, data_shape[0:-1])
    u_out = _copy_to_provided_output(u, u_out, data_shape[0:-1])
    v_out = _copy_to_provided_output(v, v_out, data_shape[0:-1])
    taup_out = _copy_to_provided_output(taup, taup_out, data_shape[0:-1])
    taum_out = _copy_to_provided_output(taum, taum_out, data_shape[0:-1])
    tau_out = _copy_to_provided_output(tau, tau_out, data_shape[0:-1])

    return PhasorResult(
        u_out,
        v_out,
        taup_out,
        taum_out,
        tau_out,
        fitted_out,
        residuals_out,
        chisq_out,
    )

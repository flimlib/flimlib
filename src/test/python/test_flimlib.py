#unit tests
import unittest
import flimlib
import numpy as np
import math
import time

samples = 128
a_in = 10.0
tau_in = 1.0
period = 0.04

# The error factor by using a plain exponential instead of integrating it over the time bin (assuming Z = 0)
error_factor = period / (tau_in * (1 - np.exp(-period/tau_in)))

tt = np.linspace(0,(samples-1)*period,samples,dtype=np.float32)
photon_count32 = a_in * np.exp(-tt/tau_in)
photon_count64 = np.asarray(photon_count32,dtype=np.float64)
photon_count2d = []
a_in1d = []
tau_in1d = []
photon_count2d_flip = []
error_factor1d = []
for i in range(samples*samples):
    a = a_in + i * .001
    tau = tau_in + i * .0001
    a_in1d += [a]
    tau_in1d += [tau]
    photon_count2d += [a * np.exp(-tt/tau)]
    photon_count2d_flip += [a * np.exp((tt-(samples-1)*period)/tau)]
    error_factor1d += [period / (tau * (1 - np.exp(-period/tau)))]

photon_count2d_flip = np.asarray(photon_count2d_flip)
trans_strided = np.flip(photon_count2d_flip[0:2])[:,0:-1:2]

linear_const = 1
photon_count_linear = linear_const * tt
    

class Test3Integral(unittest.TestCase):
    def test_output_margin(self):
        result = flimlib.GCI_triple_integral_fitting_engine(
            period, photon_count32)
        
        self.assertAlmostEqual(result.A, a_in * error_factor, 1)
        self.assertAlmostEqual(result.tau, tau_in, 1)
        self.assertAlmostEqual(result.Z, 0.0, 1)
        #self.assertTrue(result.error_code>0) error code no longer says anything about the fit
        self.assertAlmostEqual(result.fitted[0], a_in * error_factor, 1)


    def test_photon_count(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count64) #passing float 64 should be fine
        with self.assertRaises(TypeError):
            flimlib.GCI_triple_integral_fitting_engine(period, 5.7) #must be 1d array
        flimlib.GCI_triple_integral_fitting_engine(period, [1,2])

    def test_instr(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, instr=[1,2,3,4,5])
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, instr=["42"])
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, instr=list(range(300))) # no errors??? ask mark/jenu. ok we should actually just check that it's the right length
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, instr=5)
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, instr="foo")
        with self.assertRaises(TypeError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, instr={"foo": "bar"})

    def test_noise_type_input(self):
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='foo')
        with self.assertRaises(TypeError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type={"foo": "bar"})

    def test_unused_noise_types(self):
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_GAUSSIAN_FIT')
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_MLE')

    def test_noise_const(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_CONST',sig=2)
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_CONST',sig=[5])
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_CONST',sig=[1,2,3,4,5])
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_CONST',sig='foo')
    
    def test_noise_given(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_GIVEN',sig=list(range(samples)))
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_GIVEN',sig=2)
        with self.assertRaises(TypeError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_GIVEN',sig={"foo": "bar"})
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_GIVEN',sig='foo')
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_GIVEN',sig=[1,2,3,4,5])
        
    def test_noise_const_and_given(self):
        #the result should be the same if noise given is constant
        result_const = flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_CONST',sig=1.0)
        result_noise_given = flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_GIVEN',sig=np.ones(samples))
        self.assertAlmostEqual(result_const.chisq,result_noise_given.chisq)
        self.assertAlmostEqual(result_const.A,result_noise_given.A)
        self.assertAlmostEqual(result_const.Z,result_noise_given.Z)
        self.assertAlmostEqual(result_const.tau,result_noise_given.tau)

    def test_noise_poisson_data(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_POISSON_DATA')
        # this should print a warning
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_POISSON_DATA',sig=2)

    def test_noise_poisson_fit(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_POISSON_FIT')
        # this should print a warning
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count32, noise_type='NOISE_POISSON_FIT',sig=2)
    
    def test_chisq_target(self):
        result = flimlib.GCI_triple_integral_fitting_engine(period, np.flip(photon_count32), chisq_target=0.0) #impossible target causes it to run more than once and give up
        #self.assertTrue(result.error_code>1)

class TestPhasor(unittest.TestCase):
    def test_output_margin(self):
        result = flimlib.GCI_Phasor(period, photon_count32)
        self.assertEqual(result.error_code, 0)
        self.assertAlmostEqual(result.fitted[0], a_in, 1)
        self.assertAlmostEqual(result.residuals[0], 0.0, 1)


def dummy_exp_tau(x, param):
    y = param[1]*np.exp(-x / param[2])
    #y = param[0] + param[1]*np.exp(-x / param[2])
    dy_dparam = [1,0,0]
    dy_dparam[1] = np.exp(-x/param[2])
    dy_dparam[2] = param[1]*x*np.exp(-x/param[2])/(param[2]**2)
    return y, dy_dparam

def dummy_exp_tau_predicate(n_param):
    return n_param == 3

# linear fit!
def dummy_linear(x, param):
    y = param[1]*x
    #y = param[0] + param[1]*x
    dy_dparam = [1,0]
    dy_dparam[1] = x
    return y, dy_dparam

def dummy_linear_predicate(n_param):
    return n_param == 2

class TestMarquardt(unittest.TestCase):
    def test_output_margin(self):
        param_in = np.asarray([0,a_in+1,tau_in+1], dtype=np.float32) # slight offset to detect if the fitting works!
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count32, param_in)
        self.assertAlmostEqual(0.0,result.param[0],1)
        self.assertAlmostEqual(a_in,result.param[1],1)
        self.assertAlmostEqual(tau_in,result.param[2],1)

    def test_paramfree(self):
        param_in = np.asarray([0,a_in+1,tau_in+1], dtype=np.float32) # slight offset to detect if the fitting works!
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count32, 
            param_in, paramfree=[1,0,1])
        # second parameter should have been held fixed!
        self.assertAlmostEqual(a_in+1, result.param[1],1)

    def test_restraintype(self):
        with self.assertRaises(ValueError):
            flimlib.GCI_set_restrain_limits([0,0],[0],[0])
        with self.assertRaises(ValueError):
            flimlib.GCI_set_restrain_limits([[0,0],[0,0]],[0,0],[0,0])
        param_in = np.asarray([0,a_in+1,tau_in+1], dtype=np.float32)
        flimlib.GCI_set_restrain_limits([0,1,0],[0,0,0],[0,a_in-1,0])
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count32, param_in, restrain_type='ECF_RESTRAIN_USER')
        self.assertTrue(result.param[1] <= a_in-1)
    
    def test_multiexp_lambda(self):
        lambda_in = 1/tau_in # lambda is the decay rate!
        param_in = np.asarray([0,a_in+1,lambda_in+1], dtype=np.float32) # slight offset to detect if the fitting works!
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count32, param_in, fitfunc=flimlib.GCI_multiexp_lambda)

        self.assertAlmostEqual(a_in, result.param[1], 1)
        self.assertAlmostEqual(lambda_in, result.param[2],1)
    
    def test_user_defined_fitfunc(self):
        param_in = np.asarray([0,a_in+1,tau_in+1], dtype=np.float32) # slight offset to detect if the fitting works!
        fitfunc_in = flimlib.FitFunc(dummy_exp_tau, nparam_predicate=dummy_exp_tau_predicate)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count32, param_in, fitfunc=fitfunc_in)
        self.assertAlmostEqual(0.0,result.param[0],1)
        self.assertAlmostEqual(a_in,result.param[1],1)
        self.assertAlmostEqual(tau_in,result.param[2],1)

        # linear fit! (did not work with the first param being nonzero)
        param_in = np.asarray([0, linear_const+1], dtype=np.float32)
        fitfunc_in = flimlib.FitFunc(dummy_linear, nparam_predicate=dummy_linear_predicate)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count_linear, param_in, fitfunc=fitfunc_in, noise_type='NOISE_CONST',sig=1.0)
        self.assertAlmostEqual(0.0,result.param[0],1)
        self.assertAlmostEqual(linear_const,result.param[1],1)

    def test_nparam(self):
        with self.assertRaises(TypeError):
            param_in = np.asarray([a_in+1,tau_in+1], dtype=np.float32) # forgot the first parameter!
            flimlib.GCI_marquardt_fitting_engine(period, photon_count_linear, param_in)

        # any odd number of parameters greater or equal to 3 should work!
        param_in = np.asarray([0, a_in+1,tau_in+1, 1, 1], dtype=np.float32) # pass 5 parameters
        flimlib.GCI_marquardt_fitting_engine(period, photon_count_linear, param_in)
    
class TestMarquardtMany(unittest.TestCase):
    def test_output_margin(self):
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(samples*samples)], dtype=np.float32) # slight offset to detect if the fitting works!
        start = time.time_ns()
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d, param_in)
        end = time.time_ns()
        print("GCI_marquardt_fitting_engine ran on ", samples*samples, "pixels, took", (end-start)/1000000, "milliseconds")
        for i in range(samples*samples):
            self.assertAlmostEqual(result.param[i][1], a_in1d[i], 1)
            self.assertAlmostEqual(result.param[i][2], tau_in1d[i], 1)

    def test_integer_input(self):
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.int16)
        flimlib.GCI_marquardt_fitting_engine(period, np.asarray(photon_count2d[0:2],dtype=int), param_in)
    
    def test_outputs_modified(self):
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float32)
        fitted_in = np.empty((2,samples), dtype=np.float32)
        residuals_in = np.empty((2,samples), dtype=np.float32)
        chisq_in = np.empty((2,), dtype=np.float32)
        covar_in = np.empty((2,3,3), dtype=np.float32)
        alpha_in = np.empty((2,3,3), dtype=np.float32)
        erraxes_in = np.empty((2,3,3), dtype=np.float32)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, 
            fitted=fitted_in, residuals=residuals_in, chisq=chisq_in, covar=covar_in, alpha=alpha_in, erraxes=erraxes_in)
        self.assertTrue(np.all(result.param != param_in)) # param must NOT be equal
        self.assertTrue(np.all(result.fitted == fitted_in))
        self.assertTrue(np.all(result.residuals == residuals_in))
        self.assertTrue(all(result.chisq == chisq_in))
        self.assertTrue(np.all(result.covar == covar_in))
        self.assertTrue(np.all(result.alpha == alpha_in))
        self.assertTrue(np.all(result.erraxes == erraxes_in))
        # if not float 32 the inputs are still modified but a copy is made
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float64)
        fitted_in = np.empty((2,samples), dtype=np.float64)
        residuals_in = np.empty((2,samples), dtype=np.float64)
        chisq_in = np.empty((2,), dtype=np.float64)
        covar_in = np.empty((2,3,3), dtype=np.float64)
        alpha_in = np.empty((2,3,3), dtype=np.float64)
        erraxes_in = np.empty((2,3,3), dtype=np.float64)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, 
            fitted=fitted_in, residuals=residuals_in, chisq=chisq_in, covar=covar_in, alpha=alpha_in, erraxes=erraxes_in)
        self.assertTrue(np.all(result.param != param_in)) # param must NOT be equal
        self.assertTrue(np.all(result.fitted == fitted_in))
        self.assertTrue(np.all(result.residuals == residuals_in))
        self.assertTrue(all(result.chisq == chisq_in))
        self.assertTrue(np.all(result.covar == covar_in))
        self.assertTrue(np.all(result.alpha == alpha_in))
        self.assertTrue(np.all(result.erraxes == erraxes_in))
    
    def test_wrong_output_type(self):
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float32)
        with self.assertRaises(TypeError):
            fitted_in = np.empty((2,samples), dtype=int) # not compatible
            result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, fitted=fitted_in)
        with self.assertRaises(ValueError):
            fitted_in = np.empty((2,samples+3), dtype=np.float32) # too large
            result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, fitted=fitted_in)
        with self.assertRaises(ValueError):
            fitted_in = np.empty((2,samples-3), dtype=np.float32) # too small
            result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, fitted=fitted_in)

    def test_fit_start_end(self):
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float32)
        start = samples // 10
        end = samples // 10 * 9
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, fit_start=start, fit_end=end)
        self.assertTrue(result.fitted.shape == (2, samples))
        self.assertTrue(np.all(np.diff(result.fitted) < 0)) # monotonic decrease
        self.assertAlmostEqual(result.fitted[0][0], result.param[0][1], 1)

    def test_compute_flags(self):
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float32)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, 
            compute_fitted=False, compute_residuals=False, compute_chisq=False, compute_covar=False, compute_alpha=False, compute_erraxes=False)
        self.assertTrue(result.fitted == result.residuals == result.chisq == result.covar == result.alpha == result.erraxes == None)

    def test_strided(self):
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float32)
        fitted_strided = np.flip(np.empty((2,samples), dtype=np.float32)[:,0:-1:2])
        result = flimlib.GCI_marquardt_fitting_engine(period, trans_strided, param_in, fitted=fitted_strided)
        self.assertEqual(result.fitted.strides, (samples * -4, -8))

    def test_fitmask(self):
        fitted_in = np.ones((2,samples), dtype=np.float32)
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float32)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, fit_mask=[1,0], fitted=fitted_in)
        self.assertTrue(all(result.fitted[1] == 1))
        # 3-dimensional photon_count
        fitted_in = np.ones((2,2,samples), dtype=np.float32)
        param_in = np.asarray([[[0,a_in+1,tau_in+1] for i in range(2)] for i in range(2)], dtype=np.float32)
        result = flimlib.GCI_marquardt_fitting_engine(period, [photon_count2d[0:2], photon_count2d[0:2]], param_in, fit_mask=[[0,1],[1,0]], fitted=fitted_in)
        self.assertTrue(all(result.fitted[0][0] == 1) and all(result.fitted[1][1] == 1))
    
    def test_paramfree(self):
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float32)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, paramfree=[1,0,1])
        self.assertEqual(result.param[0][1], a_in+1)

    def test_size_zero_input(self):
        result = flimlib.GCI_marquardt_fitting_engine(period, [[]], [[]])
        
    def test_chisq_target(self):
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float32)
        with self.assertRaises(TypeError):
            result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, chisq_target="foo")
        
    def test_different_shapes(self):
        param_in = np.asarray([[[0,a_in+1,tau_in+1] for i in range(2)] for i in range(2)], dtype=np.float32)
        result = flimlib.GCI_marquardt_fitting_engine(period, [photon_count2d[0:2], photon_count2d[0:2]], param_in)
        self.assertEqual(result.fitted.shape, (2, 2, samples))
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count32, [0,a_in+1,tau_in+1])
        self.assertEqual(result.fitted.shape, (samples,))

    def test_restraintype(self):
        with self.assertRaises(ValueError):
            flimlib.GCI_set_restrain_limits([0,0],[0],[0])
        with self.assertRaises(ValueError):
            flimlib.GCI_set_restrain_limits([[0,0],[0,0]],[0,0],[0,0])
        param_in = np.asarray([[0,a_in+1,tau_in+1] for i in range(2)], dtype=np.float32)
        flimlib.GCI_set_restrain_limits([0,1,0],[0,0,0],[0,a_in-1,0])
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count2d[0:2], param_in, restrain_type='ECF_RESTRAIN_USER')
        self.assertTrue(all(result.param[1] <= a_in-1))
    
    def test_multiexp_lambda(self):
        lambda_in = 1/tau_in # lambda is the decay rate!
        param_in = np.asarray([0,a_in+1,lambda_in+1], dtype=np.float32)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count32, param_in, fitfunc=flimlib.GCI_multiexp_lambda)
        self.assertAlmostEqual(a_in, result.param[1], 1)
        self.assertAlmostEqual(lambda_in, result.param[2],1)
    
    def test_user_defined_fitfunc(self):
        param_in = np.asarray([0,a_in+1,tau_in+1], dtype=np.float32)
        fitfunc_in = flimlib.FitFunc(dummy_exp_tau, nparam_predicate=dummy_exp_tau_predicate)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count32, param_in, fitfunc=fitfunc_in)
        self.assertAlmostEqual(0.0,result.param[0],1)
        self.assertAlmostEqual(a_in,result.param[1],1)
        self.assertAlmostEqual(tau_in,result.param[2],1)

        # linear fit! (did not work with the first param being nonzero)
        param_in = np.asarray([0, linear_const+1], dtype=np.float32)
        fitfunc_in = flimlib.FitFunc(dummy_linear, nparam_predicate=dummy_linear_predicate)
        result = flimlib.GCI_marquardt_fitting_engine(period, photon_count_linear, param_in, fitfunc=fitfunc_in, noise_type='NOISE_CONST',sig=1.0)
        self.assertAlmostEqual(0.0,result.param[0],1)
        self.assertAlmostEqual(linear_const,result.param[1],1)

    def test_nparam(self):
        with self.assertRaises(TypeError):
            param_in = np.asarray([a_in+1,tau_in+1], dtype=np.float32) # forgot the first parameter!
            flimlib.GCI_marquardt_fitting_engine(period, photon_count32, param_in)

        # any odd number of parameters greater or equal to 3 should work!
        param_in = np.asarray([0, a_in+1,tau_in+1, 1, 1], dtype=np.float32) # pass 5 parameters
        flimlib.GCI_marquardt_fitting_engine(period, photon_count32, param_in)

class Test3IntegralMany(unittest.TestCase):
    def test_output_margin(self):
        start = time.time_ns()
        result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d)
        end = time.time_ns()
        print("GCI_triple_integral_fitting_engine ran on ", samples*samples, "pixels, took", (end-start)/1000000, "milliseconds")
        for i in range(samples*samples):
            self.assertAlmostEqual(result.tau[i], tau_in1d[i], 1)
            self.assertAlmostEqual(result.A[i], a_in1d[i] * error_factor1d[i], 1)
    
    def test_outputs_modified(self):
        fitted_in = np.empty((2,samples), dtype=np.float32)
        residuals_in = np.empty((2,samples), dtype=np.float32)
        chisq_in = np.empty((2,), dtype=np.float32)
        Z_in = np.empty((2,), dtype=np.float32)
        A_in = np.empty((2,), dtype=np.float32)
        tau_in = np.empty((2,), dtype=np.float32)
        result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2],
            fitted=fitted_in, residuals=residuals_in, chisq=chisq_in, Z=Z_in, A=A_in, tau=tau_in)
        self.assertTrue(np.all(result.fitted == fitted_in))
        self.assertTrue(np.all(result.residuals == residuals_in))
        self.assertTrue(all(result.chisq == chisq_in)) 
        self.assertTrue(np.all(result.Z == Z_in))
        self.assertTrue(np.all(result.A == A_in))
        self.assertTrue(np.all(result.tau == tau_in))
        # if not float 32 the inputs are still modified but a copy is made
        fitted_in = np.empty((2,samples), dtype=np.float64)
        residuals_in = np.empty((2,samples), dtype=np.float64)
        chisq_in = np.empty((2,), dtype=np.float64)
        Z_in = np.empty((2,), dtype=np.float64)
        A_in = np.empty((2,), dtype=np.float64)
        tau_in = np.empty((2,), dtype=np.float64)
        result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2],
            fitted=fitted_in, residuals=residuals_in, chisq=chisq_in, Z=Z_in, A=A_in, tau=tau_in)
        self.assertTrue(np.all(result.fitted == fitted_in))
        self.assertTrue(np.all(result.residuals == residuals_in))
        self.assertTrue(all(result.chisq == chisq_in)) 
        self.assertTrue(np.all(result.Z == Z_in))
        self.assertTrue(np.all(result.A == A_in))
        self.assertTrue(np.all(result.tau == tau_in))

    def test_wrong_output_type(self):
        fitted_in = np.empty((2,samples), dtype=np.float64) # float64 is compatible
        result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], fitted=fitted_in)
        self.assertTrue(np.all(result.fitted == fitted_in))
        with self.assertRaises(TypeError):
            fitted_in = np.empty((2,samples), dtype=int) # not compatible
            result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], fitted=fitted_in)
        with self.assertRaises(ValueError):
            fitted_in = np.empty((2,samples+3), dtype=np.float32) # too large
            result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], fitted=fitted_in)
        with self.assertRaises(ValueError):
            fitted_in = np.empty((2,samples-3), dtype=np.float32) # too small
            result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], fitted=fitted_in)
    
    def test_fit_start_end(self):
        start = samples // 10
        end = samples // 10 * 9
        result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], fit_start=start, fit_end=end)
        self.assertEqual(result.fitted.shape, (2, end))
        self.assertTrue(np.all(np.diff(result.fitted) < 0)) # monotonic decrease
        self.assertAlmostEqual(result.fitted[0][0], result.A[0], 1)

    def test_compute_flags(self):
        result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2],
            compute_fitted=False, compute_residuals=False, compute_chisq=False)
        self.assertTrue(result.fitted == result.residuals == result.chisq)

    def test_strided(self):
        fitted_strided = np.flip(np.empty((2,samples), dtype=np.float32)[:,0:-1:2])
        result = flimlib.GCI_triple_integral_fitting_engine(period * 2, trans_strided, fitted=fitted_strided)
        self.assertEqual(result.fitted.strides, (samples * -4, -8))

    def test_size_zero_input(self):
        result = flimlib.GCI_triple_integral_fitting_engine(period, [[]])

    def test_different_shapes(self):
        result = flimlib.GCI_triple_integral_fitting_engine(period, [photon_count2d[0:2], photon_count2d[0:2]])
        self.assertEqual(result.fitted.shape, (2, 2, samples))
        result = flimlib.GCI_triple_integral_fitting_engine(period, photon_count32)
        self.assertEqual(result.fitted.shape, (samples,))

    def test_instr(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], instr=[1,2,3,4,5])
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], instr=["42"])
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], instr=list(range(300))) # no errors??? ask mark/jenu. ok we should actually just check that it's the right length
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], instr=5)
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], instr="foo")
        with self.assertRaises(TypeError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], instr={"foo": "bar"})

    def test_noise_type_input(self):
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='foo')
        with self.assertRaises(TypeError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type={"foo": "bar"})

    def test_unused_noise_types(self):
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_GAUSSIAN_FIT')
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_MLE')

    def test_noise_const(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_CONST',sig=2)
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_CONST',sig=[5])
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_CONST',sig=[1,2,3,4,5])
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_CONST',sig='foo')
    
    def test_noise_given(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_GIVEN',sig=list(range(samples)))
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_GIVEN',sig=2)
        with self.assertRaises(TypeError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_GIVEN',sig={"foo": "bar"})
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_GIVEN',sig='foo')
        with self.assertRaises(ValueError):
            flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_GIVEN',sig=[1,2,3,4,5])
        
    def test_noise_const_and_given(self):
        #the result should be the same if noise given is constant
        result_const = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_CONST',sig=1.0)
        result_noise_given = flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_GIVEN',sig=np.ones(samples))
        self.assertTrue(all(result_const.chisq == result_noise_given.chisq))
        self.assertTrue(all(result_const.A == result_noise_given.A))
        self.assertTrue(all(result_const.Z == result_noise_given.Z))
        self.assertTrue(all(result_const.tau == result_noise_given.tau))

    def test_noise_poisson_data(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_POISSON_DATA')
        # this should print a warning
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_POISSON_DATA',sig=2)

    def test_noise_poisson_fit(self):
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_POISSON_FIT')
        # this should print a warning
        flimlib.GCI_triple_integral_fitting_engine(period, photon_count2d[0:2], noise_type='NOISE_POISSON_FIT',sig=2)

    def test_failed_fit(self):
        result = flimlib.GCI_triple_integral_fitting_engine(period, np.flip(photon_count2d[0:2]))
        self.assertTrue(np.all(np.isnan(result.Z)))
        self.assertTrue(np.all(np.isnan(result.A)))
        self.assertTrue(np.all(np.isnan(result.tau)))
        self.assertTrue(np.all(np.isnan(result.fitted)))
        self.assertTrue(np.all(np.isnan(result.residuals)))

    def test_fit_fail(self):
        result = flimlib.GCI_triple_integral_fitting_engine_many(period, np.flip(photon_count2d[0:2]))
        self.assertTrue(np.all(np.isnan(result.Z)))
        self.assertTrue(np.all(np.isnan(result.A)))
        self.assertTrue(np.all(np.isnan(result.tau)))
        self.assertTrue(np.all(np.isnan(result.fitted)))
        self.assertTrue(np.all(np.isnan(result.residuals)))

class TestPhasorMany(unittest.TestCase):
    def test_output_margin(self):
        start = time.time_ns()
        result = flimlib.GCI_Phasor(period, photon_count2d)
        end = time.time_ns()
        print("GCI_Phasor ran on ", samples*samples, "pixels, took", (end-start)/1000000, "milliseconds")
        self.assertFalse(any(result.chisq > 1)) # for this I assume that the reduced chi squared is a reliable metric
    
    def test_outputs_modified(self):
        u_in = np.empty((2,), dtype=np.float32)
        v_in = np.empty((2,), dtype=np.float32)
        taup_in = np.empty((2,), dtype=np.float32)
        taum_in = np.empty((2,), dtype=np.float32)
        tau_in = np.empty((2,), dtype=np.float32)
        fitted_in = np.empty((2,samples), dtype=np.float32)
        residuals_in = np.empty((2,samples), dtype=np.float32)
        chisq_in = np.empty((2,), dtype=np.float32)
        result = flimlib.GCI_Phasor(period, photon_count2d[0:2],
            fitted=fitted_in, residuals=residuals_in, chisq=chisq_in, u=u_in, v=v_in, taup=taup_in, taum=taum_in, tau=tau_in)
        self.assertTrue(np.all(result.v == v_in))
        self.assertTrue(np.all(result.u == u_in))
        self.assertTrue(np.all(result.taup == taup_in))
        self.assertTrue(np.all(result.taum == taum_in))
        self.assertTrue(np.all(result.tau == tau_in))
        self.assertTrue(np.all(result.fitted == fitted_in))
        self.assertTrue(np.all(result.residuals == residuals_in))
        self.assertTrue(all(result.chisq == chisq_in)) 
        # if not float 32 the inputs are still modified but a copy is made
        u_in = np.empty((2,), dtype=np.float64)
        v_in = np.empty((2,), dtype=np.float64)
        taup_in = np.empty((2,), dtype=np.float64)
        taum_in = np.empty((2,), dtype=np.float64)
        tau_in = np.empty((2,), dtype=np.float64)
        fitted_in = np.empty((2,samples), dtype=np.float64)
        residuals_in = np.empty((2,samples), dtype=np.float64)
        chisq_in = np.empty((2,), dtype=np.float64)
        result = flimlib.GCI_Phasor(period, photon_count2d[0:2],
            fitted=fitted_in, residuals=residuals_in, chisq=chisq_in, u=u_in, v=v_in, taup=taup_in, taum=taum_in, tau=tau_in)
        self.assertTrue(np.all(result.v == v_in))
        self.assertTrue(np.all(result.u == u_in))
        self.assertTrue(np.all(result.taup == taup_in))
        self.assertTrue(np.all(result.taum == taum_in))
        self.assertTrue(np.all(result.tau == tau_in))
        self.assertTrue(np.all(result.fitted == fitted_in))
        self.assertTrue(np.all(result.residuals == residuals_in))
        self.assertTrue(all(result.chisq == chisq_in)) 

    def test_wrong_output_type(self):
        fitted_in = np.empty((2,samples), dtype=np.float64) # float64 is compatible
        result = flimlib.GCI_Phasor(period, photon_count2d[0:2], fitted=fitted_in)
        self.assertTrue(np.all(result.fitted == fitted_in))
        with self.assertRaises(TypeError):
            fitted_in = np.empty((2,samples), dtype=int) # not compatible
            result = flimlib.GCI_Phasor(period, photon_count2d[0:2], fitted=fitted_in)
        with self.assertRaises(ValueError):
            fitted_in = np.empty((2,samples+3), dtype=np.float32) # too large
            result = flimlib.GCI_Phasor(period, photon_count2d[0:2], fitted=fitted_in)
        with self.assertRaises(ValueError):
            fitted_in = np.empty((2,samples-3), dtype=np.float32) # too small
            result = flimlib.GCI_Phasor(period, photon_count2d[0:2], fitted=fitted_in)

    def test_fit_start_end(self):
        start = samples // 10
        end = samples // 10 * 9
        result = flimlib.GCI_Phasor(period, photon_count2d[0:2], fit_start=start, fit_end=end)
        self.assertTrue(np.all(np.diff(result.fitted[:,start:end]) < 0)) # monotonic decrease
        self.assertTrue(result.fitted.shape == (2, end))

    def test_compute_flags(self):
        result = flimlib.GCI_Phasor(period, photon_count2d[0:2],
            compute_fitted=False, compute_residuals=False, compute_chisq=False)
        self.assertTrue(result.fitted == result.residuals == result.chisq)

    def test_strided(self):
        fitted_strided = np.flip(np.empty((2,samples), dtype=np.float32)[:,0:-1:2])
        result = flimlib.GCI_Phasor(period, trans_strided, fitted=fitted_strided)
        self.assertEqual(result.fitted.strides, (samples * -4, -8))

    def test_size_zero_input(self):
        result = flimlib.GCI_Phasor(period, [[]])
        result = flimlib.GCI_Phasor(period, [])
    
    def test_different_shapes(self):
        result = flimlib.GCI_Phasor(period, [photon_count2d[0:2], photon_count2d[0:2]])
        self.assertEqual(result.fitted.shape, (2, 2, samples))
        result = flimlib.GCI_Phasor(period, photon_count32)
        self.assertEqual(result.fitted.shape, (samples,))
        
# TODO find fixes made to the many implementation that ought to go into single
# TODO replace non-many functions with many implementation
# TODO pass error codes through on single fits?

if __name__ == '__main__':
    unittest.main() 
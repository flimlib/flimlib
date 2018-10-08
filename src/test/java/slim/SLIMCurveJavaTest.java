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

package slim;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.Random;

import org.junit.Before;
import org.junit.Test;

/**
 * Tests {@link SLIMCurve}.  Test parameters currently based 
 * off of values in files data.data and test.ini
 * 
 * @author Zach Petersen, Dasong Gao
 */
public class SLIMCurveJavaTest {
	public static final int TEST_SIZE = 100;
	public static final int SEED = 0x1226;
	final static float tolerance = 1e-5f;
	final static Random rng = new Random(SEED);
	final static int DEFAULT_RET = -100;
	final double xincd = 0.048828125; 
	final float xinc = 0.048828125f; 
	final float y[] = {
			40.000000f, 45.000000f, 34.000000f, 54.000000f, 44.000000f,
			53.000000f, 47.000000f, 56.000000f, 62.000000f, 66.000000f,
			82.000000f, 90.000000f, 108.000000f, 122.000000f, 323.000000f,
			1155.000000f, 4072.000000f, 8278.000000f, 11919.000000f, 13152.000000f,
			13071.000000f, 12654.000000f, 11946.000000f, 11299.000000f, 10859.000000f,
			10618.000000f, 10045.000000f, 9576.000000f, 9208.000000f, 9113.000000f,
			8631.000000f, 8455.000000f, 8143.000000f, 8102.000000f, 7672.000000f,
			7384.000000f, 7463.000000f, 7254.000000f, 6980.000000f, 6910.000000f,
			6411.000000f, 6355.000000f, 6083.000000f, 5894.000000f, 5880.000000f,
			5735.000000f, 5528.000000f, 5343.000000f, 5224.000000f, 4933.000000f,
			5026.000000f, 4914.000000f, 4845.000000f, 4681.000000f, 4426.000000f,
			4485.000000f, 4271.000000f, 4295.000000f, 4183.000000f, 3989.000000f,
			3904.000000f, 3854.000000f, 3801.000000f, 3600.000000f, 3595.000000f,
			3434.000000f, 3457.000000f, 3291.000000f, 3280.000000f, 3178.000000f,
			3132.000000f, 2976.000000f, 2973.000000f, 2940.000000f, 2770.000000f,
			2969.000000f, 2851.000000f, 2702.000000f, 2677.000000f, 2460.000000f,
			2536.000000f, 2528.000000f, 2347.000000f, 2382.000000f, 2380.000000f,
			2234.000000f, 2251.000000f, 2208.000000f, 2115.000000f, 2136.000000f,
			2000.000000f, 2006.000000f, 1970.000000f, 1985.000000f, 1886.000000f,
			1898.000000f, 1884.000000f, 1744.000000f, 1751.000000f, 1797.000000f,
			1702.000000f, 1637.000000f, 1547.000000f, 1526.000000f, 1570.000000f,
			1602.000000f, 1557.000000f, 1521.000000f, 1417.000000f, 1391.000000f,
			1332.000000f, 1334.000000f, 1290.000000f, 1336.000000f, 1297.000000f,
			1176.000000f, 1189.000000f, 1220.000000f, 1209.000000f, 1217.000000f,
			1140.000000f, 1079.000000f, 1059.000000f, 1074.000000f, 1061.000000f,
			1013.000000f, 1075.000000f, 1021.000000f, 1012.000000f, 940.000000f,
			982.000000f, 866.000000f, 881.000000f, 901.000000f, 883.000000f,
			893.000000f, 845.000000f, 819.000000f, 831.000000f, 758.000000f,
			794.000000f, 779.000000f, 772.000000f, 779.000000f, 791.000000f,
			729.000000f, 732.000000f, 687.000000f, 690.000000f, 698.000000f,
			661.000000f, 647.000000f, 668.000000f, 642.000000f, 619.000000f,
			629.000000f, 656.000000f, 579.000000f, 579.000000f, 600.000000f,
			563.000000f, 584.000000f, 531.000000f, 554.000000f, 526.000000f,
			484.000000f, 530.000000f, 515.000000f, 493.000000f, 502.000000f,
			479.000000f, 445.000000f, 439.000000f, 466.000000f, 431.000000f,
			423.000000f, 451.000000f, 412.000000f, 415.000000f, 393.000000f,
			404.000000f, 390.000000f, 398.000000f, 352.000000f, 394.000000f,
			376.000000f, 338.000000f, 377.000000f, 367.000000f, 355.000000f,
			352.000000f, 375.000000f, 339.000000f, 347.000000f, 316.000000f,
			295.000000f, 322.000000f, 311.000000f, 294.000000f, 304.000000f,
			264.000000f, 293.000000f, 294.000000f, 283.000000f
			};
	final double[] yd = asDoubleArr(y);
	final int fit_start = 10;
	final int fit_end = 203;
	final float instr[] = {
			0.00911042001098394393920898437500000f,
			0.03882249817252159118652343750000000f,
			0.13171100616455078125000000000000000f,
			0.25238901376724243164062500000000000f,
			0.27722999453544616699218750000000000f,
			0.18023300170898437500000000000000000f,
			0.07927619665861129760742187500000000f,
			0.03122770041227340698242187500000000f
	};  
	final int ndata = y.length;
	final double instrd[] = asDoubleArr(instr);
	final NoiseType noise = NoiseType.NOISE_POISSON_FIT;
	final double sigd[] = new double[ndata];
	final float sig[] = new float[ndata];
	final double zd[] = {0};
	final float z[] = {0};
	final double ad[]  = {1000};
	final float a[]  = {1000};
	final double taud[] = {2};
	final float tau[] = {2};
	final double fittedd[] = new double[yd.length];  
	//final float fitted[] = new float[y.length]; 
	final double chisquared[] = {0.01}; 
	final float chisquare[] = {0}; 
	int chisq_percent = 95;
	final double chisq_targetd = 237.5;
	final float chisq_target = 237.5f;
	final int nInstr = instr.length; 
	final float param0[] = { 0, 1000, 2};    // z, a, tau
	final float param1[] = { 0, 1000, 0.5f}; // z, a, lambda
	final int paramfree[] = {1, 1, 1};
	final int nparam = 3;	
	final double chisq_transd[] = {}; 
	final float chisq_trans[] = {}; 
	final int ntrans = -1;
	final int ftype = 0;
	final float chisq_delta = 0.01f;
	final double chisq_deltad = 0.01;
	final int drop_bad_transients = 0;
	final int gparam[] = {};
	int[] df = {0};
	double[] chisq_globald = {0};
	float[] chisq_global = {0};
	Float2DMatrix trans = new Float2DMatrix(new float[][]{{2, 3, 4}, {5, 6, 7}});
	Float2DMatrix param2d = new Float2DMatrix(new float[2][3]);
	Float2DMatrix covar = new Float2DMatrix(new float[nparam][nparam]);
	Float2DMatrix alpha = new Float2DMatrix(new float[nparam][nparam]);
	Float2DMatrix erraxes = new Float2DMatrix(new float[nparam][nparam]);
	float[] fitted = new float[ndata];
	float[] residuals = new float[ndata];
	RestrainType restrain = RestrainType.ECF_RESTRAIN_DEFAULT;
	
	private static double[] asDoubleArr(final float[] floatIn) {
		double[] doubleOut = new double[floatIn.length];
		for (int i = 0; i < floatIn.length; i++)
			doubleOut[i] = (float) floatIn[i];
		return doubleOut;
	}
	
	private static void compareWithLMA(float xinc, float[] y, int fit_start, int fit_end, float[] instr, 
			NoiseType noise, float[] sig, float[] param, int[] paramfree, RestrainType restrain, FitFunc func, float[] fitted, 
			float[] residuals, float[] chisquare, Float2DMatrix covar, Float2DMatrix alpha, Float2DMatrix erraxes, float chisq_target, float chisq_delta, 
			int chisq_percent,
			int actualRet, float actualParam[], Float2DMatrix actualCovar, Float2DMatrix actualAlpha, Float2DMatrix actualErraxes,
			float[] actualFitted, float[] actualResiduals, float[] actualChisquare) {
		int ret = SLIMCurve.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param, paramfree, restrain, func, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent);
		if (actualRet != DEFAULT_RET)
			assertEquals("lma value incorrect", ret, actualRet);
		// pull out matrices, those should have the same size
		int col = covar.getNcol(), row = covar.getNrow();
		float[][] cov = covar.asArray(), aCov = null;
		if (actualCovar != null)
			aCov = actualCovar.asArray();
		float[][] alph = alpha.asArray(), aAlph = null;
		if (actualAlpha != null)
			aAlph = actualAlpha.asArray();
		float[][] err = erraxes.asArray(), aErr = null;
		if (actualErraxes != null)
			aErr = actualErraxes.asArray();
		for (int i = 0; i < TEST_SIZE; i++) {
			int idx = rng.nextInt(param.length);
			if (actualParam != null)
				assertEqualsScaled("param matrix incorrect", param[idx], actualParam[idx], tolerance);
			idx = rng.nextInt(col * row - 1);
			if (aCov != null)
				assertEqualsScaled("covariance matrix incorrect", cov[idx / row][idx % col], aCov[idx / row][idx % col], tolerance);
			if (aAlph != null)
				assertEqualsScaled("alpha matrix incorrect", alph[idx / row][idx % col], aAlph[idx / row][idx % col], tolerance);
			if (aErr != null)
				assertEqualsScaled("erraxes matrix incorrect", err[idx / row][idx % col], aErr[idx / row][idx % col], tolerance);
			idx = rng.nextInt(fitted.length);
			if (actualFitted != null)
				assertEqualsScaled("fitted incorrect", fitted[idx], actualFitted[idx], tolerance);
			if (actualResiduals != null)
				assertEqualsScaled("residuals incorrect", residuals[idx], actualResiduals[idx], tolerance);
		}
		if (actualChisquare != null)
			assertEqualsScaled("Chi square incorrect", chisquare[0], actualChisquare[0], tolerance);
	}
	
	private static void assertEqualsScaled(String message, double expected, double actual, double delta) {
		assertEquals(message, expected, actual, Math.abs(expected * delta));
	}
	
	private static void assertArrayEqualsScaled(String message, float[] expecteds, float[] actuals, float delta) {
		float maxValue = expecteds[0]; 
		for(float n : expecteds)
			maxValue = n > maxValue ? n : maxValue;
		assertArrayEquals(message, expecteds, expecteds, Math.abs(maxValue * delta));
	}
	
//	private static void assertArrayEqualsScaled(String message, double[] expecteds, double[] actuals, double delta) {
//		double maxValue = expecteds[0]; 
//		for(double n : expecteds)
//			maxValue = n > maxValue ? n : maxValue;
//		assertArrayEqualsScaled(message, expecteds, actuals, Math.abs(maxValue * delta));
//	}
	
	@Before
	public void setUp() {
		for (int i = 0; i < sig.length; i++)
			sigd[i] = sig[i] = 1.0f;
	}

	/** Tests {@link SLIMCurve#fitRLD}. */
	@Test
	public void testFitRLD() {
		final int rld = SLIMCurve.RLD_fit(xincd, yd, fit_start, fit_end, instrd, 5, sigd, zd, ad, taud, fittedd, chisquared, chisq_targetd);
		assertEquals("rld incorrect", 2, rld);
		assertEqualsScaled("a incorrect", 11823.900390625, ad[0], tolerance);
		assertEqualsScaled("tau incorrect", 2.350015878677368, taud[0], tolerance);
		assertEqualsScaled("z incorrect", 105.92437744140625, zd[0], tolerance);
		assertEqualsScaled("Chi square incorrect", 78742.3359375, chisquared[0], tolerance);
	}

	/** Tests {@link SLIMCurve#fitLMA}. */
	@Test
	public void testFitLMA() {
		final double param[] = { 5.263157800072804E-5, 0.0, 1000.0, 2.0 };
		final int lma = SLIMCurve.LMA_fit(xincd, yd, fit_start, fit_end, instrd, 5, sigd, param, paramfree, fittedd, chisquared, chisq_targetd, chisquared[0]);
		assertEquals("lma incorrect", 58, lma);
		assertEqualsScaled("a incorrect", 10501.380859375, param[2], tolerance);
		assertEqualsScaled("tau incorrect", 2.8925282955169678, param[3], tolerance);
		assertEqualsScaled("z incorrect", -163.77395629882812, param[1], tolerance);
		assertEqualsScaled("Chi square incorrect", 83087.8671875, chisquared[0], tolerance);
	
	}
	
	/** Tests {@link SLIMCurve#GCI_marquardt_global_exps_instr}. */
	@Test
	public void testGCIGlobalWrapperCall() {		
		int result = DEFAULT_RET;
		Float2DMatrix fitted = new Float2DMatrix(new float[1][ndata]);
		Float2DMatrix residuals = new Float2DMatrix(new float[1][ndata]);
		result = SLIMCurve.GCI_marquardt_global_exps_instr(xinc, trans, fit_start, fit_end, instr, 
				NoiseType.swigToEnum(5), sig, FitType.FIT_GLOBAL_MULTIEXP, param2d, paramfree, restrain, chisq_delta, 
				fitted, residuals, chisq_trans, chisq_global, df, 1);
		assertTrue(result != DEFAULT_RET);
	}
	
	/** Tests {@link SLIMCurve#GCI_marquardt_global_generic_instr}. */
	@Test
	public void testGCIGlobalGenericCall() {
		int result = DEFAULT_RET;
		Float2DMatrix fitted = new Float2DMatrix(new float[1][ndata]);
		Float2DMatrix residuals = new Float2DMatrix(new float[1][ndata]);
		
		result = SLIMCurve.GCI_marquardt_global_generic_instr(xinc, trans, fit_start, fit_end,
				instr, noise, sig, param2d, paramfree, gparam, restrain, chisq_delta, FitFunc.GCI_MULTIEXP_LAMBDA, fitted, 
				residuals, chisq_trans, chisq_global, df);
		assertTrue(result != DEFAULT_RET);
	}
	
	/** Tests {@link SLIMCurve#GCI_Phasor}. */
	@Test
	public void testPhasorWrapperCall() {
		float u[] = {0}; 
		float v[] = {0};
		float taup[] = {0};
		float taum[] = {0};
		
		int ret = SLIMCurve.GCI_Phasor(xinc, y, fit_start, fit_end, z, u, v, taup, taum, tau, fitted, residuals, chisquare);
		assertEquals("phasor failed", ret, 0);
		assertEqualsScaled("z incorrect", 0.0f, z[0], tolerance);
		assertEqualsScaled("u incorrect", 0.24798244f, u[0], tolerance);
		assertEqualsScaled("v incorrect", 0.49327368f, v[0], tolerance);
		assertEqualsScaled("taup incorrect", 2.9834206f, taup[0], tolerance);
		assertEqualsScaled("taum incorrect", 2.2650633f, taum[0], tolerance);
		assertEqualsScaled("Chi square incorrect", 81201.26f, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", new float[] {293.38422f, 287.97586f, 282.66714f, 277.45636f, 272.34155f},
				Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1), tolerance);
		assertArrayEqualsScaled("residuals incorrect", new float[] {0.6157837f, 16.02414f, -18.667145f, 15.54364f, 21.658447f},
				Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1), tolerance);
		assertEqualsScaled("period incorrect", 9.423828125f, SLIMCurve.GCI_Phasor_getPeriod(), tolerance);
	}
	
	/** Tests {@link SLIMCurve#GCI_triple_integral_fitting_engine}. */
	@Test
	public void testGCI_triple_integral_fitting_engine() {
		chisquare[0] = z[0] = a[0] = tau[0] = 0;
		int ret0 = SLIMCurve.GCI_triple_integral_fitting_engine(xinc, y, fit_start, fit_end, instr, noise, sig, z, a,
				tau, fitted, residuals, chisquare, chisq_target);
		assertEquals("rld value incorrect", 2, ret0);
		assertEqualsScaled("z incorrect", 105.92438, z[0], tolerance);
		assertEqualsScaled("a incorrect", 11823.900390625, a[0], tolerance);
		assertEqualsScaled("Chi square incorrect", 78742.3359375, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1),
				new float[] {241.53043f, 236.57343f, 231.71307f, 226.94742f, 222.27469f}, tolerance);
		assertArrayEqualsScaled("residuals incorrect", Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1),
				new float[] {52.469574f, 67.426575f, 32.28694f, 66.05258f, 71.72531f}, tolerance);
		
		chisquare[0] = z[0] = a[0] = tau[0] = 0;
		int ret1 = SLIMCurve.GCI_triple_integral_fitting_engine(xinc, y, fit_start, fit_end, instr, NoiseType.NOISE_GAUSSIAN_FIT, sig, z, a,
				tau, fitted, residuals, chisquare, chisq_target);
		assertEquals("rld value incorrect", 2, ret1);
		assertEqualsScaled("z incorrect", 105.92438, z[0], tolerance);
		assertEqualsScaled("a incorrect", 11823.900390625, a[0], tolerance);
		assertEqualsScaled("Chi square incorrect", 78742.3359375, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", new float[] {241.53043f, 236.57343f, 231.71306f, 226.94742f, 222.27469f},
				Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1), tolerance);
		assertArrayEqualsScaled("residuals incorrect", new float[] {52.469574f, 67.426575f, 32.28694f, 66.05258f, 71.72531f},
				Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1), tolerance);
	}
	
	/** Tests {@link SLIMCurve#GCI_marquardt_fitting_engine}. */
	@Test
	public void testGCI_marquardt_fitting_engine() {
		int ret0 = SLIMCurve.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param0, paramfree, restrain, FitFunc.GCI_MULTIEXP_TAU, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent);
		assertEquals("lma value incorrect", 57, ret0);
		assertArrayEqualsScaled("param matrix incorrect", new float[] { -163.77396f, 10501.38f, 2.8925283f }, param0, tolerance);
		assertEqualsScaled("covariance matrix incorrect", -0.18844598531723, covar.asArray()[1][2], tolerance);
		assertEqualsScaled("alpha matrix incorrect", 19.542074, alpha.asArray()[1][2], tolerance);
		assertEqualsScaled("erraxes matrix incorrect", 1.3784504E-6, erraxes.asArray()[1][2], tolerance);
		assertEqualsScaled("Chi square incorrect", 83087.8984375, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", new float[] {231.93842f, 225.31448f, 218.80151f, 212.39758f, 206.1008f},
				Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1), tolerance);
		assertArrayEqualsScaled("residuals incorrect", new float[] {62.061584f, 78.68552f, 45.198486f, 80.60242f, 87.8992f},
				Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1), tolerance);
		
		// z, a should be approximately the same; lambda * tau should be approximately 1
		int ret1 = SLIMCurve.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param1, paramfree, restrain, FitFunc.GCI_MULTIEXP_LAMBDA, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent);
		assertEquals("lma value incorrect", 40, ret1);
		assertArrayEqualsScaled("param matrix incorrect", new float[] { -158.19781f, 10513.9375f, 0.34687048f }, param1, tolerance);
		assertEqualsScaled("covariance matrix incorrect", 0.022667186, covar.asArray()[1][2], tolerance);
		assertEqualsScaled("alpha matrix incorrect", -161.95181, alpha.asArray()[1][2], tolerance);
		assertEqualsScaled("erraxes matrix incorrect", -1.9851784E-8, erraxes.asArray()[1][2], tolerance);
		assertEqualsScaled("Chi square incorrect", 82995.6953125, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", new float[] {233.6828f, 227.10132f, 220.63046f, 214.26831f, 208.01291f},
				Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1), tolerance);
		assertArrayEqualsScaled("residuals incorrect", new float[] {60.3172f, 76.89868f, 43.369537f, 78.73169f, 85.98709f},
				Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1), tolerance);
	}
	
	/** Tests {@link SLIMCurve#GCI_EcfModelSelectionEngine}. */
	@Test
	public void testGCI_EcfModelSelectionEngine() {
		// setup arguments
		final float param0_ref[] = { 0, 1000, 2};
		final float param1_ref[] = { 0, 1000, 0.5f};
		final float chisq_diff[] = { 0 };
		final int model[] = { -1 };
		// create two models to compare
		DecayModel paramsandfits[] = { 
				new DecayModel(FitFunc.GCI_MULTIEXP_TAU,    param0, paramfree,
						restrain, fit_end + 1, chisq_target, chisq_delta, chisq_percent),
				new DecayModel(FitFunc.GCI_MULTIEXP_LAMBDA, param1, paramfree,
						restrain, fit_end + 1, chisq_target, chisq_delta, chisq_percent)
		};
		int ret = SLIMCurve.GCI_EcfModelSelectionEngine(xinc, y, fit_start, fit_end, instr, noise, sig,
				paramsandfits, chisq_diff, model);
		assertEquals("selection failed", ret, 0);
		compareWithLMA(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param0_ref, paramfree, restrain, FitFunc.GCI_MULTIEXP_TAU, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent, DEFAULT_RET, paramsandfits[0].getParams(),
				paramsandfits[0].getCovar(), paramsandfits[0].getAlpha(), paramsandfits[0].getErraxes(),
				paramsandfits[0].getFitted(), paramsandfits[0].getResiduals(),
				new float[] { paramsandfits[0].getChisq() * (fit_end - fit_start - 3) });
		compareWithLMA(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param1_ref, paramfree, restrain, FitFunc.GCI_MULTIEXP_LAMBDA, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent, DEFAULT_RET, paramsandfits[1].getParams(),
				paramsandfits[1].getCovar(), paramsandfits[1].getAlpha(), paramsandfits[1].getErraxes(),
				paramsandfits[1].getFitted(), paramsandfits[1].getResiduals(),
				new float[] { paramsandfits[1].getChisq() * (fit_end - fit_start - 3) });
		assertEqualsScaled("chisq_diff incorrect", 97.53906, chisq_diff[0], tolerance);
		assertEquals("model incorrect", 2, model[0]);
	}
	
	/** Tests builtin {@link FitFunc}s */
	@Test
	public void testBuiltinFitFunc() {
		// This test compares the behavior of GCI_MULTIEXP_LAMBDA
		// in java and c
		FitFunc lambda = FitFunc.GCI_MULTIEXP_LAMBDA;
		float[] expectedDy_dparam = new float[7];
		float expectedY = 0;
		float ex;
		float[] y = new float[1];
		for (int i = 0; i < TEST_SIZE; i++) {
			// must be odd
			int nparam = rng.nextInt(4) * 2 + 1;
			float[] param = new float[nparam];
			float x = rng.nextFloat() * 5;
			for (int j = 0; j < nparam; j++)
				param[j] = rng.nextFloat() * 5;
			float[] dy_dparam = new float[nparam];
			
			// copied from native code
			for (int j = 1; j < nparam - 1; j += 2) {
				expectedDy_dparam[j] = ex = (float) Math.exp(-param[j + 1] * x);
				ex *= param[j];
				expectedY += ex;
				expectedDy_dparam[j + 1] = -ex * x;
			}
			
			lambda.fit(x, param, y, dy_dparam);
			assertEqualsScaled("y incorrect", expectedY, y[0], tolerance);
			assertArrayEqualsScaled("expectedDy_dparam incorrect",
					expectedDy_dparam, dy_dparam, tolerance);

			y[0] = 0;
			expectedY = 0;
		}
	}
	// TODO: get rid of KEEP_ALLOC
	// TODO: delete unwanted comments
	// TODO: optimize callback to avoid allocating arrays in each single call
	/** Tests {@link FitFunc} callback */
	@Test
	public void testCustomFitFunc() {
		// This test guarantees java FitFuncs behave the same as C FitFuncs
		// outputs for reference
		final float param1_ref[] = { 0, 1000, 0.5f};
		final float fitted_ref[] = new float[ndata];
		final float residuals_ref[] = new float[ndata];
		final float chisquare_ref[] = {0}; 
		final Float2DMatrix covar_ref = new Float2DMatrix(new float[nparam][nparam]);
		final Float2DMatrix alpha_ref = new Float2DMatrix(new float[nparam][nparam]);
		final Float2DMatrix erraxes_ref = new Float2DMatrix(new float[nparam][nparam]);
		
		// copied from native code
		FitFunc customeLambda = new FitFunc() {
			public void fit(float x, float param[], float y[], float dy_dparam[]) {
				float ex;
				y[0] = 0;
				for (int i = 1; i < nparam - 1; i += 2) {
					dy_dparam[i] = ex = (float) Math.exp(-param[i + 1] * x);
					ex *= param[i];
					y[0] += ex;
					dy_dparam[i + 1] = -ex * x;
				}
			}
		};
		
		int ret = SLIMCurve.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param1, paramfree, restrain, customeLambda, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent);
		compareWithLMA(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param1_ref, paramfree, restrain, FitFunc.GCI_MULTIEXP_LAMBDA, fitted_ref, 
				residuals_ref, chisquare_ref, covar_ref, alpha_ref, erraxes_ref, chisq_target, chisq_delta, 
				chisq_percent,
				ret, param1, covar, alpha, erraxes, fitted, residuals, chisquare);
	}
	
	/** Tests {@link Float2DMatrix}.  */
	@Test
	public void testFloat2DMatrix() {
		for (int i = 0; i < TEST_SIZE; i++) {
			int row = rng.nextInt(10) + 1;
			int col = rng.nextInt(10) + 1;
			// generates random 2D array
			float[][] arr = new float[row][col];
			for (int j = 0; j < row; j++)
				for (int k = 0; k < col; k++)
					arr[j][k] = rng.nextFloat();
			// feed into C and back
			Float2DMatrix mat = new Float2DMatrix(arr);
			float[][] arrOut = mat.asArray();
			// test if C array can be freed
			mat.delete();
			assertEquals(arr.length, arrOut.length);
			for (int j = 0; j < arr.length; j++)
				assertArrayEqualsScaled("Float2DMatrix not equal", arr[j], arrOut[j], tolerance);
		}
	}
	
	/** Tests {@link Float2DMatrix}.  */
	@Test
	public void testInt2DMatrix() {
		for (int i = 0; i < TEST_SIZE; i++) {
			int row = rng.nextInt(10) + 1;
			int col = rng.nextInt(10) + 1;
			// generates random 2D array
			int[][] arr = new int[row][col];
			for (int j = 0; j < row; j++)
				for (int k = 0; k < col; k++)
					arr[j][k] = rng.nextInt();
			// feed into C and back
			Int2DMatrix mat = new Int2DMatrix(arr);
			int[][] arrOut = mat.asArray();
			// test if C array can be freed
			mat.delete();
			assertEquals(arr.length, arrOut.length);
			for (int j = 0; j < arr.length; j++)
				assertArrayEquals("Int2DMatrix not equal", arr[j], arrOut[j]);
		}
	}

}

/*
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2025 University of Oxford and Board of Regents of the
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

package flimlib;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.Random;

import org.junit.Before;
import org.junit.Test;

/**
 * Tests {@link FLIMLib}.  Test parameters currently based 
 * off of values in files data.data and test.ini
 * 
 * @author Zach Petersen
 * @author Dasong Gao
 */
public class FLIMLibJavaTest {
	public static final int TEST_SIZE = 100;
	public static final int SEED = 0x1226;
	final static float tolerance = 1e-5f;
	final static Random rng = new Random(SEED);
	final static int DEFAULT_RET = -100;
	final double xincd = 0.058628749f;//0.048828125; 
	final float xinc = 0.058628749f;//0.048828125f; 
	final float y[] = {
		// the first 38 data points are unused
/*		43, 39, 50, 46, 56, 63, 62, 74, 60, 72, 58, 47, 41, 69, 69, 58,
 *		55, 37, 55, 50, 52, 59, 51, 52, 51, 50, 53, 40, 45, 34, 54, 44,
 *		53, 47, 56, 62, 66, 82, */90, 108, 122, 323, 1155, 4072, 8278, 11919, 13152, 13071,
		12654, 11946, 11299, 10859, 10618, 10045, 9576, 9208, 9113, 8631, 8455, 8143, 8102, 7672, 7384, 7463,
		7254, 6980, 6910, 6411, 6355, 6083, 5894, 5880, 5735, 5528, 5343, 5224, 4933, 5026, 4914, 4845,
		4681, 4426, 4485, 4271, 4295, 4183, 3989, 3904, 3854, 3801, 3600, 3595, 3434, 3457, 3291, 3280,
		3178, 3132, 2976, 2973, 2940, 2770, 2969, 2851, 2702, 2677, 2460, 2536, 2528, 2347, 2382, 2380,
		2234, 2251, 2208, 2115, 2136, 2000, 2006, 1970, 1985, 1886, 1898, 1884, 1744, 1751, 1797, 1702,
		1637, 1547, 1526, 1570, 1602, 1557, 1521, 1417, 1391, 1332, 1334, 1290, 1336, 1297, 1176, 1189,
		1220, 1209, 1217, 1140, 1079, 1059, 1074, 1061, 1013, 1075, 1021, 1012, 940, 982, 866, 881,
		901, 883, 893, 845, 819, 831, 758, 794, 779, 772, 779, 791, 729, 732, 687, 690,
		698, 661, 647, 668, 642, 619, 629, 656, 579, 579, 600, 563, 584, 531, 554, 526,
		484, 530, 515, 493, 502, 479, 445, 439, 466, 431, 423, 451, 412, 415, 393, 404,
		390, 398, 352, 394, 376, 338, 377, 367, 355, 352, 375, 339, 347, 316, 295, 322,
		311, 294, 304, 264, 293, 294, 283, 278, 302, 253, 259, 252, 278, 254, 245, 246,
		242, 226, 241, 222, 198, 197, 245, 221, 228, 224, 216, 174, 166, 163, 127, 122
	};
	final double[] yd = asDoubleArr(y);
	final int fit_start = 46 - 38;
	final int fit_end = 255 - 38;
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
	final double sigd[] = null;
	final float sig[] = null;
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
	final float chisq_target = fit_end - fit_start + instr.length - 3;
	final double chisq_targetd = chisq_target;
	final float param0[] = { 160.39937f, 13382.208f, 2.4456801f };    // z, a, tau
	final float param1[] = { 160.39937f, 13382.208f, 1 / 2.4456801f}; // z, a, lambda
	final boolean paramfree[] = {true, true, true};
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
			NoiseType noise, float[] sig, float[] param, boolean[] paramfree, RestrainType restrain, FitFunc func, float[] fitted,
			float[] residuals, float[] chisquare, Float2DMatrix covar, Float2DMatrix alpha, Float2DMatrix erraxes, float chisq_target, float chisq_delta, 
			int chisq_percent,
			int actualRet, float actualParam[], Float2DMatrix actualCovar, Float2DMatrix actualAlpha, Float2DMatrix actualErraxes,
			float[] actualFitted, float[] actualResiduals, float[] actualChisquare) {
		int ret = FLIMLib.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param, paramfree, restrain, func, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent);
		if (actualRet != DEFAULT_RET)
			assertEquals("return value incorrect", ret, actualRet);
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
//		for (int i = 0; i < sig.length; i++)
//			sigd[i] = sig[i] = 0f;
	}

	/** Tests {@link FLIMLib#fitRLD}. */
	@Test
	public void testFitRLD() {
		final int rld = FLIMLib.RLD_fit(xincd, yd, fit_start, fit_end, instrd, 5, sigd, zd, ad, taud, fittedd, chisquared, chisq_targetd);
		assertEquals("rld incorrect", 2, rld);
		assertEqualsScaled("a incorrect", 12220.635, ad[0], tolerance);
		assertEqualsScaled("tau incorrect", 2.44568, taud[0], tolerance);
		assertEqualsScaled("z incorrect", 160.39937, zd[0], tolerance);
		assertEqualsScaled("Chi square incorrect", 3322.832275390625, chisquared[0], tolerance);
	}

	/** Tests {@link FLIMLib#fitLMA}. */
	@Test
	public void testFitLMA() {
		final double param[] = { 5.263157800072804E-5, 0.0, 1000.0, 2.0 };
		final int lma = FLIMLib.LMA_fit(xincd, yd, fit_start, fit_end, instrd, 5, sigd, param, paramfree, fittedd, chisquared, chisq_targetd, chisquared[0]);
		assertEquals("lma incorrect", 14, lma);
		assertEqualsScaled("a incorrect", 12726.5693359375, param[2], tolerance);
		assertEqualsScaled("tau incorrect", 2.2782721519470215, param[3], tolerance);
		assertEqualsScaled("z incorrect", 208.85049438476562, param[1], tolerance);
		assertEqualsScaled("Chi square incorrect", 2636.51171875, chisquared[0], tolerance);
	
	}
	
	/** Tests {@link FLIMLib#GCI_marquardt_global_exps_instr}. */
	@Test
	public void testGCIGlobalWrapperCall() {		
		int result = DEFAULT_RET;
		Float2DMatrix fitted = new Float2DMatrix(new float[1][ndata]);
		Float2DMatrix residuals = new Float2DMatrix(new float[1][ndata]);
		result = FLIMLib.GCI_marquardt_global_exps_instr(xinc, trans, fit_start, fit_end, instr, 
				NoiseType.swigToEnum(5), sig, FitType.FIT_GLOBAL_MULTIEXP, param2d, paramfree, restrain, chisq_delta, 
				fitted, residuals, chisq_trans, chisq_global, df, 1);
		assertTrue(result != DEFAULT_RET);
	}
	
	/** Tests {@link FLIMLib#GCI_marquardt_global_generic_instr}. */
	@Test
	public void testGCIGlobalGenericCall() {
		int result = DEFAULT_RET;
		Float2DMatrix fitted = new Float2DMatrix(new float[1][ndata]);
		Float2DMatrix residuals = new Float2DMatrix(new float[1][ndata]);
		
		result = FLIMLib.GCI_marquardt_global_generic_instr(xinc, trans, fit_start, fit_end,
				instr, noise, sig, param2d, paramfree, gparam, restrain, chisq_delta, FitFunc.GCI_MULTIEXP_LAMBDA, fitted, 
				residuals, chisq_trans, chisq_global, df);
		assertTrue(result != DEFAULT_RET);
	}
	
	/** Tests {@link FLIMLib#GCI_Phasor}. */
	@Test
	public void testPhasorWrapperCall() {
		float u[] = {0}; 
		float v[] = {0};
		float taup[] = {0};
		float taum[] = {0};
		
		int ret = FLIMLib.GCI_Phasor(xinc, y, fit_start, fit_end, z, u, v, taup, taum, tau, fitted, residuals, chisquare);
		assertEquals("phasor failed", ret, 0);
		assertEqualsScaled("z incorrect", 0.0f, z[0], tolerance);
		assertEqualsScaled("u incorrect", 0.37805783f, u[0], tolerance);
		assertEqualsScaled("v incorrect", 0.43042996f, v[0], tolerance);
		assertEqualsScaled("taup incorrect", 2.2203490f, taup[0], tolerance);
		assertEqualsScaled("taum incorrect", 2.790166f, taum[0], tolerance);
		assertEqualsScaled("Chi square incorrect", 10217.323f, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", new float[] {293.38422f, 287.97586f, 282.66714f, 277.45636f, 272.34155f},
				Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1), tolerance);
		assertArrayEqualsScaled("residuals incorrect", new float[] {0.6157837f, 16.02414f, -18.667145f, 15.54364f, 21.658447f},
				Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1), tolerance);
		assertEqualsScaled("period incorrect", 12.25340843f, FLIMLib.GCI_Phasor_getPeriod(), tolerance);
	}
	
	/** Tests {@link FLIMLib#GCI_triple_integral_fitting_engine} with/without {@code instr}. */
	@Test
	public void testGCI_triple_integral_fitting_engine() {
		chisquare[0] = z[0] = a[0] = tau[0] = 0;
		int ret0 = FLIMLib.GCI_triple_integral_fitting_engine(xinc, y, fit_start, fit_end, null, noise, sig, z, a,
				tau, fitted, residuals, chisquare, chisq_target);
		assertEquals("return value incorrect", 2, ret0);
		assertEqualsScaled("z incorrect", 160.39937, z[0], tolerance);
		assertEqualsScaled("a incorrect", 13382.208, a[0], tolerance);
		assertEqualsScaled("tau incorrect", 2.4456801, tau[0], tolerance);
		assertEqualsScaled("Chi square incorrect", 3322.8374, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1),
				new float[] {359.9866f, 358.62238f, 357.2935f, 355.999f, 0}, tolerance);
		assertArrayEqualsScaled("residuals incorrect", Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1),
				new float[] {-185.9866f, -192.62238f, -194.29349f, -228.999f, 0}, tolerance);
		
		chisquare[0] = z[0] = a[0] = tau[0] = 0;
		int ret1 = FLIMLib.GCI_triple_integral_fitting_engine(xinc, y, fit_start, fit_end, instr, noise, sig, z, a,
				tau, fitted, residuals, chisquare, chisq_target);
		assertEquals("return value incorrect", 2, ret1);
		assertEqualsScaled("z incorrect", 160.39937, z[0], tolerance);
		assertEqualsScaled("a incorrect", 12220.636, a[0], tolerance);
		assertEqualsScaled("Chi square incorrect", 3322.8323, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", new float[] {359.98657f, 358.62238f, 357.29346f, 355.99896f, 0},
				Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1), tolerance);
		assertArrayEqualsScaled("residuals incorrect", new float[] {359.98657f, 358.62238f, 357.29346f, 355.99896f, 0},
				Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1), tolerance);
	}
	
	/** Tests {@link FLIMLib#GCI_marquardt_fitting_engine}. */
	@Test
	public void testGCI_marquardt_fitting_engine() {
		int ret0 = FLIMLib.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param0, paramfree, restrain, FitFunc.GCI_MULTIEXP_TAU, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent);

		// Exact return value depends on CPU architecture.
		assertTrue("iterations greater than 10", ret0 > 10);

		assertArrayEqualsScaled("param matrix incorrect", new float[] {207.12376f, 12719.01f, 2.2813265f}, param0, tolerance);
		assertEqualsScaled("covariance matrix incorrect", -0.088296, covar.asArray()[1][2], tolerance);
		assertEqualsScaled("alpha matrix incorrect", 13.714218, alpha.asArray()[1][2], tolerance);
		assertEqualsScaled("erraxes matrix incorrect", 1.1696135E-6, erraxes.asArray()[1][2], tolerance);
		assertEqualsScaled("Chi square incorrect", 2638.129, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", new float[] {265.93097f, 264.4389f, 262.98468f, 261.56738f, 260.18604f},
				Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1), tolerance);
		assertArrayEqualsScaled("residuals incorrect", new float[] {-91.93097f, -98.4389f, -99.98468f, -134.56738f, -138.18604f},
				Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1), tolerance);
		
		// z, a should be approximately the same; lambda * tau should be approximately 1
		int ret1 = FLIMLib.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param1, paramfree, restrain, FitFunc.GCI_MULTIEXP_LAMBDA, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent);

		// Exact return value depends on CPU architecture.
		assertTrue("iterations greater than 10", ret0 > 10);

		assertArrayEqualsScaled("param matrix incorrect", new float[] {207.12376f, 12719.01f, 1 / 2.2813265f}, param1, tolerance);
		assertEqualsScaled("covariance matrix incorrect", 0.0169655140, covar.asArray()[1][2], tolerance);
		assertEqualsScaled("alpha matrix incorrect", -71.375, alpha.asArray()[1][2], tolerance);
		assertEqualsScaled("erraxes matrix incorrect", -4.318109E-08, erraxes.asArray()[1][2], tolerance);
		assertEqualsScaled("Chi square incorrect", 2638.1287, chisquare[0], tolerance);
		assertArrayEqualsScaled("fitted incorrect", new float[] {265.93094f, 264.43887f, 262.98468f, 261.56738f, 260.18604f},
				Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1), tolerance);
		assertArrayEqualsScaled("residuals incorrect", new float[] {265.93094f, 264.43887f, 262.98468f, 261.56738f, 260.18604f},
				Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1), tolerance);
	}

	/** See if the results are the same given a trivial instr. */
	@Test
	public void testInstrConsistency() {
		float[] param0 = new float[] { 160.39937f, 13382.208f, 2.4456801f };
		float[] param1 = new float[] { 160.39937f, 13382.208f, 2.4456801f };
		FLIMLib.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, null,
				noise, sig, param0, paramfree, restrain, FitFunc.GCI_MULTIEXP_TAU, fitted,
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta,
				chisq_percent);
		FLIMLib.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, new float[] {1, 0, 0, 0, 0},
				noise, sig, param1, paramfree, restrain, FitFunc.GCI_MULTIEXP_TAU, fitted,
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta,
				chisq_percent);
		assertArrayEquals(param0, param1, tolerance);
	}
	
	/** Tests {@link FLIMLib#GCI_EcfModelSelectionEngine}. */
	@Test
	public void testGCI_EcfModelSelectionEngine() {
		// setup arguments
		final float param0_ref[] = param0.clone();
		final float param1_ref[] = param1.clone();
		final float chisq_diff[] = { 0 };
		final int model[] = { -1 };
		// create two models to compare
		DecayModel paramsandfits[] = { 
				new DecayModel(FitFunc.GCI_MULTIEXP_TAU,    param0, paramfree,
						restrain, fit_end + 1, chisq_target, chisq_delta, chisq_percent),
				new DecayModel(FitFunc.GCI_MULTIEXP_LAMBDA, param1, paramfree,
						restrain, fit_end + 1, chisq_target, chisq_delta, chisq_percent)
		};
		int ret = FLIMLib.GCI_EcfModelSelectionEngine(xinc, y, fit_start, fit_end, instr, noise, sig,
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
		// they should be the same, but let's not be so strict
		assertEquals("chisq_diff incorrect", 0, chisq_diff[0], tolerance * 100);
		assertEquals("model incorrect", 1, model[0]);
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
		float y = 0;
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
			
			y = lambda.fit(x, param, dy_dparam);
			assertEqualsScaled("y incorrect", expectedY, y, tolerance);
			assertArrayEqualsScaled("expectedDy_dparam incorrect",
					expectedDy_dparam, dy_dparam, tolerance);

			y = 0;
			expectedY = 0;
		}
	}
	
	/** Tests {@link FitFunc} callback */
	@Test
	public void testCustomFitFunc() {
		// This test guarantees java FitFuncs behave the same as C FitFuncs
		// outputs for reference
		final float param1_ref[] = param1.clone();
		final float fitted_ref[] = new float[ndata];
		final float residuals_ref[] = new float[ndata];
		final float chisquare_ref[] = {0}; 
		final Float2DMatrix covar_ref = new Float2DMatrix(new float[nparam][nparam]);
		final Float2DMatrix alpha_ref = new Float2DMatrix(new float[nparam][nparam]);
		final Float2DMatrix erraxes_ref = new Float2DMatrix(new float[nparam][nparam]);
		
		// copied from native code
		FitFunc customLambda = new FitFunc() {
			public float fit(final float x, final float param[], float dy_dparam[]) {
				float ex;
				float y = 0;
				for (int i = 1; i < nparam - 1; i += 2) {
					dy_dparam[i] = ex = (float) Math.exp(-param[i + 1] * x);
					ex *= param[i];
					y += ex;
					dy_dparam[i + 1] = -ex * x;
				}
				return y;
			}
		};
		
		int ret = FLIMLib.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param1, paramfree, restrain, customLambda, fitted,
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

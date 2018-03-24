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
	final float tolerance = 1e-5f;
	final int NEVER_CALLED = -100;
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
	final float param[] = { 0, 1000, 2};
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
		assertEquals("a incorrect", 11823.900390625, ad[0], tolerance);
		assertEquals("tau incorrect", 2.350015878677368, taud[0], tolerance);
		assertEquals("z incorrect", 105.92437744140625, zd[0], tolerance);
		assertEquals("Chi square incorrect", 78742.3359375, chisquared[0], tolerance);
	}

	/** Tests {@link SLIMCurve#fitLMA}. */
	@Test
	public void testFitLMA() {
		final double param[] = { 5.263157800072804E-5, 0.0, 1000.0, 2.0 };
		final int lma = SLIMCurve.LMA_fit(xincd, yd, fit_start, fit_end, instrd, 5, sigd, param, paramfree, fittedd, chisquared, chisq_targetd, chisquared[0]);
		assertEquals("lma incorrect", 58, lma);
		assertEquals("a incorrect", 10501.380859375, param[2], tolerance);
		assertEquals("tau incorrect", 2.8925282955169678, param[3], tolerance);
		assertEquals("z incorrect", -163.77395629882812, param[1], tolerance);
		assertEquals("Chi square incorrect", 83087.8671875, chisquared[0], tolerance);
	
	}
	
	/** Tests {@link SLIMCurve#GCI_marquardt_global_exps_instr}. */
	@Test
	public void testGCIGlobalWrapperCall() {		
		int result = NEVER_CALLED;
		Float2DMatrix fitted = new Float2DMatrix(new float[1][ndata]);
		Float2DMatrix residuals = new Float2DMatrix(new float[1][ndata]);
		result = SLIMCurve.GCI_marquardt_global_exps_instr(xinc, trans, fit_start, fit_end, instr, 
				NoiseType.swigToEnum(5), sig, FitType.FIT_GLOBAL_MULTIEXP, param2d, paramfree, restrain, chisq_delta, 
				fitted, residuals, chisq_trans, chisq_global, df, 1);
		//System.out.println("result: " + result);
		assertTrue(result != NEVER_CALLED);
	}
	
	/** Tests {@link SLIMCurve#GCI_marquardt_global_generic_instr}. */
	@Test
	public void testGCIGlobalGenericCall() {
		int result = NEVER_CALLED;
		Float2DMatrix fitted = new Float2DMatrix(new float[1][ndata]);
		Float2DMatrix residuals = new Float2DMatrix(new float[1][ndata]);
		
		result = SLIMCurve.GCI_marquardt_global_generic_instr(xinc, trans, fit_start, fit_end,
				instr, noise, sig, param2d, paramfree, gparam, restrain, chisq_delta, FitFunc.GCI_MULTIEXP_LAMBDA, fitted, 
				residuals, chisq_trans, chisq_global, df);
		//System.out.println("generic result: " + result);
		assertTrue(result != NEVER_CALLED);
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
		assertEquals("z incorrect", z[0], 0.0f, tolerance);
		assertEquals("u incorrect", u[0], 0.24798244f, tolerance);
		assertEquals("v incorrect", v[0], 0.49327368f, tolerance);
		assertEquals("taup incorrect", taup[0], 2.9834206f, tolerance);
		assertEquals("taum incorrect", taum[0], 2.2650633f, tolerance);
		assertEquals("Chi square incorrect", chisquare[0], 81201.26f, tolerance);
		assertArrayEquals(Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1),
				new float[] {293.38422f, 287.97586f, 282.66714f, 277.45636f, 272.34155f}, tolerance);
		assertArrayEquals(Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1),
				new float[] {0.6157837f, 16.02414f, -18.667145f, 15.54364f, 21.658447f}, tolerance);
		assertEquals("period incorrect", SLIMCurve.GCI_Phasor_getPeriod(), 9.423828125f, tolerance);
	}
	
	/** Tests {@link SLIMCurve#GCI_triple_integral_fitting_engine}. */
	@Test
	public void testGCI_triple_integral_fitting_engine() {
		chisquare[0] = z[0] = a[0] = tau[0] = 0;
		int ret0 = SLIMCurve.GCI_triple_integral_fitting_engine(xinc, y, fit_start, fit_end, instr, noise, sig, z, a,
				tau, fitted, residuals, chisquare, chisq_target);
		assertEquals("rld value incorrect", 2, ret0);
		assertEquals("z incorrect", z[0], 105.92438, tolerance);
		assertEquals("a incorrect", a[0], 11823.900390625, tolerance);
		assertEquals("Chi square incorrect", 78742.3359375, chisquare[0], tolerance);
		assertArrayEquals("fitted incorrect", Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1),
				new float[] {241.53043f, 236.57343f, 231.71306f, 226.94742f, 222.27469f}, tolerance);
		assertArrayEquals("residuals incorrect", Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1),
				new float[] {52.469574f, 67.426575f, 32.28694f, 66.05258f, 71.72531f}, tolerance);
		
		chisquare[0] = z[0] = a[0] = tau[0] = 0;
		int ret1 = SLIMCurve.GCI_triple_integral_fitting_engine(xinc, y, fit_start, fit_end, instr, NoiseType.NOISE_GAUSSIAN_FIT, sig, z, a,
				tau, fitted, residuals, chisquare, chisq_target);
		assertEquals("rld value incorrect", 2, ret1);
		assertEquals("z incorrect", z[0], 105.92438, tolerance);
		assertEquals("a incorrect", a[0], 11823.900390625, tolerance);
		assertEquals("Chi square incorrect", 78742.3359375, chisquare[0], tolerance);
		assertArrayEquals("fitted incorrect", Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1),
				new float[] {241.53043f, 236.57343f, 231.71306f, 226.94742f, 222.27469f}, tolerance);
		assertArrayEquals("residuals incorrect", Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1),
				new float[] {52.469574f, 67.426575f, 32.28694f, 66.05258f, 71.72531f}, tolerance);
	}
	
	/** Tests {@link SLIMCurve#GCI_marquardt_fitting_engine}. */
	@Test
	public void testGCI_marquardt_fitting_engine() {
		final float param0[] = { 0, 1000, 2}; // z, a, tau
		int ret0 = SLIMCurve.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param0, paramfree, restrain, FitFunc.GCI_MULTIEXP_TAU, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent);
		assertEquals("lma value incorrect", 57, ret0);
		assertArrayEquals("param matrix incorrect", param0, new float[] { -163.77396f, 10501.38f, 2.8925283f }, tolerance);
		assertEquals("covariance matrix incorrect", covar.asArray()[1][2], -0.18844599, tolerance);
		assertEquals("alpha matrix incorrect", alpha.asArray()[1][2], 19.542074, tolerance);
		assertEquals("erraxes matrix incorrect", erraxes.asArray()[1][2], 1.3784504E-6, tolerance / 1E6);
		assertEquals("Chi square incorrect", 83087.8984375, chisquare[0], tolerance);
		assertArrayEquals("fitted incorrect", Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1),
				new float[] {231.93842f, 225.31448f, 218.80151f, 212.39758f, 206.1008f}, tolerance);
		assertArrayEquals("residuals incorrect", Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1),
				new float[] {62.061584f, 78.68552f, 45.198486f, 80.60242f, 87.8992f}, tolerance);
		
		// z, a should be approximately the same; lambda * tau should be approximately 1
		final float[] param1 = { 0, 1000, 0.5f }; // z, a, lambda
		int ret1 = SLIMCurve.GCI_marquardt_fitting_engine(xinc, y, fit_start, fit_end, instr, 
				noise, sig, param1, paramfree, restrain, FitFunc.GCI_MULTIEXP_LAMBDA, fitted, 
				residuals, chisquare, covar, alpha, erraxes, chisq_target, chisq_delta, 
				chisq_percent);
		assertEquals("lma value incorrect", 40, ret1);
		assertArrayEquals("param matrix incorrect", param1, new float[] { -158.19781f, 10513.9375f, 0.34687048f }, tolerance);
		assertEquals("covariance matrix incorrect", covar.asArray()[1][2], 0.022667186, tolerance);
		assertEquals("alpha matrix incorrect", alpha.asArray()[1][2], -161.95181, tolerance);
		assertEquals("erraxes matrix incorrect", erraxes.asArray()[1][2], -1.9851784E-8, tolerance / 1E6);
		assertEquals("Chi square incorrect", 82995.6953125, chisquare[0], tolerance);
		assertArrayEquals("fitted incorrect", Arrays.copyOfRange(fitted, fitted.length - 6, fitted.length - 1),
				new float[] {233.6828f, 227.10132f, 220.63046f, 214.26831f, 208.01291f}, tolerance);
		assertArrayEquals("residuals incorrect", Arrays.copyOfRange(residuals, residuals.length - 6, residuals.length - 1),
				new float[] {60.3172f, 76.89868f, 43.369537f, 78.73169f, 85.98709f}, tolerance);
	}
	
	/** Tests {@link SLIMCurve#GCI_EcfModelSelectionEngine}. */
	//@Test
	public void testGCI_EcfModelSelectionEngine() {
		float chisq_diff[] = { 0, 0 };
		int model[] = { -1 };
		final float param0[] = { 0, 1000, 2};
		final float param1[] = { 0, 1000, 2};
		final int paramfree0[] = {1, 1, 1};
		final int paramfree1[] = {1, 1, 1};
		Float2DMatrix covar0 = new Float2DMatrix(new float[nparam][nparam]);
		Float2DMatrix covar1 = new Float2DMatrix(new float[nparam][nparam]);
		Float2DMatrix alpha0 = new Float2DMatrix(new float[nparam][nparam]);
		Float2DMatrix alpha1 = new Float2DMatrix(new float[nparam][nparam]);
		Float2DMatrix erraxes0 = new Float2DMatrix(new float[nparam][nparam]);
		Float2DMatrix erraxes1 = new Float2DMatrix(new float[nparam][nparam]);
		float[] fitted0 = new float[ndata];
		float[] fitted1 = new float[ndata];
		float[] residuals0 = new float[ndata];
		float[] residuals1 = new float[ndata];
		DecayModelSelParamValuesAndFit paramsandfits[] = { 
				new DecayModelSelParamValuesAndFit(FitFunc.GCI_MULTIEXP_LAMBDA, param0, paramfree0, restrain, fitted0, 
						residuals0, chisq_target, chisq_delta, chisq_percent, chisquare[0], covar0, alpha0, erraxes0),
				new DecayModelSelParamValuesAndFit(FitFunc.GCI_MULTIEXP_LAMBDA, param1, paramfree1, restrain, fitted1, 
						residuals1, chisq_target, chisq_delta, chisq_percent, chisquare[0], covar1, alpha1, erraxes1)
		};
		int ret = SLIMCurve.GCI_EcfModelSelectionEngine(xinc, y, fit_start, 12, instr, noise, sig, paramsandfits, chisq_diff, model);
		System.out.println(ret);
	}
	/** Tests {@link Float2DMatrix}.  */
	@Test
	public void testFloat2DMatrix() {
		Random rng = new Random(SEED);
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
				assertArrayEquals("Float2DMatrix not equal", arr[j], arrOut[j], 1e-10f);
		}
	}
	
	/** Tests {@link Float2DMatrix}.  */
	@Test
	public void testInt2DMatrix() {
		Random rng = new Random(SEED);
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

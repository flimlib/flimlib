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

import java.util.Random;

import org.junit.Test;
//import org.scijava.nativelib.NativeLoader;

import slim.*;


/**
 * Tests {@link SLIMCurve}.  Test parameters currently based 
 * off of values in files data.data and test.ini
 * 
 * @author Zach Petersen, Dasong Gao
 */
public class SLIMCurveTest {
//static {
//        
//        try {
//            System.loadLibrary("slim-curve");
//            System.loadLibrary("slim-curve-java");
//        	//NativeLoader.loadLibrary("META-INF/lib/windows_64/slim-curve-java");
//            //NativeLoader.loadLibrary("META-INF/lib/windows_64/slim-curve");
//            
//        } catch (Exception e) {
//            //System.err.println("Native library failed to load. Exiting.\n" + e);
//            System.err.println("Cannot extract native library. Exiting.\n" + e);
//            System.exit(1);
//        }
//    }
	public static final int TEST_SIZE = 100;
	final int NEVER_CALLED = -100;
	final double xInc = 0.048828125; 
	final double y[] = {
			40.000000, 45.000000, 34.000000, 54.000000, 44.000000,
			53.000000, 47.000000, 56.000000, 62.000000, 66.000000,
			82.000000, 90.000000, 108.000000, 122.000000, 323.000000,
			1155.000000, 4072.000000, 8278.000000, 11919.000000, 13152.000000,
			13071.000000, 12654.000000, 11946.000000, 11299.000000, 10859.000000,
			10618.000000, 10045.000000, 9576.000000, 9208.000000, 9113.000000,
			8631.000000, 8455.000000, 8143.000000, 8102.000000, 7672.000000,
			7384.000000, 7463.000000, 7254.000000, 6980.000000, 6910.000000,
			6411.000000, 6355.000000, 6083.000000, 5894.000000, 5880.000000,
			5735.000000, 5528.000000, 5343.000000, 5224.000000, 4933.000000,
			5026.000000, 4914.000000, 4845.000000, 4681.000000, 4426.000000,
			4485.000000, 4271.000000, 4295.000000, 4183.000000, 3989.000000,
			3904.000000, 3854.000000, 3801.000000, 3600.000000, 3595.000000,
			3434.000000, 3457.000000, 3291.000000, 3280.000000, 3178.000000,
			3132.000000, 2976.000000, 2973.000000, 2940.000000, 2770.000000,
			2969.000000, 2851.000000, 2702.000000, 2677.000000, 2460.000000,
			2536.000000, 2528.000000, 2347.000000, 2382.000000, 2380.000000,
			2234.000000, 2251.000000, 2208.000000, 2115.000000, 2136.000000,
			2000.000000, 2006.000000, 1970.000000, 1985.000000, 1886.000000,
			1898.000000, 1884.000000, 1744.000000, 1751.000000, 1797.000000,
			1702.000000, 1637.000000, 1547.000000, 1526.000000, 1570.000000,
			1602.000000, 1557.000000, 1521.000000, 1417.000000, 1391.000000,
			1332.000000, 1334.000000, 1290.000000, 1336.000000, 1297.000000,
			1176.000000, 1189.000000, 1220.000000, 1209.000000, 1217.000000,
			1140.000000, 1079.000000, 1059.000000, 1074.000000, 1061.000000,
			1013.000000, 1075.000000, 1021.000000, 1012.000000, 940.000000,
			982.000000, 866.000000, 881.000000, 901.000000, 883.000000,
			893.000000, 845.000000, 819.000000, 831.000000, 758.000000,
			794.000000, 779.000000, 772.000000, 779.000000, 791.000000,
			729.000000, 732.000000, 687.000000, 690.000000, 698.000000,
			661.000000, 647.000000, 668.000000, 642.000000, 619.000000,
			629.000000, 656.000000, 579.000000, 579.000000, 600.000000,
			563.000000, 584.000000, 531.000000, 554.000000, 526.000000,
			484.000000, 530.000000, 515.000000, 493.000000, 502.000000,
			479.000000, 445.000000, 439.000000, 466.000000, 431.000000,
			423.000000, 451.000000, 412.000000, 415.000000, 393.000000,
			404.000000, 390.000000, 398.000000, 352.000000, 394.000000,
			376.000000, 338.000000, 377.000000, 367.000000, 355.000000,
			352.000000, 375.000000, 339.000000, 347.000000, 316.000000,
			295.000000, 322.000000, 311.000000, 294.000000, 304.000000,
			264.000000, 293.000000, 294.000000, 283.000000
			}; 
	final int fitStart = 10;
	final int fitEnd = 203;
	final double instr[] = {
			0.00911042001098394393920898437500000,
			0.03882249817252159118652343750000000,
			0.13171100616455078125000000000000000,
			0.25238901376724243164062500000000000,
			0.27722999453544616699218750000000000,
			0.18023300170898437500000000000000000,
			0.07927619665861129760742187500000000,
			0.03122770041227340698242187500000000
	};  
	final int noise = 5;
	final double sig[] = {}; 
	final double z[]  = {0};
	final double a[]  = {1000};
	final double tau[] = {2};
	final double fitted[] = new double[y.length];  
	final double chiSquare[] = {0}; 
	final double chiSquareTarget = 237.5f;
	final int chi_sq_adjust = fitEnd - fitStart - 3;
	final int nInstr = instr.length; 
	final int paramFree[] = {1, 1, 1};
	final int nParam = 2;	
	final float chisq_trans[] = {}; 
	final int ndata = 203;
	final int ntrans = -1;
	final int ftype = 0;
	final float chisq_delta = 0;
	final int drop_bad_transients = 0;
	final int gparam[] = {};
	float[] trans = {2,3,4};
	float[] param = {0, 1000, 2}; //z, a, tau
	float[] residuals = {0};
	int[] df = {0};
	float[] chisq_global = {0};
	int restrain = 1;
	String fitFunc = "GCI_MULTIEXP_LAMBDA";
	

	/** Tests {@link SLIMCurve#fitRLD}. */
	@Test
	public void testFitRLD() {
		final int rld = SLIMCurve.RLD_fit(xInc, y, fitStart, fitEnd, instr, noise, sig, z, a, tau, fitted, chiSquare, chiSquareTarget);
		System.out.println("RLD estimate = " + rld);
		System.out.println("A: " + a[0] + " Tau: " + tau[0] + " Z: " + z[0] + " X^2: " + (chiSquare[0] / chi_sq_adjust));
		int _rld = 2;
		int _a = 11823;
		int _tau = 2;
		int _z = 105;
		int _x2 = 414;
		assertEquals("rld value is not correct", _rld, rld);
		assertEquals("a value is not correct", _a, (int)a[0]);
		assertEquals("tau value is not correct", _tau, (int)tau[0]);
		assertEquals("z value is not correct", _z, (int)z[0]);
		assertEquals("Chi squared value is not correct", _x2, (int)chiSquare[0] / chi_sq_adjust);
	}

	/** Tests {@link SLIMCurve#fitLMA}. */
	@Test
	public void testFitLMA() {
		final double param[] = {chiSquare[0] / chi_sq_adjust, 0, 1000, 2}; //chi squared, z, a, tau
		final int paramFree[] = {1, 1, 1};
		final double chiSquareDelta = .01;
		final int lma = SLIMCurve.LMA_fit(xInc, y, fitStart, fitEnd, instr, noise, sig, param, paramFree, fitted, chiSquare, chiSquareTarget, chiSquareDelta);
		System.out.println("LMA estimate = " + lma);
		System.out.println("A: " + param[2] + " Tau: " + param[3] + " Z: " + param[1] + " X^2: " + (chiSquare[0] / chi_sq_adjust));
		int _lma = 58;
		int _a = 10501;
		int _tau = 2;
		int _z = -163;
		int _x2 = 437;
		assertEquals("lma value is not correct", _lma, lma);
		assertEquals("a value is not correct", _a, (int)param[2]);
		assertEquals("tau value is not correct", _tau, (int)param[3]);
		assertEquals("z value is not correct", _z, (int)param[1]);
		assertEquals("Chi squared value is not correct", _x2, (int)chiSquare[0] / chi_sq_adjust);
	
	}
	/*@Test
	public void testGCIGlobalWrapperCall() {		
		int result = NEVER_CALLED;
		result = SLIMCurve.GCI_marquardt_global_exps_instr(xInc, trans, fitStart, 
				fitEnd, instr, instr.length, noise, sig, ftype, param, paramFree, restrain, 
				chisq_delta, fitted, residuals, chisq_trans, chisq_global, null, drop_bad_transients);
		System.out.println("result: " + result);
		assertTrue(result != NEVER_CALLED);
		
	}
	@Test
	public void testGCIGlobalGenericCall() {
		int result = NEVER_CALLED;
		result = SLIMCurve.GCI_marquardt_global_generic_instr(xInc, trans, fitStart, fitEnd,
				instr, nInstr, noise, sig, param, paramFree, nParam, gparam, restrain, chisq_delta, FitFunc.GCI_MULTIEXP_LAMBDA, fitted, 
				residuals, chisq_trans, chisq_global, null);
		System.out.println("generic result: " + result);
		assertTrue(result != NEVER_CALLED);
	}
		final int nInstr = 8;
	final float sig[] = {}; 
	final int paramFree[] = {1, 1, 1};
	final int nParam = 2;	
	final float chisq_trans[] = {}; 
	final int ndata = 203;
	final int ntrans = -1;
	final int ftype = 0;
	final float chisq_delta = 0;
	final int drop_bad_transients = 0;
	final int gparam[] = {};
	float[] trans = {2,3,4};
	float[] param = {0, 1000, 2}; //z, a, tau
	float[] residuals = {0};
	int[] df = {0};
	float[] chisq_global = {0};
	int noise = 5;
	float[] fitted = {0};
	int restrain = 1;
	String fitFunc = "GCI_MULTIEXP_LAMBDA";
	
	@Test
	public void testGCIGlobalWrapperCall() {		
		int result = NEVER_CALLED;
		result = SLIMCurve.GCI_marquardt_global_exps_instr(xInc, trans, ndata, ntrans, fitStart, 
				fitEnd, instr, nInstr, noise, sig, ftype, param, paramFree, nParam, restrain, 
				chisq_delta, fitted, residuals, chisq_trans, chisq_global, df, drop_bad_transients);
		System.out.println("result: " + result);
		assertTrue(result != NEVER_CALLED);
		
	}
	@Test
	public void testGCIGlobalGenericCall() {
		int result = NEVER_CALLED;
		result = SLIMCurve.GCI_marquardt_global_generic_instr(xInc, trans, ndata, ntrans, fitStart, fitEnd,
				instr, nInstr, noise, sig, param, paramFree, nParam, gparam, restrain, chisq_delta, fitFunc, fitted, 
				residuals, chisq_trans, chisq_global, df);
		System.out.println("generic result: " + result);
		assertTrue(result != NEVER_CALLED);
	}
	
	@Test
	public void testPhasorWrapperCall() {
		float z = 0;
		float u[] = {0}; 
		float v[] = {0};
		float taup[] = {0};
		float taum[] = {0};
		float tau[] = {0};
		float fitted[] = {0};
		float residuals[] = {0};
		float chiSquare[] = {0};
		
		float result = SLIMCurve.GCI_Phasor(xInc, y, fitStart, fitEnd, z, u, v, taup, taum, tau, fitted, residuals, chiSquare);
		System.out.println("phasor result: " + result + ",  phasor period: " + cLibrary.GCI_Phasor_getPeriod() + " chisq: " + chiSquare[0]);
		assertTrue(result >= 0 || result <= -5);
	}
	*/
	
	/** Tests {@link Float2DMatrix}.  */
	@Test
	public void testFloat2DMatrix() {
		Random rng = new Random();
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
		Random rng = new Random();
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

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

package loci.slim;

import static org.junit.Assert.*;

import org.junit.Test;

/**
 * Tests {@link SLIMCurve}.
 * 
 * @author Zach Petersen
 */
public class SLIMCurveTest {
	final double xInc = .048828;
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
			0.009110,
			0.038822,
			0.131711,
			0.252389,
			0.277230,
			0.180233,
			0.079276,
			0.031228
	};  
	final int nInstr = 8;
	final int noise = 5;
	final double sig[] = {}; 
	final double z[]  = {0};
	final double a[]  = {1000};
	final double tau[] = {2};
	final double fitted[] = {0};  
	final double chiSquare[] = {0}; 
	final double chiSquareTarget = 237.5;
	final int chi_sq_adjust = fitEnd - fitStart - 3;
	
	static {
		NarSystem.loadLibrary();
	}

	/** Tests {@link SLIMCurve#fitRLD}. */
	@Test
	public void testFitRLD() {
		SLIMCurve slimCurve = new SLIMCurve();
		final int rld = slimCurve.fitRLD(xInc, y, fitStart, fitEnd, instr, nInstr, noise, sig, z, a, tau, fitted, chiSquare, chiSquareTarget);
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
		final int nParam = 3;
		final double chiSquareDelta = .01;
		SLIMCurve slimCurve = new SLIMCurve();
		final int lma = slimCurve.fitLMA(xInc, y, fitStart, fitEnd, instr, nInstr, noise, sig, param, paramFree, nParam, fitted, chiSquare, chiSquareTarget, chiSquareDelta);
		System.out.println("LMA estimate = " + lma);
		System.out.println("A: " + param[1] + " Tau: " + param[2] + " Z: " + param[0] + " X^2: " + (chiSquare[0] / chi_sq_adjust));
		int _lma = 17;
		int _a = 11795;
		int _tau = 2;
		int _z = 105;
		int _x2 = 414;
//		assertEquals("rld value is not correct", _lma, lma);
//		assertEquals("a value is not correct", _a, (int)a[0]);
//		assertEquals("tau value is not correct", _tau, (int)tau[0]);
//		assertEquals("z value is not correct", _z, (int)z[0]);
//		assertEquals("Chi squared value is not correct", _x2, (int)chiSquare[0] / chi_sq_adjust);
	
	}

}

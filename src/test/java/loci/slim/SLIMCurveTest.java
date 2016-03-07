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

import org.junit.Test;

/**
 * Tests {@link SLIMCurve}.
 * 
 * @author Zach Petersen
 */
public class SLIMCurveTest {

	static {
		NarSystem.loadLibrary();
	}

	/** Tests {@link SLIMCurve#fitRLD}. */
	@Test
	public void testFitRLD() {
//		final double xInc, final double y[], final int fitStart,
//		final int fitEnd, final double instr[], final int nInstr, final int noise,
//		final double sig[], final double z[], final double a[], final double tau[],
//		final double fitted[], final double chiSquare[],
//		final double chiSquareTarget)
	}

	/** Tests {@link SLIMCurve#fitLMA}. */
	@Test
	public void testFitLMA() {
//		final double xInc, final double y[], final int fitStart,
//		final int fitEnd, final double instr[], final int n_instr, final int noise,
//		final double sig[], final double param[], final int paramFree[],
//		final int nParam, final double fitted[], final double chiSquare[],
//		final double chiSquareTarget, final double chiSquareDelta)
	}

}

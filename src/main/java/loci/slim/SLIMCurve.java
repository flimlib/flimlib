/*
 * #%L
 * SLIM Curve package for exponential curve fitting of spectral lifetime data.
 * %%
 * Copyright (C) 2010 - 2015 University of Oxford and Board of Regents of the
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

package loci.slim;

import org.scijava.nativelib.NativeLibraryUtil;

/**
 * TODO
 * 
 * @author Aivar Grislis
 */
public class SLIMCurve {

	private static volatile boolean s_libraryLoaded = false;

	/**
	 * This supports calling the library using JNI.
	 * 
	 * @param xInc
	 * @param y
	 * @param fitStart
	 * @param fitEnd
	 * @param instr
	 * @param nInstr
	 * @param sig
	 * @param z
	 * @param a
	 * @param tau
	 * @param fitted
	 * @param chiSquare
	 * @param chiSquareTarget
	 * @return
	 */
	private native int RLD_fit(double xInc, double y[], int fitStart, int fitEnd,
		double instr[], int nInstr, int noise, double sig[], double z[],
		double a[], double tau[], double fitted[], double chiSquare[],
		double chiSquareTarget);

	/**
	 * This supports calling the library using JNI.
	 * 
	 * @param xInc
	 * @param y
	 * @param fitStart
	 * @param fitEnd
	 * @param instr
	 * @param n_instr
	 * @param sig
	 * @param param
	 * @param paramFree
	 * @param nParam
	 * @param fitted
	 * @param chiSquare
	 * @param chiSquareTarget
	 * @return
	 */
	private native int LMA_fit(double xInc, double y[], int fitStart, int fitEnd,
		double instr[], int n_instr, int noise, double sig[], double param[],
		int paramFree[], int nParam, double fitted[], double chiSquare[],
		double chiSquareTarget, double chiSquareDelta);

	public SLIMCurve() {
		// load platform-appropriate native library
		s_libraryLoaded =
			NativeLibraryUtil.loadVersionedNativeLibrary(this.getClass(),
				"slim-curve");
		System.out.println("SLIM Curve loaded " + s_libraryLoaded);
	}

	/** Does the RLD fit with native library via JNI. */
	public int fitRLD(final double xInc, final double y[], final int fitStart,
		final int fitEnd, final double instr[], final int nInstr, final int noise,
		final double sig[], final double z[], final double a[], final double tau[],
		final double fitted[], final double chiSquare[],
		final double chiSquareTarget)
	{
		final int returnValue =
			RLD_fit(xInc, y, fitStart, fitEnd, instr, nInstr, noise, sig, z, a, tau,
				fitted, chiSquare, chiSquareTarget);
		return returnValue;
	}

	/** Does the LMA fit with native library via JNI. */
	public int fitLMA(final double xInc, final double y[], final int fitStart,
		final int fitEnd, final double instr[], final int n_instr, final int noise,
		final double sig[], final double param[], final int paramFree[],
		final int nParam, final double fitted[], final double chiSquare[],
		final double chiSquareTarget, final double chiSquareDelta)
	{
		final int returnValue =
			LMA_fit(xInc, y, fitStart, fitEnd, instr, n_instr, noise, sig, param,
				paramFree, nParam, fitted, chiSquare, chiSquareTarget, chiSquareDelta);
		return returnValue;
	}

}

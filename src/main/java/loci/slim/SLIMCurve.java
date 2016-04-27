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

/**
 * TODO
 * 
 * @author Aivar Grislis
 */
public class SLIMCurve {

	static {
		NarSystem.loadLibrary();
	}
	
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
	
	public static int GCI_marquardt_global_exps_instr(float xInc, float[] trans, 
			int ndata, int ntrans, int fitStart, int fitEnd, float[] instr, int nInstr, 
			int noise, float[] sig, int ftype, float[] param, int[] paramFree, 
			int nParam, int restrain, float chisq_delta, float[] fitted, 
			float[] residuals, float[] chisq_trans, float[] chisq_global, 
			int[] df, int drop_bad_transients) {
		//handle a couple of special cases
		noise_type noiseType = getNoiseType(noise); 
		restrain_type restrainType = getRestrainType(restrain); 
		
		return cLibrary.GCI_marquardt_global_exps_instr(xInc, trans, ndata, ntrans, 
				fitStart, fitEnd, instr, nInstr, noiseType, sig, ftype, param, paramFree, nParam, 
				restrainType, chisq_delta, fitted, residuals, chisq_trans, chisq_global, df, drop_bad_transients);
		  }


	public static int GCI_marquardt_global_generic_instr(float xInc, float[] trans, 
			int ndata, int ntrans, int fitStart, int fitEnd, float[] instr, int nInstr, 
			int noise, float[] sig, float[] param, int[] paramFree, int nParam, 
			int[] gparam, int restrain, float chisq_delta, 
			String fitfunc, float[] fitted, 
			float[] residuals, float[] chisq_trans, float[] chisq_global, int[] df) {
		//handle a couple of special cases
		noise_type noiseType = getNoiseType(noise);
		restrain_type restrainType = getRestrainType(restrain);
		SWIGTYPE_p_f_float_a___float_p_float_a___float_int__void fitFunc = getFitFunc(fitfunc);
		
		return cLibrary.GCI_marquardt_global_generic_instr(xInc, trans, ndata, ntrans, fitStart, fitEnd,
				instr, nInstr, noiseType, sig, param, paramFree, nParam, gparam, restrainType, chisq_delta, fitFunc, fitted, 
				residuals, chisq_trans, chisq_global, df);
	}
	


	public static float GCI_Phasor(float xInc, float y[], int fitStart, int fitEnd, 
			float z, float[] u, float[] v, float[] taup, float[] taum, float[] tau, 
			float[] fitted, float[] residuals, float[] chiSquare) {
		
		return cLibrary.GCI_Phasor(xInc, y, fitStart, fitEnd, z, u, v, taup, taum, tau,
				fitted, residuals, chiSquare);
	}

	
	//Helper functions
	//Currently only suitable functions in Ecf.h exposed. May want to expose more.
	private static SWIGTYPE_p_f_float_a___float_p_float_a___float_int__void getFitFunc(String fitfunc) {
		if (fitfunc.equals("GCI_MULTIEXP_LAMBDA"))
			return cLibrary.GCI_MULTIEXP_LAMBDA;
		else if (fitfunc.equals("GCI_MULTIEXP_TAU"))
			return cLibrary.GCI_MULTIEXP_TAU;
		else if (fitfunc.equals("GCI_STRETCHEDEXP"))
			return cLibrary.GCI_STRETCHEDEXP;
		//default
		return cLibrary.GCI_MULTIEXP_LAMBDA;
	}
	
	private static restrain_type getRestrainType(int restrain) {
		switch (restrain) {
			case 1: return restrain_type.ECF_RESTRAIN_DEFAULT ;
			case 2: return restrain_type.ECF_RESTRAIN_USER;
			default: return restrain_type.ECF_RESTRAIN_DEFAULT;
		}
	}

	private static noise_type getNoiseType(int noise) {
		switch (noise) {
			case 1: return noise_type.NOISE_CONST;
			case 2: return noise_type.NOISE_GAUSSIAN_FIT;
			case 3: return noise_type.NOISE_GIVEN;
			case 4: return noise_type.NOISE_MLE;
			case 5: return noise_type.NOISE_POISSON_DATA;
			case 6: return noise_type.NOISE_POISSON_FIT;
			default: return noise_type.NOISE_POISSON_DATA;
		}
	}
}

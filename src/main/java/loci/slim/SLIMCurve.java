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
		//TODO: return value
		final int returnValue =
			LMA_fit(xInc, y, fitStart, fitEnd, instr, n_instr, noise, sig, param,
				paramFree, nParam, fitted, chiSquare, chiSquareTarget, chiSquareDelta);
		return returnValue;
	}
	
	//TODO: not assigning double pointers to any actual data. 
	public static int GCI_marquardt_global_exps_instr(float xInc, float tran, 
			int ndata, int ntrans, int fitStart, int fitEnd, float[] instr, int nInstr, 
			int noise, float[] sig, int ftype, float params, int[] paramFree, 
			int nParam, int restrain, float chisq_delta, float fit, 
			float res, float[] chisq_trans, float chisqGlobal, 
			int degFree, int drop_bad_transients) {
		SWIGTYPE_p_float z = cLibrary.new_floatp();
		cLibrary.floatp_assign(z, 0f);
		

		SWIGTYPE_p_p_float fitted = cLibrary.GCI_ecf_matrix(0, 0);
		SWIGTYPE_p_p_float residuals = cLibrary.GCI_ecf_matrix(0, 0);
		SWIGTYPE_p_p_float trans = cLibrary.GCI_ecf_matrix(0, 0);
		SWIGTYPE_p_p_float param = cLibrary.GCI_ecf_matrix(0, 0);
		SWIGTYPE_p_float chisq_global = cLibrary.new_floatp();
		SWIGTYPE_p_int df = cLibrary.new_intp();
		noise_type noiseType = getNoiseType(noise); 
		restrain_type restrainType = getRestrainType(restrain); 
		
		//TODO: return value
		return cLibrary.GCI_marquardt_global_exps_instr(xInc, trans, ndata, ntrans, 
				fitStart, fitEnd, instr, nInstr, noiseType, sig, ftype, param, paramFree, nParam, 
				restrainType, chisq_delta, fitted, residuals, chisq_trans, chisq_global, df, drop_bad_transients);
		  }



	//TODO: not assigning double pointers to any actual data.  Also need to fix fitFunc
	public static int GCI_marquardt_global_generic_instr(float xInc, float tran, 
			int ndata, int ntrans, int fitStart, int fitEnd, float[] instr, int nInstr, 
			int noise, float[] sig, float params, int[] paramFree, int nParam, 
			int[] gparam, int restrain, float chisq_delta, 
			int fitfunc, float fit, 
			float res, float[] chisq_trans, float chisqGlobal, int degFree) {
		SWIGTYPE_p_p_float fitted = cLibrary.GCI_ecf_matrix(0, 0);
		SWIGTYPE_p_p_float residuals = cLibrary.GCI_ecf_matrix(0, 0);
		SWIGTYPE_p_p_float trans = cLibrary.GCI_ecf_matrix(0, 0);
		SWIGTYPE_p_p_float param = cLibrary.GCI_ecf_matrix(0, 0);
		SWIGTYPE_p_float chisq_global = cLibrary.new_floatp();
		SWIGTYPE_p_int df = cLibrary.new_intp();
		noise_type noiseType = getNoiseType(noise);
		restrain_type restrainType = getRestrainType(restrain);
		
		return cLibrary.GCI_marquardt_global_generic_instr(xInc, trans, ndata, ntrans, fitStart, fitEnd,
				instr, nInstr, noiseType, sig, param, paramFree, nParam, gparam, restrainType, chisq_delta, null, fitted, 
				residuals, chisq_trans, chisq_global, df);
		  }
	
	public static float GCI_Phasor(float xInc, float y[], int fitStart, int fitEnd, 
			float _z, float _u, float _v, float _taup, float _taum, float _tau, 
			float _fitted, float _residuals, float _chiSquare) {
		SWIGTYPE_p_float z = cLibrary.new_floatp();
		SWIGTYPE_p_float u = cLibrary.new_floatp();
		SWIGTYPE_p_float v = cLibrary.new_floatp();
		SWIGTYPE_p_float taup = cLibrary.new_floatp();
		SWIGTYPE_p_float taum = cLibrary.new_floatp();
		SWIGTYPE_p_float tau = cLibrary.new_floatp();
		SWIGTYPE_p_float fitted = cLibrary.new_floatp();
		SWIGTYPE_p_float residuals = cLibrary.new_floatp();
		SWIGTYPE_p_float chiSquare = cLibrary.new_floatp();
		cLibrary.floatp_assign(z, _z);
		cLibrary.floatp_assign(u, _u);
		cLibrary.floatp_assign(v, _v);
		cLibrary.floatp_assign(taup, _taup);
		cLibrary.floatp_assign(taum, _taum);
		cLibrary.floatp_assign(fitted, _fitted);
		cLibrary.floatp_assign(residuals, _residuals);
		cLibrary.floatp_assign(chiSquare, _chiSquare);
		cLibrary.floatp_assign(z, 0f);
		
		//TODO: GCI_Phasor uses out-parameters.  However, java is pass by value for primitives.
		//So return chiSquare or do pass by reference for chiSquare
		int ret = cLibrary.GCI_Phasor(xInc, y, fitStart, fitEnd, z, u, v, taup, taum, tau, fitted, residuals, chiSquare);
		if (ret < 0 ) 
			return ret;
		return cLibrary.floatp_value(chiSquare);

	}

	
	//Helper functions
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

/*
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2022 University of Oxford and Board of Regents of the
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

/* 
Main.cpp for the stand alone executable version of FLIMLib via the C++ interface.

*/
#include <iostream>
using namespace std;

#include "flimlib.hpp"

#include "math.h"

// Test Data
float transient[] = {43, 39, 50, 46, 56, 63, 62, 74, 60, 72, 58, 47, 41, 69, 69, 58, 55, 37, 55, 50, 52, 59, 51, 52, 51, 50, 53, 40, 45, 34, 54, 44, 53, 47, 56, 62, 66, 82, 90, 108, 122, 323, 1155, 4072, 8278, 11919, 13152, 13071, 12654, 11946, 11299, 10859, 10618, 10045, 9576, 9208, 9113, 8631, 8455, 8143, 8102, 7672, 7384, 7463, 7254, 6980, 6910, 6411, 6355, 6083, 5894, 5880, 5735, 5528, 5343, 5224, 4933, 5026, 4914, 4845, 4681, 4426, 4485, 4271, 4295, 4183, 3989, 3904, 3854, 3801, 3600, 3595, 3434, 3457, 3291, 3280, 3178, 3132, 2976, 2973, 2940, 2770, 2969, 2851, 2702, 2677, 2460, 2536, 2528, 2347, 2382, 2380, 2234, 2251, 2208, 2115, 2136, 2000, 2006, 1970, 1985, 1886, 1898, 1884, 1744, 1751, 1797, 1702, 1637, 1547, 1526, 1570, 1602, 1557, 1521, 1417, 1391, 1332, 1334, 1290, 1336, 1297, 1176, 1189, 1220, 1209, 1217, 1140, 1079, 1059, 1074, 1061, 1013, 1075, 1021, 1012, 940, 982, 866, 881, 901, 883, 893, 845, 819, 831, 758, 794, 779, 772, 779, 791, 729, 732, 687, 690, 698, 661, 647, 668, 642, 619, 629, 656, 579, 579, 600, 563, 584, 531, 554, 526, 484, 530, 515, 493, 502, 479, 445, 439, 466, 431, 423, 451, 412, 415, 393, 404, 390, 398, 352, 394, 376, 338, 377, 367, 355, 352, 375, 339, 347, 316, 295, 322, 311, 294, 304, 264, 293, 294, 283, 278, 302, 253, 259, 252, 278, 254, 245, 246, 242, 226, 241, 222, 198, 197, 245, 221, 228, 224, 216, 174, 166, 163, 127, 122 };
const int ndata = 256;
int transient_rise = 38;  // near the first rise
int transient_peak = 46;  // near the peak
int transient_end = 255;
float instr[] = { 0.00911042f, 0.0388225f, 0.131711f, 0.252389f, 0.27723f, 0.180233f, 0.0792762f, 0.0312277f };
const int ninstr = 8;
float xincr = 0.058628749f;

/* Entry point
*/
int main(int argc, const char * argv[])
{
	// Return values, if required
	float fitted[ndata];
	float residuals[ndata];
	int iterations = 0;

	// Create a fitter
	FLIMLib flimLib;

	// Load up the fitter inputs
	flimLib.setupData(transient, ndata, xincr, transient_rise, transient_peak, transient_end);

	// Run the RLD fit function: Basic fit, no prompt, on data subset
	iterations = flimLib.fitRLD();
	cout << "RLD Z:" << flimLib.param[FLIMLIB_RLD_PARAM_Z]
		<< " A:"    << flimLib.param[FLIMLIB_RLD_PARAM_A] 
		<< " tau:"  << flimLib.param[FLIMLIB_RLD_PARAM_TAU]
		<< " chisq:" << flimLib.getReducedChiSq() 
		<< " iterations:" << iterations << "\n";
	if (iterations <= 0) {
		cout << "ERROR with FLIMLib.fitRLD\n";
		exit(0);
	}

	// Run the LMA fit function; all flimLib.param values will carry over into the next fit as a starting point
	iterations = flimLib.fitLMA(FLIMLIB_MONO);
	cout << "LMA Z:" << flimLib.param[FLIMLIB_MONO_PARAM_Z]
		<< " A:" << flimLib.param[FLIMLIB_MONO_PARAM_A]
		<< " tau:" << flimLib.param[FLIMLIB_MONO_PARAM_TAU]
		<< " chisq:" << flimLib.getReducedChiSq() 
		<< " iterations:" << iterations << "\n";
	if (iterations <= 0) {
		cout << "ERROR with FLIMLib.fitLMA\n";
		exit(0);
	}
	

	// Use an instrument function
	flimLib.setupIRF(instr, ninstr);

	// Run the RLD fit function: simple
	iterations = flimLib.fitRLD();
	cout << "RLDi Z:" << flimLib.param[FLIMLIB_RLD_PARAM_Z]
		<< " A:" << flimLib.param[FLIMLIB_RLD_PARAM_A]
		<< " tau:" << flimLib.param[FLIMLIB_RLD_PARAM_TAU]
		<< " chisq:" << flimLib.getReducedChiSq()
		<< " iterations:" << iterations << "\n";
	if (iterations <= 0) {
		cout << "ERROR with FLIMLib.fitRLD\n";
		exit(0);
	}

	// Assign the optional fitter Outputs to make plots if we want
	flimLib.fitted = fitted;
	flimLib.residuals = residuals;

	// Run the LMA fit function: with instr and maximum likelihood (NOISE_MLE)
	flimLib.noise_model = NOISE_MLE;
	iterations = flimLib.fitLMA(FLIMLIB_MONO);
	cout << "LMAi Z:" << flimLib.param[FLIMLIB_MONO_PARAM_Z]
		<< " A:" << flimLib.param[FLIMLIB_MONO_PARAM_A]
		<< " tau:" << flimLib.param[FLIMLIB_MONO_PARAM_TAU]
		<< " chisq:" << flimLib.getReducedChiSq()
		<< " iterations:" << iterations << "\n";
	if (iterations <= 0) {
		cout << "ERROR with FLIMLib.fitLMA\n";
		exit(0);
	}

	// Setup a resonable start for tri-exp
	flimLib.param[FLIMLIB_TRI_PARAM_A1] /= 3.0;      // divide the amplitude A1
	flimLib.param[FLIMLIB_TRI_PARAM_A2] = flimLib.param[FLIMLIB_TRI_PARAM_A1]; // let A2 = A1
	flimLib.param[FLIMLIB_TRI_PARAM_TAU2] = flimLib.param[FLIMLIB_TRI_PARAM_TAU1] / 2.0f; // let tau2 = tau1 / 2
	flimLib.param[FLIMLIB_TRI_PARAM_A3] = flimLib.param[FLIMLIB_TRI_PARAM_A1]; // let A3 = A1
	flimLib.param[FLIMLIB_TRI_PARAM_TAU3] = flimLib.param[FLIMLIB_TRI_PARAM_TAU1] / 3.0f; // let tau3 = tau1 / 2
	// all other params carry over from last fit

	// Run the LMA Bi exp
	// all starting params from last fit
	iterations = flimLib.fitLMA(FLIMLIB_BI);
	cout << "LMAi Z:" << flimLib.param[FLIMLIB_BI_PARAM_Z]
		<< " A1:" << flimLib.param[FLIMLIB_BI_PARAM_A1]
		<< " tau1:" << flimLib.param[FLIMLIB_BI_PARAM_TAU1]
		<< " A2:" << flimLib.param[FLIMLIB_BI_PARAM_A2]
		<< " tau2:" << flimLib.param[FLIMLIB_BI_PARAM_TAU2]
		<< " chisq:" << flimLib.getReducedChiSq()
		<< " iterations:" << iterations << "\n";
	if (iterations <= 0) {
		cout << "ERROR with FLIMLib.fitLMA\n";
		exit(0);
	}

	// Let's restrain Z>0
	flimLib.restrainParameter(FLIMLIB_STRETCHED_PARAM_Z, 0.0, INFINITY);

	// and fix tau3=0.5
	flimLib.fixParameter(FLIMLIB_TRI_PARAM_TAU3, 0.5);

	// Run the LMA Tri exp, fix tau3=0.5;
	iterations = flimLib.fitLMA(FLIMLIB_TRI);
	cout << "LMAi Z:" << flimLib.param[FLIMLIB_TRI_PARAM_Z]
		<< " A1:" << flimLib.param[FLIMLIB_TRI_PARAM_A1]
		<< " tau1:" << flimLib.param[FLIMLIB_TRI_PARAM_TAU1]
		<< " A2:" << flimLib.param[FLIMLIB_TRI_PARAM_A2]
		<< " tau2:" << flimLib.param[FLIMLIB_TRI_PARAM_TAU2]
		<< " A3:" << flimLib.param[FLIMLIB_TRI_PARAM_A3]
		<< " tau3:" << flimLib.param[FLIMLIB_TRI_PARAM_TAU3]
		<< " chisq:" << flimLib.getReducedChiSq()
		<< " iterations:" << iterations << "\n";
	if (iterations <= 0) {
		cout << "ERROR with FLIMLib.fitLMA\n";
		exit(0);
	}

	// Setup a resonable start for stretched-exp
	flimLib.setupIRF(NULL, 0);
	flimLib.param[FLIMLIB_STRETCHED_PARAM_Z] = 1.0;
	flimLib.param[FLIMLIB_STRETCHED_PARAM_A] = 10000.0;
	flimLib.param[FLIMLIB_STRETCHED_PARAM_TAU] = 2.0;
	flimLib.param[FLIMLIB_STRETCHED_PARAM_H] = 1.0;
	// all other params carry over from last fit

	// Clear previous fixing and restraining
	flimLib.clearFixed();
	flimLib.clearRestrained();

	// Run the LMA stretched exp
	iterations = flimLib.fitLMA(FLIMLIB_STRETCHED);
	cout << "LMAi Z:" << flimLib.param[FLIMLIB_STRETCHED_PARAM_Z]
		<< " A:" << flimLib.param[FLIMLIB_STRETCHED_PARAM_A]
		<< " tau:" << flimLib.param[FLIMLIB_STRETCHED_PARAM_TAU]
		<< " h:" << flimLib.param[FLIMLIB_STRETCHED_PARAM_H]
		<< " chisq:" << flimLib.getReducedChiSq()
		<< " iterations:" << iterations << "\n";
	if (iterations <= 0) {
		cout << "ERROR with FLIMLib.fitLMA\n";
		exit(0);
	}

	
	// Run the phasor analysis
	int error = flimLib.fitPhasor(50.0f);
	cout << "Phasor Z:" << flimLib.param[FLIMLIB_PHASOR_PARAM_Z]
		<< " tau:" << flimLib.param[FLIMLIB_PHASOR_PARAM_TAU]
		<< " u:" << flimLib.param[FLIMLIB_PHASOR_PARAM_U]
		<< " v:" << flimLib.param[FLIMLIB_PHASOR_PARAM_V]
		<< " chisq:" << flimLib.getReducedChiSq() << "\n";
	if (error != FLIMLIB_SUCCESS) {
		cout << "ERROR with FLIMLib.fitPhasor\n";
		exit(0);
	}

	return 0;
}

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

/* 
Main.c for the stand alone executable version of SLIM Curve.
Reads ini settings file that is compatible with TRI2.
Provide transient.dat and optional prompt.dat.

The settings file requires a text based prompt file, not like TRI2.
So the path/name to this must be provided as an entry in the [prompt] section, e.g.
[prompt]
file = "prompt.dat"
Also the time difference between the points must be provided in the [main] section, e.g.
[main]
x_inc = "0.048828125"

KNOWN ISSUES
Parameter fixing and restraining is not yet supported here.
This code does not read the parameter values (tau etc.) from the ini file, it always does an RLD fit for starting values. This is different to TRI2 with the phasor plot as the pre-pulse is used in that case.
With stretched exp TRI2 converts the raw tau returned by the fitting engine into an average <tau> using the Gamma function.

*/

#include "Ecf.h"
#include "GCI_Phasor.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iniparser.h"

#define FALSE 0
#define TRUE !FALSE

#define SUCCESS 0
#define BAD_USAGE -1
#define BAD_FILE -2
#define BAD_SYNTAX -3

#define MAX_DATFILE_LENGTH 4096

/* Predeclarations
*/
int fit(int fit_type, noise_type noise_model,
    float chi_sq_target, float chi_sq_delta,
    int transient_size, float *transient_values,
    int sigma_size, float *sigma_values,
    int prompt_size, float *prompt_values,
    float x_inc, int fit_start, int fit_end
    );
int parse_and_fit(dictionary *ini, const char *dat_filename);
noise_type getNoiseModel (int val);


/* A roundf() function as it does not exist in all math lib implementations, esp. Visual Studio
*/
static float my_roundf(float val)
{    
    return (float)floor((double)val + 0.5);
}

/* Entry point
*/
int main(int argc, const char * argv[])
{
    int return_value = SUCCESS;
    if (argc != 3) {
        printf("Usage: %s settings_filename transient_filename\n", argv[0]);
        return_value = BAD_USAGE;
    }
    else {
        //FILE *file = fopen(argv[1], "r");
		dictionary *ini = iniparser_load(argv[1]);

        if (NULL == ini) {
            printf("Could not open file '%s'\n", argv[1]);
            return_value = BAD_FILE;
        }
        else {
        	FILE *dat_file = fopen(argv[2], "r");
        	
        	if (NULL == dat_file) {
		        printf("Could not open file '%s'\n", argv[2]);
		        return_value = BAD_FILE;
			}
			else {
				fclose(dat_file); // close it now we have checked we can open it
		        printf("Settings: %s\nTransient Data: %s\n", argv[1], argv[2]);
		        return_value = parse_and_fit(ini, argv[2]);
				iniparser_freedict(ini);
			}
        }
    }
    return return_value;
}

/* Performs the actual fit
*/
int fit(
    int fit_type, noise_type noise_model,
    float chi_sq_target, float chi_sq_delta,
    int transient_size, float *transient_values,
    int sigma_size, float *sigma_values,
    int prompt_size, float *prompt_values,
    float x_inc, int fit_start, int fit_end
    )
{
	int return_value;
	float a, tau, z;
	float *fitted = (float *)malloc((unsigned)transient_size * sizeof(float));
	float *residuals = (float *)malloc((unsigned)transient_size * sizeof(float));
	float chi_square;
	int chi_sq_adjust;
    int n_param;
    float *params;
    int *param_free;
    int n_param_free;
    int i;
	void (*fitfunc)(float, float [], float *, float[], int) = NULL;
	restrain_type restrain = ECF_RESTRAIN_DEFAULT;
	int chi_sq_percent = 95;
    int n_data = fit_end;
	float **covar    = GCI_ecf_matrix(n_data, n_data);
	float **alpha    = GCI_ecf_matrix(n_data, n_data);
	float **err_axes = GCI_ecf_matrix(n_data, n_data);
    
    // blind initial estimates as in TRI2/SP
    a = 1000.0f;
    tau = 2.0f;
    z = 0.0f;
    
	// want non-reduced chi square target for RLD
    chi_sq_adjust = fit_end - fit_start - 3;
    
    // for RLD+LMA fits TRI2/SP adjusts as follows for initial RLD fit:
    //  fit_start becomes index of peak of transient
    //  estimated A becomes value at peak
    //  noise becomes NOISE_POISSON_FIT
    
	return_value = GCI_triple_integral_fitting_engine(
			x_inc,
			transient_values,
			fit_start,
			fit_end,
			prompt_values,
			prompt_size,
			NOISE_POISSON_FIT,
			sigma_values,
			&z,
			&a,
			&tau,
			fitted,
			residuals,
			&chi_square,
			chi_sq_target * chi_sq_adjust
			);
	
	printf("RLD return_value %d\n", return_value);
	printf("RLD estimate A %f T %f Z %f X2 %f\n", a, tau, z, chi_square / chi_sq_adjust);

	if (fit_type==0)
		goto END;
	
    // adjust single exponential estimates for multiple exponential fits
    fitfunc = GCI_multiexp_tau;
    switch (fit_type) {   // This choice of fit types is not a good thing to have! TRI2 should save something more obvious
        // single exponential
        case 1:
        case 2:
            // params are Z, A, T
            n_param = 3;
            params = (float *)malloc((size_t) n_param * sizeof(float));
            params[0] = z;
            params[1] = a;
            params[2] = tau;
            break;
		case 3:            
            // params are Z, A, T
            n_param = 3;
            params = (float *)malloc((size_t) n_param * sizeof(float));
            params[0] = z;
            params[1] = a;
            params[2] = tau;
            break;

		// double exponential
        case 4:
        case 5:
            // params are Z, A1, T1, A2, T2
            n_param = 5;
            params = (float *)malloc((size_t) n_param * sizeof(float));
            params[0] = z;
            params[1] = 0.75f * a; // values from TRI2/SP
            params[2] = tau;
            params[3] = 0.25f * a;
            params[4] = 0.6666667f * tau;
            break;
        // triple exponential
        case 6:
            // params are Z, A1, T1, A2, T2, A3, T3
            n_param = 7;
            params = (float *)malloc((size_t) n_param * sizeof(float));
            params[0] = z;
            params[1] = 0.75f * a;
            params[2] = tau;
            params[3] = 0.1666667f * a;
            params[4] = 0.6666667f * tau;
            params[5] = 0.1666667f * a;
            params[6] = 0.3333333f * tau;
            break;
        // stretched exponential
        case 11:
            // has it's own fitfunc
            fitfunc = GCI_stretchedexp;
            
            // params are Z, A, T, H
            n_param = 4;
            params = (float *)malloc((size_t) n_param * sizeof(float));
            params[0] = z;
            params[1] = a;
            params[2] = tau;
            params[3] = 1.5f;
            break;
		default:
			goto END;
    }
    
    // param_free array describes which params are free vs fixed (omits X2)
    param_free = (int *)malloc((size_t) n_param  * sizeof(int));
    n_param_free = 0;
    for (i = 0; i < n_param; ++i) {
        param_free[i] = 1;
        ++n_param_free;
    }

    chi_sq_adjust = fit_end - fit_start - n_param_free;

	if (fit_type == 3) //Phasor
	{
		float u, v, taup, taum;

		return_value = GCI_Phasor(x_inc, transient_values,
                                   fit_start,
                                   fit_end, 
                                   &(params[0]), &u, &v, &taup, &taum, &(params[2]), fitted, residuals, &chi_square);
		printf("Phasor return value %d\n", return_value);
	}
	else
	{
		return_value = GCI_marquardt_fitting_engine(
                                                x_inc,
                                                transient_values,
                                                transient_size,
                                                fit_start,
                                                fit_end,
                                                prompt_values,
                                                prompt_size,
                                                noise_model,
                                                sigma_values,
                                                params,
                                                param_free,
                                                n_param,
                                                restrain,
                                                fitfunc,
                                                fitted,
                                                residuals,
                                                &chi_square,
                                                covar,
                                                alpha,
                                                err_axes,
                                                chi_sq_target * chi_sq_adjust,
                                                chi_sq_delta,
                                                chi_sq_percent);
    
 
		printf("LMA return value %d\n", return_value);
	}
    switch (fit_type) {
        // single exponential
        case 1:
        case 2:
            // params are Z, A, T
	        printf("LMA fitted A %f T %f Z %f X2 %f\n",
                   params[1], params[2], params[0], chi_square / chi_sq_adjust);
            break;
        case 3:
            // params are Z, T
	        printf("Phasor fitted T %f Z %f X2 %f\n",
                   params[2], params[0], chi_square / chi_sq_adjust);
            break;
        case 4:
        case 5:
        // double exponential
            // params are Z, A1, T1, A2, T2
	        printf("LMA fitted A1 %f T1 %f A2 %f T2 %f Z %f X2 %f\n",
                   params[1], params[2], params[3],
                   params[4], params[0], chi_square / chi_sq_adjust);
            break;
        // triple exponential
        case 6:
            // params are Z, A1, T1, A2, T2, A3, T3
	        printf("LMA fitted A1 %f T1 %f A2 %f T2 %f A3 %f T3 %f Z %f X2 %f\n",
                   params[1], params[2], params[3],
                   params[4], params[5], params[6],
                   params[0], chi_square / chi_sq_adjust);
            break;
        // stretched exponential
        case 11:
            // params are Z, A, T, H
	        printf("LMA fitted A %f T %f H %f Z %f X2 %f\n",
                   params[1], params[2], params[3],
                   params[0], chi_square / chi_sq_adjust);
            break;
		default:;
    }
 
END: 
    free(params);
    free(fitted);
    free(residuals);
	GCI_ecf_free_matrix(covar);
	GCI_ecf_free_matrix(alpha);
	GCI_ecf_free_matrix(err_axes);
    
    return SUCCESS;
}

/* Adjusts transient array to skip time bins outside interval of interest.
 */
float *adjust_transient(float *transient_values, int transient_start_index, int transient_end_index) {
    int new_size = transient_end_index - transient_start_index + 1;
	float *new_transient = (float *)malloc((unsigned)new_size * sizeof(float));
    int i;
    
    for (i = 0; i < new_size; ++i) {
        new_transient[i] = transient_values[i + transient_start_index];
    }
    free(transient_values);
    return new_transient;
}

/* Load a *.dat file into an array of floats, where the first value is the number of values.
*/
int load_datfile (const char *filename, int *nVals, float **arrayPtr)
{
	int n, readVals, i;
	float *array;
	FILE *file=fopen(filename, "r");
	
	if (NULL == file)
		return -1;  // error opening file
		
	// Read the first value in the file which should be the number of floats to read
	readVals = fscanf(file, "%d", &n);
	
	if (readVals < 1) {
		fclose(file);
		return -2;  // error reading from file
	}

	if (n<1 || n>MAX_DATFILE_LENGTH) {
		fclose(file);
		return -3;  // max datfile length exceeded, or bad first entry
	}
	
	array = (float *)malloc(n*sizeof(float));
	for (i=0; i<n; i++) {
		readVals = fscanf(file, "%f", &(array[i]));
		if (readVals < 1)
			break;  // somehow no more values to read, should not happen
	}
		
	if (nVals!=NULL)
		*nVals=i;
		
	if (arrayPtr!=NULL)
		*arrayPtr=array;

	return 0;
}

/* Parses ".ini" file and calls SLIM Curve fitting code.
 */
int parse_and_fit(dictionary *ini, const char *dat_filename) {
    int debug = FALSE; //TRUE;
    int prompt_size = 0;
    float *prompt_values = NULL;
    int transient_size = 0;
    float *transient_values = NULL;
	float x_inc;
    int sigma_size = 0;
    float *sigma_values = NULL;
    float chi_sq_target;
    float chi_sq_delta;
    float prompt_baseline;
    float transient_start;
    float data_start;
    float transient_end;
    float prompt_delta;
    float prompt_width;
    int transient_start_index;
    int data_start_index;
    int transient_end_index;
	int fit_start;
	int fit_end;
    int fit_type;
    int noise_model;
	noise_type noise;
	int i;
    float version;
    char *s;

    s = iniparser_getstring(ini, ":filetype", NULL);
    printf("%s\n", s?s:"filetype not found.");
    version = (float)iniparser_getdouble(ini, ":version", 0.0);

    if (1.0f != version) {
        printf("Warning, unknown version %f, expecting 1.0\n", version);
    }
    
    // look for optional '[prompt]: file'
    s = iniparser_getstring(ini, "prompt:file", NULL);
    if (s) {
    	if (load_datfile(s, &prompt_size, &prompt_values)<0) {
            printf("Problem parsing prompt file.\n");
            return BAD_SYNTAX;
		}
        else if (debug) {
            printf("prompt: %d vals\n", prompt_size);
            for (i = 0; i < prompt_size; ++i) {
                printf("prompt %d is %f\n", i, prompt_values[i]);
            }
        }
    }
	
	// Load the transient
	if (load_datfile(dat_filename, &transient_size, &transient_values)<0) {
        printf("Problem parsing prompt file.\n");
        return BAD_SYNTAX;
	}
    else if (debug) {
        printf("transient: %d vals\n", transient_size);
        for (i = 0; i < transient_size; ++i) {
            printf("transient %d is %f\n", i, transient_values[i]);
        }
    }
           
    // look for required '[main]' section
    chi_sq_target = (float)iniparser_getdouble(ini, "main:chisqtarget", 1.0);
    chi_sq_delta = (float)iniparser_getdouble(ini, "main:chisqdelta", 0.000001);
	x_inc = (float)iniparser_getdouble(ini, "main:x_inc", 1.0);
    
	if (x_inc==1.0f) {
        printf("WARNING! x_inc is %f, there may not be an \'x_inc\' entry in the [main] section of the settings file.\n", x_inc);
	}
	
	if (debug) {
        printf("inc is %f\n", x_inc);
        printf("chi_sq_target %f\n", chi_sq_target);
        printf("chi_sq_delta %f\n", chi_sq_delta);
    }

    // look for required '[cursors]' section
    prompt_baseline = (float)iniparser_getdouble(ini, "cursors:pbaseline", 0.0);
    transient_start = (float)iniparser_getdouble(ini, "cursors:tstart", 1.0);
    data_start = (float)iniparser_getdouble(ini, "cursors:dstart", 1.0);
    transient_end = (float)iniparser_getdouble(ini, "cursors:tend", 10.0);
    prompt_delta = (float)iniparser_getdouble(ini, "cursors:pdelta", 0.0);
    prompt_width = (float)iniparser_getdouble(ini, "cursors:pwidth", 0.1);
    if (debug) {
        printf("prompt_baseline %f\n", prompt_baseline);
        printf("transient_start %f\n", transient_start);
        printf("data_start %f\n", data_start);
        printf("transient_end %f\n", transient_end);
        printf("prompt_delta %f\n", prompt_delta);
        printf("prompt_width %f\n", prompt_width);
    }
    
    // look for required '[fit]' section
    fit_type = iniparser_getint(ini, "fit:type", 1);
    noise_model = iniparser_getint(ini, "fit:noisemodel", 5);
	noise = getNoiseModel(noise_model);
    if (debug) {
        printf("type %d\n", fit_type);
        printf("noisemodel %d\n", noise_model);
    }

    // massage values
    transient_start_index = (int) my_roundf(transient_start / x_inc);
    data_start_index = (int) my_roundf(data_start / x_inc);
    transient_end_index = (int) my_roundf(transient_end / x_inc);
    transient_values = adjust_transient(transient_values, transient_start_index, transient_end_index);
    transient_size = transient_end_index - transient_start_index + 1;
	fit_start = data_start_index - transient_start_index;
	fit_end = transient_end_index - transient_start_index;
    if (debug) {
        printf("transient_start_index %d\n", transient_start_index);
        printf("data_start_index %d\n", data_start_index);
        printf("transient_end_index %d\n", transient_end_index);
        printf("fit_start %d\n", fit_start);
        printf("fit_end %d\n", fit_end);
    }
    
    // do the fit
    fit(
        fit_type, noise, chi_sq_target, chi_sq_delta,
        transient_size, transient_values,
        sigma_size, sigma_values,
        prompt_size, prompt_values,
        x_inc, fit_start, fit_end
    );

    // finished
    return SUCCESS;
}

noise_type getNoiseModel (int val)
{  // This is from TRI2 and must match the noise models, between TRI2 ui and SLIM Curve - This is not a good thing to have! TRI2 should save the SLIM Curve type

	switch (val)
	{
		case 1:
			return NOISE_GAUSSIAN_FIT;
			break;
		case 2:
			return NOISE_POISSON_DATA;				
			break;
		case 3:
			return NOISE_POISSON_FIT;				
			break;
		case 4:
			return NOISE_MLE;				
			break;
		default:
			return NOISE_POISSON_FIT;				
			break;
	}
}

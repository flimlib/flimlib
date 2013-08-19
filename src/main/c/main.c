/* 
This file is part of the SLIM-curve package for exponential curve fitting of spectral lifetime data.

Copyright (c) 2010-2013, Gray Institute University of Oxford & UW-Madison LOCI.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Ecf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FALSE 0
#define TRUE !FALSE
#define MAX_STRING_LENGTH 80
#define SUCCESS 0
#define BAD_USAGE -1
#define BAD_FILE -2
#define BAD_SYNTAX -3
#define UNEXPECTED_EOF -4

int fit(float x_inc, int fit_start, int fit_end, int noise, int transient_size, float *transient_values, int sigma_size, float *sigma_values, int prompt_size, float *prompt_values);
int parse_and_fit(FILE *file);
int get_section_name(FILE *file, char *section);
int validate_section(FILE *file, char *section);
char *get_string_value(FILE *file, char *tag);
int get_int_value(FILE *file, char *tag);
float get_float_value(FILE *file, char *tag);
int get_int_array_value(FILE *file, char *tag, int size, int *int_array);
int get_float_array_value(FILE *file, char *tag, int size, float *float_array);
int validate_tag(FILE *file, char *tag);
int get_next_value_as_string(FILE *file, char *buffer);
char next_non_whitespace(FILE *file);
char next_char(FILE *file);

char stringbuffer[MAX_STRING_LENGTH];

int main(int argc, const char * argv[])
{
    int return_value = SUCCESS;
    if (argc != 2) {
        printf("Usage: %s filename\n", argv[0]);
        return_value = BAD_USAGE;
    }
    else {
        FILE *file = fopen(argv[1], "r");
        if (0 == file) {
            printf("Could not open file '%s'\n", argv[1]);
            return_value = BAD_FILE;
        }
        else {
            printf("Parsing %s\n", argv[1]);
            return_value = parse_and_fit(file);
            fclose(file);
        }
    }
    return return_value;
}

int fit(
        float x_inc, int fit_start, int fit_end, int noise,
        int transient_size, float *transient_values,
		int sigma_size, float *sigma_values,
		int prompt_size, float *prompt_values
       ) {
	int return_value;
	float a, tau, z;
	float *fitted = (float *)malloc((unsigned)transient_size * sizeof(float));
	float *residuals = (float *)malloc((unsigned)transient_size * sizeof(float));
	float chi_square;
	float chi_square_target = 1.25f;
	
	return_value =  GCI_triple_integral_fitting_engine(
			x_inc,
			transient_values,
			fit_start,
			fit_end,
			prompt_values,
			prompt_size,
			noise,
			sigma_values,
			&z,
			&a,
			&tau,
			fitted,
			residuals,
			&chi_square,
			chi_square_target
			);
	
	printf("return_value %d\n", return_value);
	printf("a %f tau %f z %f\n", a, tau, z);
    return SUCCESS;
}

/* Parses ".ini" file and calls SLIM Curve fitting code.
 */
int parse_and_fit(FILE *file) {
    int debug = TRUE;
    int prompt_size = 0;
    float *prompt_values = NULL;
    int transient_size = 0;
    float *transient_values = NULL;
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
    int fit_type;
    int noise_model;
	int i;
    
    printf("%s\n", get_string_value(file, "filetype"));
    float version = get_float_value(file, "version");
    if (1.0f != version) {
        printf("Warning, unknown version %f, expecting 1.0\n", version);
    }
    if (!get_section_name(file, stringbuffer)) {
        printf("Missing '[transient]' section\n");
        return BAD_SYNTAX;
    }
    
    // look for optional '[prompt]' section
    if (0 == strcmp(stringbuffer, "prompt")) {
        prompt_size = get_int_value(file, "size");
        prompt_values = (float *) malloc((size_t) prompt_size * sizeof(float));
        if (!get_float_array_value(file, "values", prompt_size, prompt_values)) {
            printf("Problem parsing prompt values\n");
            return BAD_SYNTAX;
        }
        else if (debug) {
            printf("prompt:\n");
            for (i = 0; i < prompt_size; ++i) {
                printf("prompt %d is %f\n", i, prompt_values[i]);
            }
        }
                
        // load next section name
        get_section_name(file, stringbuffer);
    }

    // look for required '[transient]' section
    if (0 != strcmp(stringbuffer, "transient")) {
        printf("Missing '[transient]' section\n");
        return BAD_SYNTAX;
    }
    transient_size = get_int_value(file, "size");
    transient_values = (float *) malloc((size_t) transient_size * sizeof(float));
    if (!get_float_array_value(file, "values", transient_size, transient_values)) {
        printf("Problem parsing transient values\n");
        return BAD_SYNTAX;
    }
    else if (debug) {
        printf("transient:\n");
        for (i = 0; i < transient_size; ++i) {
            printf("transient %d is %f\n", i, transient_values[i]);
        }
    }
        
    if (!get_section_name(file, stringbuffer)) {
        printf("Missing '[transient]' section.\n");
        return BAD_SYNTAX;
    }
    
    // look for optional '[sigma]' section
    if (0 == strcmp(stringbuffer, "sigma")) {
        sigma_size = get_int_value(file, "size");
        sigma_values = (float *) malloc((size_t) sigma_size * sizeof(float));
        if (!get_float_array_value(file, "values", sigma_size, sigma_values)) {
            printf("Problem parsing sigma values\n");
            return BAD_SYNTAX;
        }
        else if (debug) {
            printf("sigma:\n");
            for (i = 0; i < sigma_size; ++i) {
                printf("sigma %d is %f\n", i, sigma_values[i]);
            }
        }
        
        // load next section name
        get_section_name(file, stringbuffer);
    }

    // look for required '[main]' section
    if (0 != strcmp(stringbuffer, "main")) {
        printf("Missing '[main]' section\n");
        return BAD_SYNTAX;
    }
    chi_sq_target = get_float_value(file, "chisqtarget");
    chi_sq_delta = get_float_value(file, "chisqdelta");
    if (debug) {
        printf("chi_sq_target %f\n", chi_sq_target);
        printf("chi_sq_delta %f\n", chi_sq_delta);
    }

    // look for required '[cursors]' section
    if (!validate_section(file, "cursors")) {
        return BAD_SYNTAX;
    }
    prompt_baseline = get_float_value(file, "pbaseline");
    transient_start = get_float_value(file, "tstart");
    data_start = get_float_value(file, "dstart");
    transient_end = get_float_value(file, "tend");
    prompt_delta = get_float_value(file, "pdelta");
    prompt_width = get_float_value(file, "pwidth");
    if (debug) {
        printf("prompt_baseline %f\n", prompt_baseline);
        printf("transient_start %f\n", transient_start);
        printf("data_start %f\n", data_start);
        printf("transient_end %f\n", transient_end);
        printf("prompt_delta %f\n", prompt_delta);
        printf("prompt_width %f\n", prompt_width);
    }
    
    // look for required '[fit]' section
    if (!validate_section(file, "fit")) {
        return BAD_SYNTAX;
    }
    fit_type = get_int_value(file, "type");
    noise_model = get_int_value(file, "noisemodel");
    if (debug) {
        printf("prompt_baseline %d\n", fit_type);
        printf("transient_start %d\n", noise_model);
    }
    
    // do the fit
	//TESTING
    fit(
			0.1f, 10, 60, noise_model,
			transient_size, transient_values,
			sigma_size, sigma_values,
			prompt_size, prompt_values
		);

    // finished
    return SUCCESS;
}

/* Parses for next section heading name.  Useful for
 * optional sections.
 *
 * Returns TRUE if one is present, FALSE otherwise.
 */
int get_section_name(FILE *file, char *buffer) {
    int return_value = FALSE;
    int index = 0;
    char c = next_non_whitespace(file);
    if ('[' == c) {
        while (TRUE) {
            c = next_char(file);
            if (']' == c) {
                break;
            }
            buffer[index++] = c;
        }
        return_value = TRUE;
    }
    buffer[index] = 0;
    return return_value;
}

/* Parses for required section heading.
 *
 * Returns TRUE if present, FALSE otherwise.
 */
int validate_section(FILE *file, char *section) {
    int return_value = FALSE;
    char c = next_non_whitespace(file);
    if ('[' == c) {
        while (0 != *section) {
            c = next_char(file);
            if (c != *section) {
                printf("Missing '[%s]' section\n", section);
                break;
            }
            ++section;
        }
        if (0 == *section) {
            c = next_char(file);
            if (']' == c) {
                return_value = TRUE;
            }
        }
    }
    return return_value;
}

/* Returns string value of tag or null.
 */
char *get_string_value(FILE *file, char *tag) {
    if (!validate_tag(file, tag)) {
        return NULL;
    }
    get_next_value_as_string(file, stringbuffer);
    return &stringbuffer[0];
}

/* Returns integer value of tag.
 */
int get_int_value(FILE *file, char *tag) {
    int return_value = 0;
    if (validate_tag(file, tag)) {
        get_next_value_as_string(file, stringbuffer);
        sscanf(stringbuffer, "%d", &return_value);
    }
    return return_value;
}

/* Returns float value of tag.
 */
float get_float_value(FILE *file, char *tag) {
    float return_value = 0.0f;
    if (validate_tag(file, tag)) {
        get_next_value_as_string(file, stringbuffer);
        sscanf(stringbuffer, "%f", &return_value);
        
    }
    return return_value;
}

/* Builds integer array value of tag.
 */
int get_int_array_value(FILE *file, char *tag, int size, int *int_array) {
    int return_value = FALSE;
    if (validate_tag(file, tag)) {
        int index = 0;
        while (index < size) {
            get_next_value_as_string(file, stringbuffer);
            sscanf(stringbuffer, "%d", &int_array[index++]);
        }
        return_value = TRUE;
    }
    return return_value;
}

/* Builds float array value of tag.
 */
int get_float_array_value(FILE *file, char *tag, int size, float *float_array) {
    int return_value = FALSE;
    if (validate_tag(file, tag)) {
        int index = 0;
        while (index < size) {
            get_next_value_as_string(file, stringbuffer);
            sscanf(stringbuffer, "%f", &float_array[index++]);
        }
        return_value = TRUE;
    }
    return return_value;
}

/* Parses for required tag, followed by '='.
 *
 * Returns TRUE if present, FALSE otherwise.
 */
int validate_tag(FILE *file, char *tag) {
    int first = 1;
    int c;
    while (0 != *tag) {
        if (first) {
            first = 0;
            c = next_non_whitespace(file);
        }
        else {
            c = next_char(file);
        }
        if (c != *tag) {
            printf("Missing tag %s\n", tag);
            return FALSE;
        }
        ++tag;
    }
    c = next_non_whitespace(file);
    if ('=' != c) {
        printf("Missing '=' after tag %s\n", tag);
        return FALSE;
    }
    return TRUE;
}

/* Builds string with value in buffer.
 */
int get_next_value_as_string(FILE *file, char *buffer) {
    int index = 0;
    char c = next_non_whitespace(file);
    if ('"' == c) {
        c = next_char(file);
    }
    while ('"' != c && ',' != c) {
        buffer[index++] = c;
        c = next_char(file);
    }
    buffer[index] = 0;
    return TRUE;
}

/* Returns next non-whitespace character in file.
 */
char next_non_whitespace(FILE *file) {
    char c;
    while (TRUE) {
        c = next_char(file);
        if (' ' != c && '\t' != c &&  '\n' != c && '\r' != c) {
            return c;
        }
    }
}

/* Returns next character in file, checking for EOF.
 */
char next_char(FILE *file) {
    int c = fgetc(file);
    if (EOF == c) {
        printf("Unexpected EOF\n");
        fclose(file);
        exit(UNEXPECTED_EOF);
    }
    return (char) c;
}


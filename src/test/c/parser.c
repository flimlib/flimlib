/*-
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2025 University of Oxford and Board of Regents of the
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
 *  parser.c
 *  Untitled
 *
 *  Created by Aivar Grislis on 6/2/11.
 */

#include "parser.h"
#include "../../main/c/Ecf.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Puts out an error message if 'sscanf' or 'strchr' failed.
 *
 * param count of numbers scanned
 * param type description of number scanned
 * param input string
 */
void checkError(int count, char *type, char *input) {
	if (NULL == input) {
		printf("Missing comma delimiter\n");
	}
	else if (1 != count) {
		printf("Unable to scan %s from %s\n", type, input);
	}
}

/*
 * Gets an integer from a string.
 *
 * param non-NULL input string
 * returns integer value
 */
int getIntFromString(char *input) {
	int output = 0;
	int count = sscanf(input, "%d", &output);
	checkError(count, "integer", input);
	return output;
}

/*
 * Gets an array of integers from a string.  Result needs to be freed.
 *
 * param non-NULL input string
 * param size of array
 * returns array of integers
 */
int *getInt1DFromString(char *input, int size) {
	int *output = (int *) malloc(size * sizeof(int));
	int i;
	for (i = 0; i < size; ++i) {
		int count = sscanf(input, "%d", &output[i]);
	    checkError(count, "integer", input);
		input = strchr(input, ',');
		++input;
	}
	return output;
}

/*
 * Gets a float from a string.
 *
 * param non-NULL input string
 * returns float value
 */
float getFloatFromString(char *input) {
	float output = 0.0;
	int count = sscanf(input, "%f", &output);
	checkError(count, "float", input);
	return output;
}

/*
 * Gets an array of floats from a string.  Result needs to be freed.
 *
 * param non-NULL input string
 * param size of array
 * returns array of floats
 */
float *getFloat1DFromString(char *input, int size) {
	int i;
    float *output = (float *) malloc(size * sizeof(float));
    for (i = 0; i < size; ++i) {
        int count = sscanf(input, "%f", &output[i]);
		checkError(count, "float", input);
        input = strchr(input, ',');
        ++input;
    }
    return output;
}

/*
 * Gets a matrix of floats from a string.  Result needs to be freed
 * using GCI_ecf_free_matrix().
 */
float **getFloat2DFromString(char *input, int rows, int cols) {
	int row, col;
	float **output = GCI_ecf_matrix(rows, cols);
	for (row = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col) {
			int count = sscanf(input, "%f", &output[row][col]);
			checkError(count, "float", input);
			input = strchr(input, ',');
			++input;
		}
	}
	return output;
}

/*
 * Given a cursor to a value object, checks the value is an object,
 * looks within the object for a name/value with given name,
 *
 * param cursor points to parent value object
 * param name gives name of child
 * returns cursor to value or NULL if not found
 */
json_t *getNamedChildValue(json_t *cursor, char *name) {
	json_t *returnValue = NULL;
	if (NULL != cursor && cursor->type == JSON_OBJECT) { // should be value object
		cursor = cursor->child;
		while (NULL != cursor && cursor->type == JSON_STRING) {
			if (strcmp(name, cursor->text) == 0) {
				returnValue = cursor->child; // point to value
				break;
			}
			cursor = cursor->next;
		}
	}
	return returnValue;
}

/*
 * Given a cursor to a value object, checks the value is an object,
 * looks within the object for a name/value with given name,
 * checks the value is a string, returns that string value.
 *
 * param cursor points to parent value object
 * param name gives name of child
 * returns string or NULL if not found
 */
char *getString(json_t *cursor, char *name) {
	char *returnValue = NULL;
	
	cursor = getNamedChildValue(cursor, name);
	if (NULL != cursor && JSON_STRING == cursor->type) {
		returnValue = cursor->text;
	}
	return returnValue;
}

/*
 * Finds the child with given name and returns the integer value.
 *
 * param cursor points to parent value object
 * param name gives name of child
 * param notFoundValue return value if child not present
 * returns integer value or notFoundValue if not found
 */
int getInt(json_t *cursor, char *name, int notFoundValue) {
	int returnValue = notFoundValue;
	char *string = getString(cursor, name);
	if (NULL != string) {
		returnValue = getIntFromString(string);
	}
	return returnValue;
}

/*
 * Finds the child with given name and value and returns the integer array.
 *
 * param cursor points to parent name/value object
 * param name gives name of child
 * param size returns size of array
 * returns integer array
 */
int *getInt1D(json_t *cursor, char *name, int *size) {
	int *returnValue = NULL;
	char *string;
	
	*size = 0;
	cursor = getNamedChildValue(cursor, name);
	if (NULL != cursor) {
		*size = getInt(cursor, SIZE, 0);
		string = getString(cursor, VALUE);
		if (NULL != string && *size > 0) {
			returnValue = getInt1DFromString(string, *size);
		}
	}
	return returnValue;
}

/*
 * Finds the child with given name and returns the float value.
 *
 * param cursor points to parent value object
 * param name gives name of child
 * param notFoundValue return value if child is not present
 * returns float value or notFoundValue
 */
float getFloat(json_t *cursor, char *name, float notFoundValue) {
	float returnValue = notFoundValue;
		
	char *string = getString(cursor, name);
	if (NULL != string) {
		returnValue = getFloatFromString(string);
	}
	
    return returnValue;
}

/*
 * Finds the child with given name and value and returns the float array.
 *
 * param cursor points to parent name/value object
 * param name gives name of child
 * param size returns size of array
 * returns float array or NULL if not found
 */
float *getFloat1D(json_t *cursor, char *name, int *size) {
	float *returnValue = NULL;
	char *string;
	
	*size = 0;
	cursor = getNamedChildValue(cursor, name);
	if (NULL != cursor) {
		*size = getInt(cursor, SIZE, 0);
		string = getString(cursor, VALUE);
		if (NULL != string && *size > 0) {
			returnValue = getFloat1DFromString(string, *size);
		}
		//returnValue = getFloat1DFromString(getString(cursor, "value"), *size);
	}
	return returnValue;
}

/*
 * Finds the child with given name and value and returns the float matrix.
 *
 * param cursor points to parent name/value object
 * param name gives name of child
 * param rows returns row count of matrix
 * param cols returns column count of matrix
 * returns float matrix or NULL if not found
 */
float **getFloat2D(json_t *cursor, char *name, int *rows, int *cols) {
	float **returnValue = NULL;
	char *string;
	
	*rows = 0;
	*cols = 0;
	cursor = getNamedChildValue(cursor, name);
	if (NULL != cursor) {
		*rows = getInt(cursor, ROWS, 0);
		*cols = getInt(cursor, COLS, 0);
		string = getString(cursor, VALUE);
		if (NULL != string && *rows > 0 && *cols > 0) {
			returnValue = getFloat2DFromString(string, *rows, *cols);
		}
	}
	return returnValue;
}

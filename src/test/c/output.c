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
 *  output.c
 *  Untitled
 *
 *  Created by Aivar Grislis on 6/7/11.
 *
 */

#include "output.h"
#include <assert.h>
#include <stdlib.h>

/*
 * Puts out a integer to a string.  Result needs to be freed.
 *
 * param value
 * returns string pointer
 */
char *putIntToString(int value) {
	int size;
	char *string;

    size = MAXINTSIZE;
	string = (char *) malloc(size);
// not in C88:	snprintf(string, 16, "%d", value);
	sprintf(string, "%d", value);
	return string;
}

/*
 * Puts out an array of integers to a string.  Result needs to be freed.
 *
 * param array of integers
 * param size of array
 * returns string pointer
 */
char *putInt1DToString(int *array, int size) {
    char temp[16]; //TODO ARG if I put MAXINTSIZE here won't compile NetBeans OS X
    char *string = (char *) malloc(size * 16);
    char *next = string;
    int i;
    for (i = 0; i < size; ++i) {
        if (i < size - 1) {
// not in C88:            snprintf(temp, 16, "%d, ", array[i]); //TODO ARG again MAXINTSIZE won't work here
			sprintf(temp, "%d, ", array[i]);
        }
        else {
// not in C88:            snprintf(temp, 16, "%d", array[i]);
			sprintf(temp, "%d", array[i]);
        }
        strcpy(next, temp);
        next += strlen(temp);
    }
    return string;
}


/*
 * Puts out a float to a string.  Result needs to be freed.
 *
 * param value
 * returns string pointer
 */
char *putFloatToString(float value) {
	int size;
	char *string;

    size = MAXFLOATSIZE;
	string = (char *) malloc(size);
//not in C88:	snprintf(string, 16, "%f", value);
	sprintf(string, "%f", value);
	return string;
}

/*
 * Puts out an array of floats to a string.  Result needs to be freed.
 *
 * param array of floats
 * param size of array
 * returns string pointer
 */
char *putFloat1DToString(float *array, int size) {
    char temp[16]; //TODO ARG if I put MAXFLOATSIZE here won't compile NetBeans OS X
    char *string = (char *) malloc(size * 16); //TODO ARG nor here
    char *next = string;
    int i;
    for (i = 0; i < size; ++i) {
        if (i < size - 1) {
// not in C88:            snprintf(temp, 16, "%f, ", array[i]); //TODO ARG again MAXFLOATSIZE won't work here
			sprintf(temp, "%f, ", array[i]);
        }
        else {
// not in C88:            snprintf(temp, 16, "%f", array[i]);
			sprintf(temp, "%f", array[i]);
        }
        strcpy(next, temp);
        next += strlen(temp);
    }
    return string;
}

/*
 * Puts out a matrix of floats to a string.  Result needs to be freed.
 *
 * param matrix of floats
 * param rows in matrix
 * param columns in matrix
 * returns string pointer
 */
char *putFloat2DToString(float **matrix, int rows, int cols) {
    char temp[16]; //TODO ARG if I put MAXFLOATSIZE here won't compile NetBeans OS X
	int size;
	char *string;
	char *next;
	int row, col;

    size = rows * cols * MAXFLOATSIZE;
    string = (char *) malloc(size);
    next = string;
	for (row = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col) {
			if (row == rows - 1 && col == cols - 1) {
// not in C88:                snprintf(temp, 16, "%f",matrix[row][col]);
				sprintf(temp, "%f", matrix[row][col]);
			}
			else {
// not in C88:                snprintf(temp, 16, "%f, ", matrix[row][col]);
				sprintf(temp, "%f, ", matrix[row][col]);
			}
			strcpy(next, temp);
			next += strlen(temp);													   
		}
	}
    return string;
}	

/*
 * Creates JSON name/value pair from a string.
 *
 * param name
 * param data string value
 * return pointer to JSON name/value pair
 */
json_t *putString(char *name, char *data) {
	enum json_error error;
	json_t *label, *value;
	assert(NULL != name);
	assert(NULL != data);
	label = json_new_string(name);
	value = json_new_string(data);
	error = json_insert_child(label, value);
	assert (JSON_OK == error);
	return label;
}

/*
 * Creates JSON name/value pair from an integer.
 *
 * param name
 * param value
 * return pointer to JSON name/value pair
 */
json_t *putInt(char *name, int value) {
	char *output = putIntToString(value);
	return putString(name, output);
}

/*
 * Creates JSON name/value pair from array of integers.
 *
 * param name
 * param array of integers
 * param size of array
 * return pointer to JSON name/value pair
 */
json_t *putInt1D(char *name, /*const*/ int *array, int size) { //TODO ARG why start using const here?
	enum json_error error;
	json_t *arraySize;
	json_t *arrayValue;
	json_t *object;
	json_t *label;
	
	arraySize = putInt(SIZE, size);
	arrayValue = putString(VALUE, putInt1DToString(array, size));
	
	object = json_new_object();
	error = json_insert_child(object, arraySize);
	assert(JSON_OK == error);
	error = json_insert_child(object, arrayValue);
	assert(JSON_OK == error);
	
	label = json_new_string(name);
	error = json_insert_child(label, object);
	assert(JSON_OK == error);
	
	return label;
}

/*
 * Creates JSON name/value pair from a float.
 *
 * param name
 * param value
 * return pointer to JSON name/value pair
 */
json_t *putFloat(char *name, float value) {
	char *output = putFloatToString(value);
	return putString(name, output);
}

/*
 * Creates a JSON name/value pair from array of floats.
 *
 * param name
 * param array of floats
 * param size of array
 * return pointer to JSON name/value pair
 */

json_t *putFloat1D(char *name, float *array, int size) {
	enum json_error error;
	json_t *arraySize;
	json_t *arrayValue;
	json_t *object;
	json_t *label;
	
	arraySize = putInt(SIZE, size);
	arrayValue = putString(VALUE, putFloat1DToString(array, size));
	
	object = json_new_object();
	error = json_insert_child(object, arraySize);
	assert(JSON_OK == error);
	error = json_insert_child(object, arrayValue);
	assert(JSON_OK == error);
	
	label = json_new_string(name);
	error = json_insert_child(label, object);
	assert(JSON_OK == error);
	
	return label;
}

/*
 * Creates a JSON name/value pair from matrix of floats.
 *
 * param name
 * param matrix of floats
 * param rows in matrix
 * param cols in matrix
 * return pointer to JSON name/value pair
 */

json_t *putFloat2D(char *name, float **matrix, int rows, int cols) {
	enum json_error error;
	json_t *matrixRows;
	json_t *matrixCols;
	json_t *matrixValue;
	json_t *object;
	json_t *label;
	
	matrixRows = putInt(ROWS, rows);
	matrixCols = putInt(COLS, cols);
	matrixValue = putString(VALUE, putFloat2DToString(matrix, rows, cols));
	
	object = json_new_object();
	error = json_insert_child(object, matrixRows);
	assert(JSON_OK == error);
	error = json_insert_child(object, matrixCols);
	assert(JSON_OK == error);
	error = json_insert_child(object, matrixValue);
	assert(JSON_OK == error);
	
	label = json_new_string(name);
	error = json_insert_child(label, object);
	assert(JSON_OK == error);
	
	return label;
}


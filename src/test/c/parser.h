/*
 *  parser.h
 *  Untitled
 *
 *  Created by Aivar Grislis on 6/2/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "strings.h"
#include "json.h"

#define FALSE 0
#define TRUE !FALSE

/*
 * Given a cursor to a value object, checks the value is an object,
 * looks within the object for a name/value with given name,
 *
 * param cursor points to parent value object
 * param name gives name of child
 * returns cursor to value or NULL if not found
 */
json_t *getNamedChildValue(json_t *cursor, char *name);

/*
 * Given a cursor to a value object, checks the value is an object,
 * looks within the object for a name/value with given name,
 * checks the value is a string, returns that string value.
 *
 * param cursor points to parent value object
 * param name gives name of child
 * returns string or NULL if not found
 */
char *getString(json_t *cursor, char *name);

/*
 * Finds the child with given name and returns the integer value.
 *
 * param cursor points to parent value object
 * param name gives name of child
 * param notFoundValue return value if child not present
 * returns integer value or notFoundValue if not found
 */
int getInt(json_t *cursor, char *name, int notFoundValue);

/*
 * Finds the child with given name and returns the integer array value.
 *
 * param cursor points to parent name/value object
 * param name gives name of child
 * param size returns size of array
 * returns integer array
 */
int *getInt1D(json_t *cursor, char *name, int *size);

/*
 * Finds the child with given name and returns the float value.
 *
 * param cursor points to parent value object
 * param name gives name of child
 * param notFoundValue return value if child is not present
 * returns float value or notFoundValue
 */
float getFloat(json_t *cursor, char *name, float notFoundValue);

/*
 * Finds the child with given name and returns the float array value.
 *
 * param cursor points to parent name/value object
 * param name gives name of child
 * param size returns size of array
 * returns float array or NULL if not found
 */
float *getFloat1D(json_t *cursor, char *name, int *size);

/*
 * Finds the child with given name and returns the float matrix value.
 *
 * param cursor points to parent name/value object
 * param name gives name of child
 * param rows returns row count of matrix
 * param cols returns column count of matrix
 * returns float matrix or NULL if not found
 */
float **getFloat2D(json_t *cursor, char *name, int *rows, int *cols);


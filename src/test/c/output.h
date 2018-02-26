/*
 *  output.h
 *  Untitled
 *
 *  Created by Aivar Grislis on 6/7/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "strings.h"
#include "json.h"
#include <string.h>

#define MAXFLOATSIZE 16;
#define MAXINTSIZE 16;

json_t *putFloat2D(char *name, float **matrix, int rows, int cols);
json_t *putFloat1D(char *name, float *array, int size);
json_t *putFloat(char *name, float value);
json_t *putInt1D(char *name, int *array, int size);
json_t *putInt(char *name, int value);

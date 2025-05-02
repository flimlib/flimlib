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
 *  output.h
 *  Untitled
 *
 *  Created by Aivar Grislis on 6/7/11.
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

/*-
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
/*=================================================================================*/
/* File:       bayes_DataManagement.h                                              */
/*---------------------------------------------------------------------------------*/
/* Purpose:    General purpose header file for size defines.                       */
/*---------------------------------------------------------------------------------*/
/* References: None.                                                               */
/*---------------------------------------------------------------------------------*/
/* Notes:      None.                                                               */
/*=================================================================================*/
/* Revision history:                                                               */
/*---------------------------------------------------------------------------------*/
/* Date   | Modification                                                           */
/*---------------------------------------------------------------------------------*/
/* 070906 | Creation, mrowley.                                                     */
/*---------------------------------------------------------------------------------*/
/* 080102 | Introduced double type checking functionality.                         */
/*---------------------------------------------------------------------------------*/
/*        |                                                                        */
/*=================================================================================*/


/***********************************************************************************/
/*                                                                                 */
/*                             DEFINES / SWITCHES                                  */
/*                                                                                 */
/***********************************************************************************/

#define BAYES_DM_VALID_DATA_NOT_AVAILABLE   0
#define BAYES_DM_VALID_DATA_AVAILABLE       1

#define BAYES_DM_DOUBLE_TYPE_INVALID        0
#define BAYES_DM_DOUBLE_TYPE_VALID          1

/* Moved here from the bayes version of DTYPE.h - 091028 */
#define BAYES_DM_DTYPE_NEGATIVE_INFINITY      -1
#define BAYES_DM_DTYPE_REGULAR_NUMBER          0
#define BAYES_DM_DTYPE_POSITIVE_INFINITY       1
#define BAYES_DM_DTYPE_NOT_A_NUMBER            2
#define BAYES_DM_DTYPE_UNDERFLOW               3
#define BAYES_DM_DTYPE_ZERO                    4


/***********************************************************************************/
/*                                                                                 */
/*                             FUNCTION PROTOTYPES                                 */
/*                                                                                 */
/***********************************************************************************/

float data_ComputeBinnedDataAverageArrTime(int* data,
                                           int    nbins,
										   int    fitstart,
                                           int    nphotons,
                                           float  interval);

int bayes_dm_CheckDoubleValueValid(double Value, int *Type);

void bayes_dm_CorrectInvalidDoubleValue(double *Value, int Type);

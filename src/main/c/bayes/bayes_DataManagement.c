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
/* File:       bayes_DataManagement.c                                              */
/*---------------------------------------------------------------------------------*/
/* Purpose:    Central storage for the global data required for the routines       */
/*             implementing Bayesian data analysis.                                */
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
/* 090728 | Bayes version 2.0 updates (new data management interface implemented). */
/*---------------------------------------------------------------------------------*/
/*        |                                                                        */
/*=================================================================================*/

/***********************************************************************************/
/*                                                                                 */
/*                             HEADER FILE INCLUSION                               */
/*                                                                                 */
/***********************************************************************************/

//#include "stdio.h"
//#include "string.h"
//#include "time.h"

#include "bayes_DataManagement.h"
#include "bayes_Sizes.h"

#include "matrices.h"
#include "DTYPE.h"

float data_ComputeBinnedDataAverageArrTime(int   *data,
                                           int    nbins,
										   int    fitstart,
                                           int    nphotons,
                                           float  interval)
{ // calc weighted average time of arrival
    int   bin, nphotonsbin;
    float avg;

    for (bin=fitstart, avg=0.0; bin<nbins; bin++)
    {
        nphotonsbin = data[bin];

        if (nphotonsbin)
            avg += (0.5f+(float)bin)*nphotonsbin;
    }

    return (interval*avg/(float)nbins/(float)nphotons);  // weighted average bin = avg/nPhotons, x by interval/nbins to get average time
}

/*=================================================================================*/
/* Function:  bayes_dm_CheckDoubleValueValid                                       */
/*---------------------------------------------------------------------------------*/
/* Purpose:   Checks whether a 'double variable' holds a valid value,              */
/*            i.e. determines whether underflow / underflow has occurred.          */
/*---------------------------------------------------------------------------------*/
/* Ref:       None.                                                                */
/*---------------------------------------------------------------------------------*/
/* Notes:     None.                                                                */
/*---------------------------------------------------------------------------------*/
/* Arguments: Value: the value to be checked for data integrity.                   */
/*            Type:  location to write the type of value, e.g. underflow.          */
/*---------------------------------------------------------------------------------*/
/* Return:    BAYES_DM_DOUBLE_TYPE_VALID if 'Value' is determined to be valid.     */
/*            BAYES_DM_DOUBLE_TYPE_INVALID if 'Value' is determined to be invalid. */
/*=================================================================================*/

int bayes_dm_CheckDoubleValueValid(double Value, int *Type)
{
	*Type = dtype(Value);

	if ((BAYES_DM_DTYPE_REGULAR_NUMBER == *Type) || (BAYES_DM_DTYPE_ZERO == *Type))
		return (BAYES_DM_DOUBLE_TYPE_VALID);
	else
		return (BAYES_DM_DOUBLE_TYPE_INVALID);
}


/*=================================================================================*/
/* Function:  bayes_dm_CorrectInvalidDoubleValue                                   */
/*---------------------------------------------------------------------------------*/
/* Purpose:   To overwrite a 'double variable' that has been determined not to     */
/*            hold a valid value, with a valid value.                              */
/*---------------------------------------------------------------------------------*/
/* Ref:       None.                                                                */
/*---------------------------------------------------------------------------------*/
/* Notes:     None.                                                                */
/*---------------------------------------------------------------------------------*/
/* Arguments: Value: location of the value to be overwritten.                      */
/*            Type:  the type of 'disorder' - e.g. underflow / overflow.           */
/*---------------------------------------------------------------------------------*/
/* Return:    None.                                                                */
/*=================================================================================*/

void bayes_dm_CorrectInvalidDoubleValue(double *Value, int Type)
{
	switch (Type)
	{
		case BAYES_DM_DTYPE_UNDERFLOW:
		{
//		    *Value = TINY;
            *Value = BAYES_SIZE_DOUBLE_TINY;
		//	printf("\nBAYES_DM_DTYPE_UNDERFLOW\n");
			break;
		}

		case BAYES_DM_DTYPE_NEGATIVE_INFINITY:
		{
		    *Value = MINUSINFTY;
//			printf("\nBAYES_DM_DTYPE_NEGATIVE_INFINITY\n");
			break;
		}

		case BAYES_DM_DTYPE_POSITIVE_INFINITY:
		case BAYES_DM_DTYPE_NOT_A_NUMBER:
		{
//		    *Value = BIG;
            *Value = BAYES_SIZE_DOUBLE_HUGE;

//			if (Type == BAYES_DM_DTYPE_POSITIVE_INFINITY)
//				printf("\nBAYES_DM_DTYPE_POSITIVE_INFINITY\n");
//			else
	//			printf("\nBAYES_DM_DTYPE_NOT_A_NUMBER\n");

			break;
		}

		default:
		{
		//unexpected type error message
		}
	}
}




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

int bayes_dm_CheckDoubleValueValid(double Value, int *Type);

void bayes_dm_CorrectInvalidDoubleValue(double *Value, int Type);

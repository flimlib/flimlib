/*=================================================================================*/
/* File:       bayes_Sizes.h                                                       */
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
/* 071222 | Creation, mrowley.                                                     */
/*---------------------------------------------------------------------------------*/
/*        |                                                                        */
/*=================================================================================*/

#define  TINY       1.0e-25
#define  BIG        1.0e25
#define  MINUSINFTY -1.0e25

#define BAYES_SIZE_DOUBLE_HUGE 1.0e+300
#define BAYES_SIZE_DOUBLE_TINY 1.0e-300


#define ROOTPI		1.7724538509055160272
#define PI          ROOTPI*ROOTPI

#define BAYES_INSTR_FWHM_TO_SIGMA_FACTOR 2.35482 /* sigma = FWHM/2.35482 */

#define  ROOTTWO            1.4142135623730950488
#define  ONE_OVER_ROOTTOW   1.0/ROOTTWO
#define  SMALLTIME          1.0e-2

#define BAYES_ACCURACY_2_5_PERCENT      0.025
#define BAYES_ACCURACY_5_PERCENT        0.050
#define BAYES_ACCURACY_7_5_PERCENT      0.075
#define BAYES_ACCURACY_10_PERCENT       0.100
#define BAYES_ACCURACY_15_PERCENT       0.150
#define BAYES_ACCURACY_20_PERCENT       0.200



// BAYES Error codes
#define BAYES_ERR_NO_ERROR                         0
#define BAYES_ERR_INVALID_DATA                    -1
#define BAYES_ERR_INVALID_WINDOW                  -2
#define BAYES_ERR_INVALID_MODEL                   -3
#define BAYES_ERR_FUNCTIONALITY_NOT_SUPPORTED     -4
#define BAYES_ERR_INVALID_FIXED_PARAM_VALUE       -5
#define BAYES_ERR_ALL_BAYES_PARAM_VALUES_FIXED    -6
#define BAYES_ERR_PARAM_ESTIMATION_FAILURE        -7
#define BAYES_ERR_NO_GRID_FOR_RAPID_ESTIMATION    -8
#define BAYES_ERR_MODEL_SEL_PARAM_EST_FAILURE     -9
#define BAYES_ERR_MODEL_SEL_HESSIAN_FAILURE       -10

#define BAYES_AVE_ERRS_ROUTINE_NO_ERRORS     0
#define BAYES_AVE_ERRS_ROUTINE_OUTOFRANGE  -11
#define BAYES_AVE_ERRS_ROUTINE_MP_VALS_ONLY -12   // never called
#define BAYES_AVE_ERRS_ROUTINE_ERROR        -13
#define BAYES_AVE_ERR_RAPID_INSUFFICIENT_GRID -14

#define BAYES__RESULT_USER_CANCEL         -99




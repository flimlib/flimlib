#ifndef __DTYPE__
#define __DTYPE__

#ifdef __cplusplus
    extern "C" {
#endif

/*****************************************************************************/
/* The function dtype will determine the type of double. It assumes that the */
/* double is in IEEE 754 format.                                             */
/* returns	0  reg number                                                    */
/*	1  +infinity                                                             */
/*	-1 -infinity                                                             */
/*	2  NaN                                                                   */
/*	3  underflow                                                             */
/*	4  0				                                                     */
/*****************************************************************************/
int dtype(double d);

#ifdef __cplusplus
    }
#endif

#endif

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

#include "DTYPE.h"

static const unsigned char gDoubleInfinityByteArray[sizeof(double)]
= { 0x7F, 0xF0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

static const unsigned char gDoubleNanByteArray[sizeof(double)]
= { 0x7F, 0xF8, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

static const unsigned char gDoubleExponentMask[sizeof(double)]
= { 0x7F, 0xF0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

static const unsigned char gDoubleMantissaMask[sizeof(double)]
= { 0x00, 0x0F, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };

static void CheckDouble(double number, int *mantissa, int *exponent)
{
	int i;
	unsigned char currByte;

	*mantissa = 0;
	*exponent = 1;

	// NOTE: the following code only works on x86 processors because of Endianess
	for (i = 0; i < sizeof(double); ++i)
	{
		currByte = ((unsigned char*)&number)[sizeof(double) - 1 - i];
		*exponent &= (gDoubleExponentMask[i] == (currByte & gDoubleExponentMask[i]));
		*mantissa |= (currByte & gDoubleMantissaMask[i]);
	}
}


int IsNotANumber(double number)
{
	int mantissa, exponent;

	CheckDouble(number, &mantissa, &exponent);
	return (mantissa && exponent);
}

// Returns 1 for +ive infinity, -1 for -ive infinity, and 0 otherwise
int IsInfinity(double number)
{
	int mantissa, exponent;

	CheckDouble(number, &mantissa, &exponent);
	return (!mantissa && exponent) ? (number < 0.0 ? -1 : 1) : 0;
}


/*****************************************************************************/
/* The function dtype will determine the type of double.                     */
/* returns	0  reg number                                                    */
/*	1  +infinity                                                             */
/*	-1 -infinity                                                             */
/*	2  NaN                                                                   */
/*	3  underflow   NOT IMPLEMENTED                                           */
/*	4  0				                                                     */
/*****************************************************************************/

int dtype(double d)
{
	int ret = IsInfinity(d);

	if(ret != 0){
		return ret;  // will be -1 for -Inf or +1 for Inf
	}
	else if(d==0.0){
		return 4;
	}
	else if(IsNotANumber(d)){
		return 2;
	}

	return 0;   // regular number
}

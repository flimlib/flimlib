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

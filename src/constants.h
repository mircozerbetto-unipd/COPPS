/***********************************************************************************
 * C++OPPS 2.2 - Interpretation of NMR relaxation in proteins                      *
 * Copyright (C) 2008  Mirco Zerbetto                                              * 
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or any later version.                                           *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 ***********************************************************************************
 * Author: Mirco Zerbetto                                                          *
 * Dipartimento di Scienze Chimiche - Universita' di Padova - Italy                *
 * E-mail: mirco.zerbetto@unipd.it                                                 *
 ***********************************************************************************/

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define ZERO 1.0e-15
#define STOL 1.0e-6

#define SQRT_FIVE sqrt(5.0)
#define ONE_OVER_SQRT_TWO sqrt(0.5)
#define SQRT_ONE_OVER_TWO sqrt(0.5)
#define ONE_OVER_SQRT_THREE 1.0/sqrt(3.0)
#define SQRT_ONE_OVER_FOUR sqrt(0.25)
#define ONE_OVER_SQRT_SIX 1.0/sqrt(6.0)
#define SQRT_THREE_OVER_TWO sqrt(1.5)
#define ONE_OVER_THREE 1.0/3.0
#define ONE_OVER_SIX 1.0/6.0
#define TWO_OVER_FIFTHEEN 2.0/15.0
#define ONE_OVER_TWENTYFOUR 1.0/24.0

#define gyroH 2.675198e8    /* A m N^-1 s^-1 */  // Hydrogen magnetogyric ratio
#define hbar 1.054494e-34   /* N m  s        */  // Reduced Planck constant
#define MU0_OVER_4PI 1.0e-7 /* N A^-2        */  // Vacuum magnetic permittivity

#define DEG_TO_RAD M_PI/180.0
#define RAD_TO_DEG 180.0/M_PI

#define MAX(a,b) (a > b ? a : b)

#endif

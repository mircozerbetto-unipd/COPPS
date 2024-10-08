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

/*
 ============================================================================
 Name        : prep.h
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Preprocessor directives to control program behavior
 ============================================================================
 */

#ifndef PREP_H_
#define PREP_H_

#define _LINUX_
//#define _WINDOWS_

#define WRITE_PHYSICS
//#define WRITE_BASIS
//#define WRITE_STVEC
//#define WRITE_MATRIX
//#define WRITE_LANCZOS
//#define WRITE_LANCZOS_WARNINGS
#define WRITE_RELAX
#define WRITE_ORDERPAR
//#define WRITE_DESTROY_MESSAGE
#define WRITE_ALL

//#define F77_PHASE  // This flags makes C++OPPS use the same phase as used for the dipolar interaction in the f77 copps code

#endif /*PREP_H_*/

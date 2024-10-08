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

#ifndef COPPS_TYPES_H_
#define COPPS_TYPES_H_

#include <vector>
#include <complex>

// Real

typedef std::vector <int>            ivector;
typedef std::vector <ivector>        ivectors;
typedef std::vector <double>         dvector;
typedef std::vector <dvector>        dvvector;
//typedef vector <VECTOR_double>  VDvector;
typedef std::vector <long double>    ldvector;

// Complex

typedef std::complex <ivector>            civector;
typedef std::complex <double>             dcomplex;
typedef std::complex <dvector>            cdvector;
typedef std::vector  <dcomplex>           dcvector;
typedef std::vector  <dcvector>           dcvvector;
typedef std::complex <long double>        ldcomplex;
typedef std::vector  <ldcomplex>          cldvector;
//typedef complex <CompRow_Mat_double> cmatrix;

// Other

struct potential {
	int npop;
	ldvector c20;
	ldvector c22;
	ldvector c40;
	ldvector c42;
	ldvector c44;
};
typedef std::vector <potential> pvector;

struct potential_fb2 {
	int n1Max, n2Max;
	int dim1, dim2;
	dcomplex *c;
};

#endif

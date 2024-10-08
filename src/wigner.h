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

#ifndef WIGNER_H_
#define WIGNER_H_

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <complex>

#include "constants.h"

#include "types.h"

class wigner{
	public:
		
		wigner(int);
		wigner(int,int,int,long double,long double,long double);
		~wigner();
		void setQuantumNumbers(int,int,int);
		void setAngles(long double,long double,long double);
		ldcomplex getWignerMatrix(void);
		ldcomplex getWignerMatrix(int,int,int,long double,long double,long double);
		ldcomplex getWignerMatrix_d(int,int,int,double,double,double);
		long double getReducedWignerMatrix(void);
		long double getReducedWignerMatrix(int,int,int,long double);
		long double getRank2RWM(int,int,long double);	
		
	private:
		int J,M,K;
		long double alpha,beta,gamma;
		long double *lnf;
};

#endif /*WIGNER_H_*/

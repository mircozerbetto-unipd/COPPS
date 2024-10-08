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

#ifndef BASIS_H_
#define BASIS_H_

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

#include "types.h"

class basis
{
	public:
		basis();
		~basis();
		
		void init(int,int,int,int,int,int,int,int,int,int,int,int);
		void init(int,int,int,int,int,int,int,int,int,int,int,int,int);
		void checkSymmetry(ldvector);
		void checkSymmetry(ldvector,bool*);
		void buildSpace(std::string);
		
		int getNumberOfBasisFunctions(void);
		ivector getBasisFunction(int);
		
		ivector getL1bounds(void);
		ivector getM1bounds(void);
		ivector getMmbounds(void);
		ivector getjjbounds(void);
		
		std::string toString(std::string);

		bool isOmegaVFromGeometry;
				
	private:
		
		void buildSpace_SRLS(void);
		void buildSpace_FB1(void);
		void buildSpace_FB2(void);

		int rank;
		int N1,N2;
		int L1,M1,K1;
		int L2,M2,K2;
		int jj,Mm,Mp;
		int L1min, L1max, M1min, M1max, K1min, K1max;
		int L2min, L2max, M2min, M2max, K2min, K2max;
		int jjmin, jjmax;
		int Mmmin, Mmmax, Mpmin, Mpmax;
		int nbf;
		ivectors basisIdx;
		
};

#endif /*BASIS_H_*/

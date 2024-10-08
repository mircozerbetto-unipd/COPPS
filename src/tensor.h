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

#ifndef TENSOR_H_
#define TENSOR_H_

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "constants.h"
#include "wigner.h"

#include "types.h"

class tensor
{
	public:
		tensor();
		virtual ~tensor();
		void scale(void);
		void scale(long double);
		void transform(void);
		long double iso(void);
		void setCartesianComponent(std::string,long double);
		long double getCartesianComponent(std::string);
		void setCartesianComponents(ldvector);
		ldvector getCartesianComponents(void);
		void setScale(long double);
		long double getScale(void);
		cldvector getSphericalComponents(void);
		void setAngle(std::string,long double);
		long double getAngle(std::string);
		void setAngles(ldvector);
		ldvector getAngles(void);
		
	private:

		long double dxx, dyy, dzz;
		long double Dscale;
		long double alpha, beta, gamma;
		ldvector bkp;
	    ldcomplex d00,d20,d2p1,d2m1,d2p2,d2m2;
};

#endif /*TENSOR_H_*/

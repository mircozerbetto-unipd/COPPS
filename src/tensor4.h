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

/* This class handles rotational-internal diffusion tensor for flexible
   rotetors with one degree of freedom (FB1 model for the dynamics) */

#ifndef TENSOR4_H_
#define TENSOR4_H_

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <complex>

#include "constants.h"

#include "types.h"

class tensor4
{
	public:
		tensor4();
		virtual ~tensor4();

		void setScale(long double);
		long double getScale(void);
		void scale(void);
		void scale(long double);

		//void transform(void);
		long double isoRR(void);

		void setRRComponent(std::string,long double);
		long double getRRComponent(std::string);
		void setRRComponents(ldvector);
		ldvector getRRComponents(void);

		void setRIComponent(std::string,long double);
		long double getRIComponent(std::string);
		void setRIComponents(ldvector);
		ldvector getRIComponents(void);

		void setIIComponent(long double);
		long double getIIComponent(void);

		void setComponents(ldvector);
		ldvector getComponents(void);

		long double getDpRR(void);
		long double getDmRR(void);
		ldcomplex getDpRI(void);
		ldcomplex getDmRI(void);
		
	private:

		long double dxx, dyy, dzz;	// Rotational part
		long double dii;		// Internal part
		long double dxi, dyi, dzi;	// rotational-internal part
		long double dpRR, dmRR;		// linear combinations of dxx and dyy
		ldcomplex dpRI, dmRI;		// linear combinations of dxi and dyi
		long double Dscale;		// scaling factor
		ldvector bkp;			// backup of the initial values
};

#endif /*TENSOR4_H_*/

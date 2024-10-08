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
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// NOTES                                                                           //
//                                                                                 //
// This class is intended to convert between rotation matrices and Euler angles    //
// in the ZYZ notation. Also, given three frames and two rotations the class       //
// permits to calculate the third rotation.                                        //
// In the following, E1(omega1) is the rotation matrix from a fixed frame LF to    //
// the first frame, F1. E2(omega2) is the rotation matrix tranforming grom LF      //
// to F2. Finally, E12(omega12) is the rotation matrix transforming from F1 to F2. //
// In this picture,                                                                //
// E12 = E2 E1'                                                                    //
// where the apex means the transpose.                                             //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////

#ifndef EULER_H_
#define EULER_H_

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "types.h"

class euler
{

	public:
		euler();
		virtual ~euler();
	
		void setE1(dvector);
		dvector getE1();
		void setE2(dvector);
		dvector getE2();
		void setE12(dvector);
		dvector getE12();

		dvector transpose(dvector);

		void setOmega1(dvector);
		dvector getOmega1();
		void setOmega2(dvector);
		dvector getOmega2();
		void setOmega12(dvector);
		dvector getOmega12();

		void omega1ToE1();
		void E1ToOmega1();
		void omega2ToE2();
		void E2ToOmega2();
		void omega12ToE12();
		void E12ToOmega12();
	
		void calculateE1();
		void calculateE2();
		void calculateE12();
	
		dvector mul(dvector,dvector,int,int);

	private:
		dvector E1;
		dvector E2;
		dvector E12;
		dvector omega1;
		dvector omega2;
		dvector omega12;
};

#endif

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
 ==============================================================================
 Name        : euler.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class to manage rotation matrix <-> Euler angles transformations 
 ==============================================================================
 */

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
// ********** NB: Rotation matrices are the transpose of Euler matrices ********** //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////

#include "euler.h"

/***************/
/* Constructor */
/***************/
euler::euler()
{
	E1  = dvector(9,0.0);
	E2  = dvector(9,0.0);
	E12 = dvector(9,0.0);
	omega1  = dvector(3,0.0);
	omega2  = dvector(3,0.0);
	omega12 = dvector(3,0.0);
}
/*****************/
/* Deconstructor */
/*****************/
euler::~euler()
{
}
/*****************************/
/* Set the rotation matrices */
/*****************************/
void euler::setE1(dvector e)
{
	E1 = e;
	return;
}
 void euler::setE2(dvector e)
{
	E2 = e;
	return;
}
void euler::setE12(dvector e)
{
	E12 = e;
	return;
}
/*****************************/
/* Get the rotation matrices */
/*****************************/
dvector euler::getE1()
{
	return E1;
}
dvector euler::getE2()
{
	return E2;
}
dvector euler::getE12()
{
	return E12;
}
/********************/
/* Set Euler angles */
/********************/
void euler::setOmega1(dvector o)
{
	omega1 = o;
}
void euler::setOmega2(dvector o)
{
	omega2 = o;
}
void euler::setOmega12(dvector o)
{
	omega12 = o;
}
/********************/
/* Get Euler angles */
/********************/
dvector euler::getOmega1()
{
	return omega1;
}
dvector euler::getOmega2()
{
	return omega2;
}
dvector euler::getOmega12()
{
	return omega12;
}
/***********************************************/
/* Calculate Euler angles from rotation matrix */
/***********************************************/
void euler::omega1ToE1()
{
	double ca = cos(omega1[0]), sa = sin(omega1[0]);
	double cb = cos(omega1[1]), sb = sin(omega1[1]);
	double cg = cos(omega1[2]), sg = sin(omega1[2]);
	E1[0] =  ca*cb*cg - sa*sg;
	E1[1] = -ca*cb*sg - sa*cg;	
	E1[2] =  ca*sb;
	E1[3] =  sa*cb*cg + ca*sg;
	E1[4] = -sa*cb*sg + ca*cg;
	E1[5] =  sa*sb;
	E1[6] = -sb*cg;
	E1[7] =  sb*sg;
	E1[8] =  cb;
	return;
}
void euler::omega2ToE2()
{
	double ca = cos(omega2[0]), sa = sin(omega2[0]);
	double cb = cos(omega2[1]), sb = sin(omega2[1]);
	double cg = cos(omega2[2]), sg = sin(omega2[2]);
        E2[0] =  ca*cb*cg - sa*sg;
        E2[1] = -ca*cb*sg - sa*cg;	
        E2[2] =  ca*sb;
        E2[3] =  sa*cb*cg + ca*sg;
        E2[4] = -sa*cb*sg + ca*cg;
        E2[5] =  sa*sb;
        E2[6] = -sb*cg;
        E2[7] =  sb*sg;
        E2[8] =  cb;
	return;
}
void euler::omega12ToE12()
{
	double ca = cos(omega12[0]), sa = sin(omega12[0]);
	double cb = cos(omega12[1]), sb = sin(omega12[1]);
	double cg = cos(omega12[2]), sg = sin(omega12[2]);
        E12[0] =  ca*cb*cg - sa*sg;
        E12[1] = -ca*cb*sg - sa*cg;	
        E12[2] =  ca*sb;
        E12[3] =  sa*cb*cg + ca*sg;
        E12[4] = -sa*cb*sg + ca*cg;
        E12[5] =  sa*sb;
        E12[6] = -sb*cg;
        E12[7] =  sb*sg;
        E12[8] =  cb;
	return;
}
/***********************************************/
/* Calculate rotation matrix from Euler angles */
/***********************************************/
void euler::E1ToOmega1()
{
	double sb = sqrt(1.0 - E1[8]*E1[8]);
	omega1[0] = atan2(E1[5],E1[2]);
	omega1[1] = atan2(sb,E1[8]);
	omega1[2] = atan2(E1[7],-E1[6]);
}
void euler::E2ToOmega2()
{
	double sb = sqrt(1.0 - E2[8]*E2[8]);
	omega2[0] = atan2(E2[5],E2[2]);
	omega2[1] = atan2(sb,E2[8]);
	omega2[2] = atan2(E2[7],-E2[6]);
}
void euler::E12ToOmega12()
{
	double sb = sqrt(1.0 - E12[8]*E12[8]);
	omega12[0] = atan2(E12[5],E12[2]);
	omega12[1] = atan2(sb,E12[8]);
	omega12[2] = atan2(E12[7],-E12[6]);
}
/***********************************************/
/* Obtain the third matrix given the other two */
/***********************************************/
void euler::calculateE1()
{
	E1 = mul(E12,E2,1,0);
	return;
}
void euler::calculateE2()
{
	E2 = mul(E12,E1,0,0);
	return;
}
void euler::calculateE12()
{
	E12 = mul(E2,E1,0,1);
	return;
}
/********************************/
/* Matrix matrix multiplication */
/********************************/
dvector euler::mul(dvector m1, dvector m2, int transpose1, int transpose2)
{
	dvector m3 = dvector(9,0.0);
	if (transpose1) m1 = transpose(m1);
	if (transpose2) m2 = transpose(m2);
	m3[0] = m1[0]*m2[0] + m1[1]*m2[3] + m1[2]*m2[6];
	m3[1] = m1[0]*m2[1] + m1[1]*m2[4] + m1[2]*m2[7];
	m3[2] = m1[0]*m2[2] + m1[1]*m2[5] + m1[2]*m2[8];
	m3[3] = m1[3]*m2[0] + m1[4]*m2[3] + m1[5]*m2[6];
	m3[4] = m1[3]*m2[1] + m1[4]*m2[4] + m1[5]*m2[7];
	m3[5] = m1[3]*m2[2] + m1[4]*m2[5] + m1[5]*m2[8];
	m3[6] = m1[6]*m2[0] + m1[7]*m2[3] + m1[8]*m2[6];
	m3[7] = m1[6]*m2[1] + m1[7]*m2[4] + m1[8]*m2[7];
	m3[8] = m1[6]*m2[2] + m1[7]*m2[5] + m1[8]*m2[8];
	return m3;
}
/********************/
/* Transpose matrix */
/********************/
dvector euler::transpose(dvector m)
{
	dvector mt = dvector(9,0.0);
	mt[0] = m[0];
	mt[1] = m[3];
	mt[2] = m[6];
	mt[3] = m[1];
	mt[4] = m[4];
	mt[5] = m[7];
	mt[6] = m[2];
	mt[7] = m[5];
	mt[8] = m[8];
	return mt;
}

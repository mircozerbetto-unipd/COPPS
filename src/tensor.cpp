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
 Name        : tensor.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class with methods for rank 2 tensors
 ============================================================================
 */
#include <cstdlib>

#include "tensor.h"

/***************
 * Constructor *
 ***************/

tensor::tensor()
{
	/**************************************************************************
	 * bkp vector stores the original values while dxx, dyy, dzz and the      *
	 * spherical components are all scaled by the scale factor.               *
	 * Once a new scaling is required, all cartesian and spherical components *
	 * are obtained by scaling the original values stored in bkp.             *
	 **************************************************************************/
	bkp = ldvector(3,0.0);
}

/**************
 * Destructor *
 **************/

tensor::~tensor()
{
}

/************************************************************************************************
 * Divide by a constant original cartesian components in the frame that diagonalizes the tensor *
 ************************************************************************************************/

void tensor::scale (void)
{
	dxx = bkp.at(0) / Dscale;
	dyy = bkp.at(1) / Dscale;
	dzz = bkp.at(2) / Dscale;
}

void tensor::scale (long double s)
{
	Dscale = s;
	dxx = bkp.at(0) / Dscale;
	dyy = bkp.at(1) / Dscale;
	dzz = bkp.at(2) / Dscale;
}

/****************************************************************************
 * Calculate spherical components in the frame where the tensor is diagonal *
 ****************************************************************************/

// NOTE: the euler angles OmegaO transform from OF to M2F. What we are calculating
//       is (OF)D given (M2F)D so, the transformation is
//      
//       (OF)D*(l,m) = Sum_m' D[l,m,m'](OmegaO) (M2F)D*(l,m')
//

void tensor::transform(void)
{
	wigner rot(10);
	ldcomplex mfd20, mfd22;
	mfd20 = std::complex<long double>(ONE_OVER_SQRT_SIX*(2.0*dzz-dxx-dyy),0.0);
	mfd22 = std::complex<long double>(0.5*(dxx-dyy),0.0);

	/* Calculate the complex conjugates of the spherical components */

	d00  = std::complex<long double>(-ONE_OVER_SQRT_THREE*(dxx+dyy+dzz),0.0);
	d20  = mfd20*(rot.getWignerMatrix(2, 0, 0,alpha,beta,gamma)) + mfd22*(rot.getWignerMatrix(2, 0, 2,alpha,beta,gamma) + rot.getWignerMatrix(2, 0,-2,alpha,beta,gamma));
	d2p1 = mfd20*(rot.getWignerMatrix(2, 1, 0,alpha,beta,gamma)) + mfd22*(rot.getWignerMatrix(2, 1, 2,alpha,beta,gamma) + rot.getWignerMatrix(2, 1,-2,alpha,beta,gamma));
	d2m1 = mfd20*(rot.getWignerMatrix(2,-1, 0,alpha,beta,gamma)) + mfd22*(rot.getWignerMatrix(2,-1, 2,alpha,beta,gamma) + rot.getWignerMatrix(2,-1,-2,alpha,beta,gamma));
	d2p2 = mfd20*(rot.getWignerMatrix(2, 2, 0,alpha,beta,gamma)) + mfd22*(rot.getWignerMatrix(2, 2, 2,alpha,beta,gamma) + rot.getWignerMatrix(2, 2,-2,alpha,beta,gamma));
	d2m2 = mfd20*(rot.getWignerMatrix(2,-2, 0,alpha,beta,gamma)) + mfd22*(rot.getWignerMatrix(2,-2, 2,alpha,beta,gamma) + rot.getWignerMatrix(2,-2,-2,alpha,beta,gamma));

	/* Calculate the spherical components */

	d20  = conj(d20 );
	d2p1 = conj(d2p1);
	d2m1 = conj(d2m1);
	d2p2 = conj(d2p2);
	d2m2 = conj(d2m2);

	return;				
}

/*******************************************
 * Return the isotropic part of the tensor *
 *******************************************/

long double tensor::iso (void)
{
	return ONE_OVER_THREE*(dxx+dyy+dzz); 
}

/***********************************************************************************
 * Set the cartesian components of the tensor in the frame in which it is diagonal *
 ***********************************************************************************/

void tensor::setCartesianComponent(std::string comp,long double val)
{
	if (!(comp.compare("dxx")))
		dxx = bkp[0] = val;
	else if (!(comp.compare("dyy")))
		dyy = bkp[1] = val;
	else if (!(comp.compare("dzz")))
		dzz = bkp[2] = val;
	else
		return;
	return;
}

void tensor::setCartesianComponents(ldvector comp)
{
	if (comp.size() != 3){
		std::cout << "\n\nCOOPS ERROR : cartesian components vector must be of size 3 in tensor::setCartesianComponents\n\n";
		exit(0);
	}
	dxx = bkp[0] = comp[0];
	dyy = bkp[1] = comp[1];
	dzz = bkp[2] = comp[2];
}

/****************************************************************************
 * Return the cartesan components in the frame where the tensor is diagonal *
 ****************************************************************************/

ldvector tensor::getCartesianComponents(void){
	ldvector comp;
	comp.push_back(dxx);
	comp.push_back(dyy);
	comp.push_back(dzz);
	return comp;
}

long double tensor::getCartesianComponent(std::string comp)
{
	if (!(comp.compare("dxx")))
		return dxx;
	else if (!(comp.compare("dyy")))
		return dyy;
	else if (!(comp.compare("dzz")))
		return dzz;
	else
		return -1.0; // Wrong component
}

/***********************************************
 * Set the scaling factor, but do not apply it *
 ***********************************************/

void tensor::setScale(long double s){
	Dscale = s;
}

/*****************************
 * Return the scaling factor *
 *****************************/

long double tensor::getScale(void){
	return Dscale;
}

/****************************************************** 
 * Return the spherical components ordered such that: *
 * d00 = comp[0]                                      *
 * d2k = comp[2+k+1]                                  *
 ******************************************************/

cldvector tensor::getSphericalComponents(void){
	cldvector comp;
	comp.push_back(d00);
	comp.push_back(d2m2);
	comp.push_back(d2m1);
	comp.push_back(d20);
	comp.push_back(d2p1);
	comp.push_back(d2p2);
	return comp;
}

/***********************************************************************************************
 * Set Euler angles that transform from the frame where the tensor is diagonal to target frame *
 ***********************************************************************************************/

void tensor::setAngle(std::string angle ,long double val)
{
	if (!(angle.compare("alpha"))){
		alpha = val;
	}
	else if (!(angle.compare("beta")))
		beta  = val;
	else if (!(angle.compare("gamma")))
		gamma = val;
	else
		return;
	return;
}

void tensor::setAngles(ldvector angles)
{
	if (angles.size() != 3){
		std::cout << "\n\nCOOPS ERROR : Euler angles vector must be of size 3 in tensor::setAngles\n\n";
		exit(0);
	}
	alpha = angles[0];
	beta  = angles[1];
	gamma = angles[2];
}

/********************************************************
 * Return the specified Euler angle in RAD.             *
 * If the angle is not recognized return 180.0 as error *
 ********************************************************/

long double tensor::getAngle(std::string angle)
{
	if (!(angle.compare("alpha")))
		return alpha;
	else if (!(angle.compare("beta")))
		return beta;
	else if (!(angle.compare("gamma")))
		return gamma;
	else
		return 180.0;
}

ldvector tensor::getAngles(void)
{
	ldvector angles;
	angles.push_back(alpha);
	angles.push_back(beta);
	angles.push_back(gamma);
	return angles;	
}

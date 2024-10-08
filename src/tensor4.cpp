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
 Name        : tensor4.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class with methods for diffusion tensor of FB1
 ============================================================================
 */
#include "tensor4.h"

/***************
 * Constructor *
 ***************/

tensor4::tensor4()
{
	/*******************************************
	 * bkp vector stores the original values.  *
	 *******************************************/
	bkp = ldvector(7,0.0);
}

/**************
 * Destructor *
 **************/

tensor4::~tensor4()
{
}

/***********************************************
 * Set the scaling factor, but do not apply it *
 ***********************************************/

void tensor4::setScale(long double s){
	Dscale = s;
}

/*****************************
 * Return the scaling factor *
 *****************************/

long double tensor4::getScale(void){
	return Dscale;
}

/************************************************************************************************
 * Divide by a constant original cartesian components in the frame that diagonalizes the tensor *
 ************************************************************************************************/

void tensor4::scale (void)
{
	dxx = bkp.at(0) / Dscale;
	dyy = bkp.at(1) / Dscale;
	dzz = bkp.at(2) / Dscale;
	dii = bkp.at(3) / Dscale;
	dxi = bkp.at(4) / Dscale;
	dyi = bkp.at(5) / Dscale;
	dzi = bkp.at(6) / Dscale;
}

void tensor4::scale (long double s)
{
	Dscale = s;
	dxx = bkp.at(0) / Dscale;
	dyy = bkp.at(1) / Dscale;
	dzz = bkp.at(2) / Dscale;
	dii = bkp.at(3) / Dscale;
	dxi = bkp.at(4) / Dscale;
	dyi = bkp.at(5) / Dscale;
	dzi = bkp.at(6) / Dscale;
}

/*************************************************
 * Return the isotropic value of rotational part *
 *************************************************/

long double tensor4::isoRR(void)
{
	return (dxx+dyy+dzz)*ONE_OVER_THREE;
}

/*****************************************************************
 * Set components of the rotational part of the diffusion tensor *
 *****************************************************************/

void tensor4::setRRComponent(std::string comp,long double val)
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

void tensor4::setRRComponents(ldvector comp)
{
	if (comp.size() != 3){
		std::cout << "\n\nCOOPS ERROR : vector containing rotational part of diffusion tensor must be of size 3 in tensor4::setRRComponents\n\n";
		exit(0);
	}
	dxx = bkp[0] = comp[0];
	dyy = bkp[1] = comp[1];
	dzz = bkp[2] = comp[2];
}

/******************************************************
 * Return the rotational part of the diffusion tensor *
 ******************************************************/

long double tensor4::getRRComponent(std::string comp)
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

ldvector tensor4::getRRComponents(void){
	ldvector comp;
	comp.push_back(dxx);
	comp.push_back(dyy);
	comp.push_back(dzz);
	return comp;
}

/**************************************************************************
 * Set components of the rotational-internal part of the diffusion tensor *
 **************************************************************************/

void tensor4::setRIComponent(std::string comp,long double val)
{
	if (!(comp.compare("dxi")))
		dxi = bkp[4] = val;
	else if (!(comp.compare("dyi")))
		dyi = bkp[5] = val;
	else if (!(comp.compare("dzi")))
		dzi = bkp[6] = val;
	else
		return;
	return;
}

void tensor4::setRIComponents(ldvector comp)
{
	if (comp.size() != 3){
		std::cout << "\n\nCOOPS ERROR : vector containing rotational-internal part of diffusion tensor must be of size 3 in tensor4::setRIComponents\n\n";
		exit(0);
	}
	dxi = bkp[4] = comp[0];
	dyi = bkp[5] = comp[1];
	dzi = bkp[6] = comp[2];
}

/***************************************************************
 * Return the rotational-internal part of the diffusion tensor *
 ***************************************************************/

long double tensor4::getRIComponent(std::string comp)
{
	if (!(comp.compare("dxi")))
		return dxi;
	else if (!(comp.compare("dyi")))
		return dyi;
	else if (!(comp.compare("dzi")))
		return dzi;
	else
		return -1.0; // Wrong component
}

ldvector tensor4::getRIComponents(void){
	ldvector comp;
	comp.push_back(dxi);
	comp.push_back(dyi);
	comp.push_back(dzi);
	return comp;
}

/***************************************************************
 * Set components of the internal part of the diffusion tensor *
 ***************************************************************/

void tensor4::setIIComponent(long double val)
{
	dii = bkp[3] = val;
	return;
}

/******************************************************************
 * Return components of the internal part of the diffusion tensor *
 ******************************************************************/
         
 long double tensor4::getIIComponent(void)
{
	return dii;
}

/******************************************
 * Set components of the diffusion tensor *
 ******************************************/

void tensor4::setComponents(ldvector comp)
// comp = [dxx dyy dzz dii dxi dyi dzi]
{
	if (comp.size() != 7){
		std::cout << "\n\nCOOPS ERROR : vector containing all components diffusion tensor must be of size 7 in tensor4::setComponents. The order is the folloing [dxx dyy dzz dii dxi dyi dzi]\n\n";
		exit(0);
	}
	dxx = bkp[0] = comp[0];
	dyy = bkp[1] = comp[1];
	dzz = bkp[2] = comp[2];
	dii = bkp[3] = comp[3];
	dxi = bkp[4] = comp[4];
	dyi = bkp[5] = comp[5];
	dzi = bkp[6] = comp[6];
}

/*********************************************
 * Return components of the diffusion tensor *
 *********************************************/

ldvector tensor4::getComponents(void)
{
	ldvector comp; // comp = [dxx dyy dzz dii dxi dyi dzi]
	comp.push_back(dxx);
	comp.push_back(dyy);
	comp.push_back(dzz);
	comp.push_back(dii);
	comp.push_back(dxi);
	comp.push_back(dyi);
	comp.push_back(dzi);
	return comp;
}

/********************************************************
 * Return some useful linear combinations of components *
 ********************************************************/

long double tensor4::getDpRR(void)
{
	return (dxx+dyy);
}

long double tensor4::getDmRR(void)
{
	return (dxx-dyy);
}

ldcomplex tensor4::getDpRI(void)
{
	ldcomplex dp(0.5 * dxi, -0.5 * dyi);
	return dp;
}

ldcomplex tensor4::getDmRI(void)
{
	ldcomplex dm(0.5 * dxi, 0.5 * dyi);
	return dm;
}


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
 Name        : tensor5.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class with methods for diffusion tensor of FB2
 ============================================================================
 */

#include "tensor5.h"

/***************
 * Constructor *
 ***************/

tensor5::tensor5()
{
	/*******************************************
	 * bkp vector stores the original values.  *
	 *******************************************/
	bkp = ldvector(12,0.0);
}

/**************
 * Destructor *
 **************/

tensor5::~tensor5()
{
}

/***********************************************
 * Set the scaling factor, but do not apply it *
 ***********************************************/

void tensor5::setScale(long double s){
	Dscale = s;
}

/*****************************
 * Return the scaling factor *
 *****************************/

long double tensor5::getScale(void){
	return Dscale;
}

/************************************************************************************************
 * Divide by a constant original cartesian components in the frame that diagonalizes the tensor *
 ************************************************************************************************/

void tensor5::scale (void)
{
	dxx = bkp.at(0) / Dscale; dyy = bkp.at(1)  / Dscale; dzz = bkp.at(2)  / Dscale;
	d11 = bkp.at(3) / Dscale; d22 = bkp.at(4)  / Dscale; d12 = bkp.at(5)  / Dscale;
	dx1 = bkp.at(6) / Dscale; dy1 = bkp.at(7)  / Dscale; dz1 = bkp.at(8)  / Dscale;
	dx2 = bkp.at(9) / Dscale; dy2 = bkp.at(10) / Dscale; dz2 = bkp.at(11) / Dscale;
}

void tensor5::scale (long double s)
{
	Dscale = s;
	dxx = bkp.at(0) / Dscale; dyy = bkp.at(1)  / Dscale; dzz = bkp.at(2)  / Dscale;
	d11 = bkp.at(3) / Dscale; d22 = bkp.at(4)  / Dscale; d12 = bkp.at(5)  / Dscale;
	dx1 = bkp.at(6) / Dscale; dy1 = bkp.at(7)  / Dscale; dz1 = bkp.at(8)  / Dscale;
	dx2 = bkp.at(9) / Dscale; dy2 = bkp.at(10) / Dscale; dz2 = bkp.at(11) / Dscale;
}

/*************************************************
 * Return the isotropic value of rotational part *
 *************************************************/

long double tensor5::isoRR(void)
{
	return (dxx+dyy+dzz)*ONE_OVER_THREE;
}

/*****************************************************************
 * Set components of the rotational part of the diffusion tensor *
 *****************************************************************/

void tensor5::setRRComponent(std::string comp,long double val)
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

void tensor5::setRRComponents(ldvector comp)
{
	if (comp.size() != 3){
		std::cout << "\n\nCOOPS ERROR : vector containing rotational part of diffusion tensor must be of size 3 in tensor5::setRRComponents\n\n";
		exit(0);
	}
	dxx = bkp[0] = comp[0];
	dyy = bkp[1] = comp[1];
	dzz = bkp[2] = comp[2];
}

/******************************************************
 * Return the rotational part of the diffusion tensor *
 ******************************************************/

long double tensor5::getRRComponent(std::string comp)
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

ldvector tensor5::getRRComponents(void){
	ldvector comp;
	comp.push_back(dxx);
	comp.push_back(dyy);
	comp.push_back(dzz);
	return comp;
}

/**************************************************************************
 * Set components of the rotational-internal part of the diffusion tensor *
 **************************************************************************/

void tensor5::setRIComponent(std::string comp,long double val)
{
	if (!(comp.compare("dx1")))
		dx1 = bkp[6] = val;
	else if (!(comp.compare("dy1")))
		dy1 = bkp[7] = val;
	else if (!(comp.compare("dz1")))
		dz1 = bkp[8] = val;
	else if (!(comp.compare("dx2")))
		dx2 = bkp[9] = val;
	else if (!(comp.compare("dy2")))
		dy2 = bkp[10] = val;
	else if (!(comp.compare("dz2")))
		dz2 = bkp[11] = val;
	else
		return;
	return;
}

void tensor5::setRIComponents(ldvector comp)
{
	if (comp.size() != 6){
		std::cout << "\n\nCOOPS ERROR : vector containing rotational-internal part of diffusion tensor must be of size 6 in tensor5::setRIComponents\n\n";
		exit(0);
	}
	dx1 = bkp[6]  = comp[0];
	dy1 = bkp[7]  = comp[1];
	dz1 = bkp[8]  = comp[2];
	dx2 = bkp[9]  = comp[3];
	dy2 = bkp[10] = comp[4];
	dz2 = bkp[11] = comp[5];
}

/***************************************************************
 * Return the rotational-internal part of the diffusion tensor *
 ***************************************************************/

long double tensor5::getRIComponent(std::string comp)
{
	if (!(comp.compare("dx1")))
		return dx1;
	else if (!(comp.compare("dy1")))
		return dy1;
	else if (!(comp.compare("dz1")))
		return dz1;
	else if (!(comp.compare("dx2")))
		return dx2;
	else if (!(comp.compare("dy2")))
		return dy2;
	else if (!(comp.compare("dz2")))
		return dz2;
	else
		return -1.0; // Wrong component
}

ldvector tensor5::getRIComponents(void){
	ldvector comp;
	comp.push_back(dx1);
	comp.push_back(dy1);
	comp.push_back(dz1);
	comp.push_back(dx2);
	comp.push_back(dy2);
	comp.push_back(dz2);
	return comp;
}

/***************************************************************
 * Set components of the internal part of the diffusion tensor *
 ***************************************************************/

void tensor5::setIIComponent(std::string comp, long double val)
{
	if (!(comp.compare("d11")))
		d11 = bkp[3] = val;
	else if (!(comp.compare("d22")))
		d22 = bkp[4] = val;
	else if (!(comp.compare("d12")) || !(comp.compare("d21")))
		d12 = bkp[5] = val;
	return;
}

void tensor5::setIIComponents(ldvector comp)
{
	if (comp.size() != 3){
		std::cout << "\n\nCOOPS ERROR : vector containing internal part of diffusion tensor must be of size 3 in tensor5::setIIComponents\n\n";
		exit(1);
	}
	d11 = bkp[3] = comp[0];
	d22 = bkp[4] = comp[1];
	d12 = bkp[5] = comp[2];
}

/******************************************************************
 * Return components of the internal part of the diffusion tensor *
 ******************************************************************/
         
long double tensor5::getIIComponent(std::string comp)
{
	if (!(comp.compare("d11")))
		return d11;
	else if (!(comp.compare("d22")))
		return d22;
	else if (!(comp.compare("d12")) || !(comp.compare("d21")))
		return d12;
	else
		return -1.0; // Wrong component
}

ldvector tensor5::getIIComponents(void)
{
	ldvector comp;
	comp.push_back(d11);
	comp.push_back(d22);
	comp.push_back(d12);
	return comp;
}

/******************************************
 * Set components of the diffusion tensor *
 ******************************************/

void tensor5::setComponents(ldvector comp)
// comp = [dxx dyy dzz d11 d22 d12 dx1 dy1 dz1 dx2 dy2 dz2]
{
	if (comp.size() != 12){
		std::cout << "\n\nCOOPS ERROR : vector containing all components diffusion tensor must be of size 12 in tensor5::setComponents. \n              The order is the folloing [dxx dyy dzz d11 d22 d12 dx1 dy1 dz1 dx2 dy2 dz2]\n\n";
		exit(0);
	}
	dxx = bkp[0]  = comp[0];
	dyy = bkp[1]  = comp[1];
	dzz = bkp[2]  = comp[2];
	d11 = bkp[3]  = comp[3];
	d22 = bkp[4]  = comp[4];
	d12 = bkp[5]  = comp[5];
	dx1 = bkp[6]  = comp[6];
	dy1 = bkp[7]  = comp[7];
	dz1 = bkp[8]  = comp[8];
	dx2 = bkp[9]  = comp[9];
	dy2 = bkp[10] = comp[10];
	dz2 = bkp[11] = comp[11];
}

/*********************************************
 * Return components of the diffusion tensor *
 *********************************************/

ldvector tensor5::getComponents(void)
{
	ldvector comp; // comp = [dxx dyy dzz d11 d22 d12 dx1 dy1 dz1 dx2 dy2 dz2]
	comp.push_back(dxx);
	comp.push_back(dyy);
	comp.push_back(dzz);
	comp.push_back(d11);
	comp.push_back(d22);
	comp.push_back(d12);
	comp.push_back(dx1);
	comp.push_back(dy1);
	comp.push_back(dz1);
	comp.push_back(dx2);
	comp.push_back(dy2);
	comp.push_back(dz2);
	return comp;
}

/********************************************************
 * Return some useful linear combinations of components *
 ********************************************************/

long double tensor5::getDpRR(void)
{
	return (dxx+dyy);
}

long double tensor5::getDmRR(void)
{
	return (dxx-dyy);
}

ldcomplex tensor5::getDpRI(int t)
{
	ldcomplex dp(0.0,0.0);
	if (t == 1)
		dp = ldcomplex(0.5*dx1, -0.5*dy1);
	else if (t == 2)
		dp = ldcomplex(0.5*dx2, -0.5*dy2);
	return dp;
}

ldcomplex tensor5::getDmRI(int t)
{
	ldcomplex dm(0.0,0.0);
	if (t == 1)
		dm = ldcomplex(0.5*dx1, 0.5*dy1);
	else if (t == 2)
		dm = ldcomplex(0.5*dx2, 0.5*dy2);
	return dm;
}


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
 Name        : experimental.h
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Header of experimental class
 ============================================================================
 */
#ifndef EXPERIMENTAL_H_
#define EXPERIMENTAL_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "constants.h"

#include "types.h"

#ifdef _MPI_
#include "mpi.h"
#endif

class experimental
{
	public:
		
		experimental();
		virtual ~experimental();
		void init(std::string);
		void init(std::string,int);
		void setFileName(std::string);
		bool readData(bool);
		double * getExpDataOfComp(int*);
		int getNPointsOfComp(int);
		int getNComp(void);
		std::string toString(void);
		dvector getFieldOfComp(int);
		double getAlphaV(int);
		double getBetaV(int);
		double getGammaV(int);
		int compHasPotential(int);
		potential getPotentialOfComp(int);

		bool isOmegaVFromGeometry;

		int getActiveComponent(void);
		void setActiveComponent(int);

		int getNExpData(void);
		bool hasExpValueOf(std::string);
		
	private:
		
		bool haveFile;
		bool hasT1, hasT2, hasNOE, hasCCRRT1, hasCCRRT2;
		int rank;
		int ncomp;
		int nExpData;
		int activeComponent;
		ivector nfield;
		ivector nexp;
		ivector foundPotential;
		dvector alphav, betav, gammav;
		dvvector field;
		dcvvector t1;
		dcvvector t2;
		dcvvector noe;
		dcvvector ccrrt1;
		dcvvector ccrrt2;
		pvector potentialCoeff;
		std::string fileName;
		
};

#endif /*EXPERIMENTAL_H_*/

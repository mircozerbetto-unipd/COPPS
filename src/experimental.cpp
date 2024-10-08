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
 Name        : experimental.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class to read experimental data
 ============================================================================
 */

#include <cstdlib>

#include "experimental.h"

/***************
 * Constructor *
 ***************/

experimental::experimental()
{
	haveFile = false;
}

/**************
 * Destructor *
 **************/

experimental::~experimental()
{
}

/**************
 * Initiators *
 **************/

void experimental::init(std::string f)
{
	fileName = f;
	rank = 0;
	haveFile = true;
	activeComponent = 1;
	return;
}

void experimental::init(std::string f, int r)
{
	init(f);
	rank = r;
	return;
}

/***************************************************** 
 * Set the name of file containing experimental data *
 *****************************************************/

void experimental::setFileName(std::string f)
{
	fileName = f;
	haveFile = true;
	return;
}

/**************************
 * Read experimental data *
 **************************/

bool experimental::readData(bool isFitting)
{
	if (!haveFile)
	{
		std::cout << std::endl << "ERROR in experimental::readData(void). Experimental data file was not set." << std::endl << std::endl;
		exit(0);
	}
	
	/*********************
	 * Erase all vectors *
	 *********************/
	
	nexp.clear();
	nfield.clear();
	for (unsigned int i = 0; i < field.size(); i++)
		field[i].clear();
	field.clear();
	for (unsigned int i = 0; i < t1.size(); i++)
		t1[i].clear();
	t1.clear();
	for (unsigned int i = 0; i < t2.size(); i++)
		t2[i].clear();
	t2.clear();
	for (unsigned int i = 0; i < noe.size(); i++)
		noe[i].clear();
	noe.clear();
	for (unsigned int i = 0; i < ccrrt1.size(); i++)
		ccrrt1[i].clear();
	ccrrt1.clear();
	for (unsigned int i = 0; i < ccrrt2.size(); i++)
		ccrrt2[i].clear();
	ccrrt2.clear();

	alphav.clear();
	betav.clear();
	gammav.clear();

	potentialCoeff.clear();
	foundPotential.clear();
	
	/*************
	 * Open file *
	 *************/
	
	std::fstream dataFile(fileName.c_str(),std::ios::in);

	if (!dataFile.is_open())
	{
		if (isFitting)
		{
			if (!rank)
				std::cout << std::endl << std::endl << "C++OPPS ERROR : cannot find/open exp file. " << fileName << std::endl << std::endl;
#ifdef _MPI_
				int mpi_err = MPI_Finalize();
#endif
				exit(1);
		}
		else
		{
			if (!rank)
				std::cout << std::endl << std::endl << " WARNING : _copps.exp file not found so no comparison is possible with experimental data." << fileName << std::endl << std::endl;
#ifdef _MPI_
				int mpi_err = MPI_Barrier(MPI_COMM_WORLD);
#endif
			return false;
		}
	}
	
	/*************
	 * Read file *
	 *************/
	
	ncomp = -1;
	
	hasT1     = false;
	hasT2     = false;
	hasNOE    = false;
	hasCCRRT1 = false;
	hasCCRRT2 = false;

	int position1,position2,position3,position4,position5,position6;
	char separator = ':';
	char fileLine[2048];
	double conversion;
	double re,im;
	dcomplex rel;
	dvector vzero;
	dcvector cvzero;
	std::string str, keyword, value, sigma, dimension, tmpStr;
	std::string av,bv,gv;
	double totalEuler = 0.0;
	potential voidPotential;

	voidPotential.npop = 0;
	voidPotential.c20.clear();
	voidPotential.c22.clear();
	voidPotential.c40.clear();
	voidPotential.c42.clear();
	voidPotential.c44.clear();

	while (dataFile.getline(fileLine,2048)){
		if (fileLine[0] != '#'){
			str.assign(fileLine);
			for (unsigned int i=0; i<str.length(); i++) str[i] = tolower(str[i]);
			position1 = str.find_first_of(separator,0);
			keyword = str.substr(0,position1);
			if (!keyword.compare("component"))
			{
				value = str.substr(position1+1);
				ncomp++;
				nexp.push_back(0);
				nfield.push_back(0);
				field.push_back(vzero);
				t1.push_back(cvzero);
				t2.push_back(cvzero);
				noe.push_back(cvzero);
				ccrrt1.push_back(cvzero);
				ccrrt2.push_back(cvzero);
				potentialCoeff.push_back(voidPotential);
				foundPotential.push_back(0);
			}
			if (!keyword.compare("field"))
			{
				nfield[ncomp]++;
				position2 = str.find_first_of(separator,position1+1);
				value     = str.substr(position1+1,position2-position1-1);
				dimension = str.substr(position2+1);
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				field[ncomp].push_back(conversion * atof(value.c_str()));
			}
			if (!keyword.compare("t1"))
			{
				position2 = str.find_first_of(separator,position1+1);
				position3 = str.find_first_of(separator,position2+1);
				value     = str.substr(position1+1,position2-position1-1);
				sigma     = str.substr(position2+1,position3-position2-1);
				dimension = str.substr(position3+1);
				conversion = 1.0;
				if (!(dimension.compare("ms")))
					conversion = 1.0e-3;
				else if (!(dimension.compare("us")))
					conversion = 1.0e-6;
				re = conversion * atof(value.c_str());
				im = conversion * atof(sigma.c_str());
				rel = dcomplex(re,im);
				t1[ncomp].push_back(rel);
				hasT1 = true;
			}
			if (!keyword.compare("t2"))
			{
				position2 = str.find_first_of(separator,position1+1);
				position3 = str.find_first_of(separator,position2+1);
				value     = str.substr(position1+1,position2-position1-1);
				sigma     = str.substr(position2+1,position3-position2-1);
				dimension = str.substr(position3+1);
				conversion = 1.0;
				if (!(dimension.compare("ms")))
					conversion = 1.0e-3;
				else if (!(dimension.compare("us")))
					conversion = 1.0e-6;
				re = conversion * atof(value.c_str());
				im = conversion * atof(sigma.c_str());
				rel = dcomplex(re,im);
				t2[ncomp].push_back(rel);
				hasT2 = true;
			}
			if (!keyword.compare("noe"))
			{
				position2 = str.find_first_of(separator,position1+1);
				value     = str.substr(position1+1,position2-position1-1);
				sigma     = str.substr(position2+1);
				re = atof(value.c_str());
				im = atof(sigma.c_str());
				rel = dcomplex(re,im);
				noe[ncomp].push_back(rel);
				hasNOE = true;
			}
                        if (!keyword.compare("ccrrt1"))
                        {
                                position2 = str.find_first_of(separator,position1+1);
                                position3 = str.find_first_of(separator,position2+1);
                                value     = str.substr(position1+1,position2-position1-1);
                                sigma     = str.substr(position2+1,position3-position2-1);
                                dimension = str.substr(position3+1);
                                rel.real(atof(value.c_str()));
                                rel.imag(atof(sigma.c_str()));
                                ccrrt1[ncomp].push_back(rel);
				hasCCRRT1 = true;
                        }
                        if (!keyword.compare("ccrrt2"))
                        {
                                position2 = str.find_first_of(separator,position1+1);
                                position3 = str.find_first_of(separator,position2+1);
                                value     = str.substr(position1+1,position2-position1-1);
                                sigma     = str.substr(position2+1,position3-position2-1);
                                dimension = str.substr(position3+1);
                                rel.real(atof(value.c_str()));
                                rel.imag(atof(sigma.c_str()));
                                ccrrt2[ncomp].push_back(rel);
				hasCCRRT2 = true;
                        }

			//===================================================================================================//
			//                         WHAT FOLLOWS PERTAINS TO SRLS CALCULATIONS ONLY                           //
			//===================================================================================================//
			if (!keyword.compare("euler_angles"))
			{
				position2 = str.find_first_of(separator,position1+1);
				position3 = str.find_first_of(separator,position2+1);
				position4 = str.find_first_of(separator,position3+1);
				av        = str.substr(position1+1,position2-position1-1);
				bv        = str.substr(position2+1,position3-position2-1);
				gv        = str.substr(position3+1,position4-position3-1);
				dimension = str.substr(position4+1);
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				alphav.push_back(conversion * atof(av.c_str()));
				betav.push_back(conversion * atof(bv.c_str()));
				gammav.push_back(conversion * atof(gv.c_str()));
				totalEuler += fabs(alphav.at(ncomp)) + fabs(betav.at(ncomp)) + fabs(gammav.at(ncomp));
			}
			if (!keyword.compare("potential"))
			{
				foundPotential[ncomp] = 1;

				position2 = str.find_first_of(separator,position1+1);
				position3 = str.find_first_of(separator,position2+1);
				position4 = str.find_first_of(separator,position3+1);
				position5 = str.find_first_of(separator,position4+1);
				position6 = str.find_first_of(separator,position5+1);

	                        //======================================================//
	                        //======================================================//
	                        //==                 NOTE ON SIGNS                    ==//
	                        //======================================================//
	                        //== The sign of the coefficients is inverted with    ==//
	                        //== respect to the input values because the code is  ==//
	                        //== implemented with a different notation to the one ==//
	                        //== that is found in the papers of Polimeno, Freed   ==//
	                        //== and Meirovitch. In particular, their definition  ==//
        	                //== of the potential is the following:               ==//
	                        //==                                                  ==//
	                        //== U = -c20 D(2,0,0) - c22 [D(2,0,2) + D(2,0,-2)]   ==//
	                        //==                                                  ==//
	                        //== while C++OPPS uses the following notation:       ==//
	                        //==                                                  ==//
	                        //== U = c20 D(2,0,0) + c22 [D(2,0,2) + D(2,0,-2)]    ==//
	                        //==                                                  ==//
	                        //== Input and Output of C++OPPS mantain the notation ==//
	                        //== found in the past leterature, so that it is      ==//
	                        //== easier to compare results with past simulations. ==//
        	                //======================================================//
	                        //======================================================//

				tmpStr = str.substr(position1+1,position2-position1-1);
				potentialCoeff[ncomp].c20.push_back(-atof(tmpStr.c_str()));
				tmpStr = str.substr(position2+1,position3-position2-1);
				potentialCoeff[ncomp].c22.push_back(-atof(tmpStr.c_str()));
				tmpStr = str.substr(position3+1,position4-position3-1);
				potentialCoeff[ncomp].c40.push_back(-atof(tmpStr.c_str()));
				tmpStr = str.substr(position4+1,position5-position4-1);
				potentialCoeff[ncomp].c42.push_back(-atof(tmpStr.c_str()));
				tmpStr = str.substr(position5+1,position6-position5-1);
				potentialCoeff[ncomp].c44.push_back(-atof(tmpStr.c_str()));

				potentialCoeff[ncomp].npop = potentialCoeff[ncomp].c20.size();
			}
			
		}
	}
	
	ncomp++;

	if (totalEuler > ZERO) isOmegaVFromGeometry = true; else isOmegaVFromGeometry = false;

	nExpData = 0;
	if (hasT1)     nExpData++;
	if (hasT2)     nExpData++;
	if (hasNOE)    nExpData++;
	if (hasCCRRT1) nExpData++;
	if (hasCCRRT2) nExpData++;
	
	return true;
}

/************************************************************************************************************************************
 * Return the experimental data for a given component or for all following the scheme                                               *
 * v[c] = T1(B1) sigmaT1(B1) T2(B1) sigmaT2(B1) NOE(B1) sigmaNOE(B1) ... T1(Bn) sigmaT1(Bn) T2(Bn) sigmaT2(Bn) NOE(Bn) sigmaNOE(Bn) * 
 * c = IN = number of component / OUT = size of v                                                                                   *
 ************************************************************************************************************************************/

double * experimental::getExpDataOfComp(int* c)
{
	int ic = c[0];
	if (ic > ncomp)
	{
		std::cout << std::endl << "ERROR in experimental::getDataOfComp(int c). c = " << ic << " excedes n components = " << ncomp << std::endl << std::endl;
		exit(0);
	}
	if (ic <= 0)
		{
			std::cout << std::endl << "ERROR in experimental::getDataOfComp(int c). c = " << ic << " must be grater than 0 " << std::endl << std::endl;
			exit(0);
		}
	ic--;
	double *v = new double[(nExpData*2)*field[ic].size()];
	int pos = 0;
	int j;
	for (unsigned int i = 0; i < field[ic].size(); i++)
	{
		j = 0;
		if (hasT1)
		{
			v[pos+j]   = t1[ic][i].real();
			v[pos+j+1] = t1[ic][i].imag();
			j += 2;
		}
		if (hasT2)
		{
			v[pos+j]   = t2[ic][i].real();
			v[pos+j+1] = t2[ic][i].imag();
			j += 2;
		}
		if (hasNOE)
		{
			v[pos+j]   = noe[ic][i].real();
			v[pos+j+1] = noe[ic][i].imag();
			j += 2;
		}
		if (hasCCRRT1)
		{
			v[pos+j]   = ccrrt1[ic][i].real();
			v[pos+j+1] = ccrrt1[ic][i].imag();
			j += 2;
		}
		if (hasCCRRT2)
		{
			v[pos+j]   = ccrrt2[ic][i].real();
			v[pos+j+1] = ccrrt2[ic][i].imag();
			j += 2;
		}

		pos += 2*nExpData;

	}

	*c = (2*nExpData) * field[ic].size();
	return v;
}

/************************************************************
 * Return the number of experimental points for component c *
 ************************************************************/

int experimental::getNPointsOfComp(int c)
{
	if (c > ncomp)
	{
		std::cout << std::endl << "ERROR in experimental::getDataOfComp(int c). c = " << c << " excedes n components = " << ncomp << std::endl << std::endl;
		exit(0);
	}
	if (c <= 0)
		{
			std::cout << std::endl << "ERROR in experimental::getDataOfComp(int c). c = " << c << " must be grater than 0 " << std::endl << std::endl;
			exit(0);
		}
	c--;
	return 3*field[c].size();
}

/***********************************
 * Return the number of components *
 ***********************************/

int experimental::getNComp(void)
{
	return ncomp;
}

/***************************************
 * Return the alphaV of i-th component *
 ***************************************/

double experimental::getAlphaV(int i)
{
	return alphav.at(i);
}

/**************************************
 * Return the betaV of i-th component *
 **************************************/

double experimental::getBetaV(int i)
{
	return betav.at(i);
}

/***************************************
 * Return the gammaaV of i-th component *
 ***************************************/

double experimental::getGammaV(int i)
{
	return gammav.at(i);
}

/***************************************************
 * Return the list of fields for a given component *
 ***************************************************/

dvector experimental::getFieldOfComp(int c)
{
	if (c > ncomp)
	{
		std::cout << std::endl << "ERROR in experimental::getDataOfComp(int c). c = " << c << " excedes n components = " << ncomp << std::endl << std::endl;
		exit(0);
	}
	if (c <= 0)
		{
			std::cout << std::endl << "ERROR in experimental::getDataOfComp(int c). c = " << c << " must be greater than 0 " << std::endl << std::endl;
			exit(0);
		}
	c--;
	return field[c];
}

/***************************
 * Potential related stuff *
 ***************************/

int experimental::compHasPotential(int c)
{
	if (c > ncomp)
	{
		std::cout << std::endl << "ERROR in experimental::compHasPotential(int c). c = " << (c+1) << " excedes n components = " << ncomp << std::endl << std::endl;
		exit(0);
	}
	if (c <= 0)
		{
			std::cout << std::endl << "ERROR in experimental::getDataOfComp(int c). c = " << c << " must be greater than 0 " << std::endl << std::endl;
			exit(0);
		}
	c--;
	return foundPotential[c];
}

potential experimental::getPotentialOfComp(int c)
{
	if (c > ncomp)
	{
		std::cout << std::endl << "ERROR in experimental::compHasPotential(int c). c = " << (c+1) << " excedes n components = " << ncomp << std::endl << std::endl;
		exit(0);
	}
	if (c <= 0)
		{
			std::cout << std::endl << "ERROR in experimental::getDataOfComp(int c). c = " << c << " must be greater than 0 " << std::endl << std::endl;
			exit(0);
		}
	c--;
	return potentialCoeff[c];
}

/****************************************************/
/* Set / get the component number that is being fit */
/****************************************************/

int experimental::getActiveComponent(void)
{
	return activeComponent;
}

void experimental::setActiveComponent(int ac)
{
	activeComponent = ac;
	return;
}

/************************************************/
/* Get the number of experimental data measured */
/************************************************/

int experimental::getNExpData(void)
{
	return nExpData;
}

/*******************************************************************/
/* Return if an experimental observable is present in the exp file */
/*******************************************************************/

bool experimental::hasExpValueOf(std::string s)
{
	if (!(s.compare("t1")) || !(s.compare("T1")))
		return hasT1;
	else if (!(s.compare("t2")) || !(s.compare("T2")))
		return hasT2;
	else if (!(s.compare("noe")) || !(s.compare("NOE")))
		return hasNOE;
	else if (!(s.compare("ccrrt1")) || !(s.compare("CCRRT1")))
		return hasCCRRT1;
	else if (!(s.compare("ccrrt2")) || !(s.compare("CCRRT2")))
		return hasCCRRT2;
	else
		return false;
}

/*******************************
 * Write out experimental data *
 *******************************/

std::string experimental::toString(void)
{
	std::ostringstream ostr;
	
	ostr << "TASK " << rank << ": Experimental data:" << std::endl;
	for (int i = 0; i < ncomp; i++)
	{
		ostr << "TASK " << rank << ": --------------------------------" << std::endl;
		ostr << "TASK " << rank << ": Component " << i+1 << std::endl;
		ostr << "TASK " << rank << ": --------------------------------" << std::endl;
		for (unsigned int j = 0; j < field[i].size(); j++)
		{
			ostr << "TASK " << rank << ": - Field = " << std::fixed << std::setprecision(2) << field[i][j]*1.0e-6 << " MHz" << std::endl;
			if (hasT1)
				ostr << "TASK " << rank << ":   T1      = " << std::fixed << std::setprecision(2) << t1[i][j].real()*1.0e3 << " +/- " << t1[i][j].imag()*1.0e3 << " ms" << std::endl;
			if (hasT2)
				ostr << "TASK " << rank << ":   T2      = " << std::fixed << std::setprecision(2) << t2[i][j].real()*1.0e3 << " +/- " << t2[i][j].imag()*1.0e3 << " ms" << std::endl;
			if (hasNOE)
				ostr << "TASK " << rank << ":   NOE     = " << std::fixed << std::setprecision(4) << noe[i][j].real() << " +/- " << noe[i][j].imag() << std::endl;
			if (hasCCRRT1)
				ostr << "TASK " << rank << ":   CCRRT1  = " << std::fixed << std::setprecision(4) << ccrrt1[i][j].real() << " +/- " << ccrrt1[i][j].imag() << std::endl;
			if (hasCCRRT2)
				ostr << "TASK " << rank << ":   CCRRT2  = " << std::fixed << std::setprecision(4) << ccrrt2[i][j].real() << " +/- " << ccrrt2[i][j].imag() << std::endl;
		}
		ostr << "TASK " << rank << ":   Euler angles = (" << std::fixed << std::setprecision(2) << alphav.at(i)*RAD_TO_DEG << ", " << betav.at(i)*RAD_TO_DEG << ", " << gammav.at(i)*RAD_TO_DEG << ") deg" << std::endl;
		if (foundPotential[i])
		{
			ostr << "TASK " << rank << ":   List of potential coefficients:" << std::endl;
			for (int j = 0; j < potentialCoeff[i].npop; j++)
			{
				ostr << "TASK " << rank << ":     " << (j+1) << ") ";
				ostr << std::fixed << std::setprecision(2) << "c20 = " << potentialCoeff[i].c20.at(j) << "; ";
				ostr << std::fixed << std::setprecision(2) << "c22 = " << potentialCoeff[i].c22.at(j) << "; ";
				ostr << std::fixed << std::setprecision(2) << "c40 = " << potentialCoeff[i].c40.at(j) << "; ";
				ostr << std::fixed << std::setprecision(2) << "c42 = " << potentialCoeff[i].c42.at(j) << "; ";
				ostr << std::fixed << std::setprecision(2) << "c44 = " << potentialCoeff[i].c44.at(j);
				ostr << std::endl;
			}
		}
	}
	ostr << "TASK " << rank << ": --------------------------------" << std::endl;
	
	return ostr.str();
}


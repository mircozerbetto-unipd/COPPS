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
 Name        : physics.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class to read user defined physical data
 ============================================================================
 */
#include "physics.h"
#include <string>
#include <fstream>
#include <sstream>

#ifdef __EXPCOPPS__
#undef __EXPCOPPS__
#endif

bool cdy1, cdz1, cdy2, cdz2;
bool ratio22, ratio42, ratio44;
bool acfFound;

// Constructor
physics::physics()
{
	// set some defaults for 15N-1H {Damberg et al. JACS (2005) 127, 1995-2005}
	gyromag = -2.7116e7;
	bondLength = 1.015e-10;
	deltaCSA = -169.0;
	// Set some defaults for the TS-X calculation
	population = 1.0;
	jumpFrequency = 0.0;
	
}

// Deconstructor
physics::~physics()
{
#ifdef WRITE_DESTROY_MESSAGE
	std::cout << "TASK " << rank << ": Cleared physics object" << std::endl;
#endif
}

// Initiators
void physics::init(void)
{
	rank = 0;
	std::cout << "TASK " << rank << ": Created physics object" << std::endl;
}

void physics::init(int r)
{
	rank = r;
	std::cout << "TASK " << rank << ": Created physics object" << std::endl;
}

// Get rank of task in MPI calculation
int physics::getRank(void)
{
	return rank;
}

// Set rank of task in MPI calculation
void physics::setRank(int r)
{
	rank = r;
	return;
}

// Return the DXX component of protein D in the frame where it is diagonal
long double physics::getProteinDxx(void)
{
	return proteinD.getCartesianComponent("dxx");
}

// Modify the DXX component of protein D in the frame where it is diagonal
void physics::setProteinDxx(long double val)
{
	proteinD.setCartesianComponent("dxx",val);
	return;
}

// Return the DYY component of protein D in the frame where it is diagonal
long double physics::getProteinDyy(void)
{
	return proteinD.getCartesianComponent("dyy");
}

// Modify the DYY component of protein D in the frame where it is diagonal
void physics::setProteinDyy(long double val)
{
	proteinD.setCartesianComponent("dyy",val);
	return;
}

// Return the DZZ component of protein D in the frame where it is diagonal
long double physics::getProteinDzz(void)
{
	return proteinD.getCartesianComponent("dzz");
}

// Modify the DZZ component of protein D in the frame where it is diagonal
void physics::setProteinDzz(long double val)
{
	proteinD.setCartesianComponent("dzz",val);
	return;
}

// Scale the DXX, DYY and DZZ components of protein D
void physics::scaleProteinD(long double val)
{
	proteinD.scale(val);
}

// Make spherical components of protein D
void physics::transformProteinD(void)
{
	proteinD.transform();
}

// Return D(l,m) spherical component of protein D in potential frame
ldcomplex physics::getProteinDlm(int l, int m)
{
	cldvector sphcomp = proteinD.getSphericalComponents();
	switch (l){
		case 0:
			return sphcomp[0];
		case 2:
			return sphcomp[2+m+1];
	}

	return std::complex<long double>(0.0,0.0);
}

// Return the DXX component of probe D in the frame where it is diagonal
long double physics::getProbeDxx(void)
{
	return probeD.getCartesianComponent("dxx");
}

// Modify the DXX component of probe D in the frame where it is diagonal
void physics::setProbeDxx(long double val)
{
	probeD.setCartesianComponent("dxx",val);
	return;
}

// Return the DYY component of probe D in the frame where it is diagonal
long double physics::getProbeDyy(void)
{
	return probeD.getCartesianComponent("dyy");
}

// Modify the DYY component of probe D in the frame where it is diagonal
void physics::setProbeDyy(long double val)
{
	probeD.setCartesianComponent("dyy",val);
	return;
}

// Return the DZZ component of probe D in the frame where it is diagonal
long double physics::getProbeDzz(void)
{
	return probeD.getCartesianComponent("dzz");
}

// Modify the DZZ component of probe D in the frame where it is diagonal
void physics::setProbeDzz(long double val)
{
	probeD.setCartesianComponent("dzz",val);
	return;
}

// Scale the DXX, DYY and DZZ components of probe D
void physics::scaleProbeD(long double val)
{
	probeD.scale(val);
}

// Make spherical components of probe D
void physics::transformProbeD(void)
{
	probeD.transform();
}

// Return D(l,m) spherical component of probe D in oriented frame
ldcomplex physics::getProbeDlm(int l, int m)
{
	cldvector sphcomp = probeD.getSphericalComponents();
	switch (l){
		case 0:
			return sphcomp[0];
		case 2:
			return sphcomp[2+m+1];
	}

	return std::complex<long double>(0.0,0.0);
}

/*******************/
/* Read input file */
/*******************/

void physics::readInputFile(std::string inFileName){
	
	bool toFit = false, toFit2 = false;;
	int position1,position2,position3;
	char fileLine[2048];
	char separator;
	dcomplex cval;
	long double conversion;
	std::string str,parameter,value,dimension,toFitStr;

	outputString.str("");
	
	/************************************/
	/* INITIALIZATION OF SOME VARIABLES */
        /************************************/

	/* For SRLS  we consider that the potential have the fixed form
	 * U = c20d200(b) + c22[d202(b)+d20-2(b)]cos(2g) + c40d400(b) +c42[d402(b)+d40-2(b)]cos(2g) +c44[d404(b)+d40-4(b)]cos(4g)
	 * with dlmk = reduced Wigner matrix, b = beta and g = gamma.
	 * However the calculation of matrix elements has been implemented for the general case:
	 * Sum(n=0 to N) Sum(m=-n to n) cnm*Dn0m(Omega)
	 * with the assumption cnm = (-)^m cn-m*
	 * so the input will be adapted to the implementation
	 */
	n_coeff = 5;
	coeff = new long double[15];
	for (int i=0; i<15; i++) coeff[i] = 0.0;
	potentialFromExpFile = false;

        Rexchange = 0.0;
      	
	separator = ':';
	fitting = false;
	nfit = nfit2 = 0;
	nstep = 0;
	probeLmax = 0;
	fitMethod = "minpack";
	field.clear();
	fb1_coeff.clear();
	fb2_coeff.n1Max = fb2_coeff.n2Max = 0;
	fb2_coeff.dim1  = fb2_coeff.dim2  = 0;
	fb2_coeff.c = NULL;
	
	cdy1 = false; cdz1 = false;
	cdy2 = false; cdz2 = false;

	ratio22 = false;
	ratio42 = false;
	ratio44 = false;
	
	fitTol1 = 1.0e-5; // xtol = ftol = gtol
	fitTol2 = 1.0e-2; // epsfcn

	Dmin = 1.0e20;
	Krect = 0.0;
	
	dynModel = "srls";

	constrain_OmegaV_OmegaD = false;
	constrain_OmegaD_OmegaV = false;

	for (int i = 0; i < NPARS; i++) fitpar[i] = fitpar2[i] = false;

	order_parameters_tolerance = 0.01;

	nHydrogens = 1;

	/**********************************************
	 * boolean operators to check input integrity *
	 **********************************************/

	bool dynModelSetByUser = false;

	bool fieldFound = false;

	bool LmaxFound = false;
	bool NmaxFound = false;
	bool N1maxFound = false;
	bool N2maxFound = false;
	
	bool proteinDxxFound = false, proteinDyyFound = false, proteinDzzFound = false;
	bool proteinAlphaFound = false, proteinBetaFound = false, proteinGammaFound = false;

	bool probeDxxFound = false, probeDyyFound = false, probeDzzFound = false;
	bool probeAlphaFound = false, probeBetaFound = false, probeGammaFound = false;

	bool fb1DxxFound = false, fb1DyyFound = false, fb1DzzFound = false;
	bool fb1DxiFound = false, fb1DyiFound = false, fb1DziFound = false;
	bool fb1DiiFound = false;

	bool fb2DxxFound = false, fb2DyyFound = false, fb2DzzFound = false;
	bool fb2D11Found = false, fb2D22Found = false, fb2D12Found = false;
	bool fb2Dx1Found = false, fb2Dy1Found = false, fb2Dz1Found = false;
	bool fb2Dx2Found = false, fb2Dy2Found = false, fb2Dz2Found = false;
	
	bool alphaDipFound  = false, betaDipFound  = false, gammaDipFound  = false;	
	bool alphaDip2Found = false, betaDip2Found = false, gammaDip2Found = false;	
	bool alphaCsaFound  = false, betaCsaFound  = false, gammaCsaFound  = false;

	bool c20Found = false, c22Found = false, c40Found = false, c42Found = false, c44Found = false;

	bool rexFound = false;

	bool nucleusFound = false;

	bool rectFound = false;

	bool tolFound = false;

	bool cubatureParamsFound = false;

	bool consistentD1yConstraint = true;
	bool consistentD1zConstraint = true;
	bool consistentD2yConstraint = true;
	bool consistentD2zConstraint = true;

	acfFound = false;
	bool spdFound = false;

	bool orderParametersToleranceFound = false;

	bool nHydrogensFound = false;

	/* TS-X checks */
	
	bool populationFound = false;
	bool jumpFrequencyFound = false;
	bool hchSigmaFound = false;
	bool alphaDipState1Found = false, betaDipState1Found = false, gammaDipState1Found = false;
	bool alphaCsaState1Found = false, betaCsaState1Found = false, gammaCsaState1Found = false;
	bool alphaDip2State1Found = false, betaDip2State1Found = false, gammaDip2State1Found = false;
	bool alphaDipState2Found = false, betaDipState2Found = false, gammaDipState2Found = false;
	bool alphaCsaState2Found = false, betaCsaState2Found = false, gammaCsaState2Found = false;
	bool alphaDip2State2Found = false, betaDip2State2Found = false, gammaDip2State2Found = false;

	/* 3S-DCM checks */
	int _3SFB_OmegaDip1Site1Found = 0;
	int _3SFB_OmegaDip1Site2Found = 0;
	int _3SFB_OmegaDip1Site3Found = 0;
	int _3SFB_OmegaDip2Site1Found = 0;
	int _3SFB_OmegaDip2Site2Found = 0;
	int _3SFB_OmegaDip2Site3Found = 0;
	int _3SFB_OmegaCSASite1Found  = 0;
	int _3SFB_OmegaCSASite2Found  = 0;
	int _3SFB_OmegaCSASite3Found  = 0;
	int _3SFB_populationFound     = 0;
	int _3SFB_jumpFrequencyFound  = 0;

	/****************************/
	/* START READING INPUT FILE */
	/****************************/

	double dtemp;

#ifdef WRITE_ALL	
	std::cout << "TASK " << rank << ": reading input file " << inFileName << std::endl;
#endif

	std::fstream inFile(inFileName.c_str(),std::ios::in);

	if(!inFile.is_open())
	{
		if (!rank)
			std::cout << std::endl << std::endl << "C++OPPS ERROR : cannot find/open input file " << inFileName << std::endl << std::endl;
#ifdef _MPI_
		int mpi_err = MPI_Finalize();
#endif
		exit(1);
	}

	while (inFile.getline(fileLine,2048)){
		if (fileLine[0] != '#' &&  fileLine[0] != ' '){
			
			str.assign(fileLine);
					
			for (unsigned int i=0; i<str.length(); i++) str[i] = tolower(str[i]);

			position1 = str.find_first_of(separator,0);
			position2 = str.find_first_of(separator,position1+1);
			position3 = str.find_first_of(separator,position2+1);
			parameter = str.substr(0,position1);
			value     = str.substr(position1+1,position2-position1-1);
			dimension = str.substr(position2+1,position3-position2-1);
			toFitStr  = str.substr(position3+1);
			toFit     = (!toFitStr.compare("fit")  ? true : false);
			toFit2    = (!toFitStr.compare("fit2") ? true : false);

			fitProteinDiffWithConstraints = false;
			fitProbeDiffWithConstraints = false;

			/****************************
			 * Scan for dynamical model *
			 ****************************/

			if (!(parameter.compare("dynamics"))){
				dynModel = value;
				dynModelSetByUser = true;
				if (!rank)
					outputString << "TASK 0: " << "Model for dynamics = " << dynModel << std::endl;
				isFB1 = !(dynModel.compare("fb1")) || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb")); 
			}

			/*****************************
			 * Diffusion tensor for SRLS *
			 *****************************/

			// Scan for protein D
			else if (!(parameter.compare("protein_dxx"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				proteinD.setCartesianComponent("dxx", conversion * atof(value.c_str()));
				fitpar[1] = !cdy1 && toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[1] = !cdy1 && toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				consistentD1yConstraint = toFit || toFit2;
				consistentD1zConstraint = consistentD1zConstraint && (toFit || toFit2);
				proteinDxxFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "Protein Dxx = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
			}
			else if (!(parameter.compare("protein_dyy"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				proteinD.setCartesianComponent("dyy", conversion * atof(value.c_str()));
				fitpar[0] = toFit || cdz1;
				nfit += (toFit ? 1 : 0);
				fitpar2[0] = toFit2 || cdz1;
				nfit2 += (toFit2 ? 1 : 0);
				if (!(toFitStr).compare("constrain"))
				{
					cdy1 = true;
					fitpar[1] = false; 
					fitpar[0] = true;
				}
				else if (!(toFitStr).compare("constrain2"))
				{
					cdy1 = true;
					fitpar2[1] = false; 
					fitpar2[0] = true;
				}
				else
					cdy1 = false;
				consistentD1zConstraint = consistentD1zConstraint && ((toFit || toFit2) || cdy1);
				proteinDyyFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "Protein Dyy = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
			}
			else if (!(parameter.compare("protein_dzz"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				proteinD.setCartesianComponent("dzz", conversion * atof(value.c_str()));
				fitpar[2] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[2] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				if (!((toFitStr).compare("constrain"))) 
				{ 
					cdz1 = true; 
					fitpar[2] = false; 
					fitpar[0] = true;
				}
				else if (!((toFitStr).compare("constrain2"))) 
				{ 
					cdz1 = true; 
					fitpar2[2] = false; 
					fitpar2[0] = true;
				}
				else
					cdz1 = false;
				proteinDzzFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "Protein Dzz = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
			}
			
			// Scan for Omega_V Euler angels
			else if (!(parameter.compare("protein-potential_alpha"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				proteinD.setAngle("alpha", conversion * atof(value.c_str()));
				if (!(toFitStr).compare("constrain"))
				{
					toFit = true;
					constrain_OmegaD_OmegaV = true;
				}
				else if (!(toFitStr).compare("constrain2"))
				{
					toFit2 = true;
					constrain_OmegaD_OmegaV = true;
				}
				fitpar[6] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[6] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				proteinAlphaFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "TASK 0:" <<  "alpha(VF->M1F) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("protein-potential_beta"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				proteinD.setAngle("beta", conversion * atof(value.c_str()));
				if (!(toFitStr).compare("constrain"))
				{
					toFit = true;
					constrain_OmegaD_OmegaV = true;
				}
				else if (!(toFitStr).compare("constrain2"))
				{
					toFit2 = true;
					constrain_OmegaD_OmegaV = true;
				}
				fitpar[7] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[7] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				proteinBetaFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "beta(VF->M1F) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("protein-potential_gamma"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				proteinD.setAngle("gamma", conversion * atof(value.c_str()));
				if (!(toFitStr).compare("constrain"))
				{
					toFit = true;
					constrain_OmegaD_OmegaV = true;
				}
				else if (!(toFitStr).compare("constrain2"))
				{
					toFit2 = true;
					constrain_OmegaD_OmegaV = true;
				}
				fitpar[8] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[8] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				proteinGammaFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "gamma(VF->M1F) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}

			// Scan for probe D
			else if (!(parameter.compare("probe_dxx"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				probeD.setCartesianComponent("dxx", conversion * atof(value.c_str()));
				fitpar[4] = !cdy2 && toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[4] = !cdy2 && toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				consistentD2yConstraint = toFit || toFit2;
				consistentD2zConstraint = consistentD2zConstraint && (toFit || toFit2);
				probeDxxFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "Probe Dxx = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
			}
			else if (!(parameter.compare("probe_dyy"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				probeD.setCartesianComponent("dyy", conversion * atof(value.c_str()));
				fitpar[3] = toFit || cdz2;
				nfit += (toFit ? 1 : 0);
				fitpar2[3] = toFit2 || cdz2;
				nfit2 += (toFit2 ? 1 : 0);
				if (!((toFitStr).compare("constrain"))) 
				{
					cdy2 = true;
					fitpar[4] = false; 
					fitpar[3] = true; 
				}
				else if (!((toFitStr).compare("constrain2"))) 
				{
					cdy2 = true;
					fitpar2[4] = false; 
					fitpar2[3] = true; 
				}
				else
					cdy2 = false;
				consistentD2zConstraint = consistentD2zConstraint && ((toFit || toFit2) || cdy2);
				probeDyyFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "Probe Dyy = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
			}
			else if (!(parameter.compare("probe_dzz"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				probeD.setCartesianComponent("dzz", conversion * atof(value.c_str()));
				fitpar[5] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[5] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				if (!((toFitStr).compare("constrain")))
				{ 
					cdz2 = true;
					fitpar[5] = false; 
					fitpar[3] = true; 
				}
				else if (!((toFitStr).compare("constrain2")))
				{ 
					cdz2 = true;
					fitpar2[5] = false; 
					fitpar2[3] = true; 
				}
				else
					cdz2 = false;
				probeDzzFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "Probe Dzz = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
			}
			
			// Scan for Omega_O Euler angels
			else if (!(parameter.compare("probe-oriented_alpha"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				probeD.setAngle("alpha", conversion * atof(value.c_str()));
				fitpar[9] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[9] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				probeAlphaFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "alpha(OF->M2F) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("probe-oriented_beta"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				probeD.setAngle("beta", conversion * atof(value.c_str()));
				fitpar[10] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[10] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				probeBetaFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "beta(OF->M2F) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("probe-oriented_gamma"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				probeD.setAngle("gamma", conversion * atof(value.c_str()));
				fitpar[11] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[11] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				probeGammaFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "gamma(OF->M2F) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}

			/****************************
			 * Diffusion tensor for FBn *
			 ****************************/

			else if (!(parameter.compare("dxx"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (isFB1)
				{
					fb1_diften.setRRComponent("dxx", conversion * atof(value.c_str()));
					fb1DxxFound = true;
					if (!rank) outputString << "TASK 0: " << "Dxx = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				else if (!dynModel.compare("fb2"))
				{
					fb2_diften.setRRComponent("dxx", conversion * atof(value.c_str()));
					fb2DxxFound = true;
					if (!rank) outputString << "TASK 0: " << "Dxx = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				fitpar[25] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[25] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("dyy"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (isFB1)
				{
					fb1_diften.setRRComponent("dyy", conversion * atof(value.c_str()));
					fb1DyyFound = true;
					if (!rank) outputString << "TASK 0: " << "Dyy = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				else if (!dynModel.compare("fb2"))
				{
					fb2_diften.setRRComponent("dyy", conversion * atof(value.c_str()));
					fb2DyyFound = true;
					if (!rank) outputString << "TASK 0: " << "Dyy = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("dzz"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (isFB1)
				{
					fb1_diften.setRRComponent("dzz", conversion * atof(value.c_str()));
					fb1DzzFound = true;;
					if (!rank) outputString << "TASK 0: " << "Dzz = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				else if (!dynModel.compare("fb2"))
				{
					fb2_diften.setRRComponent("dzz", conversion * atof(value.c_str()));
					fb2DzzFound = true;;
					if (!rank) outputString << "TASK 0: " << "Dzz = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("d11"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (isFB1)
				{
					fb1_diften.setIIComponent(conversion * atof(value.c_str()));
					fb1DiiFound = true;
					if (!rank) outputString << "TASK 0: " << "D11 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				else if (!dynModel.compare("fb2"))
				{
					fb2_diften.setIIComponent("d11",conversion * atof(value.c_str()));
					fb2D11Found = true;
					if (!rank) outputString << "TASK 0: " << "D11 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("d22"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (!dynModel.compare("fb2"))
				{
					fb2_diften.setIIComponent("d22",conversion * atof(value.c_str()));
					fb2D22Found = true;
					if (!rank) outputString << "TASK 0: " << "D22 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("d12")) || !(parameter.compare("d21"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (!dynModel.compare("fb2"))
				{
					fb2_diften.setIIComponent("d12",conversion * atof(value.c_str()));
					fb2D12Found = true;
					if (!rank) outputString << "TASK 0: " << "D12 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("dx1"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (isFB1)
				{
					fb1_diften.setRIComponent("dxi", conversion * atof(value.c_str()));
					fb1DxiFound = true;
					if (!rank) outputString << "TASK 0: " << "Dx1 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				else if (!dynModel.compare("fb2"))
				{
					fb2_diften.setRIComponent("dx1", conversion * atof(value.c_str()));
					fb2Dx1Found = true;
					if (!rank) outputString << "TASK 0: " << "Dx1 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("dy1"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (isFB1)
				{
					fb1_diften.setRIComponent("dyi", conversion * atof(value.c_str()));
					fb1DyiFound = true;
					if (!rank) outputString << "TASK 0: " << "Dy1 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				else if (!dynModel.compare("fb2"))
				{
					fb2_diften.setRIComponent("dy1", conversion * atof(value.c_str()));
					fb2Dy1Found = true;
					if (!rank) outputString << "TASK 0: " << "Dy1 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("dz1"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (isFB1)
				{
					fb1_diften.setRIComponent("dzi", conversion * atof(value.c_str()));
					fb1DziFound = true;
					if (!rank) outputString << "TASK 0: " << "Dz1 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				else if (!dynModel.compare("fb2"))
				{
					fb2_diften.setRIComponent("dz1", conversion * atof(value.c_str()));
					fb2Dz1Found = true;
					if (!rank) outputString << "TASK 0: " << "Dz1 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("dx2"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (!dynModel.compare("fb2"))
				{
					fb2_diften.setRIComponent("dx2", conversion * atof(value.c_str()));
					fb2Dx2Found = true;
					if (!rank) outputString << "TASK 0: " << "Dx2 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("dy2"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (!dynModel.compare("fb2"))
				{
					fb2_diften.setRIComponent("dy2", conversion * atof(value.c_str()));
					fb2Dy2Found = true;
					if (!rank) outputString << "TASK 0: " << "Dy2 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
			else if (!(parameter.compare("dz2"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				if (!dynModel.compare("fb2"))
				{
					fb2_diften.setRIComponent("dz2", conversion * atof(value.c_str()));
					fb2Dz2Found = true;
					if (!rank) outputString << "TASK 0: " << "Dz2 = " << conversion * atof(value.c_str()) << " Hz" << std::endl;
				}
				// TODO
				//fitpar[1] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[1] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
			}
		
			/************************************
			 * Euler angles of magnetic tensors *
			 ************************************/

			// Scan for Omega_D Euler angels (MF -> DF rotation)
			else if (!(parameter.compare("dipolar_alpha"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omega_D[0] = omega_D_bkp[0] = conversion * atof(value.c_str());
				if (!(toFitStr).compare("constrain"))
				{
					toFit = true;
					constrain_OmegaV_OmegaD = true;
				}
				else if (!(toFitStr).compare("constrain2"))
				{
					toFit2 = true;
					constrain_OmegaV_OmegaD = true;
				}
				fitpar[12] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[12] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				alphaDipFound = true;
				if (!rank) outputString << "TASK 0: " << "alpha(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_beta"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omega_D[1] = omega_D_bkp[1] = conversion * atof(value.c_str());
				if (!(toFitStr).compare("constrain"))
				{
					toFit = true;
					constrain_OmegaV_OmegaD = true;
				}
				else if (!(toFitStr).compare("constrain2"))
				{
					toFit2 = true;
					constrain_OmegaV_OmegaD = true;
				}
				fitpar[13] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[13] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				betaDipFound = true;
				if (!rank) outputString << "TASK 0: " << "beta(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_gamma"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omega_D[2] = omega_D_bkp[2] = conversion * atof(value.c_str());
				if (!(toFitStr).compare("constrain"))
				{
					toFit = true;
					constrain_OmegaV_OmegaD = true;
				}
				else if (!(toFitStr).compare("constrain2"))
				{
					toFit2 = true;
					constrain_OmegaV_OmegaD = true;
				}
				fitpar2[14] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				gammaDipFound = true;
				if (!rank) outputString << "TASK 0: " << "gamma(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}

			// Scan for Omega_D2 Euler angels (DF -> D2F rotation)
			else if (!(parameter.compare("dipolar_alpha2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omega_D2[0] = conversion * atof(value.c_str());
				fitpar[26] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[26] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				alphaDip2Found = true;
				if (!rank) outputString << "TASK 0: " << "alpha(DF->D2F) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_beta2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omega_D2[1] = conversion * atof(value.c_str());
				fitpar[27] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[27] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				betaDip2Found = true;
				if (!rank) outputString << "TASK 0: " << "beta(DF->D2F) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_gamma2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omega_D2[2] = conversion * atof(value.c_str());
				fitpar[28] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[28] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				gammaDip2Found = true;
				if (!rank) outputString << "TASK 0: " << "gamma(DF->D2F) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}

			/* TS-X DIPOLAR ANGLES */ // NOTE: NO FITTING IMPLEMENTED FOR THESE ANGLES AT THE MOMENT
			else if (!(parameter.compare("dipolar_alpha_1"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip_1[0] = conversion * atof(value.c_str());
				//fitpar[12] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[12] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				alphaDipState1Found = true;
				if (!rank) outputString << "TASK 0: " << "alpha_1(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_beta_1"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip_1[1] = conversion * atof(value.c_str());
				//fitpar[13] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[13] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				betaDipState1Found = true;
				if (!rank) outputString << "TASK 0: " << "beta_1(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_gamma_1"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip_1[2] = conversion * atof(value.c_str());
				//fitpar[14] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[14] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				gammaDipState1Found = true;
				if (!rank) outputString << "TASK 0: " << "gamma_1(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_alpha_2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip_2[0] = conversion * atof(value.c_str());
				//fitpar[12] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[12] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				alphaDipState2Found = true;
				if (!rank) outputString << "TASK 0: " << "alpha_2(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_beta_2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip_2[1] = conversion * atof(value.c_str());
				//fitpar[13] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[13] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				betaDipState2Found = true;
				if (!rank) outputString << "TASK 0: " << "beta_2(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_gamma_2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip_2[2] = conversion * atof(value.c_str());
				//fitpar[14] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[14] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				gammaDipState2Found = true;
				if (!rank) outputString << "TASK 0: " << "gamma_2(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			
			else if (!(parameter.compare("dipolar_alpha2_1"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip2_1[0] = conversion * atof(value.c_str());
				//fitpar[12] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[12] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				alphaDip2State1Found = true;
				if (!rank) outputString << "TASK 0: " << "alpha2_2(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_beta2_1"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip2_1[1] = conversion * atof(value.c_str());
///////////////////////////////////////////////////
// ONLY TEMPORARY
				fitpar[27] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[27] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
///////////////////////////////////////////////////
				betaDip2State1Found = true;
				if (!rank) outputString << "TASK 0: " << "beta2_1(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_gamma2_1"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip2_1[2] = conversion * atof(value.c_str());
				//fitpar[14] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[14] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				gammaDip2State1Found = true;
				if (!rank) outputString << "TASK 0: " << "gamma2_1(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_alpha2_2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip2_2[0] = conversion * atof(value.c_str());
				//fitpar[12] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[12] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				alphaDip2State2Found = true;
				if (!rank) outputString << "TASK 0: " << "alpha2_2(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_beta2_2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip2_2[1] = conversion * atof(value.c_str());
				//fitpar[13] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[13] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				betaDip2State2Found = true;
				if (!rank) outputString << "TASK 0: " << "beta2_2(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("dipolar_gamma2_2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaDip2_2[2] = conversion * atof(value.c_str());
				//fitpar[14] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[14] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				gammaDip2State2Found = true;
				if (!rank) outputString << "TASK 0: " << "gamma2_2(MF->DF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}

			// Scan for Omega_CSA Euler angels (DF -> CF rotation)
			else if (!(parameter.compare("csa_alpha"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omega_CSA[0] = conversion * atof(value.c_str());
				fitpar[15] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[15] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				alphaCsaFound = true;
				if (!rank) outputString << "TASK 0: " << "alpha(DF->CF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("csa_beta"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omega_CSA[1] = conversion * atof(value.c_str());
				fitpar[16] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[16] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				betaCsaFound = true;
				if (!rank) outputString << "TASK 0: " << "beta(DF->CF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("csa_gamma"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omega_CSA[2] = conversion * atof(value.c_str());
				fitpar[17] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[17] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				gammaCsaFound = true;
				if (!rank) outputString << "TASK 0: " << "gamma(DF->CF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}

			/* TS-X CSA ANGLES */ // NOTE: NO FITTING IMPLEMENTED FOR THESE ANGLES AT THE MOMENT
			else if (!(parameter.compare("csa_alpha_1"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaCsa_1[0] = conversion * atof(value.c_str());
				//fitpar[12] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[12] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				alphaCsaState1Found = true;
				if (!rank) outputString << "TASK 0: " << "csa_1(DF->CF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("csa_beta_1"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaCsa_1[1] = conversion * atof(value.c_str());
				//fitpar[13] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[13] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				betaCsaState1Found = true;
				if (!rank) outputString << "TASK 0: " << "beta_1(DF->CF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("csa_gamma_1"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaCsa_1[2] = conversion * atof(value.c_str());
				//fitpar[14] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[14] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				gammaCsaState1Found = true;
				if (!rank) outputString << "TASK 0: " << "gamma_1(DF->CF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("csa_alpha_2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaCsa_2[0] = conversion * atof(value.c_str());
				//fitpar[12] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[12] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				alphaCsaState2Found = true;
				if (!rank) outputString << "TASK 0: " << "alpha_2(DF->CF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("csa_beta_2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaCsa_2[1] = conversion * atof(value.c_str());
				//fitpar[13] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[13] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				betaCsaState2Found = true;
				if (!rank) outputString << "TASK 0: " << "beta_2(DF->CF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}
			else if (!(parameter.compare("csa_gamma_2"))){
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				omegaCsa_2[2] = conversion * atof(value.c_str());
				//fitpar[14] = toFit;
				//nfit += (toFit ? 1 : 0);
				//fitpar2[14] = toFit2;
				//nfit2 += (toFit2 ? 1 : 0);
				gammaCsaState2Found = true;
				if (!rank) outputString << "TASK 0: " << "gamma_2(DF->CF) = " << conversion * atof(value.c_str()) << " rad" << std::endl;
			}

			/*****************************************/
			/* TS-SRLS POPULATION AND JUMP RATE DATA */
			/*****************************************/

			else if (!(parameter.compare("population")))
			{
				population = atof(value.c_str());
				fitpar[29] = toFit;
				nfit += toFit ? 1 : 0;
				fitpar2[29] = toFit2;
				nfit2 += toFit2 ? 1 : 0;
				populationFound = true;
				if (!rank && (!(dynModel.compare("ts-srls")) || !(dynModel.compare("ts-fb1")))) outputString << "TASK 0: " << "Populations of the two states: P1 = " << population <<", P2 = " << 1.0-population << std::endl;
			}

			else if (!(parameter.compare("jump_frequency")))
			{
				jumpFrequency = atof(value.c_str());
				fitpar[30] = toFit;
				nfit += toFit ? 1 : 0;
				fitpar2[30] = toFit2;
				nfit2 += toFit2 ? 1 : 0;
				jumpFrequencyFound = true;
				if (!rank && (!(dynModel.compare("ts-srls")) || !(dynModel.compare("ts-fb1")))) outputString << "TASK 0: " << "Jump frequency: " << jumpFrequency << " Hz" << std::endl;
			}

			else if (!(parameter.compare("hch_sigma")))
			{
				hchSigma = atof(value.c_str());
				fitpar[31] = toFit;
				nfit += toFit ? 1 : 0;
				fitpar2[31] = toFit2;
				nfit2 += toFit2 ? 1 : 0;
				hchSigmaFound = true;
				if (!rank && (!(dynModel.compare("ts-srls")) || !(dynModel.compare("ts-fb1")))) outputString << "TASK 0: " << "Sigma HCH libration = " << hchSigma <<" deg" << std::endl;
			}

			/********************/
			/* 3S-FB PARAMETERS */
			/********************/

			/* DIP 1 ANGLES */
			else if (!(parameter.compare("3sfb_alpha-dip1-site1")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip1[0][0] = dtemp;
				_3SFB_OmegaDip1Site1Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Alpha Dip1 Site 1 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_alpha-dip1-site2")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip1[0][1] = dtemp;
				_3SFB_OmegaDip1Site2Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Alpha Dip1 Site 2 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_alpha-dip1-site3")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip1[0][2] = dtemp;
				_3SFB_OmegaDip1Site3Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Alpha Dip1 Site 3 = " << dtemp <<" deg" << std::endl;
			}

			else if (!(parameter.compare("3sfb_beta-dip1-site1")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip1[1][0] = dtemp;
				_3SFB_OmegaDip1Site1Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Beta Dip1 Site 1 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_beta-dip1-site2")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip1[1][1] = dtemp;
				_3SFB_OmegaDip1Site2Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Beta Dip1 Site 2 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_beta-dip1-site3")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip1[1][2] = dtemp;
				_3SFB_OmegaDip1Site3Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Beta Dip1 Site 3 = " << dtemp <<" deg" << std::endl;
			}

			else if (!(parameter.compare("3sfb_gamma-dip1-site1")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip1[2][0] = dtemp;
				_3SFB_OmegaDip1Site1Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Gamma Dip1 Site 1 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_gamma-dip1-site2")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip1[2][1] = dtemp;
				_3SFB_OmegaDip1Site2Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Gamma Dip1 Site 2 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_gamma-dip1-site3")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip1[2][2] = dtemp;
				_3SFB_OmegaDip1Site3Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Gamma Dip1 Site 3 = " << dtemp <<" deg" << std::endl;
			}

			/* DIP 2 ANGLES */
			else if (!(parameter.compare("3sfb_alpha-dip2-site1")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip2[0][0] = dtemp;
				_3SFB_OmegaDip2Site1Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Alpha Dip2 Site 1 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_alpha-dip2-site2")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip2[0][1] = dtemp;
				_3SFB_OmegaDip2Site2Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Alpha Dip2 Site 2 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_alpha-dip2-site3")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip2[0][2] = dtemp;
				_3SFB_OmegaDip2Site3Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Alpha Dip2 Site 3 = " << dtemp <<" deg" << std::endl;
			}

			else if (!(parameter.compare("3sfb_beta-dip2-site1")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip2[1][0] = dtemp;
				_3SFB_OmegaDip2Site1Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Beta Dip2 Site 1 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_beta-dip2-site2")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip2[1][1] = dtemp;
				_3SFB_OmegaDip2Site2Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Beta Dip2 Site 2 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_beta-dip2-site3")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip2[1][2] = dtemp;
				_3SFB_OmegaDip2Site3Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Beta Dip2 Site 3 = " << dtemp <<" deg" << std::endl;
			}

			else if (!(parameter.compare("3sfb_gamma-dip2-site1")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip2[2][0] = dtemp;
				_3SFB_OmegaDip2Site1Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Gamma Dip2 Site 1 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_gamma-dip2-site2")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip2[2][1] = dtemp;
				_3SFB_OmegaDip2Site2Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Gamma Dip2 Site 2 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_gamma-dip2-site3")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaDip2[2][2] = dtemp;
				_3SFB_OmegaDip2Site3Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Gamma Dip2 Site 3 = " << dtemp <<" deg" << std::endl;
			}

			/* CSA ANGLES */
			else if (!(parameter.compare("3sfb_alpha-csa-site1")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaCSA[0][0] = dtemp;
				_3SFB_OmegaCSASite1Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Alpha CSA Site 1 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_alpha-csa-site2")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaCSA[0][1] = dtemp;
				_3SFB_OmegaCSASite2Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Alpha CSA Site 2 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_alpha-csa-site3")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaCSA[0][2] = dtemp;
				_3SFB_OmegaCSASite3Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Alpha CSA Site 3 = " << dtemp <<" deg" << std::endl;
			}

			else if (!(parameter.compare("3sfb_beta-csa-site1")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaCSA[1][0] = dtemp;
				_3SFB_OmegaCSASite1Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Beta CSA Site 1 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_beta-csa-site2")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaCSA[1][1] = dtemp;
				_3SFB_OmegaCSASite2Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Beta CSA Site 2 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_beta-csa-site3")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaCSA[1][2] = dtemp;
				_3SFB_OmegaCSASite3Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Beta CSA Site 3 = " << dtemp <<" deg" << std::endl;
			}

			else if (!(parameter.compare("3sfb_gamma-csa-site1")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaCSA[2][0] = dtemp;
				_3SFB_OmegaCSASite1Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Gamma CSA Site 1 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_gamma-csa-site2")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaCSA[2][1] = dtemp;
				_3SFB_OmegaCSASite2Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Gamma CSA Site 2 = " << dtemp <<" deg" << std::endl;
			}
			else if (!(parameter.compare("3sfb_gamma-csa-site3")))
			{
				conversion = 1.0;
				if (!(dimension.compare("deg")))
					conversion = DEG_TO_RAD;
				dtemp = conversion * atof(value.c_str());
				_3SFB_OmegaCSA[2][2] = dtemp;
				_3SFB_OmegaCSASite3Found += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Gamma CSA Site 3 = " << dtemp <<" deg" << std::endl;
			}

			/* BOLTZMANN POPULATIONS */
			else if (!(parameter.compare("3sfb_pop_site1")))
			{
				dtemp = atof(value.c_str());
				_3SFB_population[0] = dtemp;
				_3SFB_sqPopulation[0] = sqrt(dtemp);
				_3SFB_populationFound += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Boltzmann population of Site 1 = " << dtemp << std::endl;
			}
			else if (!(parameter.compare("3sfb_pop_site2")))
			{
				dtemp = atof(value.c_str());
				_3SFB_population[1] = dtemp;
				_3SFB_sqPopulation[1] = sqrt(dtemp);
				_3SFB_populationFound += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Boltzmann population of Site 2 = " << dtemp << std::endl;
			}
			else if (!(parameter.compare("3sfb_pop_site3")))
			{
				dtemp = atof(value.c_str());
				_3SFB_population[2] = dtemp;
				_3SFB_sqPopulation[2] = sqrt(dtemp);
				_3SFB_populationFound += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "Boltzmann population of Site 3 = " << dtemp << std::endl;
			}

			/* JUMP FREQUENCIES */
			else if (!(parameter.compare("3sfb_w_1-2")))
			{
				dtemp = atof(value.c_str());
				_3SFB_jumpFrequency[0] = dtemp;
				_3SFB_jumpFrequencyFound += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "w1->2 = " << dtemp << " Hz" << std::endl;
			}
			else if (!(parameter.compare("3sfb_w_1-3")))
			{
				dtemp = atof(value.c_str());
				_3SFB_jumpFrequency[1] = dtemp;
				_3SFB_jumpFrequencyFound += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "w1->3 = " << dtemp << " Hz" << std::endl;
			}
			else if (!(parameter.compare("3sfb_w_2-3")))
			{
				dtemp = atof(value.c_str());
				_3SFB_jumpFrequency[2] = dtemp;
				_3SFB_jumpFrequencyFound += 1;
				if (!rank && !(dynModel.compare("3s-fb"))) outputString << "TASK 0: " << "w2->3 = " << dtemp << " Hz" << std::endl;
			}

			/**********************************
			 * Coefficients of SRLS potential *
			 **********************************/

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

			// Scan for potential
			else if(!(parameter.compare("c20"))){
				conversion = 1.0;
				coeff[(2*(2+1))/2+0] = -1.0 * conversion * atof(value.c_str());
				fitpar[18] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[18] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				c20Found = true;
				if (!potentialFromExpFile)
                                        potentialFromExpFile = !((toFitStr).compare("exp"));
				if ( ( !(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")) ) && !rank)
				{
					outputString << "TASK 0: " << "c20 = ";
					if (!potentialFromExpFile)
						outputString << -conversion * atof(value.c_str()) << std::endl;
					else
						outputString << "from _copps.exp file" << std::endl;
				}
			}
			else if(!(parameter.compare("c22"))){
				conversion = 1.0;
				coeff[(2*(2+1))/2+2] = -1.0 * conversion * atof(value.c_str());
				if (!(toFitStr.compare("ratio")))
				{
					toFit = true;
					ratio22 = true;
				}
				else if (!(toFitStr.compare("ratio2")))
				{
					toFit2 = true;
					ratio22 = true;
				}
				fitpar[19] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[19] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				c22Found = true;
				if (!potentialFromExpFile)
                                        potentialFromExpFile = !((toFitStr).compare("exp"));
				if ( ( !(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")) ) && !rank)
				{
					outputString << "TASK 0: " << "c22 = ";
					if (!potentialFromExpFile)
						outputString << -conversion * atof(value.c_str()) << std::endl;
					else
						outputString << "from _copps.exp file" << std::endl;
				}
			}
			else if(!(parameter.compare("c40"))){
				conversion = 1.0;
				coeff[(4*(4+1))/2+0] = -1.0 * conversion * atof(value.c_str());
				fitpar[20] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[20] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				c40Found = true;
				if (!potentialFromExpFile)
                                        potentialFromExpFile = !((toFitStr).compare("exp"));
				if ( ( !(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")) ) && !rank)
				{
					outputString << "TASK 0: " << "c40 = ";
					if (!potentialFromExpFile)
						outputString << -conversion * atof(value.c_str()) << std::endl;
					else
						outputString << "from _copps.exp file" << std::endl;
				}
			}
			else if(!(parameter.compare("c42"))){
				conversion = 1.0;
				coeff[(4*(4+1))/2+2] = -1.0 * conversion * atof(value.c_str());
				if (!(toFitStr.compare("ratio")))
				{
					toFit = true;
					ratio42 = true;
				}
				else if (!(toFitStr.compare("ratio2")))
				{
					toFit2 = true;
					ratio42 = true;
				}
				fitpar[21] = toFit;
				nfit += (toFit ? 1 : 0);
				fitpar2[21] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				c42Found = true;
				if (!potentialFromExpFile)
                                        potentialFromExpFile = !((toFitStr).compare("exp"));
				if ( ( !(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")) ) && !rank)
				{
					outputString << "TASK 0: " << "c42 = ";
					if (!potentialFromExpFile)
						outputString << -conversion * atof(value.c_str()) << std::endl;
					else
						outputString << "from _copps.exp file" << std::endl;
				}
			}
			else if(!(parameter.compare("c44"))){
				conversion = 1.0;
				coeff[(4*(4+1))/2+4] = -1.0 * conversion * atof(value.c_str());
				if (!(toFitStr.compare("ratio")))
				{
					toFit = true;
					ratio44 = true;
				}
				else if (!(toFitStr.compare("ratio2")))
				{
					toFit2 = true;
					ratio44 = true;
				}
				fitpar2[22] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				c44Found = true;
				if (!potentialFromExpFile)
                                        potentialFromExpFile = !((toFitStr).compare("exp"));
				if ( ( !(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")) ) && !rank)
				{
					outputString << "TASK 0: " << "c44 = ";
					if (!potentialFromExpFile)
						outputString << -conversion * atof(value.c_str()) << std::endl;
					else
						outputString << "from _copps.exp file" << std::endl;
				}
			}

			/*******************************************
			 * Coefficients for FB1 internal potential *
			 *******************************************/

			else if (!(parameter.compare("potential_coefficient_t1"))){
				cval = dcomplex(atof(value.c_str()), atof(dimension.c_str()));
				if (isFB1)
				{
					fb1_coeff.push_back(cval);
					if (!rank) outputString << "TASK 0: " << "potential coeff t1 =  (" << cval.real() << ", " << cval.imag() << ")" << std::endl;
				}
			}

			/*******************************************
			 * Coefficients for FB2 internal potential *
			 *******************************************/

			else if (!(parameter.compare("potential_coefficients_t2")))
			{
				if (!(dynModel.compare("fb2")))
				{
					fb2_coeff.n1Max = atoi(value.c_str());
					fb2_coeff.dim1  = 2 * fb2_coeff.n1Max + 1;
					fb2_coeff.n2Max = atoi(dimension.c_str());
					fb2_coeff.dim2  = 2 * fb2_coeff.n2Max + 1;
					fb2_coeff.c     = (dcomplex *)calloc(fb2_coeff.dim1 * fb2_coeff.dim2, sizeof(dcomplex));
					if (!rank) outputString << "TASK 0: " << "FB2 potential matrix (" << fb2_coeff.dim1 << " x " << fb2_coeff.dim2 << "):" << std::endl;
					double cr, ci;
					for (int it1 = 0; it1 < fb2_coeff.dim1; it1++)
					{
						std::stringstream ls (std::stringstream::in | std::stringstream::out);
						inFile.getline(fileLine,2048);
						ls.write(fileLine,2048);
						if (!rank) outputString << "TASK 0: ";
						for (int it2 = 0; it2 < fb2_coeff.dim2; it2++)
						{
							ls >> cr >> ci;
							fb2_coeff.c[it1 * fb2_coeff.dim2 + it2] = dcomplex (cr,ci);
							if (!rank) outputString << fb2_coeff.c[it1 * fb2_coeff.dim2 + it2] << "\t";
						}
						if (!rank) outputString << std::endl;
					}
				}
			}

			/********************************/
			/* Max Eval in cubature library */
			/********************************/

			else if (!(parameter.compare("cubatureparams")))
			{
				maxEval = atoi(value.c_str());
				cubatureParamsFound = true;
				relErr  = atof(dimension.c_str());
				if (!rank && !(dynModel.compare("fb2"))) outputString << "TASK 0: maxEval in 2D integation routine: " << maxEval << std::endl << "TASK 0: relErr in 2D integration routine: " << relErr << std::endl;
			}
			
                        /**********************
			 * Scan for R exhange *
			 **********************/

			else if(!(parameter.compare("rexchange"))){
				conversion = 1.0;
				Rexchange = conversion * atof(value.c_str());
				fitpar[23] = toFit;
				nfit += (toFit ? 1 : 0);
				rexFound = true;
				fitpar2[23] = toFit2;
				nfit2 += (toFit2 ? 1 : 0);
				if (!rank) outputString << "TASK 0: " << "Rex =  " << conversion * atof(value.c_str()) << " Hz" << std::endl;
			}

			/****************************
			 * Hydrodynamics parameters *
			 ****************************/

			// Scan for temperature
			else if (!(parameter.compare("temperature"))){
				conversion = 0.0;
				if (!(dimension.compare("c")))
					conversion = 273.15;
				temperature = conversion + atof(value.c_str());
				if (!rank) outputString << "TASK 0: " << "Temperature =  " << temperature << " K" << std::endl;
			}
			// Scan for viscosity
			else if (!(parameter.compare("viscosity"))){
				conversion = 1.0;
				if (!(dimension.compare("p")))
					conversion = 1.0e-1;
				else if (!(dimension.compare("cp")))
					conversion = 1.0e-3;
				else if (!(dimension.compare("mpas")))
					conversion = 1.0e-3;
				else if (!(dimension.compare("upas")))
					conversion = 1.0e-6;
				viscosity = conversion * atof(value.c_str());
				if (!rank) outputString << "TASK 0: " << "Viscosity =  " << viscosity << " Pa s" << std::endl;
			}
			
			// Scan for effective radius
			else if (!(parameter.compare("effective_radius"))){
				conversion = 1.0;
				if (!(dimension.compare("m")))
					conversion = 1.0e10;
				else if (!(dimension.compare("nm")))
					conversion = 1.0e1;
				else if (!(dimension.compare("pm")))
					conversion = 1.0e-2;
				effective_radius = conversion * atof(value.c_str());
				if (!rank) outputString << "TASK 0: " << "Effective radius =  " << effective_radius << " A" << std::endl;
			}
			
			/***************************
			 * Fitting related options *
			 ***************************/

			// Scan for fitting flag
			else if (!(parameter.compare("fitting"))){
				int v = atoi(value.c_str());
				if (v) fitting = true;
				else   fitting = false;
				if (fitting && !rank)
					outputString << "TASK 0: " << "Fitting procedure required  " << std::endl;
				else if (!fitting && !rank)
					outputString << "TASK 0: " << "No fitting procedure required  " << std::endl;
			}
			
			// Scan for fitting method
			else if (!(parameter.compare("fitting_method"))){
				fitMethod.assign(value);
				if (fitting && !rank)
					outputString << "TASK 0: " << "Fitting routine =  " << fitMethod << std::endl;
			}

			// Scan for fit Tolerance (bshomer, 5 Feb 09)
			// *** Changed to accept two values (mz, 9 Nov 09)
                        else if (!(parameter.compare("fit_tolerance")))
                        {
                                fitTol1 = atof(value.c_str());
				fitTol2 = atof(dimension.c_str());
				tolFound = true;
				if (fitting && !rank)
					outputString << "TASK 0: " << "Fitting tolerances = " << fitTol1 << ", " << fitTol2 << std::endl;
                        }

			// Scan for rectification parameters
                        else if (!(parameter.compare("rectification")))
                        {
                                Dmin  = atof(value.c_str());
				Krect = atof(dimension.c_str());
				rectFound = true;
				if (fitting && !rank)
				{
					outputString << "TASK 0: " << "Diffusion tensor chi-square rectification with:" << std::endl;
					outputString << "TASK 0: " << "  Dmin = " << Dmin << " Hz" << std::endl;
					outputString << "TASK 0: " << "  Krec = " << Krect << std::endl;
				}
                        }

			// Scan for oder parameters stop tolerance
			else if (!(parameter.compare("order_parameters_tolerance")))
			{
				order_parameters_tolerance = atof(value.c_str());
				orderParametersToleranceFound = true;
				if (fitting && !rank)
					outputString << "TASK 0: " << "Order parameters tolerance = " << order_parameters_tolerance << std::endl;
			}
			
			/*********************************
			 * Truncation parameter for SRLS *
			 *********************************/

			// Scan for Lmax of probe
			else if (!(parameter.compare("lmax"))){
				probeLmax = atoi(value.c_str());
				LmaxFound = true;
				if (!(dynModel.compare("srls")) && !rank)
					outputString << "TASK 0: " << "Lmax = " << probeLmax << std::endl;
			}

			/********************************
			 * Truncation parameter for FBn *
			 ********************************/

			// Scan for Nmax for internal rotation
			else if (!(parameter.compare("nmax"))){
				Nmax = atoi(value.c_str());
				NmaxFound = true;
				if (isFB1 && !rank)
					outputString << "TASK 0: " << "Nmax = " << Nmax << std::endl;
			}
			
			else if (!(parameter.compare("n1max"))){
				N1max = atoi(value.c_str());
				N1maxFound = true;
				if (isFB1 && !rank)
				{
					Nmax = atoi(value.c_str());
					NmaxFound = true;
					outputString << "TASK 0: " << "Nmax = " << Nmax << std::endl;
				}
				if (!(dynModel.compare("fb2")) && !rank)
					outputString << "TASK 0: " << "N1max = " << N1max << std::endl;
			}
			else if (!(parameter.compare("n2max"))){
				N2max = atoi(value.c_str());
				N2maxFound = true;
				if (!(dynModel.compare("fb2")) && !rank)
					outputString << "TASK 0: " << "N2max = " << N2max << std::endl;
			}
			/**************************
			 * Scan for Lanczos steps *
			 **************************/

			else if(!(parameter.compare("lanczos_max_n_step"))){
				nstep = atoi(value.c_str());
				if (!rank) outputString << "TASK 0: " << "Max Lanczos steps = " << nstep << std::endl;
			}
			
			/**************************************************
			 * Scan for fileds (for single-shot calculations) *
			 **************************************************/

			else if (!(parameter.compare("field"))){
				conversion = 1.0;
				if (!(dimension.compare("khz")))
					conversion = 1.0e3;
				else if (!(dimension.compare("mhz")))
					conversion = 1.0e6;
				else if (!(dimension.compare("ghz")))
					conversion = 1.0e9;
				field.push_back(conversion * atof(value.c_str()));
				fieldFound = true;
				if (!fitting && !rank)
					outputString << "TASK 0: " << "Frequency  = " << field.back()*1.0e-6 << " MHz"  << std::endl;
			}

			/*************************
			 * Scan for nucleus type *
			 *************************/

			else if (!(parameter.compare("nucleus")))
			{
				if (!value.compare("n15"))
					gyromag = -2.7116e7;
				else if (!value.compare("c13"))
					gyromag = 6.7283e7;
				bondLength = atof(dimension.c_str()) * 1.0e-10;
				deltaCSA = atof(toFitStr.c_str());
				nucleusFound = true;
				if (!rank)
				{
					outputString << "TASK 0: " << "Nucleus = " << value << std::endl;
					outputString << "TASK 0: " << "Magnetogyric ratio = " << gyromag << " MHz/T" << std::endl;
					outputString << "TASK 0: " << "Bondlength = " << bondLength << " A" << std::endl;
					outputString << "TASK 0: " << "CSA = " << deltaCSA << " ppm" << std::endl;
				}
			}

			/***********************************************************/
			/* Scan for probe type, i.e. how many H atoms are attached */
			/***********************************************************/
			// NOTE: actual implementation is for 1 (e.g. CH, NH) or 2 (e.g. CH2)

			else if (!(parameter.compare("hydrogens")))
			{
				nHydrogens = atoi(value.c_str());
				if (nHydrogens < 1 || nHydrogens > 2)
				{
					if (!rank) std::cout << std::endl << std::endl << "ERROR: admitted number of hydrogen atoms is 1 or 2" << std::endl << std::endl;
#ifdef _MPI_
					MPI_Finalize();
#endif
					exit(1);
				}
				nHydrogensFound = true;
				if (!rank) outputString << "TASK 0: " << "Number of hydrogen atoms = " << nHydrogens << std::endl;
			}

			/*******************************/
			/* Scan for ACF output control */
			/*******************************/

			else if (!(parameter.compare("acf")))
			{
				Nt = atoi(value.c_str());
				ti = atof(dimension.c_str());
				tf = atof(toFitStr.c_str());
				acfFound = true;
			}

			/*******************************/
			/* Scan for SPD output control */
			/*******************************/

			else if (!(parameter.compare("spd")))
			{
				Nw = atoi(value.c_str());
				wi = atof(dimension.c_str());
				wf = atof(toFitStr.c_str());
				spdFound = true;
			}

			/************************/
			/* Scan for smallj read */
			/************************/

			else if (!(parameter.compare("smallj")))
				if (!(value.compare("read"))) readSmallJs = true;
		}
	}

	inFile.close();

	/****************************
	 * Check integrity of input *
	 ****************************/

	// WARNINGS

	int warningNumber = 0;
	std::ostringstream warningString;
	warningString << std::endl;
	warningString << "====================================================================================================" << std::endl;
	warningString << "=                                C++OPPS WARNINGS - READ CAREFULLY                                 =" << std::endl;
	warningString << "====================================================================================================" << std::endl;

	if (!dynModelSetByUser &&  !rank)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": model for dynamics not set, defaulting to SRLS." << std::endl;
	}

	if (!(dynModel.compare("srls")) && !proteinAlphaFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": alpha_V not found, defaulting to 0.0 deg." << std::endl;
		proteinD.setAngle("alpha", 0.0);
		fitpar[6] = fitpar2[6] = false;
	}
	if (!(dynModel.compare("srls")) && !proteinBetaFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": beta_V not found, defaulting to 0.0 deg." << std::endl;
		proteinD.setAngle("beta", 0.0);
		fitpar[7] = fitpar2[7] = false;
	}
	if (!(dynModel.compare("srls")) && !proteinGammaFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": gamma_V not found, defaulting to 0.0 deg." << std::endl;
		proteinD.setAngle("gamma", 0.0);
		fitpar[8] = fitpar2[8] = false;
	}

	if (!(dynModel.compare("srls")) && !probeAlphaFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": alpha_O not found, defaulting to 0.0 deg." << std::endl;
		probeD.setAngle("alpha", 0.0);
		fitpar[9] = fitpar2[9] = false;
	}
	if (!(dynModel.compare("srls")) && !probeBetaFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": beta_O not found, defaulting to 0.0 deg." << std::endl;
		probeD.setAngle("beta", 0.0);
		fitpar[10] = fitpar2[10] = false;
	}
	if (!(dynModel.compare("srls")) && !probeGammaFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": gamma_V not found, defaulting to 0.0 deg." << std::endl;
		probeD.setAngle("gamma", 0.0);
		fitpar[11] = fitpar2[11] = false;
	}

	if (!alphaDipFound && (dynModel.compare("3s-fb")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": dipolar_alpha not found, defaulting to 0.0 deg." << std::endl;
		omega_D[0] = 0.0;
		fitpar[12] = fitpar2[12] = false;
	}	
	if (!betaDipFound && (dynModel.compare("3s-fb")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": dipolar_beta not found, defaulting to 0.0 deg." << std::endl;
		omega_D[1] = 0.0;
		fitpar[13] = fitpar2[13] = false;
	}	
	if (!gammaDipFound && (dynModel.compare("3s-fb")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": dipolar_gamma not found, defaulting to 0.0 deg." << std::endl;
		omega_D[2] = 0.0;
		fitpar[14] = fitpar2[14] = false;
	}	

	if (!alphaDip2Found && (dynModel.compare("3s-fb")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": dipolar_alpha2 not found, defaulting to 0.0 deg." << std::endl;
		omega_D2[0] = 0.0;
		////fitpar[12] = fitpar2[12] = false;	
	}	
	if (!betaDip2Found && (dynModel.compare("3s-fb")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": dipolar_beta2 not found, defaulting to 0.0 deg." << std::endl;
		omega_D2[1] = 0.0;
		///fitpar[13] = fitpar2[13] = false;
	}	
	if (!gammaDip2Found && (dynModel.compare("3s-fb")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": dipolar_gamma2 not found, defaulting to 0.0 deg." << std::endl;
		omega_D2[2] = 0.0;
		////fitpar[14] = fitpar2[14] = false;
	}

	if (!alphaCsaFound && (dynModel.compare("3s-fb")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": csa_alpha not found, defaulting to 0.0 deg." << std::endl;
		omega_CSA[0] = 0.0;
		fitpar[15] = fitpar2[15] = false;
	}	
	if (!betaCsaFound && (dynModel.compare("3s-fb")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": csa_beta not found, defaulting to 0.0 deg." << std::endl;
		omega_CSA[1] = 0.0;
		fitpar[16] = fitpar2[16] = false;
	}	
	if (!gammaCsaFound && (dynModel.compare("3s-fb")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": csa_gamma not found, defaulting to 0.0 deg." << std::endl;
		omega_CSA[2] = 0.0;
		fitpar[17] = fitpar2[17] = false;
	}	

	if (!(dynModel.compare("srls")) && !c20Found)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": c20 potential coefficient not found, defaulting to 0.0." << std::endl;
		coeff[(2*(2+1))/2+0] = 0.0;
                fitpar[18] = fitpar2[18] = false;
	}
	if (!(dynModel.compare("srls")) && !c22Found)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": c22 potential coefficient not found, defaulting to 0.0." << std::endl;
		coeff[(2*(2+1))/2+2] = 0.0;
                fitpar[19] = fitpar2[19] = false;
	}
	if (!(dynModel.compare("srls")) && !c40Found)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": c40 potential coefficient not found, defaulting to 0.0." << std::endl;
		coeff[(4*(4+1))/2+0] = 0.0;
                fitpar[20] = fitpar2[20] = false;
	}
	if (!(dynModel.compare("srls")) && !c42Found)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": c42 potential coefficient not found, defaulting to 0.0." << std::endl;
		coeff[(4*(4+1))/2+2] = 0.0;
                fitpar[21] = fitpar2[21] = false;
	}
	if (!(dynModel.compare("srls")) && !c44Found)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": c44 potential coefficient not found, defaulting to 0.0." << std::endl;
		coeff[(4*(4+1))/2+4] = 0.0;
                fitpar[22] = fitpar2[22] = false;
	}

	if (isFB1)
	{
		n_coeff = fb1_coeff.size();
		if (n_coeff <= 0)
		{
			warningNumber++;
			warningString << "WARNING #" << warningNumber << ": no internal potential: motion may be coupled only through diffusion tensor" << std::endl;
		}
	}

	if (!(dynModel.compare("fb2")))
	{
		if (fb2_coeff.dim1 <= 0)
		{
			warningNumber++;
			warningString << "WARNING #" << warningNumber << ": no internal potential set for torsion 1" << std::endl;
		}
		if (fb2_coeff.dim2 <= 0)
		{
			warningNumber++;
			warningString << "WARNING #" << warningNumber << ": no internal potential set for torsion 2" << std::endl;
		}
	}

	if (!rexFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": no exchange rate found, defaulting to 0.0 Hz." << std::endl;
		Rexchange = 0.0;
		fitpar[23] = fitpar2[23] = false;
	}

	if (!nucleusFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": no spin probe defined, defaulting to 15N - 1H." << std::endl;
	}
	
	if (!fieldFound && !fitting)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": no fields found in one-shot calculation input. Setting a default value to 100.0 MHz" << std::endl;
		field.push_back(1.0e8);
	}
	
	if((!LmaxFound || probeLmax <= 0) && !(dynModel.compare("srls")))
	{
		long double dmin = std::min(proteinD.iso(),probeD.iso());
		if (dmin > 1.0e9) probeLmax = (fitting ? 8: 4);
		else if (dmin > 1.0e8) probeLmax = (fitting ? 16: 8);
		else if (dmin > 1.0e7) probeLmax = (fitting ? 20: 10);
		else probeLmax = (fitting ? 30: 25);
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": no probeLmax found. Value has been tentatively set to " << probeLmax << "." << std::endl;
	}
	
	if (nstep <= 0)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": nstep <= 0 for Lanczos. Value will be chosen as the 15% of basis dimensions at runtime" << std::endl;
	}
	
	if (!tolFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": no fit tolerances set. Using the defaul value of 1.0e-5 for xtol, ftol, gtol and 1.0e-2 for epsfcn (see the MINPACK LMDIF.F file as reference)" << std::endl;
	}

	if (!cubatureParamsFound && !(dynModel.compare("fb2")))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber <<": no cubature routine parameters found. Defaulting to maxEval = 5000 and relErr = 1.0e-8. It is recommended to check if this is enough for precision of integrals" << std::endl;
		maxEval = 10000;
		relErr  = 1.0e-10;
	}

	if (!rectFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": no rectification keyword: fitting without X^2 rectification for diffusion tensor" << std::endl;
	}

	if (!orderParametersToleranceFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": no order_parameters_tolerance found. Dafaulting to 1%" << std::endl;
	}
	
	if (!nHydrogensFound)
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": no Hydrogens found. Dafaulting to 1" << std::endl;
	}

	if ((!(dynModel.compare("ts-srls")) || !(dynModel.compare("ts-fb1"))) && (!alphaDip2State1Found || !betaDip2State1Found || !gammaDip2State1Found))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": One or more angles for the dipolar2 interaction in state 1 were not found. Setting to 0.0." << std::endl;
		omegaDip2_1[0] = 0.0; omegaDip2_1[1] = 0.0; omegaDip2_1[2] = 0.0;
	}
	if ((!(dynModel.compare("ts-srls")) || !(dynModel.compare("ts-fb1"))) && (!alphaDip2State2Found || !betaDip2State2Found || !gammaDip2State2Found))
	{
		warningNumber++;
		warningString << "WARNING #" << warningNumber << ": One or more angles for the dipolar2 interaction in state 2 were not found. Setting to 0.0." << std::endl;
		omegaDip2_2[0] = 0.0; omegaDip2_2[1] = 0.0; omegaDip2_2[2] = 0.0;
	}

	warningString << "====================================================================================================" << std::endl;

	// ERRORS

	int errorNumber = 0;
	std::ostringstream errorString;
	
	errorString << std::endl;
	errorString << "====================================================================================================" << std::endl;
	errorString << "=                                           C++OPPS ERRORS                                         =" << std::endl;
	errorString << "====================================================================================================" << std::endl;

	/* SRLS diffusion tensors */

	if (!(dynModel.compare("srls")) && !proteinDxxFound)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": protein_Dxx not found." << std::endl;
	}
	if (!(dynModel.compare("srls")) && !proteinDyyFound)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": protein_Dyy not found." << std::endl;
	}
	if (!(dynModel.compare("srls")) && !proteinDzzFound)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": protein_Dzz not found." << std::endl;
	}

	if (!(dynModel.compare("srls")) && !probeDxxFound)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": probe_Dxx not found." << std::endl;
	}
	if (!(dynModel.compare("srls")) && !probeDyyFound)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": probe_Dyy not found." << std::endl;
	}
	if (!(dynModel.compare("srls")) && !probeDzzFound)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": probe_Dzz not found." << std::endl;
	}

	/* TS-X specific data */
	if (!dynModel.compare("ts-srls") || !dynModel.compare("ts-fb1"))
	{
		if (!hchSigmaFound)
		{
#ifdef __EXPCOPPS__
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": hch angle sigma is required as input." << std::endl;
#else
			hchSigma = 1.0e-13;
			hchSigmaFound = true;
#endif
		}
		if (!populationFound)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": population keyword not found. It is mandatory for TS-SRLS model." << std::endl;
		}
		else if (population < 0.0 || population > 1.0)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": population = " << population << " out of range [0.0, 1.0]." << std::endl;
		}

		if (!jumpFrequencyFound)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": jump_frequency keyword not found. It is mandatory for TS-SRLS model." << std::endl;
		}

		if (!alphaDipState1Found || !betaDipState1Found || !gammaDipState1Found)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": One or more angles for the dipolar interaction in state 1 were not found." << std::endl;
		}
		if (!alphaCsaState1Found || !betaCsaState1Found || !gammaCsaState1Found)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": One or more angles for the CSA interaction in state 1 were not found." << std::endl;
		}
		if (!alphaDipState2Found || !betaDipState2Found || !gammaDipState2Found)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": One or more angles for the dipolar interaction in state 2 were not found." << std::endl;
		}
		if (!alphaCsaState2Found || !betaCsaState2Found || !gammaCsaState2Found)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": One or more angles for the CSA interaction in state 2 were not found." << std::endl;
		}
		if (readSmallJs)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": Interpolation of pre-calculated small j's not yet implemented for TS-SRLS model." << std::endl;
		}
		if (acfFound || spdFound)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": Output of full correlation functions and spectral densities not yet implemented for TS-SRLS model." << std::endl;
		}
	}

	/* FBn diffusion tensor */

	if ( (isFB1 && !fb1DxxFound) || (!(dynModel.compare("fb2")) && !fb2DxxFound) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": dxx not found." << std::endl;
	}
	if ( (isFB1 && !fb1DyyFound) || (!(dynModel.compare("fb2")) && !fb2DyyFound) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": dyy not found." << std::endl;
	}
	if ( (isFB1 && !fb1DzzFound) || (!(dynModel.compare("fb2")) && !fb2DzzFound) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": dzz not found." << std::endl;
	}
	if ( (isFB1 && !fb1DiiFound) || (!(dynModel.compare("fb2")) && !fb2D11Found) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": d11 not found." << std::endl;
	}
	if ( (!(dynModel.compare("fb2")) && !fb2D22Found) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": d22 not found." << std::endl;
	}
	if ( (!(dynModel.compare("fb2")) && !fb2D12Found) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": d12 not found." << std::endl;
	}
	if ( (isFB1 && !fb1DxiFound) || (!(dynModel.compare("fb2")) && !fb2Dx1Found) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": dx1 not found." << std::endl;
	}
	if ( (isFB1 && !fb1DyiFound) || (!(dynModel.compare("fb2")) && !fb2Dy1Found) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": dy1 not found." << std::endl;
	}
	if ( (isFB1 && !fb1DziFound) || (!(dynModel.compare("fb2")) && !fb2Dz1Found) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": dz1 not found." << std::endl;
	}
	if ( (!(dynModel.compare("fb2")) && !fb2Dx2Found) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": dx2 not found." << std::endl;
	}
	if ( (!(dynModel.compare("fb2")) && !fb2Dy2Found) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": dy2 not found." << std::endl;
	}
	if ( (!(dynModel.compare("fb2")) && !fb2Dz2Found) )
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": dz2 not found." << std::endl;
	}

	/* Nmax */

	if (!NmaxFound && isFB1)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": Nmax not found." << std::endl;
	}

	if (!N1maxFound && !(dynModel.compare("fb2")))
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": N1max not found." << std::endl;
	}

	if (!N2maxFound && !(dynModel.compare("fb2")))
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": N2max not found." << std::endl;
	}

	/* 3S-FB checks */
	if (!(dynModel.compare("3s-fb")))
	{
		if (_3SFB_OmegaDip1Site1Found < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Euler angles for Dip1 interaction in site 1 are missing." << std::endl;
		}
		if (_3SFB_OmegaDip1Site2Found < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Euler angles for Dip1 interaction in site 2 are missing." << std::endl;
		}
		if (_3SFB_OmegaDip1Site3Found < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Euler angles for Dip1 interaction in site 3 are missing." << std::endl;
		}

		if (_3SFB_OmegaDip2Site1Found < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Euler angles for Dip2 interaction in site 1 are missing." << std::endl;
		}
		if (_3SFB_OmegaDip2Site2Found < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Euler angles for Dip2 interaction in site 2 are missing." << std::endl;
		}
		if (_3SFB_OmegaDip2Site3Found < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Euler angles for Dip2 interaction in site 3 are missing." << std::endl;
		}

		if (_3SFB_OmegaCSASite1Found < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Euler angles for CSA interaction in site 1 are missing." << std::endl;
		}
		if (_3SFB_OmegaCSASite2Found < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Euler angles for CSA interaction in site 2 are missing." << std::endl;
		}
		if (_3SFB_OmegaCSASite3Found < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Euler angles for CSA interaction in site 3 are missing." << std::endl;
		}

		if (_3SFB_populationFound < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more Boltzmann populations of the sites are missing." << std::endl;
		}
		if (fabs(_3SFB_population[0] + _3SFB_population[1] + _3SFB_population[2] - 1.0) > 0.01)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": the Boltzmann populations do not sum to 1." << std::endl;
		}
		if (fabs(_3SFB_population[0] + _3SFB_population[1] + _3SFB_population[2] - 1.0) > 1.0e-10)
		{
			warningNumber++;
			warningString << "WARNING #" << warningNumber << ": the Boltzmann populations sum to 1 within 1\% difference. Renormalizing...." << std::endl;

			double s = 1.0 / (_3SFB_population[0] + _3SFB_population[1] + _3SFB_population[2]);

			_3SFB_population[0] *= s;
			_3SFB_population[1] *= s;
			_3SFB_population[2] *= s;

			_3SFB_sqPopulation[0] = sqrt(_3SFB_population[0]);
			_3SFB_sqPopulation[1] = sqrt(_3SFB_population[1]);
			_3SFB_sqPopulation[2] = sqrt(_3SFB_population[2]);
			warningString << "          " <<  ": the normalized Boltzmann populations are: " << _3SFB_population[0] << ", " << _3SFB_population[1] << ", " << _3SFB_population[2] << std::endl;
		}
		

		if (_3SFB_jumpFrequencyFound < 3)
		{
			errorNumber++;
			errorString << "ERROR #" << errorNumber << ": one or more jumping frequencies are missing." << std::endl;
		}

	}

	/* Fitting constraints among parameters */

	if (( ((fitpar[12] || fitpar[13] || fitpar[14]) && (fitpar[6] || fitpar[7] || fitpar[8])) || ((fitpar2[12] || fitpar2[13] || fitpar2[14]) && (fitpar2[6] || fitpar2[7] || fitpar2[8])) ) && (constrain_OmegaD_OmegaV || constrain_OmegaV_OmegaD))
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": cannot fit both omegaV and omegaD with constraints." << std::endl;
	}

	if (cdy1 && !consistentD1yConstraint)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": setting [ protein_dxx:FIX & protein_dyy:CONSTRAIN ] is not allowed." << std::endl;
	}

	if (cdz1 && !consistentD1zConstraint)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": setting [ (protein_dxx:FIX | protein_dyy:FIX) & protein_dzz:CONSTRAIN ] is not allowed." << std::endl;
	}

	if (cdy2 && !consistentD2yConstraint)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": setting [ probe_dxx:FIX & probe_dyy:CONSTRAIN ] is not allowed." << std::endl;
	}

	if (cdz2 && !consistentD2zConstraint)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": setting [ (probe_dxx:FIX | probe_dyy:FIX) & probe_dzz:CONSTRAIN ] is not allowed." << std::endl;
	}

	/* Further output control */

	if ((acfFound &&  !spdFound) | (!acfFound &&  spdFound))
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber << ": ACF keyword requires SPD and vice-versa." << std::endl;
	}
	if (readSmallJs && !spdFound)
	{
		errorNumber++;
		errorString << "ERROR #" << errorNumber <<": ACF and SPD keywords are required wheh readind small j's. Also, their ";
		errorString << "paramters must correspond to those emplyed in the generation of the small j files" << std::endl;
	}

	errorString << "====================================================================================================" << std::endl;
	
	if (!rank)
	{
		if (warningNumber > 0)
			std::cout << warningString.str();
		if (errorNumber > 0)
		{
			std::cout << errorString.str();
			exit(0);
		}
	}

	/******************************** 
	 * Scales the diffusion tensors *
	 ********************************/
	
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
		// Scale the diffusion tensors by D2YY
		scale = probeD.getCartesianComponent("dyy");
		proteinD.setScale(scale);
		proteinD.scale();
		probeD.setScale(scale);
		probeD.scale();

		// Scales the jumping frequency by D2YY
		if (!(dynModel.compare("ts-srls"))) jumpFrequency /= scale;

		/********************************************************** 
		 * Calculates the spherical components of protein D in VF *
		 **********************************************************/
		
		proteinD.transform();
		
		/******************************************************** 
		 * Calculates the spherical components of probe D in OF *
		 ********************************************************/
		
		probeD.transform();
	}

	else if (isFB1)
	{
		scale = fb1_diften.getRRComponent("dzz");
		fb1_diften.setScale(scale);
		fb1_diften.scale();

		// Scales the jumping frequency by Dzz
		if (!(dynModel.compare("ts-fb1"))) jumpFrequency /= scale;

		// Scales the jump frequencies in 3S-FB model
		if (!(dynModel.compare("3s-fb")))
		{
			_3SFB_jumpFrequency[0] /= scale;
			_3SFB_jumpFrequency[1] /= scale;
			_3SFB_jumpFrequency[2] /= scale;
		}
	}

	else if (!(dynModel.compare("fb2")))
	{
		scale = fb2_diften.getRRComponent("dzz");
		fb2_diften.setScale(scale);
		fb2_diften.scale();
	}

	Dmin /= (double)scale;
	
	/*********************************************************
	 * Decide which spectral densities have to be calculated *
	 *********************************************************/

	spectralDensities.clear();

	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
		{
		bool isPotentialAxial = !potentialFromExpFile && (fabs(coeff[5]) <= ZERO && (!fitpar[19] || !fitpar2[19])) && (fabs(coeff[12]) <= ZERO &&  (!fitpar[21] || !fitpar2[21])) && (fabs(coeff[14]) <= ZERO && (!fitpar[22] || !fitpar2[22]));
	
		
		if (isPotentialAxial)
		{
			spectralDensities.push_back(12);
			spectralDensities.push_back(18);
			spectralDensities.push_back(24);
		}
		else
		{
			spectralDensities.push_back(12);
			spectralDensities.push_back(18);
			spectralDensities.push_back(24);
			spectralDensities.push_back(4);
			spectralDensities.push_back(8);
			spectralDensities.push_back(14);
			spectralDensities.push_back(9);
			spectralDensities.push_back(13);
			spectralDensities.push_back(19);
		}
	}

	else if ( isFB1 || !(dynModel.compare("fb2")) )
	{
		spectralDensities.push_back(12);
		spectralDensities.push_back(18);
		spectralDensities.push_back(24);
		spectralDensities.push_back(4);
		spectralDensities.push_back(8);
		spectralDensities.push_back(14);
		spectralDensities.push_back(9);
		spectralDensities.push_back(13);
		spectralDensities.push_back(19);
	}

	/****************************************/
	/* Create the _tw.par file if requested */
	/****************************************/

	if (acfFound)
	{
		OUT_CONTROL = 1;
		std::fstream file;
		position1 = inFileName.find_last_of("_");
		std::string fileName = inFileName.substr(0,position1) + "_tw.par";
		std::cout << "TASK " << rank << ": Creating " << fileName << " file for ACF/SPD functions" << std::endl;
		file.open(fileName.c_str(),std::ios::out);
		if (file.is_open())
		{
			file << Nt << std::endl << ti << std::endl << tf << std::endl;
			file << Nw << std::endl << wi << std::endl << wf << std::endl;
			file.close();
		}
		else
		{
			std::cout << std::endl << "TASK " << rank << ": ERROR: cannot create " << fileName << " file" << std::endl << std::endl;
			exit(1);
		}
	}

	/*******************************************************/
	/* Generates the structure with potential coefficients */
	/*******************************************************/

	Ucoeff.npop = 1;
	Ucoeff.c20.clear(); Ucoeff.c20.push_back(coeff[3]);
	Ucoeff.c22.clear(); Ucoeff.c22.push_back(coeff[5]);
	Ucoeff.c40.clear(); Ucoeff.c40.push_back(coeff[10]);
	Ucoeff.c42.clear(); Ucoeff.c42.push_back(coeff[12]);
	Ucoeff.c44.clear(); Ucoeff.c44.push_back(coeff[14]);

	/***********************************************************/
	/* Determins if potential coefficients are going to be fit */
	/***********************************************************/

	potentialFit = false;
	if (!dynModel.compare("srls") && ( (fitpar[18] || fitpar[19] || fitpar[20] || fitpar[21] || fitpar[22]) || (fitpar2[18] || fitpar2[19] || fitpar2[20] || fitpar2[21] || fitpar2[22]) )&& !potentialFromExpFile)
		potentialFit = true;

	return;
}

/*************************
 *  Write out parameters *
 *************************/

std::string physics::toString(void){

	std::ostringstream ostr;
	
	ostr << "TASK " << rank << ": ********************************************" << std::endl;
	ostr << outputString.str();

	if (readSmallJs)
	{
		ostr << "TASK " << rank << ": ********************************************" << std::endl;
		ostr << "TASK " << rank << ": Small j's will be red from a previous" << std::endl;
		ostr << "TASK " << rank << ": calculation ran with the ACF/SPD keywords" << std::endl;
	}


	if (acfFound)
	{
		ostr << "TASK " << rank << ": ********************************************" << std::endl;
		ostr << "TASK " << rank << ":   Autocorrelation functions will be calculated with Nt =  " << Nt << ", ti = " << ti << ", tf = " << tf << std::endl;
		ostr << "TASK " << rank << ":   Spectral densities will be calculated with Nw =  " << Nt << ", wi = " << ti << ", wf = " << tf << std::endl;
	}

	
	int i = 0;
	if (fitting)
	{
		ostr << "TASK " << rank << ": ********************************************" << std::endl;
		ostr << "TASK " << rank << ": Parameters to fit:" << std::endl;
		if ( (fitpar[0] && fitpar[1] && fitpar[2]) || (fitpar2[0] && fitpar2[1] && fitpar2[2]) )
		{
			ostr << "TASK " << rank << ": " << ++i << ") Protein Dxx" << std::endl;
			ostr << "TASK " << rank << ": " << ++i << ") Protein Dyy" << std::endl;
			ostr << "TASK " << rank << ": " << ++i << ") Protein Dzz" << std::endl;
		}
		else if (cdy1 && !cdz1)
		{
			ostr << "TASK " << rank << ": " << ++i << ") Protein Dxx = Dyy" << std::endl;
			if (fitpar[2] || fitpar2[2])
				ostr << "TASK " << rank << ": " << ++i << ") Protein Dzz" << std::endl;

		}
		else if (cdz1)
			ostr << "TASK " << rank << ": " << ++i << ") Trace of protein D" << std::endl;
		else
		{
			if (fitpar[1] || fitpar2[1])
				ostr << "TASK " << rank << ": " << ++i << ") Protein Dxx" << std::endl;
			if (fitpar[0] || fitpar2[0])
				ostr << "TASK " << rank << ": " << ++i << ") Protein Dyy" << std::endl;
			if (fitpar[2] || fitpar2[2])
				ostr << "TASK " << rank << ": " << ++i << ") Protein Dzz" << std::endl;
		}
		if ( (fitpar[3] && fitpar[4] && fitpar[5]) || (fitpar2[3] && fitpar2[4] && fitpar2[5]) )
		{
			ostr << "TASK " << rank << ": " << ++i << ") Probe Dxx" << std::endl;
			ostr << "TASK " << rank << ": " << ++i << ") Probe Dyy" << std::endl;
			ostr << "TASK " << rank << ": " << ++i << ") Probe Dzz" << std::endl;
		}
		else if (cdy2 && !cdz2)
		{
			ostr << "TASK " << rank << ": " << ++i << ") Probe Dxx = Dyy" << std::endl;
			if (fitpar[5] || fitpar2[5])
				ostr << "TASK " << rank << ": " << ++i << ") Probe Dzz" << std::endl;

		}
		else if (cdz2)
			ostr << "TASK " << rank << ": " << ++i << ") Trace of probe D" << std::endl;
		else
		{
			if (fitpar[4] || fitpar2[4])
				ostr << "TASK " << rank << ": " << ++i << ") Probe Dxx" << std::endl;
			if (fitpar[3] || fitpar2[3])
				ostr << "TASK " << rank << ": " << ++i << ") Probe Dyy" << std::endl;
			if (fitpar[5] || fitpar2[5])
				ostr << "TASK " << rank << ": " << ++i << ") Probe Dzz" << std::endl;
		}
		if (fitpar[6]  || fitpar2[6])  ostr << "TASK " << rank << ": " << ++i << ") Alpha M1F -> VF"<< std::endl;
		if (fitpar[7]  || fitpar2[7])  ostr << "TASK " << rank << ": " << ++i << ") Beta  M1F -> VF"<< std::endl;
		if (fitpar[8]  || fitpar2[8])  ostr << "TASK " << rank << ": " << ++i << ") Gamma M1F -> VF"<< std::endl;
		if (fitpar[9]  || fitpar2[9])  ostr << "TASK " << rank << ": " << ++i << ") Alpha M2F -> OF"<< std::endl;
		if (fitpar[10] || fitpar2[10]) ostr << "TASK " << rank << ": " << ++i << ") Beta  M2F -> OF"<< std::endl;
		if (fitpar[11] || fitpar2[11]) ostr << "TASK " << rank << ": " << ++i << ") Gamma M2F -> OF"<< std::endl;
		if (fitpar[12] || fitpar2[12]) ostr << "TASK " << rank << ": " << ++i << ") Alpha OF -> DF"<< std::endl;
		if (fitpar[13] || fitpar2[13]) ostr << "TASK " << rank << ": " << ++i << ") Beta  OF -> DF"<< std::endl;
		if (fitpar[14] || fitpar2[14]) ostr << "TASK " << rank << ": " << ++i << ") Gamma OF -> DF"<< std::endl;
		if (fitpar[15] || fitpar2[15]) ostr << "TASK " << rank << ": " << ++i << ") Alpha DF -> CF"<< std::endl;
		if (fitpar[16] || fitpar2[16]) ostr << "TASK " << rank << ": " << ++i << ") Beta  DF -> CF"<< std::endl;
		if (fitpar[17] || fitpar2[17]) ostr << "TASK " << rank << ": " << ++i << ") Gamma DF -> CF"<< std::endl;
		if (fitpar[26] || fitpar2[26]) ostr << "TASK " << rank << ": " << ++i << ") Alpha DF -> D2F"<< std::endl;
		if (fitpar[27] || fitpar2[27]) ostr << "TASK " << rank << ": " << ++i << ") Beta  DF -> D2F"<< std::endl;
		if (fitpar[28] || fitpar2[28]) ostr << "TASK " << rank << ": " << ++i << ") Gamma DF -> D2F"<< std::endl;
		if (fitpar[18] || fitpar2[18]) ostr << "TASK " << rank << ": " << ++i << ") c20"<< std::endl;
		if (fitpar[19] || fitpar2[19])
		{
			if (ratio22)
				ostr << "TASK " << rank << ": " << ++i << ") c22/c20"<< std::endl;
			else
				ostr << "TASK " << rank << ": " << ++i << ") c22"<< std::endl;
		}
		if (fitpar[20] || fitpar2[20]) ostr << "TASK " << rank << ": " << ++i << ") c40"<< std::endl;
		if (fitpar[21] || fitpar2[21])
		{
			if (ratio42)
				ostr << "TASK " << rank << ": " << ++i << ") c42/c40"<< std::endl;
			else
				ostr << "TASK " << rank << ": " << ++i << ") c42"<< std::endl;
		}
		if (fitpar[22] || fitpar2[22])
		{
			if (ratio44)
				ostr << "TASK " << rank << ": " << ++i << ") c44/c40"<< std::endl;
			else
				ostr << "TASK " << rank << ": " << ++i << ") c44"<< std::endl;
		}
		if (fitpar[23] || fitpar2[23]) ostr << "TASK " << rank << ": " << ++i << ") R exchange"<< std::endl;
		if (fitpar[25] || fitpar2[25]) ostr << "TASK " << rank << ": " << ++i << ") D correction"<< std::endl;
		if (fitpar[29] || fitpar2[29]) ostr << "TASK " << rank << ": " << ++i << ") State 1 population" << std::endl;
		if (fitpar[30] || fitpar2[30]) ostr << "TASK " << rank << ": " << ++i << ") Jump frequency" << std::endl;
		if (fitpar[31] || fitpar2[31]) ostr << "TASK " << rank << ": " << ++i << ") hch sigma" << std::endl;
	}

	ostr << "TASK " << rank << ": ********************************************" << std::endl;

	return ostr.str();
}

/*************************************************
 * Obtain a specific parameter in the list below *
 *************************************************/

long double physics::getValueOf(std::string par){
	if (!(par.compare("protein_dxx")))
		return proteinD.getCartesianComponent("dxx");
	else if (!(par.compare("protein_dyy")))
		return proteinD.getCartesianComponent("dyy");
	else if (!(par.compare("protein_dzz")))
		return proteinD.getCartesianComponent("dzz");
	else if (!(par.compare("protein_alpha")))
		return proteinD.getAngle("alpha");
	else if (!(par.compare("protein_beta")))
		return proteinD.getAngle("beta");
	else if (!(par.compare("protein_gamma")))
		return proteinD.getAngle("gamma");
	else if (!(par.compare("probe_dxx")))
		return probeD.getCartesianComponent("dxx");
	else if (!(par.compare("probe_dyy")))
		return probeD.getCartesianComponent("dyy");
	else if (!(par.compare("probe_dzz")))
		return probeD.getCartesianComponent("dzz");
	else if (!(par.compare("probe_alpha")))
		return probeD.getAngle("alpha");
	else if (!(par.compare("probe_beta")))
		return probeD.getAngle("beta");
	else if (!(par.compare("probe_gamma")))
		return probeD.getAngle("gamma");
	else if (!(par.compare("c20")))
		return coeff[3];
	else if (!(par.compare("c22")))
		return coeff[5];
	else if (!(par.compare("c40")))
		return coeff[10];
	else if (!(par.compare("c42")))
		return coeff[12];
	else if (!(par.compare("c44")))
		return coeff[14];
	else if (!(par.compare("dipolar_alpha")))
		return omega_D[0];
	else if (!(par.compare("dipolar_beta")))
		return omega_D[1];
	else if (!(par.compare("dipolar_gamma")))
		return omega_D[2];
	else if (!(par.compare("dipolar_alpha2")))
		return omega_D2[0];
	else if (!(par.compare("dipolar_beta2")))
		return omega_D2[1];
	else if (!(par.compare("dipolar_gamma2")))
		return omega_D2[2];
	else if (!(par.compare("csa_alpha")))
		return omega_CSA[0];
	else if (!(par.compare("csa_beta")))
		return omega_CSA[1];
	else if (!(par.compare("csa_gamma")))
		return omega_CSA[2];
	else if (!(par.compare("Rexchange")))
		return Rexchange;
	else if (!(par.compare("fit_tolerance_1")))
		return fitTol1;
	else if (!(par.compare("fit_tolerance_2")))
		return fitTol2;
	// FB1 parameters
	else if (!(par.compare("dxx")))
		return fb1_diften.getRRComponent("dxx");
	else if (!(par.compare("dyy")))
		return fb1_diften.getRRComponent("dyy");
	else if (!(par.compare("dzz")))
		return fb1_diften.getRRComponent("dzz");
	else if (!(par.compare("dxi")))
		return fb1_diften.getRIComponent("dxi");
	else if (!(par.compare("dyi")))
		return fb1_diften.getRIComponent("dyi");
	else if (!(par.compare("dzi")))
		return fb1_diften.getRIComponent("dzi");
	else if (!(par.compare("dii")))
		return fb1_diften.getIIComponent();
	
	return 0.0;

}

/*****************************************
 * Set the value of a specific parameter *
 *****************************************/

void physics::setValueOf(std::string par, long double val){
	if (!(par.compare("protein_dxx")))
		proteinD.setCartesianComponent("dxx",val);
	else if (!(par.compare("protein_dyy")))
		proteinD.setCartesianComponent("dyy",val);
	else if (!(par.compare("protein_dzz")))
		proteinD.setCartesianComponent("dzz",val);
	else if (!(par.compare("protein_alpha")))
		proteinD.setAngle("alpha",val);
	else if (!(par.compare("protein_beta")))
		proteinD.setAngle("beta",val);
	else if (!(par.compare("protein_gamma")))
		proteinD.setAngle("gamma",val);
	else if (!(par.compare("probe_dxx")))
		probeD.setCartesianComponent("dxx",val);
	else if (!(par.compare("probe_dyy")))
		probeD.setCartesianComponent("dyy",val);
	else if (!(par.compare("probe_dzz")))
		probeD.setCartesianComponent("dzz",val);
	else if (!(par.compare("probe_alpha")))
		probeD.setAngle("alpha",val);
	else if (!(par.compare("probe_beta")))
		probeD.setAngle("beta",val);
	else if (!(par.compare("probe_gamma")))
		probeD.setAngle("gamma",val);
	else if (!(par.compare("c20")))
		coeff[3] = val;
	else if (!(par.compare("c22")))
		coeff[5] = val;
	else if (!(par.compare("c40")))
		coeff[10] = val;
	else if (!(par.compare("c42")))
		coeff[12] = val;
	else if (!(par.compare("c44")))
		coeff[14] = val;
	else if (!(par.compare("dipolar_alpha")))
		omega_D[0] = val;
	else if (!(par.compare("dipolar_beta")))
		omega_D[1] = val;
	else if (!(par.compare("dipolar_gamma")))
		omega_D[2] = val;
	else if (!(par.compare("dipolar_alpha2")))
		omega_D2[0] = val;
	else if (!(par.compare("dipolar_beta2")))
		omega_D2[1] = val;
	else if (!(par.compare("dipolar_gamma2")))
		omega_D2[2] = val;
	else if (!(par.compare("csa_alpha")))
		omega_CSA[0] = val;
	else if (!(par.compare("csa_beta")))
		omega_CSA[1] = val;
	else if (!(par.compare("csa_gamma")))
		omega_CSA[2] = val;
	else if (!(par.compare("Rexchange")))
		Rexchange = val;
	else if (!(par.compare("fit_tolerance_1")))
		fitTol1 = val;
	else if (!(par.compare("fit_tolerance_2")))
		fitTol2 = val;
	// FB1 parameters
	else if (!(par.compare("dxx")))
		fb1_diften.setRRComponent("dxx",val);
	else if (!(par.compare("dyy")))
		fb1_diften.setRRComponent("dyy",val);
	else if (!(par.compare("dzz")))
		fb1_diften.setRRComponent("dzz",val);
	else if (!(par.compare("dxi")))
		fb1_diften.setRIComponent("dxi",val);
	else if (!(par.compare("dyi")))
		fb1_diften.setRIComponent("dyi",val);
	else if (!(par.compare("dzi")))
		fb1_diften.setRIComponent("dzi",val);
	else if (!(par.compare("dii")))
		fb1_diften.setIIComponent(val);
	return;
}

/***************************************************/
/* Get/set coefficients from 'potential' structure */
/***************************************************/

potential physics::getPotentialCoefficients(void)
{
	Ucoeff.npop = 1;
	Ucoeff.c20.clear(); Ucoeff.c20.push_back(coeff[3]);
	Ucoeff.c22.clear(); Ucoeff.c22.push_back(coeff[5]);
	Ucoeff.c40.clear(); Ucoeff.c40.push_back(coeff[10]);
	Ucoeff.c42.clear(); Ucoeff.c42.push_back(coeff[12]);
	Ucoeff.c44.clear(); Ucoeff.c44.push_back(coeff[14]);
	return Ucoeff;
}

void physics::setPotentialCoefficients(potential u)
{
	Ucoeff    = u;
	coeff[3]  = u.c20.at(0);
	coeff[5]  = u.c22.at(0);
	coeff[10] = u.c40.at(0);
	coeff[12] = u.c42.at(0);
	coeff[14] = u.c44.at(0);
	return;
}

/******************************
 * Return data of not 1H atom *
 ******************************/

dvector physics::getNuclearData(void)
{
	dvector v;
	v.push_back(gyromag);
	v.push_back(bondLength);
	v.push_back(deltaCSA);
	return v;
}

/************************************************ 
 * Return all spherical components of protein D *
 ************************************************/

cldvector physics::getProteinDSph(void)
{
	return proteinD.getSphericalComponents();
}

/**********************************************
 * Return all spherical components of probe D *
 **********************************************/

cldvector physics::getProbeDSph(void)
{
	return probeD.getSphericalComponents();
}

/****************************************************************************
 * Return the number of coefficients of the expansion of coupling potential *
 ****************************************************************************/

int physics::getNCoeff(void)
{
	return n_coeff;
}

/************************************************************************************
 * Return a vector containing the coefficients of the potential following the rule  *
 * c(l,m) = l*(l+1)/2 + m,    for m >=0                                             *
 * if m < 0 the coefficient is equal to c(l,-m)*                                    *
 ************************************************************************************/

ldvector physics::getCoeff(void)
{
	int size = 15;
	ldvector c;
	for (int i = 0; i < size; i++)
		c.push_back(coeff[i]);
	return c;
}

/*********************************************************
 * Return the Euler angles that troansform from OF to DF *
 *********************************************************/

ldvector physics::getOmegaDipolar()
{
	ldvector od;
	od.push_back(omega_D[0]);
	od.push_back(omega_D[1]);
	od.push_back(omega_D[2]);
	return od;
}

ldvector physics::getOmegaDipolar2()
{
	ldvector od;
	od.push_back(omega_D2[0]);
	od.push_back(omega_D2[1]);
	od.push_back(omega_D2[2]);
	return od;
}

/*********************************************************
 * Return the Euler angles that troansform from DF to CF *
 *********************************************************/

ldvector physics::getOmegaCSA(void)
{
	ldvector ocsa;
	ocsa.push_back(omega_CSA[0]);
	ocsa.push_back(omega_CSA[1]);
	ocsa.push_back(omega_CSA[2]);
	return ocsa;	
}

/**********************************************************************
 * Return the Lmax for the expansion of the rotational space of probe *
 **********************************************************************/

int physics::getProbeLMax(void)
{
	return probeLmax;
}

/***************************
 * Return the fitting flag *
 ***************************/

bool physics::getFitting(void)
{
	return fitting;
}

/*****************************
 * Return the fitting method *
 *****************************/

std::string physics::getFittingMethod(void)
{
	return fitMethod;
}

/*****************************************************************************
 * Return a vector indicating which spectral densities have to be calculated *
 *****************************************************************************/

ivector physics::getSpectralDensities(void)
{
	return spectralDensities;
}

/**********************************************
 * Return the scaling constant of frequencies *
 **********************************************/

long double physics::getScale(void)
{
	return scale;
}

/************************************************
 * Return the number of parameters to be fitted *
 ************************************************/

int physics::getNFit(void)
{
	return nfit;
}

int physics::getNFit2(void)
{
	return nfit2;
}

/************************************************************************************************
 * Return the ordered list of booleans indicating if every parameter have, or not, to be fitted *
 ************************************************************************************************/

bool* physics::getFitPar(void)
{
	bool* f = fitpar;
	return f;
}

bool* physics::getFitPar2(void)
{
	bool* f = fitpar2;
	return f;
}

/*********************************
 * Return the vector with fields *
 *********************************/

dvector physics::getField(void)
{
	return field;
}

/************************************
 * Return the max nstep for Lanczos *
 ************************************/

int physics::getLanczosNStep(void)
{
	return nstep;
}

/********************************************
 * Return the string of the dynamics model *
 *******************************************/

std::string physics::getDynamicsModel(void)
{
	return dynModel;
}

/***********************
 * FB1 related methods *
 ***********************/

int physics::getNmax(void)
{
	return Nmax;
}

int physics::getFB1NCoeff(void)
{
	return fb1_coeff.size();
}

dcvector physics::getFB1Coeff(void)
{
	return fb1_coeff;
}

dvector physics::getFB1Diften(void)
{
	dvector d;
	d.push_back((double)fb1_diften.getRRComponent("dxx"));
	d.push_back((double)fb1_diften.getRRComponent("dyy"));
	d.push_back((double)fb1_diften.getRRComponent("dzz"));
	d.push_back((double)fb1_diften.getIIComponent());
	d.push_back((double)fb1_diften.getRIComponent("dxi"));
	d.push_back((double)fb1_diften.getRIComponent("dyi"));
	d.push_back((double)fb1_diften.getRIComponent("dzi"));
	return d;
}

void physics::scaleFB1diften(long double s)
{
	fb1_diften.scale(s);
}

/***********************/
/* FB2 related methods */
/***********************/

int physics::getNmax(int t)
{
	switch (t)
	{
		case 1:
			return N1max;
		case 2:
			return N2max;
	}
	return -1;
}
potential_fb2 physics::getFB2Coeff()
{
	return fb2_coeff;
}

int physics::getFB2CoeffDim(int t)
{
	int z = 0;
	switch (t)
	{
		case 1:
			return fb2_coeff.dim1;
		case 2:
			return fb2_coeff.dim2;
	}
	return z;
}

dcomplex physics::getFB2CoeffAt(int n1, int n2)
{
	dcomplex c = dcomplex(0.0, 0.0);
	if (n1 >= -fb2_coeff.n1Max && n1 <= fb2_coeff.n1Max && n2 >= -fb2_coeff.n2Max && n2 <= fb2_coeff.n2Max)
		c = fb2_coeff.c[(n1 + fb2_coeff.n1Max) * fb2_coeff.dim2 + (n2 + fb2_coeff.n2Max)];
	return c;
}

dvector physics::getFB2Diften(void)
{
	dvector d;
	d.push_back((double)fb2_diften.getRRComponent("dxx"));
	d.push_back((double)fb2_diften.getRRComponent("dyy"));
	d.push_back((double)fb2_diften.getRRComponent("dzz"));
	d.push_back((double)fb2_diften.getIIComponent("d11"));
	d.push_back((double)fb2_diften.getIIComponent("d22"));
	d.push_back((double)fb2_diften.getIIComponent("d12"));
	d.push_back((double)fb2_diften.getRIComponent("dx1"));
	d.push_back((double)fb2_diften.getRIComponent("dy1"));
	d.push_back((double)fb2_diften.getRIComponent("dz1"));
	d.push_back((double)fb2_diften.getRIComponent("dx2"));
	d.push_back((double)fb2_diften.getRIComponent("dy2"));
	d.push_back((double)fb2_diften.getRIComponent("dz2"));
	return d;
}

void physics::scaleFB2diften(long double s)
{
	fb2_diften.scale(s);
}
/**********************************************
 * Return the values of OmegaD before fitting *
 **********************************************/
long double physics::getBkpDipolarAngle(std::string a)
{
	if (!(a.compare("alpha")))
		return omega_D_bkp[0];
	else if (!(a.compare("beta")))
		return omega_D_bkp[1];
	else if (!(a.compare("gamma")))
		return omega_D_bkp[2];
	else
	{
		std::cout << std::endl <<  "ERROR: " << a << " not recognized in physics::getBkpDipolarAngle. Choose among alpha, beta and gamma." << std::endl << std::endl;
		exit(1);
	}
	return 0.0;
}
/***************************************************/
/* Return booleans for OmegaV - OmegaD constraints */
/***************************************************/
bool physics::getConstrainVD()
{
	return constrain_OmegaV_OmegaD;
}
bool physics::getConstrainDV()
{
	return constrain_OmegaD_OmegaV;
}

/********************************/
/* Return info on rectification */
/********************************/
double physics::getDmin()
{
	return Dmin;
}
double physics::getKrect()
{
	return Krect;
}

/**************************************/
/* Order parameters tolerance methods */
/**************************************/
double physics::getOrderParametersTolerance()
{
	return order_parameters_tolerance;
}
void physics::setOrderParametersTolerance(double d)
{
	order_parameters_tolerance = d;
	return;
}

/*************************************************/
/* Get the number of hydrogen atoms in the probe */
/*************************************************/
int physics::getNHydrogens()
{
	return nHydrogens;
}

/****************************************************************************/
/* Return the number of points for the discretization of spectral densities */
/****************************************************************************/
int physics::getNw()
{
	return Nw;
}

/*************************************************************/
/* Flag indicating if potential is read from _copps.exp file */
/*************************************************************/
bool physics::getPotentialFromExpFile(void)
{
        return potentialFromExpFile;
}

bool physics::isPotentialFit(void)
{
	return potentialFit;
}

/******************************/
/* Return cubature parameters */
/******************************/
void physics::getCubatureParams(int *me, double *re)
{
	*me = maxEval;
	*re = relErr;
	return;
}

/*******************/
/* TS-SRLS methods */
/*******************/

dvector physics::getOmegaDipolarState(int s)
{
	dvector v(3,0.0);
	switch (s)
	{
		case 1:
		{
			v.at(0) = omegaDip_1[0];
			v.at(1) = omegaDip_1[1];
			v.at(2) = omegaDip_1[2];
			break;
		}
		case 2:
		{
			v.at(0) = omegaDip_2[0];
			v.at(1) = omegaDip_2[1];
			v.at(2) = omegaDip_2[2];
			break;
		}
	}
	return v;
}

dvector physics::getOmegaCsaState(int s)
{
	dvector v(3,0.0);
	switch (s)
	{
		case 1:
		{
			v.at(0) = omegaCsa_1[0];
			v.at(1) = omegaCsa_1[1];
			v.at(2) = omegaCsa_1[2];
			break;
		}
		case 2:
		{
			v.at(0) = omegaCsa_2[0];
			v.at(1) = omegaCsa_2[1];
			v.at(2) = omegaCsa_2[2];
			break;
		}
	}
	return v;
}

dvector physics::getOmegaDipolar2State(int s)
{
	dvector v(3,0.0);
	switch (s)
	{
		case 1:
		{
			v.at(0) = omegaDip2_1[0];
			v.at(1) = omegaDip2_1[1];
			v.at(2) = omegaDip2_1[2];
			break;
		}
		case 2:
		{
			v.at(0) = omegaDip2_2[0];
			v.at(1) = omegaDip2_2[1];
			v.at(2) = omegaDip2_2[2];
			break;
		}
	}
	return v;
}

double physics::getPopulation(void)
{
	return population;
}
void physics::setPopulation(double p)
{
	population = p; return;
}

double physics::getJumpFrequency(void)
{
	return jumpFrequency;
}
void physics::setJumpFrequency(double j)
{
	jumpFrequency = j; return;
}

/* 3S-FB */

dvector physics::_3SFB_getOmegaDipolar1Site(int s) // site is numbered 1-based
{
	s = s-1;
	dvector v;
	v.push_back(_3SFB_OmegaDip1[s][0]);
	v.push_back(_3SFB_OmegaDip1[s][1]);
	v.push_back(_3SFB_OmegaDip1[s][2]);
	return v;
}

dvector physics::_3SFB_getOmegaDipolar2Site(int s)
{
	s = s-1;
	dvector v;
	v.push_back(_3SFB_OmegaDip2[s][0]);
	v.push_back(_3SFB_OmegaDip2[s][1]);
	v.push_back(_3SFB_OmegaDip2[s][2]);
	return v;
}

dvector physics::_3SFB_getOmegaCSASite(int s)
{
	s = s-1;
	dvector v;
	v.push_back(_3SFB_OmegaCSA[s][0]);
	v.push_back(_3SFB_OmegaCSA[s][1]);
	v.push_back(_3SFB_OmegaCSA[s][2]);
	return v;
}

double  physics::_3SFB_getPopulation(int s)
{
	return _3SFB_population[s];
}

double  physics::_3SFB_getSqPopulation(int s)
{
	return _3SFB_sqPopulation[s];
}

dvector physics::_3SFB_getSqPopulations(void)
{
	dvector v;
	v.push_back(_3SFB_sqPopulation[0]);
	v.push_back(_3SFB_sqPopulation[1]);
	v.push_back(_3SFB_sqPopulation[2]);
	return v;
}

void    physics::_3SFB_setPopulation(int s, double p)
{
	_3SFB_population[s] = p;
	_3SFB_sqPopulation[s] = p;
	return;
}

double  physics::_3SFB_getJumpFrequency(int j)
{
	return _3SFB_jumpFrequency[j];
}

dvector physics::_3SFB_getJumpFrequencies(void)
{
	dvector v;
	v.push_back(_3SFB_jumpFrequency[0]);
	v.push_back(_3SFB_jumpFrequency[1]);
	v.push_back(_3SFB_jumpFrequency[2]);
	return v;
}

void    physics::_3SFB_setJumpFrequency(int j, double w)
{
	_3SFB_jumpFrequency[j] = w;
	return;
}

///////////////////////////////////////////
// ONLY TEMPORARY
void physics::setDB2(double d)
{
	omegaDip2_1[1] = omegaDip2_2[1] = d;
	return;
}

double physics::getHchSigma(void)
{
	return hchSigma;
}
void physics::setHchSigma(double s)
{
	hchSigma = s;
	return;
}
///////////////////////////////////////////

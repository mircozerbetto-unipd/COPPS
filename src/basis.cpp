/***********************************************************************************
 * C++OPPS 2.2 - Interpretation of NMR relaxation in proteins			   *
 * Copyright (C) 2008  Mirco Zerbetto						   * 
 * 										   *
 * This program is free software; you can redistribute it and/or	 	   *
 * modify it under the terms of the GNU General Public License			   *
 * as published by the Free Software Foundation; either version 2		   *
 * of the License, or any later version.					   *
 *										   *
 * This program is distributed in the hope that it will be useful,		   *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of		   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		   *
 * GNU General Public License for more details.					   *
 * 										   *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 ***********************************************************************************
 * Author: Mirco Zerbetto							   *
 * Dipartimento di Scienze Chimiche - Universita' di Padova - Italy		   *
 * E-mail: mirco.zerbetto@unipd.it						   *
 ***********************************************************************************/

/*
 ============================================================================
 Name        : basis.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class to determine, also with symmetries, the basis indexes
 ============================================================================
 */

#include "basis.h"

// Constructor
basis::basis()
{
}

// Destructor
basis::~basis()
{
#ifdef WRITE_DESTROY_MESSAGE
	std::cout << "TASK " << rank << ": Cleared basis object" << std::endl;
#endif
}

// Initiator: set the boundaries to the quantum numbers
void basis::init(int L1min, int L1max, int M1min, int M1max, int K1min, int K1max,
		int L2min, int L2max, int M2min, int M2max, int K2min, int K2max)
{

	rank = 0;
	
	std::cout << "TASK " << rank << "Created basis object" << std::endl;
	
	this->L1min = L1min;
	this->L1max = L1max;
	this->M1min = M1min;
	this->M1max = M1max;
	this->K1min = K1min;
	this->K1max = K1max;
	
	this->L2min = L2min;
	this->L2max = L2max;
	this->M2min = M2min;
	this->M2max = M2max;
	this->K2min = K2min;
	this->K2max = K2max;
	
	this->Mmmin = K1min-M2max;
	this->Mmmax = K1max-M2min;
	
	this->Mpmin = K1min+M2min;
	this->Mpmax = K1max+M2max;
	
	jjmin = -1;
	jjmax =  1;
	
	nbf = 0;
	
	return;

}

void basis::init(int L1min, int L1max, int M1min, int M1max, int K1min, int K1max,
		int L2min, int L2max, int M2min, int M2max, int K2min, int K2max, int r)
{
	
	rank = r;
	
	this->L1min = L1min;
	this->L1max = L1max;
	this->M1min = M1min;
	this->M1max = M1max;
	this->K1min = K1min;
	this->K1max = K1max;
	
	this->L2min = L2min;
	this->L2max = L2max;
	this->M2min = M2min;
	this->M2max = M2max;
	this->K2min = K2min;
	this->K2max = K2max;
	
	this->Mmmin = K1min-M2max;
	this->Mmmax = K1max-M2min;
	
	this->Mpmin = K1min+M2min;
	this->Mpmax = K1max+M2max;
	
	jjmin = -1;
	jjmax =  1;
	
	nbf = 0;
	
	return;

}

// Sets limits on quantum numbers basing on symmetry

void basis::checkSymmetry(ldvector data){
	
	/* The entries of vectord data are
	 * data[0] = alpha_V
	 * data[1] = beta_V
	 * data[2] = gamma_V
	 * data[3] = D1xx - D1yy
	 * data[4] = abs[im{D1(2,1)}] + abs[im{D1(2,2)}] + abs[im{D2(2,1)}] + abs[im{D2(2,2)}]
	 * data[5] = abs(c22) + abs(c42) + abs(c44)
	 */
	
	
	/* 
	 * There is no global orienting potential
	 * so L1 and M1 are constrained to the values of the physical observable.
	 * Here we study only observables that transform as D20K, so we fix
	 * L1 = 2 and M1 = 0
	 */
	L1max = 2;
	L1min = 2;
	M1min = 0;
	M1max = 0;
	
	/*
	 * In case of no potential tilt and axially simmetryc diffusion tensor
	 * of the protein, Mm is diagonal and only basis functions with Mm = 0
	 * contribute to the spectral density
	 */
	
	if (data.size() < 4)
		std::cout << "TASK " << rank << ": WARNING : wrong number of arguments in basis::checkSymmetry (found " << data.size() << " instead of 4). Cannot determine Mm symmetry." << std::endl;
	else
	if  ((fabs(data[0])<1.0e-13 && fabs(data[1])<1.0e-13 && fabs(data[2])<1.0e-13 && fabs(data[3])<1.0e-13))
		Mmmin = Mmmax = 0;
	else  // TO CHECK
	{
		Mmmin = -2;
		Mmmax =  2;
	}
	
	if (data.size() < 5)
		std::cout << "TASK " << rank << ": WARNING : wrong number of arguments in basis::checkSymmetry (found " << data.size() << " instead of 5). Cannot determine j symmetry." << std::endl;
	else
	if  (data[4]<1.0e-13)
		;//jjmin = jjmax = 1;

	if (data.size() < 6)
		std::cout << "TASK " << rank << ": WARNING : wrong number of arguments in basis::checkSymmetry (found " << data.size() << " instead of 6). Cannot determine K2 symmetry." << std::endl;
	else
	if  (data[5]<1.0e-13)
		K2min = -(K2max = 2);

	return;
	
}

extern bool cdy1, cdz1, cdy2, cdz2;

void basis::checkSymmetry(ldvector data, bool* fit){
	
	/* The entries of vectord data are
	 * data[0] = alpha_V
	 * data[1] = beta_V
	 * data[2] = gamma_V
	 * data[3] = D1xx - D1yy
	 * data[4] = abs[im{D1(2,1)}] + abs[im{D1(2,2)}] + abs[im{D2(2,1)}] + abs[im{D2(2,2)}]
	 * data[5] = abs(c22) + abs(c42) + abs(c44)
	 */
	
	
	/* 
	 * There is no global orienting potential
	 * so L1 and M1 are constrained to the values of the physical observable.
	 * Here we study only observables that transform as D20K, so we fix
	 * L1 = 2 and M1 = 0
	 */
	L1max = 2;
	L1min = 2;
	M1min = 0;
	M1max = 0;
	
	/*
	 * In case of no potential tilt and axially simmetryc diffusion tensor
	 * of the protein, Mm is diagonal and only basis functions with Mm = 0
	 * contribute to the spectral density
	 */
	
	bool checkMmSym = !fit[6] && !fit[7] && !fit[8] && ((!fit[0] && !fit[1]) || (fit[0] && cdy1) | cdz1) && !isOmegaVFromGeometry;
	bool checkjjSym = !fit[6] && !fit[7] && !fit[8] && !fit[9] && !fit[10] && !fit[11];
	bool checkK2Sym = !fit[19] && !fit[21] && !fit[22];
	
	if (data.size() < 4)
		std::cout << "TASK " << rank << ": WARNING : wrong number of arguments in basis::checkSymmetry (found " << data.size() << " instead of 4). Cannot determine Mm symmetry." << std::endl;
	else
	if  (fabs(data[0])<1.0e-13 && fabs(data[1])<1.0e-13 && fabs(data[2])<1.0e-13 && fabs(data[3])<1.0e-13 && checkMmSym)
		Mmmin = Mmmax = 0;
	else  // TO CHECK
	{
		Mmmin = -2;
		Mmmax =  2;
	}
	
	if (data.size() < 5)
		std::cout << "TASK " << rank << ": WARNING : wrong number of arguments in basis::checkSymmetry (found " << data.size() << " instead of 5). Cannot determine j symmetry." << std::endl;
	else
	if  (data[4]<1.0e-13 && checkjjSym)
		;//jjmin = jjmax = 1;

	if (data.size() < 6)
		std::cout << "TASK " << rank << ": WARNING : wrong number of arguments in basis::checkSymmetry (found " << data.size() << " instead of 6). Cannot determine K2 symmetry." << std::endl;
	else
	if  (data[5]<1.0e-13 && checkK2Sym)
		K2min = -(K2max = 2);

	return;
	
}

/******************************************
 * Handles call to basis indexes building *
 ******************************************/

void basis::buildSpace(std::string dynModel)
{
	std::cout << "TASK " << rank << ": Building basis indexes..." << std::endl;
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
		 buildSpace_SRLS();
	else if (!(dynModel.compare("fb1")) || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb")))
		buildSpace_FB1();
	else if (!(dynModel.compare("fb2")))
		buildSpace_FB2();
	std::cout << "TASK " << rank << ": ...indexs created with success" << std::endl;
}

/**************************************
 * Build basis indexes for SRLS model *
 **************************************/

void basis::buildSpace_SRLS()
{
	
	//K2min = -(K2max = 10);

	std::cout << "TASK " << rank << ": Building rotational space..." << std::endl;
	
	ivector idx(9);
	nbf = 0;
	basisIdx.clear();
	int K2low, jMin, jMax;
	
	for (L1 = L1min; L1 <= L1max; L1++){
		for (M1 = std::max(M1min,-L1); M1 <= std::min(M1max,L1); M1++){
			for (K1 = std::max(K1min,0); K1 <= std::min(K1max,L1); K1++){
				for (L2 = L2min; L2 <= L2max; L2++){
					for (M2 = std::max(M2min,-L2); M2 <= std::min(M2max,L2); M2++){
						K2low = (K1 > 0 ? -L2 : (M2 < 0) ? 1 : 0);
						for (K2 = std::max(K2min,K2low); K2 <= std::min(K2max,L2); K2++){
							jMin = ((K1 == 0 && K2 == 0 && M2 == 0) ? (!((L1+L2)%2) ? 1 : -1) : -1);
							jMax = ((K1 == 0 && K2 == 0 && M2 == 0) ? jMin : 1);
							for (jj = jMin; jj <= jMax; jj+=2){
								idx[0] = L1; idx[1] = M1; idx[2] = K1;
								idx[3] = L2; idx[4] = M2; idx[5] = K2;
								idx[6] = jj; idx[7] = K1-M2; idx[8] = K1+M2;
								if (idx[7] >= Mmmin && idx[7] <= Mmmax && idx[8] >= Mpmin && idx[8] <= Mpmax && idx[6] >= jjmin && idx[6] <= jjmax){
									basisIdx.push_back(idx);
									nbf++;
								}
							}
						}
					}
				}
			}
		}
	}

	std::cout << "TASK " << rank << ": ... rotational space created with success" << std::endl;
	
	return;
}

/*************************************
 * Build basis indexes for FB1 model *
 *************************************/

void basis::buildSpace_FB1()
{
	std::cout << "TASK " << rank << ": Building space..." << std::endl;
	ivector idx(5);
	nbf = 0;
	basisIdx.clear();
	int Klow, jMin, jMax;

	for (N1 = -K2max; N1 <= K2max; N1++){
		for (L1 = L1min; L1 <= L1max; L1++){
			for (M1 = std::max(M1min,-L1); M1 <= std::min(M1max,L1); M1++){
				Klow = (N1 < 0 ? 1 : 0);
				for (K1 = std::max(K1min,Klow); K1 <= std::min(K1max,L1); K1++){

					if (!K1 && !N1) jMin = jMax = ((L1%2)==0 ? 1 : -1);
					else jMax = -(jMin = -1);
					for (jj = jMin; jj <= jMax; jj+=2){
						idx[0] = N1;
						idx[1] = L1;
						idx[2] = M1;
						idx[3] = K1;
						idx[4] = jj;
						basisIdx.push_back(idx);
						nbf++;
					}
				}
			}
		}
	}
	std::cout << "TASK " << rank << ": ... space created with success" << std::endl;
	return;
}

/*************************************
 * Build basis indexes for FB2 model *
 *************************************/

void basis::buildSpace_FB2()
{
	std::cout << "TASK " << rank << ": Building space..." << std::endl;
	ivector idx(6);
	nbf = 0;
	basisIdx.clear();
	int Klow, jMin, jMax;
	
	for (N1 = -K2min; N1 <= K2min; N1++){
		for (N2 = -K2max; N2 <= K2max; N2++){
			for (L1 = L1min; L1 <= L1max; L1++){
				for (M1 = std::max(M1min,-L1); M1 <= std::min(M1max,L1); M1++){
					Klow = (((N1 < 0) || (N1 == 0 && N2 < 0)) ? 1 : 0);
					for (K1 = std::max(K1min,Klow); K1 <= std::min(K1max,L1); K1++){

						if (!K1 && !N1 && !N2) jMin = jMax = ((L1%2)==0 ? 1 : -1);
						else jMax = -(jMin = -1);
						for (jj = jMin; jj <= jMax; jj+=2){
							idx[0] = N1;
							idx[1] = N2;
							idx[2] = L1;
							idx[3] = M1;
							idx[4] = K1;
							idx[5] = jj;
							basisIdx.push_back(idx);
							nbf++;
						}
					}
				}
			}
		}
	}
	std::cout << "TASK " << rank << ": ... space created with success" << std::endl;
	return;
}

/*****************************************/
/* Returns the number of basis functions */
/*****************************************/

int basis::getNumberOfBasisFunctions()
{
	return nbf;
}

ivector basis::getBasisFunction(int n)
{
	if (n >= 0 && n < nbf)
		return basisIdx[n];
	else{
		std::cout << "TASK " << rank << ": COOPS ERROR. Integer parameter of basis::getBasisFunction must be in range [0 <= n < nbf]. Found n = " <<  n << ", while nbf = " << nbf << ".";
		exit(0);
	}
}

// Returns L1min and L1max
ivector basis::getL1bounds(void)
{
	ivector iv(2);
	iv[0] = L1min;
	iv[1] = L1max;
	return iv;
}

// Returns M1min and M1max
ivector basis::getM1bounds(void)
{
	ivector iv(2);
	iv[0] = M1min;
	iv[1] = M1max;
	return iv;
}

// Returns Mmmin and Mmmax
ivector basis::getMmbounds(void)
{
	ivector iv(2);
	iv[0] = Mmmin;
	iv[1] = Mmmax;
	return iv;
}

// Returns Mmmin and Mmmax
ivector basis::getjjbounds(void)
{
	ivector iv(2);
	iv[0] = jjmin;
	iv[1] = jjmax;
	return iv;
}

// Writes to std out all the rotational space
std::string basis::toString(std::string dynModel)
{
	std::ostringstream ostr;
	
	/**************/
	/* SRLS MODEL */
	/**************/

	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
		if (basisIdx.size() > 0){
			ivector idx;
			ostr << "TASK " << rank << ": ------------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << "             Basis functions indexes" << std::endl;
			ostr << "TASK " << rank << ": ------------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << ": \tL1\tM1\tK1\tL2\tM2\tK2\tj\tMm\tMp" << std::endl;
			for (unsigned int i=0; i<basisIdx.size(); i++){
				idx = basisIdx[i];
				ostr << "TASK " << rank << ": " << i+1 << "\t" << idx[0] << "\t" << idx[1] << "\t" << idx[2] << "\t" << idx[3] << "\t" << idx[4] << "\t" << idx[5] << "\t" << idx[6] << "\t" << idx[7] << "\t" << idx[8] << std::endl;
				}
			ostr << "TASK " << rank << ": ------------------------------------------------------------------" << std::endl;
			return ostr.str();
		}

		ostr << "TASK " << rank << ": Rotational space was not builded yet." << std::endl;
	}

	/*************/
	/* FB1  MOEL */
	/*************/

	else if (!(dynModel.compare("fb1")) || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb")))
	{
		if (basisIdx.size() > 0){
			ivector idx;
			ostr << "TASK " << rank << ": ------------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << "             Basis functions indexes" << std::endl;
			ostr << "TASK " << rank << ": ------------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << ": \tN\tL\tM\tK\tj" << std::endl;
			for (unsigned int i=0; i<basisIdx.size(); i++){
				idx = basisIdx[i];
				ostr << "TASK " << rank << ": " << i+1 << "\t" << idx[0] << "\t" << idx[1] << "\t" << idx[2] << "\t" << idx[3] << "\t" << idx[4] << std::endl;
			}
			ostr << "TASK " << rank << ": ------------------------------------------------------------------" << std::endl;
			return ostr.str();
		}

		ostr << "TASK " << rank << ": Rotational space was not builded yet." << std::endl;
	}

        /*************/
        /* FB2  MOEL */
        /*************/

	else if (!(dynModel.compare("fb2")))
	{
		if (basisIdx.size() > 0){
			ivector idx;
			ostr << "TASK " << rank << ": ------------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << "             Basis functions indexes" << std::endl;
			ostr << "TASK " << rank << ": ------------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << ": \tN1\tN2\tL\tM\tK\tj" << std::endl;
			for (unsigned int i=0; i<basisIdx.size(); i++){
				idx = basisIdx[i];
				ostr << "TASK " << rank << ": " << i+1 << "\t" << idx[0] << "\t" << idx[1] << "\t" << idx[2] << "\t" << idx[3] << "\t" << idx[4] << "\t" << idx[5] << std::endl;
			}
			ostr << "TASK " << rank << ": ------------------------------------------------------------------" << std::endl;
			return ostr.str();
		}
		
		ostr << "TASK " << rank << ": Rotational space was not builded yet." << std::endl;
	}


	return ostr.str();
	
}

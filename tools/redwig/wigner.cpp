/***********************************************************************************
 * C++OPPS 1.0 - Interpretation of NMR relaxation in proteins                      *
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
 Name        : wigner.cpp
 Author      : Mirco Zerbetto
 Version     : 1.0
 Copyright   : 2008 Mirco Zerbetto
 Description : Class with methods fo calculate both
               Wigner matrices and reduced Wigner matrices
 ============================================================================
 */
#include "wigner.h"

// Constructors
wigner::wigner(int RJ)
{	
	J = 0;
	M = 0;
	K = 0;
	alpha = 0.0;
	beta = 0.0;
	gamma = 0.0;
	lnf = new long double[6*RJ+1];
	lnf[0]=0.0;
	long double j1=1.0;
	for (int i=1; i<=6*RJ; i++) {
		lnf[i] = lnf[i-1] + log(j1);
		j1 += 1.0;
	}
}

wigner::wigner(int J, int M, int K, long double alpha, long double beta, long double gamma)
{
	this->J = J;
	this->M = M;
	this->K = K;
	this->alpha = alpha;
	this->beta = beta;
	this->gamma = gamma;
	
	lnf = new long double[6*J+1];
	lnf[0]=0.0;
	long double j1=1.0;
	for (int i=1; i<=6*J; i++) {
		lnf[i] = lnf[i-1] + log(j1);
		j1 += 1.0;
	}

}

//Destructor
wigner::~wigner()
{
}

// Set J,M,K
void wigner::setQuantumNumbers(int J, int M, int K)
{
	this->J = J;
	this->M = M;
	this->K = K;
	return;
}

// Set Euler angles
void wigner::setAngles(long double alpha, long double beta, long double gamma)
{
	this->alpha = alpha;
	this->beta = beta;
	this->gamma = gamma;
	return;
}

ldcomplex wigner::getWignerMatrix(){
	ldcomplex A, ris;
	long double C, mult, halfBeta;
	int i, imin, imax;

	A = complex<long double>(cos(alpha*(long double)M+gamma*(long double)K),-sin(alpha*(long double)M+gamma*(long double)K));
	  
	/************************************************************/
	/* Calculation of reduced Wigner matrix is based on formula */
	/* at page 22 of "Angular Momentum", D. M. Brink e G.       */
	/* R. Satchler. "k" is "n" in the book and "i" is "t".      */
	/************************************************************/

	imin = (0 > (M-K) ? 0 : (M-K));
	imax = ((J+M) < (J-K) ? (J+M) : (J-K));
	mult = (imin%2 == 0 ? 1.0 : -1.0);
	C = 0.0;
	halfBeta = 0.5*beta;

	for(i=imin;i<=imax;i++) {
		C += mult*(long double)(pow(cos(halfBeta),2*J+M-K-2*i)*pow(sin(halfBeta),2*i+K-M))*exp(-(lnf[J+M-i]+lnf[J-K-i]+lnf[i]+lnf[i+K-M]));
	    mult *= -1.0;
	}
	C *= exp(.5*(lnf[J+M]+lnf[J-M]+lnf[J+K]+lnf[J-K]));

	ris = C*A;

	return ris;
	
}

ldcomplex wigner::getWignerMatrix(int J, int M, int K, long double alpha, long double beta, long double gamma){

	this->J = J;
	this->M = M;
	this->K = K;
	this->alpha = alpha;
	this->beta = beta;
	this->gamma = gamma;
	
	ldcomplex A, ris;
	long double C, mult, halfBeta;
	int i, imin, imax;
	
	A = complex<long double>(cos(alpha*(long double)M+gamma*(long double)K),-sin(alpha*(long double)M+gamma*(long double)K));
	  
	/************************************************************/
	/* Calculation of reduced Wigner matrix is based on formula */
	/* at page 22 of "Angular Momentum", D. M. Brink e G.       */
	/* R. Satchler. "k" is "n" in the book and "i" is "t".      */
	/************************************************************/

	imin = (0 > (M-K) ? 0 : (M-K));
	imax = ((J+M) < (J-K) ? (J+M) : (J-K));
	mult = (imin%2 == 0 ? 1.0 : -1.0);
	C = 0.0;
	halfBeta = 0.5*beta;
	
	for(i=imin;i<=imax;i++) {
		C += mult*(long double)(pow(cos(halfBeta),2*J+M-K-2*i)*pow(sin(halfBeta),2*i+K-M))*exp(-(lnf[J+M-i]+lnf[J-K-i]+lnf[i]+lnf[i+K-M]));
	    	mult *= -1.0;
	}	
	C *= exp(.5*(lnf[J+M]+lnf[J-M]+lnf[J+K]+lnf[J-K]));

	ris = C*A;

	return ris;
  
}

long double wigner::getReducedWignerMatrix()
{
	long double bkpAlpha, bkpGamma;
	bkpAlpha = alpha;
	bkpGamma = gamma;
	ldcomplex d = this->getWignerMatrix(J,M,K,0.0,beta,0.0);
	alpha = bkpAlpha;
	gamma = bkpGamma;
	return d.real();
}

long double wigner::getReducedWignerMatrix(int J, int M, int K, long double beta)
{
	long double bkpAlpha, bkpGamma;
	bkpAlpha = alpha;
	bkpGamma = gamma;
	ldcomplex d = this->getWignerMatrix(J,M,K,0.0,beta,0.0);
	alpha = bkpAlpha;
	gamma = bkpGamma;
	return d.real();
}

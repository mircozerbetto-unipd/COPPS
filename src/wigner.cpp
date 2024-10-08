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
 Name        : wigner.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
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

// Return the Wigner matrix
ldcomplex wigner::getWignerMatrix()
{

	ldcomplex phaseAlpha, phaseGamma, bigD;
	long double smalld;
	
	phaseAlpha.real( cos((long double)M * alpha));
	phaseAlpha.imag(-sin((long double)M * alpha));
	phaseGamma.real( cos((long double)K * gamma));
	phaseGamma.imag(-sin((long double)K * gamma));

	smalld = getReducedWignerMatrix();

	bigD = phaseAlpha * smalld * phaseGamma;

	return bigD;

}

ldcomplex wigner::getWignerMatrix(int J, int M, int K, long double alpha, long double beta, long double gamma)
{

	this->J = J;
	this->M = M;
	this->K = K;
	this->alpha = alpha;
	this->beta  = beta;
	this->gamma = gamma;
	
	ldcomplex bigD;
	long double smalld;
	
	ldcomplex phaseAlpha(cos((long double)M * alpha), -sin((long double)M * alpha));
	ldcomplex phaseGamma(cos((long double)K * gamma), -sin((long double)K * gamma));

	smalld = getReducedWignerMatrix();

	bigD = phaseAlpha * smalld * phaseGamma;

	return bigD;
  
}

ldcomplex wigner::getWignerMatrix_d(int J, int M, int K, double alpha, double beta, double gamma)
{

	this->J = J;
	this->M = M;
	this->K = K;
	this->alpha = (long double)alpha;
	this->beta  = (long double)beta;
	this->gamma = (long double)gamma;
	
	ldcomplex bigD;
	long double smalld;
	
	ldcomplex phaseAlpha(cos((long double)M * alpha), -sin((long double)M * alpha));
	ldcomplex phaseGamma(cos((long double)K * gamma), -sin((long double)K * gamma));

	smalld = getReducedWignerMatrix();

	bigD = phaseAlpha * smalld * phaseGamma;

	return bigD;
  
}

long double wigner::getReducedWignerMatrix()
{
	long double dfun = getReducedWignerMatrix(this->J,this->M,this->K,this->beta);
	return dfun;
}

long double wigner::getReducedWignerMatrix(int JJ, int MM, int KK, long double beta)
{

	// Some analytical solutions //

	if (JJ == 0) return 1.0;
	if (JJ == 2) return getRank2RWM(MM,KK,beta);

	// Explicit solution //

	long double sb2=sin(0.5*beta);
	long double cb2=cos(0.5*beta);

	long double lognum = 0.5 * (lnf[JJ+MM] + lnf[JJ-MM] + lnf[JJ+KK] + lnf[JJ-KK]);
        long double dfun = 0.0;
	long double logden;
	long double mul;

	for (int n = 0; n <= 2*JJ; n++)
	{
		if (n >= 0 && (JJ-MM-n) >= 0 && (JJ-KK-n) >= 0 && (MM+KK+n) >= 0)
		{
			logden = lnf[n] + lnf[JJ-MM-n] + lnf[JJ-KK-n] + lnf[MM+KK+n];
			mul = ( !((n+JJ-KK)%2) ? 1.0 : -1.0);
			dfun += mul * exp(lognum-logden) * pow(cb2,MM+KK+2*n) * pow(sb2,2*JJ-MM-KK-2*n);
		}
	}

	return dfun;
}

long double wigner::getRank2RWM(int MM, int KK, long double b)
{
	int idx = (MM+2)*5+(KK+2);
	long double cb, sb, mul;
	switch(idx)
	{
		case 0: case 24:
		{
			cb = cos(0.5*b);
			return cb*cb*cb*cb;
		}
		case 1: case 5: case 19: case 23:
		{
			cb = cos(0.5*b);
			sb = sin(b);
			mul = (MM > KK ? -1.0 : 1.0);
			return mul*cb*cb*sb;
		}
		case 2: case 10: case 14: case 22:
		{
			sb = sin(b);
			return 0.5*SQRT_THREE_OVER_TWO*sb*sb;
		}
		case 3: case 9: case 15: case 21:
		{
			cb = sin(0.5*b);
			sb = sin(b);
			mul = (MM > KK ? -1.0 : 1.0);
			return mul*cb*cb*sb;
		}
		case 4: case 20:
		{
			sb = sin(0.5*b);
			return sb*sb*sb*sb;
		}
		case 6: case 18:
		{
			cb = cos(b);
			sb = sin(0.5*b);
			return cb*cb - sb*sb;
		}
		case 7: case 11: case 13: case 17:
		{
			sb = sin(2.0*b);
			mul = (MM > KK ? -1.0 : 1.0);
			return mul*0.5*SQRT_THREE_OVER_TWO*sb;
		}
		case 8: case 16:
		{
			cb = cos(0.5*b);
			sb = cos(b);
			return cb*cb - sb*sb;
		}
		case 12:
		{
			cb = cos(b);
			return 0.5*(3.0*cb*cb-1.0);
		}
		default:
		{
			std::cout << std::endl << std::endl << "ERROR : wrong value for M and/or K in wigner::getRank2RWM" << std::endl << std::endl;
			exit(1);
		}
	}
}

// OLD CODE //

//	/************************************************************/
//	/* Calculation of reduced Wigner matrix is based on formula */
//	/* at page 22 of "Angular Momentum", D. M. Brink e G.       */
//	/* R. Satchler. "k" is "n" in the book and "i" is "t".      */
//	/************************************************************/
//
//	imin = (0 > (M-K) ? 0 : (M-K));
//	imax = ((J+M) < (J-K) ? (J+M) : (J-K));
//	mult = (imin%2 == 0 ? 1.0 : -1.0);
//	C = 0.0;
//	halfBeta = 0.5*beta;
//	
//	for(i=imin;i<=imax;i++) {
//		C += mult*(long double)(pow(cos(halfBeta),2*J+M-K-2*i)*pow(sin(halfBeta),2*i+K-M))*exp(-(lnf[J+M-i]+lnf[J-K-i]+lnf[i]+lnf[i+K-M]));
//		mult *= -1.0;
//	}	
//	C *= exp(.5*(lnf[J+M]+lnf[J-M]+lnf[J+K]+lnf[J-K]));
//
//	ris = C*A;

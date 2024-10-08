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

#ifndef STVEC_H_
#define STVEC_H_

/* Standard C++ headers */
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <math.h>

/* coops headers */
#include "constants.h"
#include "physics.h"
#include "basis.h"
#include "wigner.h"
#include "prep.h"

/* MPI */
#ifdef _MPI_
#include "mpi.h"
#endif

/* SparseLib headers */
#include "ilupre_double.h"
#include "compcol_double.h"
#include "iohb_double.h"
#include "spblas.h"

/* Special functions */
extern "C"{
#include "cquadpak.h"
#include "bessel.h"
}

/* cubature library */
#include "cubature.h"

#include "types.h"

class stvec
{
	public:
		stvec();
		virtual ~stvec();
		
		void init(physics*);
		void init(physics*,int);
		void update(void);
		void updatePotential(void);
		
		void projectOnBasis(basis*,int,int);
		void scatterProjections(bool,int);
		VECTOR_double* getProjections(void);
		VECTOR_double getProjections(int);
		VECTOR_double getProjections(unsigned int);
		
		void setL1Constrain(int);
		int getL1Constrain(void);
		void setM1Constrain(int);
		int getM1Constrain(void);
		void setjjConstrain(int);
		int getjjConstrain(void);
		void setMmConstrain(int);
		int getMmConstrain(void);
		
		dvector getNpmKK1(void);
		double getNpmKK1(int);
		double getNpmKK1(unsigned int);
		dvector calculateOrdersParams(void);
		std::string toString(basis*);
		
		bool hasBeenInit(void);
				
	private:
		
		int rank;
		int nrows, globalRow, nbf;
		int KK,KK1;
		int L1constrain, M1constrain, Mmconstrain, jjconstrain;
		static int nTors;
		ivector sd;
		bool isL1constrained, isM1constrained, isMmconstrained, isjjconstrained;
		dvector npmkk1, normpm;
		physics *phy;
		static double delta;
		static wigner w;
		VECTOR_double *projections;

                static int nMult, nMult2, tors;
                static int L1,mK1,mK,L2,K1,K2;
                static long double c20,c22,c40,c42,c44;
                static dcvector coef_int_pot;
                static potential_fb2 ufb2;

		// FBN methods

		void projectOnBasis_FB1(basis*,int,int);
		void projectOnBasis_FB2(basis*,int,int);
		static double internal_integrand_real(double);
		static double internal_integrand_imag(double);
		static double fb1_peq(double);
		static double fb2_peq(const double *, int);
		static double peq_ratio(double, double, double, double);
		static int fb2_integrand_real(unsigned int, const double *, void *, unsigned int, double *);
		static int fb2_integrand_imag(unsigned int, const double *, void *, unsigned int, double *);
		
		static double fb2_integrand_real_j (double *, int *);
		static double fb2_integrand_imag_j (double *, int *);

		// SRLS methods

		void projectOnBasis_SRLS(basis*,int,int);
		inline double integrate(int,int,int,int,int,int,int);
		static double integrand_no_bessel(double);
		static double integrand_one_bessel_u2(double);
		static double integrand_one_bessel_u4(double);
		static double integrand_inf_bessel(double);
		static double integrandOrderParams(double);
		
		bool initiated;
		
		std::string dynModel; // string indicating the model for the dynamics

};

#endif /*STVEC_H_*/

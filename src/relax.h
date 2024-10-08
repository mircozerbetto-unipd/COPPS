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

#ifndef RELAX_H_
#define RELAX_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>

#include "constants.h"
#include "stvec.h"
#include "lanczos.h"
#include "wigner.h"
#include "spline.h"

extern "C"{
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
}


#include "types.h"

class relax
{
	public:
		
		relax();
		virtual ~relax();

		void init (stvec*,lanczos*,ldvector,ldvector,ldvector,dvector,double,ivector,int);
		void init (stvec*,lanczos*,ldvector,ldvector,ldvector,dvector,double,ivector,int,int);
		void init (stvec*,lanczos*,dvector,dvector,dvector,dvector,dvector,dvector,double,double,dvector,double,ivector,int,int);
		void update(ldvector,ldvector,ldvector,dvector,double,double);
		void update_tssrls (dvector, dvector, dvector, dvector, dvector, dvector, double, double, dvector);
		void makeSymmJ(int);
		void makeSmallJ(int);
		void makeBigJ(int);
		void makeBigJ_tssrls(int);
		dvector makeT1T2NOE(int);
		dvvector getT1T2NOE(bool);
		dvector getField(void);
                void setRex(double);
                void setRex(long double);
		
		void setNuclearData(dvector);

		void readSmallJs(int);
		void interpSmallJ(int);

		void setHchSigma(double);

		/* 3S-FB */
		void init (stvec*, lanczos*, dvector, dvector, dvector, dvector, dvector, dvector, dvector, dvector, dvector, dvector, dvector, dvector, double, ivector, int, int);
		void buildXFunctions(void);
		void makeSmallJ_3SFB(int);
		void makeBigJ_3SFB(int);
		/* END OF 3S-FB */

	private:
		
		int rank;
		int nH;
		ivector sd;
		dvector field;
		double larmorFreqHmin;
		double freqScale;
		double spectrumFactor;
                double rex;
		
		double gyroN, rNH, deltaCSA;
		
		ldvector omegaDip, omegaCsa, omegaDip2;
		dvector larmorFreqH, larmorFreqN, ww;
		dvector *freq;
		dcvector *symmj, *symmjB;
		dvector *smallj, *smalljAllFreq, *smalljB;
		dvector bigJdip, bigJcsa, bigJdip2, bigKdip;
		stvec* v;
		lanczos* lcz;

		bool isTsSrlsCalculation, isTsFb1Calculation, is3SFBCalculation;
		dvector omegaDip_1, omegaCsa_1, omegaDip2_1;
		dvector omegaDip_2, omegaCsa_2, omegaDip2_2;
		dvector omegaDip_3, omegaCsa_3, omegaDip2_3;
		double population, jumpFrequency;
		dvector A_dip_ac1, A_dip_ac2, A_dip_cc12, A_csa_ac;
		dvector B_dip_ac1, B_dip_ac2, B_dip_cc12, B_csa_ac;
		
		bool smallJdone;
		bool bigJdone;
		
		void makeLarmor(double);
		void buildBigJsWeights(void);
		ldcomplex averaged_reduced_d2m0(int, long double);

		double hchSigma;

		/* 3S-FB */
		double lambda[3];
		double K[9];
		dvector sqPopulation, jumpFrequencies;
		dcomplex *X_Dip1, *X_Dip2, *X_CSA;
		dcvector *symmj_1, *symmj_2, *symmj_3;
		dvector *smallj_1, *smallj_2, *smallj_3;
};

#endif /*RELAX_H_*/

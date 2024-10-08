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

#ifndef LANCZOS_H_
#define LANCZOS_H_

#include <cstdlib>
#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "prep.h"
#ifdef _MPI_
#include "mpi.h"
#endif

#include "constants.h"
#include "stvec.h"
#include "matrix.h"

#include "spblas.h"
#include "mvvd.h"
#include "comprow_double.h"

#include "types.h"

class lanczos
{
	public:
		
		lanczos();
		virtual ~lanczos();
		void init(int,stvec*,matrix*,ivector);
		void init(int,stvec*,matrix*,ivector,int,int);
		void update(void);
		
		void runLaczos(void);
		dcvector calculateSpectralDensity(dvector,dvector,dvector,double,bool);
		dcvector* calculateSpectralDensities(dvector,double,bool,bool,ldvector,ldvector,ldvector,int);
		
		dvvector* getAlphaPtr(void);
		dvvector* getBetaPtr(void);
		
		dvvector tqli(unsigned int,dvector,dvector);
		dvvector corrFunc(dvector,dvector,int,double,double);
		dvvector specDens(dvector,dvector,int,double,double);
		
		void setScale(double);

		std::string toString(bool,bool,int);
		
	private:
		
		int rank;
		int nproc;
		int nstep;
		int nstep0;
		int nbf;
		int jnumber;
		
		ivector sd;
		
		dvvector alpha;
		dvvector beta;
		
		stvec *v;
		matrix *mat;
		
		void smvm(CompRow_Mat_double*,VECTOR_double*,VECTOR_double*,VECTOR_double*);
		double ddotu(VECTOR_double*,VECTOR_double*);
		void daxpy(VECTOR_double*,VECTOR_double*,double);
		void dscsw(VECTOR_double*,VECTOR_double*,double);
		
		bool lastRun;
		bool oldRun;

		double scale;
};

#endif /*LANCZOS_H_*/

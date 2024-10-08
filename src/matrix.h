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

#ifndef MATRIX_H_
#define MATRIX_H_

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <complex>
#include <vector>
#include <iomanip>
#include <math.h>

#include "constants.h"
#include "basis.h"
#include "physics.h"
#include "wigner.h"
#include "s3j.h"

#include "comprow_double.h"
#include "mvvd.h"
#include "spblas.h"

#include "types.h"

#define MAKEJA
#define MAKEJB
#define MAKEJC
#define MAKEJD
#define MAKEFA_OR_MAKEFB
#define MAKEFA
#define MAKEFB

#define CMENO(x,y) sqrt((x)*((x)+1)-(y)*((y)-1))
#define CPIU(x,y)  sqrt((x)*((x)+1)-(y)*((y)+1))

class matrix
{
public:
	
	matrix();
	virtual ~matrix();
	
	void init(basis*,physics*,s3j*);
	void init(basis*,physics*,s3j*,int,int,int,int,int);
	void setupCalculation(void);
	void update(void);
	
	void storeFaFb(void);
	void storeFIIcoeff(int);
	void setDiften(void);
	void buildMatrix(void);
	dcomplex buildElement_SRLS(int,int,int,double*,double*,int*,int*);
	dcomplex buildElement_FB1(int,int,int,double*,double*,int*,int*);
	dcomplex buildElement_FB2(int,int,int,double*,double*,int*,int*);
	
	int getNrows(void);
	void setNrows(int);
	int getGlobalRow(void);
	void setGlobalRow(int);
	int getRealNZ(void);
	int getImagNZ(void);
	
	CompRow_Mat_double* getAptr(void);
	VECTOR_double* getDiagPtr(void);
	
	inline double cplm(int,int);
	inline double cmlm(int,int);
	
	std::string toString(int,int);
	
private:
	
	bool isInit;
	int rank;
	int nrows;
	int global_row;
	int ncols;
	int global_col;
	int realNZ;
	int imagNZ;
	int mulim, jlim;
	int oldi;
	int iL1, iK1, iM1, iL2, iK2, iM2, ijj, iN11, iN21; // global row indexes
	double dxx, dyy, dzz, dpr, dmr, d11, d22, d12, dx1, dy1, dz1, dx2, dy2, dz2; dcomplex dp1, dm1, dp2, dm2; // diffusion coefficients for FBn
	dcomplex probeD[6], proteinD[6]; // diffusion tensors SRLS
	dvector L2Mult;
	dvector *fa;
	dvector *fb;
	dvector cp,cm;
	dcvector f11; // array of coefficients of potential dependent part of internal operator in FB1 model
	potential_fb2 ufb2;
	basis *bas;
	physics *phy;
	s3j *trj;
	CompRow_Mat_double A;
	VECTOR_double diag;
	
	void grow(CompRow_Mat_double*,int,dvector,ivector);
        double getCP(int, int);
        double getCM(int, int);
	
	string dynModel; // info about model of dynamics
	
};

#endif /*MATRIX_H_*/

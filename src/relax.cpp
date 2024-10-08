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
 Name        : relax.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class to evaluate relaxation times T1, T2 and NOE
 ============================================================================
 */
#include "relax.h"

extern int OUT_CONTROL;
extern std::string path, project;

/***************
 * Constructor *
 ***************/

relax::relax()
{
	isTsSrlsCalculation = false;
	isTsFb1Calculation = false;
	is3SFBCalculation = false;
        // set some defaults for 15N-1H {Damberg et al. JACS (2005) 127, 1995-2005}
        gyroN = -2.7116e7;
        rNH = 1.015e-10;
        deltaCSA = -169.0;
}

/*****************
 * Deconstructor *
 *****************/

relax::~relax()
{
#ifdef WRITE_DESTROY_MESSAGE
	std::cout << "TASK " << rank << ": Cleared relax object" << std::endl;
#endif
}

/**************
 * Initiators *
 **************/

void relax::init(stvec* s, lanczos* l, ldvector dipEul, ldvector dipEul2, ldvector csaEul, dvector b, double scale, ivector js, int nHydrogens)
{
	rank = 0;
	v = s;
	lcz = l;
	omegaDip = dipEul;
	omegaDip2 = dipEul2;
	omegaCsa = csaEul;
	field = b;
	sd = js;
	nH = nHydrogens;
	freq = new dvector[field.size()];
	freqScale = 1.0/scale;
	spectrumFactor = 1.0;//6.0*1.0e8*11.5e-9;
	makeLarmor(freqScale);
	smallJdone = bigJdone = false;
        rex = 0.0;
	return;
}

void relax::init(stvec* s, lanczos* l, ldvector dipEul, ldvector dipEul2, ldvector csaEul, dvector b, double scale, ivector js, int nHydrogens, int r)
{
	rank = r;
	v = s;
	lcz = l;
	omegaDip = dipEul;
	omegaDip2 = dipEul2;
	omegaCsa = csaEul;
	field = b;
	sd = js;
	nH = nHydrogens;
	freq = new dvector[field.size()];
	freqScale = 1.0/scale;
	spectrumFactor = 1.0;//6.0*1.0e8*11.5e-9;
	makeLarmor(freqScale);
	smallJdone = bigJdone = false;
        rex = 0.0;
	return;
}

void relax::init (stvec* s, lanczos* l, dvector od_1, dvector oc_1, dvector od2_1, dvector od_2, dvector oc_2, dvector od2_2, double pop, double jf, dvector b, double scale , ivector js, int nHydrogens, int r)
{
	isTsSrlsCalculation = true;
	isTsFb1Calculation = true;
	rank = r;
	v = s;
	lcz = l;
	omegaDip_1  = od_1;
	omegaDip2_1 = od2_1;
	omegaCsa_1  = oc_1;
	omegaDip_2  = od_2;
	omegaDip2_2 = od2_2;
	omegaCsa_2  = oc_2;
	population = pop;
	jumpFrequency = jf;
	field = b;
	sd = js;
	nH = nHydrogens;
	freq = new dvector[field.size()];
	freqScale = 1.0/scale;
	spectrumFactor = 1.0;//6.0*1.0e8*11.5e-9;
	makeLarmor(freqScale);
	smallJdone = bigJdone = false;
        rex = 0.0;
	return;
}

/* 3S-FB */
void relax::init (stvec* s, lanczos* l, dvector od1_1, dvector oc_1, dvector od2_1, dvector od1_2, dvector oc_2, dvector od2_2, dvector od1_3, dvector oc_3, dvector od2_3, dvector pop, dvector jf, dvector b, double scale , ivector js, int nHydrogens, int r)
{
	is3SFBCalculation = true;
	rank = r;
	v = s;
	lcz = l;
	omegaDip_1  = od1_1;
	omegaDip2_1 = od2_1;
	omegaCsa_1  = oc_1;
	omegaDip_2  = od1_2;
	omegaDip2_2 = od2_2;
	omegaCsa_2  = oc_2;
	omegaDip_3  = od1_3;
	omegaDip2_3 = od2_3;
	omegaCsa_3  = oc_3;

	sqPopulation = pop;
	jumpFrequencies = jf;
	field = b;
	sd = js;
	nH = nHydrogens;
	freq = new dvector[field.size()];
	freqScale = 1.0/scale;
	spectrumFactor = 1.0;//6.0*1.0e8*11.5e-9;
	makeLarmor(freqScale);
	smallJdone = bigJdone = false;
        rex = 0.0;

	/* Build eigenvalues and eigenvectors of the rates matrix */
	K[0] = jf[0] * pop[1] * pop[1] + jf[1] * pop[2] * pop[2];
	K[1] = -jf[0] * pop[0] * pop[1];
	K[2] = -jf[1] * pop[0] * pop[2];
	K[3] = K[1];
	K[4] = jf[0] * pop[0] * pop[0] + jf[2] * pop[2] * pop[2];
	K[5] = -jf[2] * pop[1] * pop[2];
	K[6] = K[2];
	K[7] = K[5];
	K[8] = jf[1] * pop[0] * pop[0] + jf[2] * pop[1] * pop[1];

	for (int i = 0; i < 9; ++i)
		std::cout <<" K:: " << K[i] << std::endl;

	lambda[0]=lambda[1]=(lambda[2]=0.0);
	int info;
	int lwo=1000;
	double *work;
	work=(double *)calloc(lwo,sizeof(double));
	int i_index=3;
	char UPLO='U';
	char JOBZ='V';
	dsyev_(&JOBZ,&UPLO,&i_index,K,&i_index,lambda,work,&lwo,&info);

	return;
}

void relax::buildXFunctions(void)
{
	X_Dip1 = new dcomplex[3 * 5];
	X_Dip2 = new dcomplex[3 * 5];
	X_CSA  = new dcomplex[3 * 5];

	dcomplex *dk0csa, *dk0dip, *dk0dip2, cK;
	ldcomplex dlmk;
	wigner w(10);

	cK.imag(0.0);
	
	dk0csa = new dcomplex[5];
	dk0csa[0] = complex<double>(0.0,0.0);
	dk0csa[1] = complex<double>(0.0,0.0);
	dk0csa[2] = complex<double>(0.0,0.0);
	dk0csa[3] = complex<double>(0.0,0.0);
	dk0csa[4] = complex<double>(0.0,0.0);

	dk0dip = new dcomplex[5];
	dk0dip[0] = complex<double>(0.0,0.0);
	dk0dip[1] = complex<double>(0.0,0.0);
	dk0dip[2] = complex<double>(0.0,0.0);
	dk0dip[3] = complex<double>(0.0,0.0);
	dk0dip[4] = complex<double>(0.0,0.0);

	if (nH == 2)
	{
		dk0dip2 = new dcomplex[5];
		dk0dip2[0] = complex<double>(0.0,0.0);
		dk0dip2[1] = complex<double>(0.0,0.0);
		dk0dip2[2] = complex<double>(0.0,0.0);
		dk0dip2[3] = complex<double>(0.0,0.0);
		dk0dip2[4] = complex<double>(0.0,0.0);
	}

	/**********/
	/* SITE 1 */
	/**********/

	/* Calculate weights for CSA interaction */
	for (int i = -2; i <= 2; i++)
	{
		dlmk = w.getWignerMatrix(2,i,0,omegaCsa_1[0],omegaCsa_1[1],omegaCsa_1[2]);
		dk0csa[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
		dk0csa[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
		dk0csa[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
		dk0csa[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
		dk0csa[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
	}
	
	/* Calculate weights for Dipolar interaction */
	
	dk0dip[0] = (dcomplex)w.getWignerMatrix(2,-2,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);
	dk0dip[1] = (dcomplex)w.getWignerMatrix(2,-1,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);
	dk0dip[2] = (dcomplex)w.getWignerMatrix(2, 0,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);
	dk0dip[3] = (dcomplex)w.getWignerMatrix(2, 1,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);
	dk0dip[4] = (dcomplex)w.getWignerMatrix(2, 2,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);

	if (nH == 2) 
	{
		for (int i = -2; i <= 2; i++)
		{
			dlmk = w.getWignerMatrix(2,i,0,omegaDip2_1[0],omegaDip2_1[1],omegaDip2_1[2]);
			// Apply vibrational averaging
			//dlmk = averaged_reduced_d2m0(i,omegaDip2[1]);
			dk0dip2[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
			dk0dip2[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
			dk0dip2[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
			dk0dip2[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
			dk0dip2[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
		}
	}

	/* Calculate X of site 1 */
	X_Dip1[0 * 3 + 0] = complex<double>(0.0,0.0);
	X_Dip1[0 * 3 + 1] = complex<double>(0.0,0.0);
	X_Dip1[0 * 3 + 2] = complex<double>(0.0,0.0);
	X_Dip1[0 * 3 + 3] = complex<double>(0.0,0.0);
	X_Dip1[0 * 3 + 4] = complex<double>(0.0,0.0);
	for (int j = 0; j < 3; ++j)
	{
		cK.real(K[j * 3 + 0]);
		for (int k = -2; k <= 2; ++k)
			X_Dip1[0 * 3 + k + 2] += cK * dk0dip[k+2] * sqPopulation[j];
	}

	if (nH == 2) 
	{
		X_Dip2[0 * 3 + 0] = complex<double>(0.0,0.0);
		X_Dip2[0 * 3 + 1] = complex<double>(0.0,0.0);
		X_Dip2[0 * 3 + 2] = complex<double>(0.0,0.0);
		X_Dip2[0 * 3 + 3] = complex<double>(0.0,0.0);
		X_Dip2[0 * 3 + 4] = complex<double>(0.0,0.0);
		for (int j = 0; j < 3; ++j)
		{
			cK.real(K[j * 3 + 0]);
			for (int k = -2; k <= 2; ++k)
				X_Dip2[0 * 3 + k + 2] += cK * dk0dip2[k+2] * sqPopulation[j];
		}
	}

	X_CSA[0 * 3 + 0]  = complex<double>(0.0,0.0);
	X_CSA[0 * 3 + 1]  = complex<double>(0.0,0.0);
	X_CSA[0 * 3 + 2]  = complex<double>(0.0,0.0);
	X_CSA[0 * 3 + 3]  = complex<double>(0.0,0.0);
	X_CSA[0 * 3 + 4]  = complex<double>(0.0,0.0);
	for (int j = 0; j < 3; ++j)
	{
		cK.real(K[j * 3 + 0]);
		for (int k = -2; k <= 2; ++k)
			X_CSA[0 * 3 + k + 2] += cK * dk0csa[k+2] * sqPopulation[j];
	}

	/**********/
	/* SITE 2 */
	/**********/

	/* Calculate weights for CSA interaction */
	for (int i = -2; i <= 2; i++)
	{
		dlmk = w.getWignerMatrix(2,i,0,omegaCsa_2[0],omegaCsa_2[1],omegaCsa_2[2]);
		dk0csa[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		dk0csa[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		dk0csa[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		dk0csa[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		dk0csa[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
	}
	
	/* Calculate weights for Dipolar interaction */
	
	dk0dip[0] = (dcomplex)w.getWignerMatrix(2,-2,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);
	dk0dip[1] = (dcomplex)w.getWignerMatrix(2,-1,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);
	dk0dip[2] = (dcomplex)w.getWignerMatrix(2, 0,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);
	dk0dip[3] = (dcomplex)w.getWignerMatrix(2, 1,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);
	dk0dip[4] = (dcomplex)w.getWignerMatrix(2, 2,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);

	if (nH == 2) 
	{
		for (int i = -2; i <= 2; i++)
		{
			dlmk = w.getWignerMatrix(2,i,0,omegaDip2_2[0],omegaDip2_2[1],omegaDip2_2[2]);
			// Apply vibrational averaging
			//dlmk = averaged_reduced_d2m0(i,omegaDip2[1]);
			dk0dip2[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
			dk0dip2[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
			dk0dip2[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
			dk0dip2[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
			dk0dip2[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		}
	}

	/* Calculate X of site 2 */
	X_Dip1[1 * 3 + 0] = complex<double>(0.0,0.0);
	X_Dip1[1 * 3 + 1] = complex<double>(0.0,0.0);
	X_Dip1[1 * 3 + 2] = complex<double>(0.0,0.0);
	X_Dip1[1 * 3 + 3] = complex<double>(0.0,0.0);
	X_Dip1[1 * 3 + 4] = complex<double>(0.0,0.0);
	for (int j = 0; j < 3; ++j)
	{
		cK.real(K[j * 3 + 1]);
		for (int k = -2; k <= 2; ++k)
			X_Dip1[1 * 3 + k + 2] += cK * dk0dip[k+2] * sqPopulation[j];
	}

	if (nH == 2)
	{
		X_Dip2[1 * 3 + 0] = complex<double>(0.0,0.0);
		X_Dip2[1 * 3 + 1] = complex<double>(0.0,0.0);
		X_Dip2[1 * 3 + 2] = complex<double>(0.0,0.0);
		X_Dip2[1 * 3 + 3] = complex<double>(0.0,0.0);
		X_Dip2[1 * 3 + 4] = complex<double>(0.0,0.0);
		for (int j = 0; j < 3; ++j)
		{
			cK.real(K[j * 3 + 1]);
			for (int k = -2; k <= 2; ++k)
				X_Dip2[1 * 3 + k + 2] += cK * dk0dip2[k+2] * sqPopulation[j];
		}
	}

	X_CSA[1 * 3 + 0]  = complex<double>(0.0,0.0);
	X_CSA[1 * 3 + 1]  = complex<double>(0.0,0.0);
	X_CSA[1 * 3 + 2]  = complex<double>(0.0,0.0);
	X_CSA[1 * 3 + 3]  = complex<double>(0.0,0.0);
	X_CSA[1 * 3 + 4]  = complex<double>(0.0,0.0);
	for (int j = 0; j < 3; ++j)
	{
		cK.real(K[j * 3 + 1]);
		for (int k = -2; k <= 2; ++k)
			X_CSA[1 * 3 + k + 2] += cK * dk0csa[k+2] * sqPopulation[j];
	}

	/**********/
	/* SITE 3 */
	/**********/

	/* Calculate weights for CSA interaction */
	for (int i = -2; i <= 2; i++)
	{
		dlmk = w.getWignerMatrix(2,i,0,omegaCsa_3[0],omegaCsa_3[1],omegaCsa_3[2]);
		dk0csa[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]) * dlmk);
		dk0csa[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]) * dlmk);
		dk0csa[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]) * dlmk);
		dk0csa[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]) * dlmk);
		dk0csa[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]) * dlmk);
	}
	
	/* Calculate weights for Dipolar interaction */
	
	dk0dip[0] = (dcomplex)w.getWignerMatrix(2,-2,0,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]);
	dk0dip[1] = (dcomplex)w.getWignerMatrix(2,-1,0,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]);
	dk0dip[2] = (dcomplex)w.getWignerMatrix(2, 0,0,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]);
	dk0dip[3] = (dcomplex)w.getWignerMatrix(2, 1,0,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]);
	dk0dip[4] = (dcomplex)w.getWignerMatrix(2, 2,0,omegaDip_3[0],omegaDip_3[1],omegaDip_3[2]);

	if (nH == 2) 
	{
		for (int i = -2; i <= 2; i++)
		{
			dlmk = w.getWignerMatrix(2,i,0,omegaDip2_3[0],omegaDip2_3[1],omegaDip2_3[2]);
			// Apply vibrational averaging
			//dlmk = averaged_reduced_d2m0(i,omegaDip2[1]);
			dk0dip2[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip_3[0],omegaDip_3[1],omegaDip_2[2]) * dlmk);
			dk0dip2[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip_3[0],omegaDip_3[1],omegaDip_2[2]) * dlmk);
			dk0dip2[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip_3[0],omegaDip_3[1],omegaDip_2[2]) * dlmk);
			dk0dip2[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip_3[0],omegaDip_3[1],omegaDip_2[2]) * dlmk);
			dk0dip2[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip_3[0],omegaDip_3[1],omegaDip_2[2]) * dlmk);
		}
	}

	/* Calculate X of site 3 */
	X_Dip1[2 * 3 + 0] = complex<double>(0.0,0.0);
	X_Dip1[2 * 3 + 1] = complex<double>(0.0,0.0);
	X_Dip1[2 * 3 + 2] = complex<double>(0.0,0.0);
	X_Dip1[2 * 3 + 3] = complex<double>(0.0,0.0);
	X_Dip1[2 * 3 + 4] = complex<double>(0.0,0.0);
	for (int j = 0; j < 3; ++j)
	{
		cK.real(K[j * 3 + 2]);
		for (int k = -2; k <= 2; ++k)
			X_Dip1[2 * 3 + k + 2] += cK * dk0dip[k+2] * sqPopulation[j];
	}

	if (nH == 2) 
	{
		X_Dip2[2 * 3 + 0] = complex<double>(0.0,0.0);
		X_Dip2[2 * 3 + 1] = complex<double>(0.0,0.0);
		X_Dip2[2 * 3 + 2] = complex<double>(0.0,0.0);
		X_Dip2[2 * 3 + 3] = complex<double>(0.0,0.0);
		X_Dip2[2 * 3 + 4] = complex<double>(0.0,0.0);
		for (int j = 0; j < 3; ++j)
		{
			cK.real(K[j * 3 + 2]);
			for (int k = -2; k <= 2; ++k)
				X_Dip2[2 * 3 + k + 2] += cK * dk0dip2[k+2] * sqPopulation[j];
		}
	}

	X_CSA[2 * 3 + 0]  = complex<double>(0.0,0.0);
	X_CSA[2 * 3 + 1]  = complex<double>(0.0,0.0);
	X_CSA[2 * 3 + 2]  = complex<double>(0.0,0.0);
	X_CSA[2 * 3 + 3]  = complex<double>(0.0,0.0);
	X_CSA[2 * 3 + 4]  = complex<double>(0.0,0.0);
	for (int j = 0; j < 3; ++j)
	{
		cK.real(K[j * 3 + 2]);
		for (int k = -2; k <= 2; ++k)
			X_CSA[2 * 3 + k + 2] += cK * dk0csa[k+2] * sqPopulation[j];
	}

	return;
}

/* END OF 3S-FB */

/******************************************
 * Method to update object while fitting *
 *****************************************/

void relax::update(ldvector dipEul, ldvector dipEul2, ldvector csaEul, dvector b, double P, double F)
{
	omegaDip = dipEul;
	omegaDip2 = dipEul2;
	omegaCsa = csaEul;
	field = b;
	freq = new dvector[field.size()];
	smallJdone = bigJdone = false;
	makeLarmor(freqScale);
        rex = 0.0;
	population = P;
	jumpFrequency = F;
}

void relax::update_tssrls (dvector od_1, dvector oc_1, dvector od2_1, dvector od_2, dvector oc_2, dvector od2_2, double pop, double jf, dvector b)
{
	omegaDip_1  = od_1;
	omegaDip2_1 = od2_1;
	omegaCsa_1  = oc_1;
	omegaDip_2  = od_2;
	omegaDip2_2 = od2_2;
	omegaCsa_2  = oc_2;
	population = pop;
	jumpFrequency = jf;
	field = b;
	freq = new dvector[field.size()];
	smallJdone = bigJdone = false;
	makeLarmor(freqScale);
        rex = 0.0;
}

/********************
 * Set Nuclear Data *
 ********************/

void relax::setNuclearData(dvector nuc)
{
	gyroN = nuc.at(0);
	rNH = nuc.at(1);
	deltaCSA = nuc.at(2);
	freq = new dvector[field.size()];
	makeLarmor(freqScale);
	return;
}

/******************
 * Set R exchange *
 ******************/

void relax::setRex(double r)
{
  rex = r;
  return;
}
void relax::setRex(long double r)
{
  rex = (double)r;
  return;
}


/**************************************************************************
 * Calculate the Larmor frequencies for every value of the magnetic field *
 **************************************************************************/

//#define PRINT_LARMOR

void relax::makeLarmor(double scale)
{
	larmorFreqH.clear();
	larmorFreqN.clear();
#ifdef PRINT_LARMOR
	std::cout << "*********************************" << std::endl;
#endif
	for (unsigned int i = 0; i < field.size(); i++)
	{
		larmorFreqH.push_back(field[i] * 2.0*M_PI * scale * spectrumFactor);
		larmorFreqN.push_back(gyroN * larmorFreqH.back() / gyroH);
		freq[i].push_back(0.0);
		freq[i].push_back(larmorFreqH.back());
		freq[i].push_back(-(larmorFreqN.back()));
		freq[i].push_back((larmorFreqH.back()-larmorFreqN.back()));
		freq[i].push_back((larmorFreqH.back()+larmorFreqN.back()));
#ifdef PRINT_LARMOR
		std::cout << freq[i][0] << " "  << freq[i][1] << " " << freq[i][2] << " " << freq[i][3] << " " << freq[i][4] << std::endl;
#endif
	}
#ifdef PRINT_LARMOR
	std::cout << "*********************************" << std::endl;
#endif
	larmorFreqHmin = larmorFreqH[0];
	return;
}

/**********************************************
 * Calculate symmetrized j's at a given field *
 **********************************************/

void relax::makeSymmJ(int ifield)
{
	if (!isTsSrlsCalculation && !isTsFb1Calculation && !is3SFBCalculation)
	{
		if (!ifield && OUT_CONTROL == 1 && !rank)
			symmj = lcz->calculateSpectralDensities(freq[ifield],0.0,false,true,omegaDip,omegaDip2,omegaCsa,nH);
		else
			symmj = lcz->calculateSpectralDensities(freq[ifield],0.0,false,false,omegaDip,omegaDip2,omegaCsa,nH);
	}
	else if (is3SFBCalculation)
	{
		ldvector vzero (3, 0.0);
		std::cout << "LAMBDA = " << lambda[0] << "; " << lambda[1] << "; " << lambda[2] << std::endl;
		symmj_1 = lcz->calculateSpectralDensities(freq[ifield],lambda[0],false,OUT_CONTROL,vzero,vzero,vzero,nH);
		symmj_2 = lcz->calculateSpectralDensities(freq[ifield],lambda[1],false,OUT_CONTROL,vzero,vzero,vzero,nH);
		symmj_3 = lcz->calculateSpectralDensities(freq[ifield],lambda[2],false,OUT_CONTROL,vzero,vzero,vzero,nH);

	}
	else // NOTE: the output of correlation functions and spectral densities is not implemented for the TS-X model
	{
		ldvector vzero (3, 0.0);
		symmj  = lcz->calculateSpectralDensities(freq[ifield],0.0,false,false,vzero,vzero,vzero,nH);
		symmjB = lcz->calculateSpectralDensities(freq[ifield],jumpFrequency,false,false,vzero,vzero,vzero,nH);
	}
	return;
}

/****************************************************************************************
 * Calculate real part of small j's as combinations of symmetrized j's at a given field *
 ****************************************************************************************/

//#define PRINT_SMALLJ

void relax::makeSmallJ(int ifield)
{
	int ikk1, k, k1, ikk, ik1k1, z, z1;
	unsigned int y;
	smallj = new dvector[sd.size()];
	for (unsigned int i = 0; i < 2*sd.size(); i+=2)
	{
		y = i>>1;
		smallj[y].clear();
		ikk1 = sd[y];
		k = ikk1/5;  ikk   = (abs(k-2)+2)*5+(abs(k-2)+2);
		k1 = ikk1%5; ik1k1 = (abs(k1-2)+2)*5+(abs(k1-2)+2);
#ifdef PRINT_SMALLJ
		std::cout << (k-2) << " " << (k1-2) << " " << ikk << " " << ik1k1 << std::endl;
#endif
		for (unsigned int j = 0; j < sd.size(); j++)
		{
			if (sd[j] == ikk) z = 2*j;
			if (sd[j] == ik1k1) z1 = 2*j;
		}
#ifdef PRINT_SMALLJ
		std::cout<<"====================================="<<std::endl;
#endif
		for (unsigned int j = 0; j < freq[ifield].size(); j++)
		{
			smallj[y].push_back(4.0*(v->getNpmKK1(i)*symmj[i][j].real() + v->getNpmKK1(i+1)*symmj[i+1][j].real()));
			smallj[y][j] -= (v->getNpmKK1(z)*symmj[z][j].real() + v->getNpmKK1(z+1)*symmj[z+1][j].real());
			smallj[y][j] -= (v->getNpmKK1(z1)*symmj[z1][j].real() + v->getNpmKK1(z1+1)*symmj[z1+1][j].real());
			smallj[y][j] *= 0.125;
#ifdef PRINT_SMALLJ
			std::cout << j+1 << " "<<(k-2)<<" " << (k1-2)<< " " <<z<<" " <<z1<<  " " <<smallj[y][j] << std::endl;
#endif
		}
#ifdef PRINT_SMALLJ
		std::cout<<"====================================="<<std::endl;
#endif
	}

	/**********************************************/
	/* Small j's with boradening in TS-X model */
	/**********************************************/

#ifdef PRINT_SMALLJ
	std::cout << std::endl << "Spectral densities with broadening given by jump frequency" << std::endl;
#endif
	if (isTsSrlsCalculation || isTsFb1Calculation)
	{
		smalljB = new dvector[sd.size()];
		for (unsigned int i = 0; i < 2*sd.size(); i+=2)
		{
			y = i>>1;
			smalljB[y].clear();
			ikk1 = sd[y];
			k = ikk1/5;  ikk   = (abs(k-2)+2)*5+(abs(k-2)+2);
			k1 = ikk1%5; ik1k1 = (abs(k1-2)+2)*5+(abs(k1-2)+2);
#ifdef PRINT_SMALLJ
			std::cout << (k-2) << " " << (k1-2) << " " << ikk << " " << ik1k1 << std::endl;
#endif
			for (unsigned int j = 0; j < sd.size(); j++)
			{
				if (sd[j] == ikk) z = 2*j;
				if (sd[j] == ik1k1) z1 = 2*j;
			}
#ifdef PRINT_SMALLJ
			std::cout<<"====================================="<<std::endl;
#endif
			for (unsigned int j = 0; j < freq[ifield].size(); j++)
			{
				smalljB[y].push_back(4.0*(v->getNpmKK1(i)*symmjB[i][j].real() + v->getNpmKK1(i+1)*symmjB[i+1][j].real()));
				smalljB[y][j] -= (v->getNpmKK1(z)*symmjB[z][j].real() + v->getNpmKK1(z+1)*symmjB[z+1][j].real());
				smalljB[y][j] -= (v->getNpmKK1(z1)*symmjB[z1][j].real() + v->getNpmKK1(z1+1)*symmjB[z1+1][j].real());
				smalljB[y][j] *= 0.125;
#ifdef PRINT_SMALLJ
				std::cout << j+1 << " "<<(k-2)<<" " << (k1-2)<< " " <<z<<" " <<z1<<  " " <<smalljB[y][j] << std::endl;
#endif
			}
#ifdef PRINT_SMALLJ
			std::cout<<"====================================="<<std::endl;
#endif
		}
	}
	return;
}

/* 3S-FB */
void relax::makeSmallJ_3SFB(int ifield)
{
	int ikk1, k, k1, ikk, ik1k1, z, z1;
	unsigned int y;

	/* SITE 1 */

	smallj_1 = new dvector[sd.size()];
	for (unsigned int i = 0; i < 2*sd.size(); i+=2)
	{
		y = i>>1;
		smallj_1[y].clear();
		ikk1 = sd[y];
		k = ikk1/5;  ikk   = (abs(k-2)+2)*5+(abs(k-2)+2);
		k1 = ikk1%5; ik1k1 = (abs(k1-2)+2)*5+(abs(k1-2)+2);
#ifdef PRINT_SMALLJ
		std::cout << (k-2) << " " << (k1-2) << " " << ikk << " " << ik1k1 << std::endl;
#endif
		for (unsigned int j = 0; j < sd.size(); j++)
		{
			if (sd[j] == ikk) z = 2*j;
			if (sd[j] == ik1k1) z1 = 2*j;
		}
#ifdef PRINT_SMALLJ
		std::cout<<"====================================="<<std::endl;
#endif
		for (unsigned int j = 0; j < freq[ifield].size(); j++)
		{
			smallj_1[y].push_back(4.0*(v->getNpmKK1(i)*symmj_1[i][j].real() + v->getNpmKK1(i+1)*symmj_1[i+1][j].real()));
			smallj_1[y][j] -= (v->getNpmKK1(z)*symmj_1[z][j].real() + v->getNpmKK1(z+1)*symmj_1[z+1][j].real());
			smallj_1[y][j] -= (v->getNpmKK1(z1)*symmj_1[z1][j].real() + v->getNpmKK1(z1+1)*symmj_1[z1+1][j].real());
			smallj_1[y][j] *= 0.125;
#ifdef PRINT_SMALLJ
			std::cout << "SMALL j 1: " <<  j+1 << " "<<(k-2)<<" " << (k1-2)<< " " <<z<<" " <<z1<<  " " <<smallj_1[y][j] << std::endl;
#endif
		}
#ifdef PRINT_SMALLJ
		std::cout<<"====================================="<<std::endl;
#endif
	}

	/* SITE 2 */

	smallj_2 = new dvector[sd.size()];
	for (unsigned int i = 0; i < 2*sd.size(); i+=2)
	{
		y = i>>1;
		smallj_2[y].clear();
		ikk1 = sd[y];
		k = ikk1/5;  ikk   = (abs(k-2)+2)*5+(abs(k-2)+2);
		k1 = ikk1%5; ik1k1 = (abs(k1-2)+2)*5+(abs(k1-2)+2);
#ifdef PRINT_SMALLJ
		std::cout << (k-2) << " " << (k1-2) << " " << ikk << " " << ik1k1 << std::endl;
#endif
		for (unsigned int j = 0; j < sd.size(); j++)
		{
			if (sd[j] == ikk) z = 2*j;
			if (sd[j] == ik1k1) z1 = 2*j;
		}
#ifdef PRINT_SMALLJ
		std::cout<<"====================================="<<std::endl;
#endif
		for (unsigned int j = 0; j < freq[ifield].size(); j++)
		{
			smallj_2[y].push_back(4.0*(v->getNpmKK1(i)*symmj_2[i][j].real() + v->getNpmKK1(i+1)*symmj_2[i+1][j].real()));
			smallj_2[y][j] -= (v->getNpmKK1(z)*symmj_2[z][j].real() + v->getNpmKK1(z+1)*symmj_2[z+1][j].real());
			smallj_2[y][j] -= (v->getNpmKK1(z1)*symmj_2[z1][j].real() + v->getNpmKK1(z1+1)*symmj_2[z1+1][j].real());
			smallj_2[y][j] *= 0.125;
#ifdef PRINT_SMALLJ
			std::cout << j+1 << " "<<(k-2)<<" " << (k1-2)<< " " <<z<<" " <<z1<<  " " <<smallj_2[y][j] << std::endl;
#endif
		}
#ifdef PRINT_SMALLJ
		std::cout<<"====================================="<<std::endl;
#endif
	}

	/* SITE 3 */

	smallj_3 = new dvector[sd.size()];
	for (unsigned int i = 0; i < 2*sd.size(); i+=2)
	{
		y = i>>1;
		smallj_3[y].clear();
		ikk1 = sd[y];
		k = ikk1/5;  ikk   = (abs(k-2)+2)*5+(abs(k-2)+2);
		k1 = ikk1%5; ik1k1 = (abs(k1-2)+2)*5+(abs(k1-2)+2);
#ifdef PRINT_SMALLJ
		std::cout << (k-2) << " " << (k1-2) << " " << ikk << " " << ik1k1 << std::endl;
#endif
		for (unsigned int j = 0; j < sd.size(); j++)
		{
			if (sd[j] == ikk) z = 2*j;
			if (sd[j] == ik1k1) z1 = 2*j;
		}
#ifdef PRINT_SMALLJ
		std::cout<<"====================================="<<std::endl;
#endif
		for (unsigned int j = 0; j < freq[ifield].size(); j++)
		{
			smallj_3[y].push_back(4.0*(v->getNpmKK1(i)*symmj_3[i][j].real() + v->getNpmKK1(i+1)*symmj_3[i+1][j].real()));
			smallj_3[y][j] -= (v->getNpmKK1(z)*symmj_3[z][j].real() + v->getNpmKK1(z+1)*symmj_3[z+1][j].real());
			smallj_3[y][j] -= (v->getNpmKK1(z1)*symmj_3[z1][j].real() + v->getNpmKK1(z1+1)*symmj_3[z1+1][j].real());
			smallj_3[y][j] *= 0.125;
#ifdef PRINT_SMALLJ
			std::cout << j+1 << " "<<(k-2)<<" " << (k1-2)<< " " <<z<<" " <<z1<<  " " <<smallj_3[y][j] << std::endl;
#endif
		}
#ifdef PRINT_SMALLJ
		std::cout<<"====================================="<<std::endl;
#endif
	}
	return;
}

/*****************************************************************************************
 * Calculate real part of big J's for both dipolar and csa interactions at a given field *
 *****************************************************************************************/

#define PHASE(a) (a%2 == 0 ? 1.0 : -1.0)

extern "C"{
#include "cquadpak.h"
}

double peqBeta(double x);
double avgD200(double x);
double avgD210(double x);
double avgD220(double x);

double hch;
double sigma;

void relax::makeBigJ(int ifield)
{
	bigJdip.clear();
	bigJcsa.clear();
	bigJdip2.clear();
	bigKdip.clear();
	dcomplex *dk0csa, *dk0dip, *dk0dip2;
	ldcomplex dlmk;
	wigner w(10);
	
	/*****************************************
	 * Calculate weights for CSA interaction *
	 *****************************************/
	
	dk0csa = new dcomplex[5];
	dk0csa[0] = complex<double>(0.0,0.0);
	dk0csa[1] = complex<double>(0.0,0.0);
	dk0csa[2] = complex<double>(0.0,0.0);
	dk0csa[3] = complex<double>(0.0,0.0);
	dk0csa[4] = complex<double>(0.0,0.0);
	for (int i = -2; i <= 2; i++)
	{
		dlmk = w.getWignerMatrix(2,i,0,omegaCsa[0],omegaCsa[1],omegaCsa[2]);
		dk0csa[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
		dk0csa[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
		dk0csa[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
		dk0csa[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
		dk0csa[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
	}
	
	/*********************************************
	 * Calculate weights for Dipolar interaction *
	 *********************************************/
	
	dk0dip = new dcomplex[5];
	dk0dip[0] = complex<double>(0.0,0.0);
	dk0dip[1] = complex<double>(0.0,0.0);
	dk0dip[2] = complex<double>(0.0,0.0);
	dk0dip[3] = complex<double>(0.0,0.0);
	dk0dip[4] = complex<double>(0.0,0.0);
	dk0dip[0] = (dcomplex)w.getWignerMatrix(2,-2,0,omegaDip[0],omegaDip[1],omegaDip[2]);
	dk0dip[1] = (dcomplex)w.getWignerMatrix(2,-1,0,omegaDip[0],omegaDip[1],omegaDip[2]);
	dk0dip[2] = (dcomplex)w.getWignerMatrix(2, 0,0,omegaDip[0],omegaDip[1],omegaDip[2]);
	dk0dip[3] = (dcomplex)w.getWignerMatrix(2, 1,0,omegaDip[0],omegaDip[1],omegaDip[2]);
	dk0dip[4] = (dcomplex)w.getWignerMatrix(2, 2,0,omegaDip[0],omegaDip[1],omegaDip[2]);

	if (nH == 2) 
	{
		dk0dip2 = new dcomplex[5];
		dk0dip2[0] = complex<double>(0.0,0.0);
		dk0dip2[1] = complex<double>(0.0,0.0);
		dk0dip2[2] = complex<double>(0.0,0.0);
		dk0dip2[3] = complex<double>(0.0,0.0);
		dk0dip2[4] = complex<double>(0.0,0.0);
		for (int i = -2; i <= 2; i++)
		{
			dlmk = w.getWignerMatrix(2,i,0,omegaDip2[0],omegaDip2[1],omegaDip2[2]);
			// Apply vibrational averaging
			//dlmk = averaged_reduced_d2m0(i,omegaDip2[1]);
			dk0dip2[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
			dk0dip2[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
			dk0dip2[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
			dk0dip2[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
			dk0dip2[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
		}
	}

	/********************************
	 * Calculate spectral densities *
	 ********************************/

	double DD[5][5];
	double DD2[5][5];
	double KK[5][5];
	double CC[5][5];
	dcomplex ctmp;

	for (int k1 = -2; k1 <= 2; k1++)
	{
		for (int k2 = -2; k2 <= 2; k2++)
		{
			/* dipolar-dipolar factors */
#ifdef F77_PHASE
			ctmp = dk0dip[-k1+2] * dk0dip[k2+2];             // F77 copps phase
#else
			ctmp = PHASE(k1) * dk0dip[-k1+2] * dk0dip[k2+2];
#endif
			DD[k1+2][k2+2] = ctmp.real();
			/* csa-csa factors */
			ctmp = PHASE(k1) * dk0csa[-k1+2] * dk0csa[k2+2];
			CC[k1+2][k2+2] = ctmp.real();
			if (nH == 2)
			{
				/* dipolar2-dipolar2 factors */
				ctmp = PHASE(k1) * dk0dip2[-k1+2] * dk0dip2[k2+2];
				DD2[k1+2][k2+2] = ctmp.real();
				/* dipolar1-dipolar2 factors */
				ctmp = PHASE(k1) * dk0dip[-k1+2] * dk0dip2[k2+2];
				KK[k1+2][k2+2] = ctmp.real();
			}
		}
	}


	for (unsigned int i = 0; i < freq[ifield].size(); i++)
	{
		
		bigJdip.push_back(0.0);
		bigJcsa.push_back(0.0);
		bigJdip2.push_back(0.0);
		bigKdip.push_back(0.0);
		
		/***************
		 * K = K' part *
		 ***************/

		/* (0, 0) */
		bigJdip[i] += smallj[0][i] * DD[2][2];
		bigJcsa[i] += smallj[0][i] * CC[2][2];
		if (nH == 2)
		{
			bigJdip2[i] += smallj[0][i] * DD2[2][2];
			bigKdip[i]  += smallj[0][i] * KK[2][2];
		}
		/* (-1,-1) + (1,1) */
		bigJdip[i] += smallj[1][i] * (DD[1][1] + DD[3][3]);
		bigJcsa[i] += smallj[1][i] * (CC[1][1] + CC[3][3]);
		if (nH == 2)
		{
			bigJdip2[i] += smallj[1][i] * (DD2[1][1] + DD2[3][3]);
			bigKdip[i]  += smallj[1][i] * (KK[1][1]  + KK[3][3]);
		}
		/* (-2,-2) + (2,2) */
		bigJdip[i] += smallj[2][i] * (DD[0][0] + DD[4][4]);
		bigJcsa[i] += smallj[2][i] * (CC[0][0] + CC[4][4]);
		if (nH == 2)
		{
			bigJdip2[i] += smallj[2][i] * (DD2[0][0] + DD2[4][4]);
			bigKdip[i]  += smallj[2][i] * (KK[0][0]  + KK[4][4]);
		}

		/* DEBUGGING */
		/*	
		std::cout << "SMALL J 0 (" << i+1<<") = " << smallj[0][i]<<std::endl;
		std::cout << "SMALL J 1 (" << i+1<<") = " << smallj[1][i]<<std::endl;
		std::cout << "SMALL J 2 (" << i+1<<") = " << smallj[2][i]<<std::endl;
		*/

		/****************
		 * K != K' part *
		 ****************/

		//     j(K,K') = (-)^(K-K') j(-K,-K')    // CORRECT?
		
		if (sd.size() > 3)
		{
			/* (-2, 2) */
			bigJdip[i] += 2.0*smallj[3][i] * DD[0][4];
			bigJcsa[i] += 2.0*smallj[3][i] * CC[0][4];
			if (nH == 2)
			{
				bigJdip2[i] += 2.0*smallj[3][i] * DD2[0][4];
				bigKdip[i]  += 2.0*smallj[3][i] * KK[0][4];
			}
			/* (-1, 1) */
			bigJdip[i] += 2.0*smallj[4][i] * DD[1][3];
			bigJcsa[i] += 2.0*smallj[4][i] * CC[1][3];
			if (nH == 2)
			{
				bigJdip2[i] += 2.0*smallj[4][i] * DD2[1][3];
				bigKdip[i]  += 2.0*smallj[4][i] * KK[1][3];
			}
			/* (-2, 0) + (0, 2) */
			bigJdip[i] += 2.0*smallj[5][i] * (DD[0][2] + DD[2][4]);
			bigJcsa[i] += 2.0*smallj[5][i] * (CC[0][2] + CC[2][4]);
			if (nH == 2)
			{
				bigJdip2[i] += 2.0*smallj[5][i] * (DD2[0][2] + DD2[2][4]);
				bigKdip[i]  += 2.0*smallj[5][i] * (KK[0][2]  + KK[2][4]);
			}

			/* DEBUGGING */
			/*
			std::cout << "SMALL J 3 (" << i+1<<") = " << smallj[3][i]<<std::endl;
			std::cout << "SMALL J 4 (" << i+1<<") = " << smallj[4][i]<<std::endl;
			std::cout << "SMALL J 5 (" << i+1<<") = " << smallj[5][i]<<std::endl;
			*/

			if (sd.size() > 6)
			{
				/* (-1, 2) + (-2, 1) */
				bigJdip[i] += 2.0*smallj[6][i] * (DD[1][4] - DD[0][3]);
				bigJcsa[i] += 2.0*smallj[6][i] * (CC[1][4] - CC[0][3]);
				if (nH == 2)
				{
					bigJdip2[i] += 2.0*smallj[6][i] * (DD2[1][4] - DD2[0][3]);
					bigKdip[i]  += 2.0*smallj[6][i] * (KK[1][4]  - KK[0][3]);
				}
				/* (0, 1) + (-1, 0) */
				bigJdip[i] += 2.0*smallj[7][i] * (DD[2][3] - DD[1][2]);
				bigJcsa[i] += 2.0*smallj[7][i] * (CC[2][3] - CC[1][2]);
				if (nH == 2)
				{
					bigJdip2[i] += 2.0*smallj[7][i] * (DD2[2][3] - DD2[1][2]);
					bigKdip[i]  += 2.0*smallj[7][i] * (KK[2][3]  - KK[1][2]);
				}
				/* (1, 2) + (-2, -1) */
				bigJdip[i] += 2.0*smallj[8][i] * (DD[3][4] - DD[0][1]);
				bigJcsa[i] += 2.0*smallj[8][i] * (CC[3][4] - CC[0][1]);
				if (nH == 2)
				{
					bigJdip2[i] += 2.0*smallj[8][i] * (DD2[3][4] - DD2[0][1]);
					bigKdip[i]  += 2.0*smallj[8][i] * (KK[3][4]  - KK[0][1]);
				}

				/* DEBUGGING */
				/*
				std::cout << "SMALL J 6 (" << i+1<<") = " << smallj[6][i]<<std::endl;
				std::cout << "SMALL J 7 (" << i+1<<") = " << smallj[7][i]<<std::endl;
				std::cout << "SMALL J 8 (" << i+1<<") = " << smallj[8][i]<<std::endl;
				*/

			}
		}
	}

	return;
}

// Methods needed by the previous routine to calculate vibrational average

#define SQ_3_OVER_2 sqrt(1.5)
double peqBeta(double x)
{
	double xe = hch;
	double s2 = sigma * sigma;
        return exp(-0.5 * (x - xe) * (x - xe) / s2);
}

double avgD200(double x)
{
	return sin(x) * (0.5 * (3. * cos(x) * cos(x) - 1.0)) * peqBeta(x);
}
double avgD210(double x)
{
	return sin(x) * (-SQ_3_OVER_2 * sin(x) * cos(x)) * peqBeta(x);
}
double avgD220(double x)
{
	return sin(x) * (0.5 * SQ_3_OVER_2 * sin(x) * sin(x)) * peqBeta(x);
}

/******************************************/
/* Create weights for the TS-X big J's */
/******************************************/

void relax::buildBigJsWeights(void)
{

	double rho1, rho2;
	ldcomplex dlmk;
	dcomplex dk0csa_1[5], dk0dip_1[5], dk0dip2_1[5];
	dcomplex dk0csa_2[5], dk0dip_2[5], dk0dip2_2[5];

	// Init Wigner matrices calculator

	wigner w(10);

	// Init arrays of weights

	A_dip_ac1  = dvector (5 * 5, 0.0);
	A_dip_ac2  = dvector (5 * 5, 0.0);
	A_dip_cc12 = dvector (5 * 5, 0.0);
	A_csa_ac   = dvector (5 * 5, 0.0);

	B_dip_ac1  = dvector (5 * 5, 0.0);
	B_dip_ac2  = dvector (5 * 5, 0.0);
	B_dip_cc12 = dvector (5 * 5, 0.0);
	B_csa_ac   = dvector (5 * 5, 0.0);
	
	// Calculate Wigner matrices for CSA interaction
	
	dk0csa_1[0] = complex<double>(0.0,0.0);
	dk0csa_1[1] = complex<double>(0.0,0.0);
	dk0csa_1[2] = complex<double>(0.0,0.0);
	dk0csa_1[3] = complex<double>(0.0,0.0);
	dk0csa_1[4] = complex<double>(0.0,0.0);

	dk0csa_2[0] = complex<double>(0.0,0.0);
	dk0csa_2[1] = complex<double>(0.0,0.0);
	dk0csa_2[2] = complex<double>(0.0,0.0);
	dk0csa_2[3] = complex<double>(0.0,0.0);
	dk0csa_2[4] = complex<double>(0.0,0.0);

	for (int i = -2; i <= 2; i++)
	{
		dlmk = w.getWignerMatrix_d(2,i,0,omegaCsa_1[0],omegaCsa_1[1],omegaCsa_1[2]);
		dk0csa_1[0] += (dcomplex)(w.getWignerMatrix_d(2,-2,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
		dk0csa_1[1] += (dcomplex)(w.getWignerMatrix_d(2,-1,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
		dk0csa_1[2] += (dcomplex)(w.getWignerMatrix_d(2, 0,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
		dk0csa_1[3] += (dcomplex)(w.getWignerMatrix_d(2, 1,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
		dk0csa_1[4] += (dcomplex)(w.getWignerMatrix_d(2, 2,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);

		dlmk = w.getWignerMatrix(2,i,0,omegaCsa_2[0],omegaCsa_2[1],omegaCsa_2[2]);
		dk0csa_2[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		dk0csa_2[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		dk0csa_2[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		dk0csa_2[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		dk0csa_2[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
	}

	// Calculate Wigner matrices for Dipolar interactions
	
	dk0dip_1[0] = (dcomplex)w.getWignerMatrix(2,-2,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);
	dk0dip_1[1] = (dcomplex)w.getWignerMatrix(2,-1,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);
	dk0dip_1[2] = (dcomplex)w.getWignerMatrix(2, 0,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);
	dk0dip_1[3] = (dcomplex)w.getWignerMatrix(2, 1,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);
	dk0dip_1[4] = (dcomplex)w.getWignerMatrix(2, 2,0,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]);

	dk0dip_2[0] = (dcomplex)w.getWignerMatrix(2,-2,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);
	dk0dip_2[1] = (dcomplex)w.getWignerMatrix(2,-1,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);
	dk0dip_2[2] = (dcomplex)w.getWignerMatrix(2, 0,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);
	dk0dip_2[3] = (dcomplex)w.getWignerMatrix(2, 1,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);
	dk0dip_2[4] = (dcomplex)w.getWignerMatrix(2, 2,0,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]);

	if (nH == 2) 
	{
		dk0dip2_1[0] = complex<double>(0.0,0.0);
		dk0dip2_1[1] = complex<double>(0.0,0.0);
		dk0dip2_1[2] = complex<double>(0.0,0.0);
		dk0dip2_1[3] = complex<double>(0.0,0.0);
		dk0dip2_1[4] = complex<double>(0.0,0.0);

		dk0dip2_2[0] = complex<double>(0.0,0.0);
		dk0dip2_2[1] = complex<double>(0.0,0.0);
		dk0dip2_2[2] = complex<double>(0.0,0.0);
		dk0dip2_2[3] = complex<double>(0.0,0.0);
		dk0dip2_2[4] = complex<double>(0.0,0.0);

		for (int i = -2; i <= 2; i++)
		{
			dlmk = w.getWignerMatrix(2,i,0,omegaDip2_1[0],omegaDip2_1[1],omegaDip2_1[2]);
			// Apply vibrational averaging
			//dlmk = averaged_reduced_d2m0(i,omegaDip2_1[1]);
			dk0dip2_1[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
			dk0dip2_1[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
			dk0dip2_1[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
			dk0dip2_1[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);
			dk0dip2_1[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip_1[0],omegaDip_1[1],omegaDip_1[2]) * dlmk);

			dlmk = w.getWignerMatrix(2,i,0,omegaDip2_2[0],omegaDip2_2[1],omegaDip2_2[2]);
			// Apply vibrational averaging
			//dlmk = averaged_reduced_d2m0(i,omegaDip2_2[1]);
			dk0dip2_2[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
			dk0dip2_2[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
			dk0dip2_2[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
			dk0dip2_2[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
			dk0dip2_2[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip_2[0],omegaDip_2[1],omegaDip_2[2]) * dlmk);
		}
	}

	// Calculate weights for the magnetic interactions

	rho1 = population;
	rho2 = 1.0 - population;

	for (int k1 = 0; k1 <= 4; k1++)
	{
		for (int k2 = 0; k2 <= 4; k2++)
		{
			// 1. CSA-CSA acf
			A_csa_ac[k1 * 5 + k2]  = rho1 * rho1 * (conj(dk0csa_1[k1]) * dk0csa_1[k2]).real();
			A_csa_ac[k1 * 5 + k2] += rho2 * rho2 * (conj(dk0csa_2[k1]) * dk0csa_2[k2]).real();
			A_csa_ac[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0csa_1[k1]) * dk0csa_2[k2]).real();
			A_csa_ac[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0csa_2[k1]) * dk0csa_1[k2]).real();

			B_csa_ac[k1 * 5 + k2]  = rho1 * rho2 * (conj(dk0csa_1[k1]) * dk0csa_1[k2]).real();
			B_csa_ac[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0csa_2[k1]) * dk0csa_2[k2]).real();
			B_csa_ac[k1 * 5 + k2] -= rho1 * rho2 * (conj(dk0csa_1[k1]) * dk0csa_2[k2]).real();
			B_csa_ac[k1 * 5 + k2] -= rho1 * rho2 * (conj(dk0csa_2[k1]) * dk0csa_1[k2]).real();

			// 2. Dip-Dip acf 1st hydrogen
			A_dip_ac1[k1 * 5 + k2]  = rho1 * rho1 * (conj(dk0dip_1[k1]) * dk0dip_1[k2]).real();
			A_dip_ac1[k1 * 5 + k2] += rho2 * rho2 * (conj(dk0dip_2[k1]) * dk0dip_2[k2]).real();
			A_dip_ac1[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0dip_1[k1]) * dk0dip_2[k2]).real();
			A_dip_ac1[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0dip_2[k1]) * dk0dip_1[k2]).real();

			B_dip_ac1[k1 * 5 + k2]  = rho1 * rho2 * (conj(dk0dip_1[k1]) * dk0dip_1[k2]).real();
			B_dip_ac1[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0dip_2[k1]) * dk0dip_2[k2]).real();
			B_dip_ac1[k1 * 5 + k2] -= rho1 * rho2 * (conj(dk0dip_1[k1]) * dk0dip_2[k2]).real();
			B_dip_ac1[k1 * 5 + k2] -= rho1 * rho2 * (conj(dk0dip_2[k1]) * dk0dip_1[k2]).real();

			if (nH == 2)
			{
				// 2. Dip-Dip acf 2nd hydrogen
				A_dip_ac2[k1 * 5 + k2]  = rho1 * rho1 * (conj(dk0dip2_1[k1]) * dk0dip2_1[k2]).real();
				A_dip_ac2[k1 * 5 + k2] += rho2 * rho2 * (conj(dk0dip2_2[k1]) * dk0dip2_2[k2]).real();
				A_dip_ac2[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0dip2_1[k1]) * dk0dip2_2[k2]).real();
				A_dip_ac2[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0dip2_2[k1]) * dk0dip2_1[k2]).real();

				B_dip_ac2[k1 * 5 + k2]  = rho1 * rho2 * (conj(dk0dip2_1[k1]) * dk0dip2_1[k2]).real();
				B_dip_ac2[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0dip2_2[k1]) * dk0dip2_2[k2]).real();
				B_dip_ac2[k1 * 5 + k2] -= rho1 * rho2 * (conj(dk0dip2_1[k1]) * dk0dip2_2[k2]).real();
				B_dip_ac2[k1 * 5 + k2] -= rho1 * rho2 * (conj(dk0dip2_2[k1]) * dk0dip2_1[k2]).real();

				// 2. Dip-Dip ccf hydrogens 1 - 2
				A_dip_cc12[k1 * 5 + k2]  = rho1 * rho1 * (conj(dk0dip_1[k1]) * dk0dip2_1[k2]).real();
				A_dip_cc12[k1 * 5 + k2] += rho2 * rho2 * (conj(dk0dip_2[k1]) * dk0dip2_2[k2]).real();
				A_dip_cc12[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0dip_1[k1]) * dk0dip2_2[k2]).real();
				A_dip_cc12[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0dip_2[k1]) * dk0dip2_1[k2]).real();

				B_dip_cc12[k1 * 5 + k2]  = rho1 * rho2 * (conj(dk0dip_1[k1]) * dk0dip2_1[k2]).real();
				B_dip_cc12[k1 * 5 + k2] += rho1 * rho2 * (conj(dk0dip_2[k1]) * dk0dip2_2[k2]).real();
				B_dip_cc12[k1 * 5 + k2] -= rho1 * rho2 * (conj(dk0dip_1[k1]) * dk0dip2_2[k2]).real();
				B_dip_cc12[k1 * 5 + k2] -= rho1 * rho2 * (conj(dk0dip_2[k1]) * dk0dip2_1[k2]).real();
			}
		}
	}

}

/************************************/
/* Create big J's for TS-X model */
/************************************/

void relax::makeBigJ_tssrls(int ifield)
{
	bigJdip.clear();
	bigJcsa.clear();
	bigJdip2.clear();
	bigKdip.clear();

	for (unsigned int i = 0; i < freq[ifield].size(); i++)
	{
		
		bigJdip.push_back(0.0);
		bigJcsa.push_back(0.0);
		bigJdip2.push_back(0.0);
		bigKdip.push_back(0.0);
		
		/***************
		 * K = K' part *
		 ***************/

		/* (0, 0) */
		bigJdip[i] += smallj[0][i] * A_dip_ac1[12] + smalljB[0][i] * B_dip_ac1[12];
		bigJcsa[i] += smallj[0][i] * A_csa_ac[12]  + smalljB[0][i] * B_csa_ac[12];
		if (nH == 2)
		{
			bigJdip2[i] += smallj[0][i] * A_dip_ac2[12]  + smalljB[0][i] * B_dip_ac2[12];
			bigKdip[i]  += smallj[0][i] * A_dip_cc12[12] + smalljB[0][i] * B_dip_cc12[12];
		}
		/* (-1,-1) + (1,1) */
		bigJdip[i] += smallj[1][i] * (A_dip_ac1[6] + A_dip_ac1[18]) + smalljB[1][i] * (B_dip_ac1[6] + B_dip_ac1[18]);
		bigJcsa[i] += smallj[1][i] * (A_csa_ac[6] +  A_csa_ac[18])  + smalljB[1][i] * (B_csa_ac[6] +  B_csa_ac[18]);
		if (nH == 2)
		{
			bigJdip2[i] += smallj[1][i] * (A_dip_ac2[6]  + A_dip_ac2[18])  + smalljB[1][i] * (B_dip_ac2[6]  + B_dip_ac2[18]);
			bigKdip[i]  += smallj[1][i] * (A_dip_cc12[6] + A_dip_cc12[18]) + smalljB[1][i] * (B_dip_cc12[6] + B_dip_cc12[18]);
		}
		/* (-2,-2) + (2,2) */
		bigJdip[i] += smallj[2][i] * (A_dip_ac1[0] + A_dip_ac1[24]) + smalljB[2][i] * (B_dip_ac1[0] + B_dip_ac1[24]);
		bigJcsa[i] += smallj[2][i] * (A_csa_ac[0]  + A_csa_ac[24])  + smalljB[2][i] * (B_csa_ac[0]  + B_csa_ac[24]);
		if (nH == 2)
		{
			bigJdip2[i] += smallj[2][i] * (A_dip_ac2[0]  + A_dip_ac2[24])  + smalljB[2][i] * (B_dip_ac2[0]  + B_dip_ac2[24]);
			bigKdip[i]  += smallj[2][i] * (A_dip_cc12[0] + A_dip_cc12[24]) + smalljB[2][i] * (B_dip_cc12[0] + B_dip_cc12[24]);
		}

		/* DEBUGGING */
		/*	
		std::cout << "SMALL J 0 (" << i+1<<") = " << smallj[0][i]<<std::endl;
		std::cout << "SMALL J 1 (" << i+1<<") = " << smallj[1][i]<<std::endl;
		std::cout << "SMALL J 2 (" << i+1<<") = " << smallj[2][i]<<std::endl;
		*/

		/****************
		 * K != K' part *
		 ****************/

		//     j(K,K') = (-)^(K-K') j(-K,-K')    // CORRECT?
		
		if (sd.size() > 3)
		{
			/* (-2, 2) */
			bigJdip[i] += 2.0 * smallj[3][i] * A_dip_ac1[4] + 2.0 * smalljB[3][i] * B_dip_ac1[4];
			bigJcsa[i] += 2.0 * smallj[3][i] * A_csa_ac[4]  + 2.0 * smalljB[3][i] * B_csa_ac[4];
			if (nH == 2)
			{
				bigJdip2[i] += 2.0 * smallj[3][i] * A_dip_ac2[4]  + 2.0 * smalljB[3][i] * B_dip_ac2[4];
				bigKdip[i]  += 2.0 * smallj[3][i] * A_dip_cc12[4] + 2.0 * smalljB[3][i] * B_dip_cc12[4];
			}
			/* (-1, 1) */
			bigJdip[i] += 2.0 * smallj[4][i] * A_dip_ac1[8] + 2.0 * smalljB[4][i] * B_dip_ac1[8];
			bigJcsa[i] += 2.0 * smallj[4][i] * A_csa_ac[8]  + 2.0 * smalljB[4][i] * B_csa_ac[8];
			if (nH == 2)
			{
				bigJdip2[i] += 2.0 * smallj[4][i] * A_dip_ac2[8]  + 2.0 * smalljB[4][i] * B_dip_ac2[8];
				bigKdip[i]  += 2.0 * smallj[4][i] * A_dip_cc12[8] + 2.0 * smalljB[4][i] * B_dip_cc12[8];
			}
			/* (-2, 0) + (0, 2) */
			bigJdip[i] += 2.0 * smallj[5][i] * (A_dip_ac1[2] + A_dip_ac1[14]) + 2.0 * smalljB[5][i] * (B_dip_ac1[2] + B_dip_ac1[14]);
			bigJcsa[i] += 2.0 * smallj[5][i] * (A_csa_ac[2]  + A_csa_ac[14])  + 2.0 * smalljB[5][i] * (B_csa_ac[2]  + B_csa_ac[14]);
			if (nH == 2)
			{
				bigJdip2[i] += 2.0 * smallj[5][i] * (A_dip_ac2[2]  + A_dip_ac2[14])  + 2.0 * smalljB[5][i] * (B_dip_ac2[2]  + B_dip_ac2[14]);
				bigKdip[i]  += 2.0 * smallj[5][i] * (A_dip_cc12[2] + A_dip_cc12[14]) + 2.0 * smalljB[5][i] * (B_dip_cc12[2] + B_dip_cc12[14]);
			}

			/* DEBUGGING */
			/*
			std::cout << "SMALL J 3 (" << i+1<<") = " << smallj[3][i]<<std::endl;
			std::cout << "SMALL J 4 (" << i+1<<") = " << smallj[4][i]<<std::endl;
			std::cout << "SMALL J 5 (" << i+1<<") = " << smallj[5][i]<<std::endl;
			*/

			if (sd.size() > 6)
			{
				/* (-1, 2) + (-2, 1) */
				bigJdip[i] += 2.0 * smallj[6][i] * (A_dip_ac1[9] - A_dip_ac1[3]) + 2.0 * smalljB[6][i] * (B_dip_ac1[9] - B_dip_ac1[3]);
				bigJcsa[i] += 2.0 * smallj[6][i] * (A_csa_ac[9]  - A_csa_ac[3])  + 2.0 * smalljB[6][i] * (B_csa_ac[9]  - B_csa_ac[3]);
				if (nH == 2)
				{
					bigJdip2[i] += 2.0 * smallj[6][i] * (A_dip_ac2[9]  - A_dip_ac2[3])  + 2.0 * smalljB[6][i] * (B_dip_ac2[9]  - B_dip_ac2[3]);
					bigKdip[i]  += 2.0 * smallj[6][i] * (A_dip_cc12[9] - A_dip_cc12[3]) + 2.0 * smalljB[6][i] * (B_dip_cc12[9] - B_dip_cc12[3]);
				}
				/* (0, 1) + (-1, 0) */
				bigJdip[i] += 2.0 * smallj[7][i] * (A_dip_ac1[13] - A_dip_ac1[7]) + 2.0 * smalljB[7][i] * (B_dip_ac1[13] - B_dip_ac1[7]);
				bigJcsa[i] += 2.0 * smallj[7][i] * (A_csa_ac[13]  - A_csa_ac[7])  + 2.0 * smalljB[7][i] * (B_csa_ac[13]  - B_csa_ac[7]);
				if (nH == 2)
				{
					bigJdip2[i] += 2.0 * smallj[7][i] * (A_dip_ac2[13]  - A_dip_ac2[7])  + 2.0 * smalljB[7][i] * (B_dip_ac2[13]  - B_dip_ac2[7]);
					bigKdip[i]  += 2.0 * smallj[7][i] * (A_dip_cc12[13] - A_dip_cc12[7]) + 2.0 * smalljB[7][i] * (B_dip_cc12[13] - B_dip_cc12[7]);;
				}
				/* (1, 2) + (-2, -1) */
				bigJdip[i] += 2.0 * smallj[8][i] * (A_dip_ac1[19] - A_dip_ac1[1]) + 2.0 * smalljB[8][i] * (B_dip_ac1[19] - B_dip_ac1[1]);
				bigJcsa[i] += 2.0 * smallj[8][i] * (A_csa_ac[19]  - A_csa_ac[1])  + 2.0 * smalljB[8][i] * (B_csa_ac[19]  - B_csa_ac[1]);
				if (nH == 2)
				{
					bigJdip2[i] += 2.0 * smallj[8][i] * (A_dip_ac2[19]  - A_dip_ac2[1])  + 2.0 * smalljB[8][i] * (B_dip_ac2[19]  - B_dip_ac2[1]);
					bigKdip[i]  += 2.0 * smallj[8][i] * (A_dip_cc12[19] - A_dip_cc12[1]) + 2.0 * smalljB[8][i] * (B_dip_cc12[19] - B_dip_cc12[1]);;
				}

				/* DEBUGGING */
				/*
				std::cout << "SMALL J 6 (" << i+1<<") = " << smallj[6][i]<<std::endl;
				std::cout << "SMALL J 7 (" << i+1<<") = " << smallj[7][i]<<std::endl;
				std::cout << "SMALL J 8 (" << i+1<<") = " << smallj[8][i]<<std::endl;
				*/

			}
		}
	}

	return;
}

/**********************************/
/* Create big J's for 3S-FB model */
/**********************************/

void relax::makeBigJ_3SFB(int ifield)
{
	bigJdip.clear();
	bigJcsa.clear();
	bigJdip2.clear();
	bigKdip.clear();

	for (unsigned int i = 0; i < freq[ifield].size(); i++)
	{
		
		bigJdip.push_back(0.0);
		bigJcsa.push_back(0.0);
		bigJdip2.push_back(0.0);
		bigKdip.push_back(0.0);
		
		/***************
		 * K = K' part *
		 ***************/

		/* (0, 0) */
		bigJdip[i] += smallj_1[0][i] * (conj(X_Dip1[0 * 3 + 2]) * X_Dip1[0 * 3 + 2]).real() + smallj_2[0][i] * (conj(X_Dip1[1 * 3 + 2]) * X_Dip1[1 * 3 + 2]).real() + smallj_3[0][i] * (conj(X_Dip1[2 * 3 + 2]) * X_Dip1[2 * 3 + 2]).real();
		bigJcsa[i] += smallj_1[0][i] * (conj(X_CSA[0 * 3 + 2]) * X_CSA[0 * 3 + 2]).real() + smallj_2[0][i] * (conj(X_CSA[1 * 3 + 2]) * X_CSA[1 * 3 + 2]).real() + smallj_3[0][i] * (conj(X_CSA[2 * 3 + 2]) * X_CSA[2 * 3 + 2]).real();

		if (nH == 2)
		{
			bigJdip2[i] += smallj_1[0][i] * (conj(X_Dip2[0 * 3 + 2]) * X_Dip2[0 * 3 + 2]).real() + smallj_2[0][i] * (conj(X_Dip2[1 * 3 + 2]) * X_Dip2[1 * 3 + 2]).real() + smallj_3[0][i] * (conj(X_Dip2[2 * 3 + 2]) * X_Dip2[2 * 3 + 2]).real();
			bigKdip[i]  += smallj_1[0][i] * (conj(X_Dip1[0 * 3 + 2]) * X_Dip2[0 * 3 + 2]).real() + smallj_2[0][i] * (conj(X_Dip1[1 * 3 + 2]) * X_Dip2[1 * 3 + 2]).real() + smallj_3[0][i] * (conj(X_Dip1[2 * 3 + 2]) * X_Dip2[2 * 3 + 2]).real();;
		}
		/* (-1,-1) + (1,1) */
		bigJdip[i] += smallj_1[1][i] * ((conj(X_Dip1[0 * 3 + 1]) * X_Dip1[0 * 3 + 1]).real() + (conj(X_Dip1[0 * 3 + 3]) * X_Dip1[0 * 3 + 3]).real());
		bigJdip[i] += smallj_2[1][i] * ((conj(X_Dip1[1 * 3 + 1]) * X_Dip1[1 * 3 + 1]).real() + (conj(X_Dip1[1 * 3 + 3]) * X_Dip1[1 * 3 + 3]).real());
		bigJdip[i] += smallj_3[1][i] * ((conj(X_Dip1[2 * 3 + 1]) * X_Dip1[2 * 3 + 1]).real() + (conj(X_Dip1[2 * 3 + 3]) * X_Dip1[2 * 3 + 3]).real());

		bigJcsa[i] += smallj_1[1][i] * ((conj(X_CSA[0 * 3 + 1]) * X_CSA[0 * 3 + 1]).real() + (conj(X_CSA[0 * 3 + 3]) * X_CSA[0 * 3 + 3]).real());
		bigJcsa[i] += smallj_2[1][i] * ((conj(X_CSA[1 * 3 + 1]) * X_CSA[1 * 3 + 1]).real() + (conj(X_CSA[1 * 3 + 3]) * X_CSA[1 * 3 + 3]).real());
		bigJcsa[i] += smallj_3[1][i] * ((conj(X_CSA[2 * 3 + 1]) * X_CSA[2 * 3 + 1]).real() + (conj(X_CSA[2 * 3 + 3]) * X_CSA[2 * 3 + 3]).real());

		if (nH == 2)
		{
			bigJdip2[i] += smallj_1[1][i] * ((conj(X_Dip2[0 * 3 + 1]) * X_Dip2[0 * 3 + 1]).real() + (conj(X_Dip2[0 * 3 + 3]) * X_Dip2[0 * 3 + 3]).real());
			bigJdip2[i] += smallj_2[1][i] * ((conj(X_Dip2[1 * 3 + 1]) * X_Dip2[1 * 3 + 1]).real() + (conj(X_Dip2[1 * 3 + 3]) * X_Dip2[1 * 3 + 3]).real());
			bigJdip2[i] += smallj_3[1][i] * ((conj(X_Dip2[2 * 3 + 1]) * X_Dip2[2 * 3 + 1]).real() + (conj(X_Dip2[2 * 3 + 3]) * X_Dip2[2 * 3 + 3]).real());

			bigKdip[i]  += smallj_1[1][i] * ((conj(X_Dip1[0 * 3 + 1]) * X_Dip2[0 * 3 + 1]).real() + (conj(X_Dip1[0 * 3 + 3]) * X_Dip2[0 * 3 + 3]).real());
			bigKdip[i]  += smallj_2[1][i] * ((conj(X_Dip1[1 * 3 + 1]) * X_Dip2[1 * 3 + 1]).real() + (conj(X_Dip1[1 * 3 + 3]) * X_Dip2[1 * 3 + 3]).real());
			bigKdip[i]  += smallj_3[1][i] * ((conj(X_Dip1[2 * 3 + 1]) * X_Dip2[2 * 3 + 1]).real() + (conj(X_Dip1[2 * 3 + 3]) * X_Dip2[2 * 3 + 3]).real());
		}
		/* (-2,-2) + (2,2) */
		bigJdip[i] += smallj_1[2][i] * ((conj(X_Dip1[0 * 3 + 0]) * X_Dip1[0 * 3 + 0]).real() + (conj(X_Dip1[0 * 3 + 4]) * X_Dip1[0 * 3 + 4]).real());
		bigJdip[i] += smallj_2[2][i] * ((conj(X_Dip1[1 * 3 + 0]) * X_Dip1[1 * 3 + 0]).real() + (conj(X_Dip1[1 * 3 + 4]) * X_Dip1[1 * 3 + 4]).real());
		bigJdip[i] += smallj_3[2][i] * ((conj(X_Dip1[2 * 3 + 0]) * X_Dip1[2 * 3 + 0]).real() + (conj(X_Dip1[2 * 3 + 4]) * X_Dip1[2 * 3 + 4]).real());

		bigJcsa[i] += smallj_1[2][i] * ((conj(X_CSA[0 * 3 + 0]) * X_CSA[0 * 3 + 0]).real() + (conj(X_CSA[0 * 3 + 4]) * X_CSA[0 * 3 + 4]).real());
		bigJcsa[i] += smallj_2[2][i] * ((conj(X_CSA[1 * 3 + 0]) * X_CSA[1 * 3 + 0]).real() + (conj(X_CSA[1 * 3 + 4]) * X_CSA[1 * 3 + 4]).real());
		bigJcsa[i] += smallj_3[2][i] * ((conj(X_CSA[2 * 3 + 0]) * X_CSA[2 * 3 + 0]).real() + (conj(X_CSA[2 * 3 + 4]) * X_CSA[2 * 3 + 4]).real());

		if (nH == 2)
		{
			bigJdip2[i] += smallj_1[2][i] * ((conj(X_Dip2[0 * 3 + 0]) * X_Dip2[0 * 3 + 0]).real() + (conj(X_Dip2[0 * 3 + 4]) * X_Dip2[0 * 3 + 4]).real());
			bigJdip2[i] += smallj_2[2][i] * ((conj(X_Dip2[1 * 3 + 0]) * X_Dip2[1 * 3 + 0]).real() + (conj(X_Dip2[1 * 3 + 4]) * X_Dip2[1 * 3 + 4]).real());
			bigJdip2[i] += smallj_3[2][i] * ((conj(X_Dip2[2 * 3 + 0]) * X_Dip2[2 * 3 + 0]).real() + (conj(X_Dip2[2 * 3 + 4]) * X_Dip2[2 * 3 + 4]).real());

			bigKdip[i]  += smallj_1[2][i] * ((conj(X_Dip1[0 * 3 + 0]) * X_Dip2[0 * 3 + 0]).real() + (conj(X_Dip1[0 * 3 + 4]) * X_Dip2[0 * 3 + 4]).real());
			bigKdip[i]  += smallj_2[2][i] * ((conj(X_Dip1[1 * 3 + 0]) * X_Dip2[1 * 3 + 0]).real() + (conj(X_Dip1[1 * 3 + 4]) * X_Dip2[1 * 3 + 4]).real());
			bigKdip[i]  += smallj_3[2][i] * ((conj(X_Dip1[2 * 3 + 0]) * X_Dip2[2 * 3 + 0]).real() + (conj(X_Dip1[2 * 3 + 4]) * X_Dip2[2 * 3 + 4]).real());
		}

		/* DEBUGGING */
		/*	
		std::cout << "SMALL J 0 (" << i+1<<") = " << smallj[0][i]<<std::endl;
		std::cout << "SMALL J 1 (" << i+1<<") = " << smallj[1][i]<<std::endl;
		std::cout << "SMALL J 2 (" << i+1<<") = " << smallj[2][i]<<std::endl;
		*/

		/****************
		 * K != K' part *
		 ****************/

		//     j(K,K') = (-)^(K-K') j(-K,-K')    // CORRECT?
		
		if (sd.size() > 3)
		{
			/* (-2, 2) */
			bigJdip[i] += 2.0 * smallj_1[3][i] * (conj(X_Dip1[0 * 3 + 0]) * X_Dip1[0 * 3 + 4]).real();
			bigJdip[i] += 2.0 * smallj_2[3][i] * (conj(X_Dip1[1 * 3 + 0]) * X_Dip1[1 * 3 + 4]).real();
			bigJdip[i] += 2.0 * smallj_3[3][i] * (conj(X_Dip1[2 * 3 + 0]) * X_Dip1[2 * 3 + 4]).real();

			bigJcsa[i] += 2.0 * smallj_1[3][i] * (conj(X_CSA[0 * 3 + 0]) * X_CSA[0 * 3 + 4]).real();
			bigJcsa[i] += 2.0 * smallj_2[3][i] * (conj(X_CSA[1 * 3 + 0]) * X_CSA[1 * 3 + 4]).real();
			bigJcsa[i] += 2.0 * smallj_3[3][i] * (conj(X_CSA[2 * 3 + 0]) * X_CSA[2 * 3 + 4]).real();

			if (nH == 2)
			{
				bigJdip2[i] += 2.0 * smallj_1[3][i] * (conj(X_Dip2[0 * 3 + 0]) * X_Dip2[0 * 3 + 4]).real();
				bigJdip2[i] += 2.0 * smallj_2[3][i] * (conj(X_Dip2[1 * 3 + 0]) * X_Dip2[1 * 3 + 4]).real();
				bigJdip2[i] += 2.0 * smallj_3[3][i] * (conj(X_Dip2[2 * 3 + 0]) * X_Dip2[2 * 3 + 4]).real();

				bigKdip[i]  += 2.0 * smallj_1[3][i] * (conj(X_Dip1[0 * 3 + 0]) * X_Dip2[0 * 3 + 4]).real();
				bigKdip[i]  += 2.0 * smallj_2[3][i] * (conj(X_Dip1[1 * 3 + 0]) * X_Dip2[1 * 3 + 4]).real();
				bigKdip[i]  += 2.0 * smallj_3[3][i] * (conj(X_Dip1[2 * 3 + 0]) * X_Dip2[2 * 3 + 4]).real();
			}
			/* (-1, 1) */
			bigJdip[i] += 2.0 * smallj_1[4][i] * (conj(X_Dip1[0 * 3 + 1]) * X_Dip1[0 * 3 + 3]).real();
			bigJdip[i] += 2.0 * smallj_2[4][i] * (conj(X_Dip1[1 * 3 + 1]) * X_Dip1[1 * 3 + 3]).real();
			bigJdip[i] += 2.0 * smallj_3[4][i] * (conj(X_Dip1[2 * 3 + 1]) * X_Dip1[2 * 3 + 3]).real();

			bigJcsa[i] += 2.0 * smallj_1[4][i] * (conj(X_CSA[0 * 3 + 1]) * X_CSA[0 * 3 + 3]).real();
			bigJcsa[i] += 2.0 * smallj_2[4][i] * (conj(X_CSA[1 * 3 + 1]) * X_CSA[1 * 3 + 3]).real();
			bigJcsa[i] += 2.0 * smallj_3[4][i] * (conj(X_CSA[2 * 3 + 1]) * X_CSA[2 * 3 + 3]).real();

			if (nH == 2)
			{
				bigJdip2[i] += 2.0 * smallj_1[4][i] * (conj(X_Dip2[0 * 3 + 1]) * X_Dip2[0 * 3 + 3]).real();
				bigJdip2[i] += 2.0 * smallj_2[4][i] * (conj(X_Dip2[1 * 3 + 1]) * X_Dip2[1 * 3 + 3]).real();
				bigJdip2[i] += 2.0 * smallj_3[4][i] * (conj(X_Dip2[2 * 3 + 1]) * X_Dip2[2 * 3 + 3]).real();

				bigKdip[i]  += 2.0 * smallj_1[4][i] * (conj(X_Dip1[0 * 3 + 1]) * X_Dip2[0 * 3 + 3]).real();
				bigKdip[i]  += 2.0 * smallj_2[4][i] * (conj(X_Dip1[1 * 3 + 1]) * X_Dip2[1 * 3 + 3]).real();
				bigKdip[i]  += 2.0 * smallj_3[4][i] * (conj(X_Dip1[2 * 3 + 1]) * X_Dip2[2 * 3 + 3]).real();
			}
			/* (-2, 0) + (0, 2) */

			bigJdip[i] += 2.0 * smallj_1[5][i] * ((conj(X_Dip1[0 * 3 + 0]) * X_Dip1[0 * 3 + 2]).real() + (conj(X_Dip1[0 * 3 + 2]) * X_Dip1[0 * 3 + 4]).real());
			bigJdip[i] += 2.0 * smallj_2[5][i] * ((conj(X_Dip1[1 * 3 + 0]) * X_Dip1[1 * 3 + 2]).real() + (conj(X_Dip1[1 * 3 + 2]) * X_Dip1[1 * 3 + 4]).real());
			bigJdip[i] += 2.0 * smallj_3[5][i] * ((conj(X_Dip1[2 * 3 + 0]) * X_Dip1[2 * 3 + 2]).real() + (conj(X_Dip1[2 * 3 + 2]) * X_Dip1[2 * 3 + 4]).real());

			bigJcsa[i] += 2.0 * smallj_1[5][i] * ((conj(X_CSA[0 * 3 + 0]) * X_CSA[0 * 3 + 2]).real() + (conj(X_CSA[0 * 3 + 2]) * X_CSA[0 * 3 + 4]).real());
			bigJcsa[i] += 2.0 * smallj_2[5][i] * ((conj(X_CSA[1 * 3 + 0]) * X_CSA[1 * 3 + 2]).real() + (conj(X_CSA[1 * 3 + 2]) * X_CSA[1 * 3 + 4]).real());
			bigJcsa[i] += 2.0 * smallj_3[5][i] * ((conj(X_CSA[2 * 3 + 0]) * X_CSA[2 * 3 + 2]).real() + (conj(X_CSA[2 * 3 + 2]) * X_CSA[2 * 3 + 4]).real());

			if (nH == 2)
			{
				bigJdip2[i] += 2.0 * smallj_1[5][i] * ((conj(X_Dip2[0 * 3 + 0]) * X_Dip2[0 * 3 + 2]).real() + (conj(X_Dip2[0 * 3 + 2]) * X_Dip2[0 * 3 + 4]).real());
				bigJdip2[i] += 2.0 * smallj_2[5][i] * ((conj(X_Dip2[1 * 3 + 0]) * X_Dip2[1 * 3 + 2]).real() + (conj(X_Dip2[1 * 3 + 2]) * X_Dip2[1 * 3 + 4]).real());
				bigJdip2[i] += 2.0 * smallj_3[5][i] * ((conj(X_Dip2[2 * 3 + 0]) * X_Dip2[2 * 3 + 2]).real() + (conj(X_Dip2[2 * 3 + 2]) * X_Dip2[2 * 3 + 4]).real());

				bigKdip[i]  += 2.0 * smallj_1[5][i] * ((conj(X_Dip1[0 * 3 + 0]) * X_Dip2[0 * 3 + 2]).real() + (conj(X_Dip1[0 * 3 + 2]) * X_Dip2[0 * 3 + 4]).real());
				bigKdip[i]  += 2.0 * smallj_2[5][i] * ((conj(X_Dip1[1 * 3 + 0]) * X_Dip2[1 * 3 + 2]).real() + (conj(X_Dip1[1 * 3 + 2]) * X_Dip2[1 * 3 + 4]).real());
				bigKdip[i]  += 2.0 * smallj_3[5][i] * ((conj(X_Dip1[2 * 3 + 0]) * X_Dip2[2 * 3 + 2]).real() + (conj(X_Dip1[2 * 3 + 2]) * X_Dip2[2 * 3 + 4]).real());
			}

			/* DEBUGGING */
			/*
			std::cout << "SMALL J 3 (" << i+1<<") = " << smallj[3][i]<<std::endl;
			std::cout << "SMALL J 4 (" << i+1<<") = " << smallj[4][i]<<std::endl;
			std::cout << "SMALL J 5 (" << i+1<<") = " << smallj[5][i]<<std::endl;
			*/

			if (sd.size() > 6)
			{
				/* (-1, 2) + (-2, 1) */
				bigJdip[i] += 2.0 * smallj_1[6][i] * ((conj(X_Dip1[0 * 3 + 1]) * X_Dip1[0 * 3 + 4]).real() + (conj(X_Dip1[0 * 3 + 0]) * X_Dip1[0 * 3 + 3]).real());
				bigJdip[i] += 2.0 * smallj_2[6][i] * ((conj(X_Dip1[1 * 3 + 1]) * X_Dip1[1 * 3 + 4]).real() + (conj(X_Dip1[1 * 3 + 0]) * X_Dip1[1 * 3 + 3]).real());
				bigJdip[i] += 2.0 * smallj_3[6][i] * ((conj(X_Dip1[2 * 3 + 1]) * X_Dip1[2 * 3 + 4]).real() + (conj(X_Dip1[2 * 3 + 0]) * X_Dip1[2 * 3 + 3]).real());

				bigJcsa[i] += 2.0 * smallj_1[6][i] * ((conj(X_CSA[0 * 3 + 1]) * X_CSA[0 * 3 + 4]).real() + (conj(X_CSA[0 * 3 + 0]) * X_CSA[0 * 3 + 3]).real());
				bigJcsa[i] += 2.0 * smallj_2[6][i] * ((conj(X_CSA[1 * 3 + 1]) * X_CSA[1 * 3 + 4]).real() + (conj(X_CSA[1 * 3 + 0]) * X_CSA[1 * 3 + 3]).real());
				bigJcsa[i] += 2.0 * smallj_3[6][i] * ((conj(X_CSA[2 * 3 + 1]) * X_CSA[2 * 3 + 4]).real() + (conj(X_CSA[2 * 3 + 0]) * X_CSA[2 * 3 + 3]).real());

				if (nH == 2)
				{
					bigJdip2[i] += 2.0 * smallj_1[6][i] * ((conj(X_Dip2[0 * 3 + 1]) * X_Dip2[0 * 3 + 4]).real() + (conj(X_Dip2[0 * 3 + 0]) * X_Dip2[0 * 3 + 3]).real());
					bigJdip2[i] += 2.0 * smallj_2[6][i] * ((conj(X_Dip2[1 * 3 + 1]) * X_Dip2[1 * 3 + 4]).real() + (conj(X_Dip2[1 * 3 + 0]) * X_Dip2[1 * 3 + 3]).real());
					bigJdip2[i] += 2.0 * smallj_3[6][i] * ((conj(X_Dip2[2 * 3 + 1]) * X_Dip2[2 * 3 + 4]).real() + (conj(X_Dip2[2 * 3 + 0]) * X_Dip2[2 * 3 + 3]).real());

					bigKdip[i]  += 2.0 * smallj_1[6][i] * ((conj(X_Dip1[0 * 3 + 1]) * X_Dip2[0 * 3 + 4]).real() + (conj(X_Dip1[0 * 3 + 0]) * X_Dip2[0 * 3 + 3]).real());
					bigKdip[i]  += 2.0 * smallj_2[6][i] * ((conj(X_Dip1[1 * 3 + 1]) * X_Dip2[1 * 3 + 4]).real() + (conj(X_Dip1[1 * 3 + 0]) * X_Dip2[1 * 3 + 3]).real());
					bigKdip[i]  += 2.0 * smallj_3[6][i] * ((conj(X_Dip1[2 * 3 + 1]) * X_Dip2[2 * 3 + 4]).real() + (conj(X_Dip1[2 * 3 + 0]) * X_Dip2[2 * 3 + 3]).real());
				}
				/* (0, 1) + (-1, 0) */
				bigJdip[i] += 2.0 * smallj_1[7][i] * ((conj(X_Dip1[0 * 3 + 2]) * X_Dip1[0 * 3 + 3]).real() + (conj(X_Dip1[0 * 3 + 1]) * X_Dip1[0 * 3 + 2]).real());
				bigJdip[i] += 2.0 * smallj_2[7][i] * ((conj(X_Dip1[1 * 3 + 2]) * X_Dip1[1 * 3 + 3]).real() + (conj(X_Dip1[1 * 3 + 1]) * X_Dip1[1 * 3 + 2]).real());
				bigJdip[i] += 2.0 * smallj_3[7][i] * ((conj(X_Dip1[2 * 3 + 2]) * X_Dip1[2 * 3 + 3]).real() + (conj(X_Dip1[2 * 3 + 1]) * X_Dip1[2 * 3 + 2]).real());

				bigJcsa[i] += 2.0 * smallj_1[7][i] * ((conj(X_CSA[0 * 3 + 2]) * X_CSA[0 * 3 + 3]).real() + (conj(X_CSA[0 * 3 + 1]) * X_CSA[0 * 3 + 2]).real());
				bigJcsa[i] += 2.0 * smallj_2[7][i] * ((conj(X_CSA[1 * 3 + 2]) * X_CSA[1 * 3 + 3]).real() + (conj(X_CSA[1 * 3 + 1]) * X_CSA[1 * 3 + 2]).real());
				bigJcsa[i] += 2.0 * smallj_3[7][i] * ((conj(X_CSA[2 * 3 + 2]) * X_CSA[2 * 3 + 3]).real() + (conj(X_CSA[2 * 3 + 1]) * X_CSA[2 * 3 + 2]).real());

				if (nH == 2)
				{
					bigJdip2[i] += 2.0 * smallj_1[7][i] * ((conj(X_Dip2[0 * 3 + 2]) * X_Dip2[0 * 3 + 3]).real() + (conj(X_Dip2[0 * 3 + 1]) * X_Dip2[0 * 3 + 2]).real());
					bigJdip2[i] += 2.0 * smallj_2[7][i] * ((conj(X_Dip2[1 * 3 + 2]) * X_Dip2[1 * 3 + 3]).real() + (conj(X_Dip2[1 * 3 + 1]) * X_Dip2[1 * 3 + 2]).real());
					bigJdip2[i] += 2.0 * smallj_3[7][i] * ((conj(X_Dip2[2 * 3 + 2]) * X_Dip2[2 * 3 + 3]).real() + (conj(X_Dip2[2 * 3 + 1]) * X_Dip2[2 * 3 + 2]).real());

					bigKdip[i]  += 2.0 * smallj_1[7][i] * ((conj(X_Dip1[0 * 3 + 2]) * X_Dip2[0 * 3 + 3]).real() + (conj(X_Dip1[0 * 3 + 1]) * X_Dip2[0 * 3 + 2]).real());
					bigKdip[i]  += 2.0 * smallj_2[7][i] * ((conj(X_Dip1[1 * 3 + 2]) * X_Dip2[1 * 3 + 3]).real() + (conj(X_Dip1[1 * 3 + 1]) * X_Dip2[1 * 3 + 2]).real());
					bigKdip[i]  += 2.0 * smallj_3[7][i] * ((conj(X_Dip1[2 * 3 + 2]) * X_Dip2[2 * 3 + 3]).real() + (conj(X_Dip1[2 * 3 + 1]) * X_Dip2[2 * 3 + 2]).real());
				}
				/* (1, 2) + (-2, -1) */
				bigJdip[i] += 2.0 * smallj_1[8][i] * ((conj(X_Dip1[0 * 3 + 3]) * X_Dip1[0 * 3 + 4]).real() + (conj(X_Dip1[0 * 3 + 0]) * X_Dip1[0 * 3 + 1]).real());
				bigJdip[i] += 2.0 * smallj_2[8][i] * ((conj(X_Dip1[1 * 3 + 3]) * X_Dip1[1 * 3 + 4]).real() + (conj(X_Dip1[1 * 3 + 0]) * X_Dip1[1 * 3 + 1]).real());
				bigJdip[i] += 2.0 * smallj_3[8][i] * ((conj(X_Dip1[2 * 3 + 3]) * X_Dip1[2 * 3 + 4]).real() + (conj(X_Dip1[2 * 3 + 0]) * X_Dip1[2 * 3 + 1]).real());

				bigJcsa[i] += 2.0 * smallj_1[8][i] * ((conj(X_CSA[0 * 3 + 3]) * X_CSA[0 * 3 + 4]).real() + (conj(X_CSA[0 * 3 + 0]) * X_CSA[0 * 3 + 1]).real());
				bigJcsa[i] += 2.0 * smallj_2[8][i] * ((conj(X_CSA[1 * 3 + 3]) * X_CSA[1 * 3 + 4]).real() + (conj(X_CSA[1 * 3 + 0]) * X_CSA[1 * 3 + 1]).real());
				bigJcsa[i] += 2.0 * smallj_3[8][i] * ((conj(X_CSA[2 * 3 + 3]) * X_CSA[2 * 3 + 4]).real() + (conj(X_CSA[2 * 3 + 0]) * X_CSA[2 * 3 + 1]).real());

				if (nH == 2)
				{
					bigJdip2[i] += 2.0 * smallj_1[8][i] * ((conj(X_Dip2[0 * 3 + 3]) * X_Dip2[0 * 3 + 4]).real() + (conj(X_Dip2[0 * 3 + 0]) * X_Dip2[0 * 3 + 1]).real());
					bigJdip2[i] += 2.0 * smallj_2[8][i] * ((conj(X_Dip2[1 * 3 + 3]) * X_Dip2[1 * 3 + 4]).real() + (conj(X_Dip2[1 * 3 + 0]) * X_Dip2[1 * 3 + 1]).real());
					bigJdip2[i] += 2.0 * smallj_3[8][i] * ((conj(X_Dip2[2 * 3 + 3]) * X_Dip2[2 * 3 + 4]).real() + (conj(X_Dip2[2 * 3 + 0]) * X_Dip2[2 * 3 + 1]).real());

					bigKdip[i]  += 2.0 * smallj_1[8][i] * ((conj(X_Dip1[0 * 3 + 3]) * X_Dip2[0 * 3 + 4]).real() + (conj(X_Dip1[0 * 3 + 0]) * X_Dip2[0 * 3 + 1]).real());
					bigKdip[i]  += 2.0 * smallj_2[8][i] * ((conj(X_Dip1[1 * 3 + 3]) * X_Dip2[1 * 3 + 4]).real() + (conj(X_Dip1[1 * 3 + 0]) * X_Dip2[1 * 3 + 1]).real());
					bigKdip[i]  += 2.0 * smallj_3[8][i] * ((conj(X_Dip1[2 * 3 + 3]) * X_Dip2[2 * 3 + 4]).real() + (conj(X_Dip1[2 * 3 + 0]) * X_Dip2[2 * 3 + 1]).real());
				}

				/* DEBUGGING */
				/*
				std::cout << "SMALL J 6 (" << i+1<<") = " << smallj[6][i]<<std::endl;
				std::cout << "SMALL J 7 (" << i+1<<") = " << smallj[7][i]<<std::endl;
				std::cout << "SMALL J 8 (" << i+1<<") = " << smallj[8][i]<<std::endl;
				*/

			}
		}
	}

	return;
}


/*********************************************
 * Calculate T1, T2 and NOE at a given field *
 *********************************************/

dvector relax::makeT1T2NOE(int ifield)
{
	double dipFactor, csaFactor;
	dvector relaxationTimes;
	
	dipFactor =  MU0_OVER_4PI*gyroH*gyroN*hbar/(rNH*rNH*rNH);
	dipFactor *= dipFactor;
	dipFactor *= 0.1;
	dipFactor *= freqScale*freqScale;
	dipFactor *= spectrumFactor;
	
	csaFactor = larmorFreqN[ifield]*larmorFreqN[ifield]*deltaCSA*deltaCSA*TWO_OVER_FIFTHEEN*1.0e-12;
	csaFactor /= spectrumFactor;
	
        double exchangeCorr = (larmorFreqH[ifield]/larmorFreqHmin);
        exchangeCorr *= exchangeCorr;
        exchangeCorr *= rex*freqScale;


	std::cout<<std::endl<<"Jdip(0) = " << bigJdip[0] << std::endl;

	double T1 = dipFactor*(bigJdip[3]+3.0*bigJdip[2]+6.0*bigJdip[4]) + csaFactor*bigJcsa[2];
        if (nH ==  2) T1 += dipFactor*(bigJdip2[3]+3.0*bigJdip2[2]+6.0*bigJdip2[4]);
	T1 = 1.0/T1;
	relaxationTimes.push_back(T1*freqScale);

	double T2 = 0.5*dipFactor*(4.0*bigJdip[0]+bigJdip[3]+3.0*bigJdip[2]+6.0*bigJdip[1]+6.0*bigJdip[4])+ONE_OVER_SIX*csaFactor*(3.0*bigJcsa[2]+4.0*bigJcsa[0]) + exchangeCorr;
	if (nH == 2) T2 += 0.5*dipFactor*(4.0*bigJdip2[0]+bigJdip2[3]+3.0*bigJdip2[2]+6.0*bigJdip2[1]+6.0*bigJdip2[4]);
	T2 = 1.0/T2;
	relaxationTimes.push_back(T2*freqScale);

	double NOE = 0.0;
	if (nH == 1)
		NOE = 1.0 + (gyroH/gyroN) * dipFactor * (6.0*bigJdip[4]-bigJdip[3]) * T1;
	else if (nH == 2)
		NOE = 1.0 + (gyroH/gyroN) * dipFactor* ( (6.0*bigJdip[4]-bigJdip[3]) + (6.0*bigJdip2[4]-bigJdip2[3]) ) * T1;
	relaxationTimes.push_back(NOE);
	
	if (nH == 2)
	{
		double dipFactor2 = 10.0*dipFactor/freqScale;
		double CCRRT1 = 0.6 * dipFactor2*bigKdip[2];
		relaxationTimes.push_back(CCRRT1);
		double CCRRT2 = 0.3 * dipFactor2*((4.0/3.0)*bigKdip[0] + bigKdip[2]);
		relaxationTimes.push_back(CCRRT2);
	}

	return relaxationTimes;	
}

 /******************************************
  * Return the T1 T2 and NOE at all fields *
  ******************************************/

dvvector relax::getT1T2NOE(bool readSmallJs)
{
	dvvector relaxationTimes(field.size());

	if (is3SFBCalculation)
	{
		buildXFunctions();
	}
	else if (isTsSrlsCalculation || isTsFb1Calculation)
	{
		buildBigJsWeights();
	}

	for (unsigned int i = 0; i < field.size(); i++)
	{
		/* Build small j's */
		if (!readSmallJs)
		{
			makeSymmJ(i);
			if (is3SFBCalculation)
				makeSmallJ_3SFB(i);
			else
				makeSmallJ(i);
		}
		else
			interpSmallJ(i);

		/* Build big J's */
		if (isTsSrlsCalculation || isTsFb1Calculation)
			makeBigJ_tssrls(i);
		else if (is3SFBCalculation)
			makeBigJ_3SFB(i);
		else
			makeBigJ(i);

		/* Calculate relaxation times */
		relaxationTimes[i] = makeT1T2NOE(i);
	}
	return relaxationTimes;
}

/***********************
 * Return field vector *
 ***********************/

dvector relax::getField(void)
{
        return field;
}

/******************************************************/
/* Read small j spectral densites from a previous run */
/******************************************************/

void relax::readSmallJs(int NP)
{
	int y;
	int jnumber = sd.size();
	char smjn[2];
	double tmp1, tmp2;
	std::string separator;
#ifdef _LINUX_
	separator = "/";
#else
	separator = "\\";
#endif

	std::string fileName;
	fstream file;

	/*************************/
	/* Read stored small j's */
	/*************************/

	ww = dvector(NP,0.0);
	smalljAllFreq = new dvector[jnumber];

	for (int i = 0; i < 2*jnumber; i+=2) {
		y = i>>1; 
		smalljAllFreq[y].clear();
		sprintf(smjn,"%d",y);
		fileName = path + separator + project + "_smallj_" + smjn + ".dat";
		file.open(fileName.c_str(),ios::in);
		smalljAllFreq[y] = dvector(NP,0.0);
		for (int n = 0; n < NP; n++)
		{
			file >> tmp1 >> tmp2;
		        smalljAllFreq[y][n] = tmp2;
			if (i == 0) ww[n] = tmp1;
		}
		file.close();
	}
	return;
}

/*************************************************************************/
/* Evaluate small j's at relevant frequencies using spline interpolation */
/*************************************************************************/

void relax::interpSmallJ(int ifield)
{
	int y;
	double tmp;
	smallj = new dvector[sd.size()];
        for (unsigned int i = 0; i < sd.size(); i++)
        {
                smallj[i].clear();
		spline spl;
		spl.init(rank);
		spl.setExperimentalPoints(ww,smalljAllFreq[i]);
		spl.makeSpline();
                for (unsigned int j = 0; j < freq[ifield].size(); j++)
                {
			// NOTE: ABSOLUTE VALUE IN THE CALL TO SPLINE INTERPOLATION 
			//       IS BECAUSE FREQ CAN BE NEGATIVE, BUT RESULTS ARE
			//       INDEPENDENT ON THE SIGN
			tmp = spl.splint(fabs(freq[ifield].at(j)));
                        smallj[i].push_back(tmp);
                }
	}
	return;
}

/*******************************************************************************/
/* Method to set the sigma of the Gaussian distribution of the H-C-H libration */
/*******************************************************************************/
void relax::setHchSigma(double s)
{
	hchSigma = sigma = s;
	return;
}

/*****************************************************************/
/* Method to perform averaging over libratio of H-C-H bond angle */
/*****************************************************************/
ldcomplex relax::averaged_reduced_d2m0(int m, long double beta12)
{
	int nv = 0, er = 0;
	double ae = 0.0;
	long double norm;
	hch = (double)beta12;
	sigma = hchSigma;
	ldcomplex dlmk;
	norm = 1./(long double)dqag(peqBeta,0.0,M_PI,0.0,1e-10,6,&ae,&nv,&er);
	switch(m)
	{
		case -2:
			dlmk.real(norm * (long double)dqag(avgD220,0.0,M_PI,0.0,1e-10,6,&ae,&nv,&er));
			break;
		case -1:
			dlmk.real(-norm * (long double)dqag(avgD210,0.0,M_PI,0.0,1e-10,6,&ae,&nv,&er));
			break;
		case 0:
			dlmk.real(norm * (long double)dqag(avgD200,0.0,M_PI,0.0,1e-10,6,&ae,&nv,&er));
			break;
		case 1:
			dlmk.real(norm * (long double)dqag(avgD210,0.0,M_PI,0.0,1e-10,6,&ae,&nv,&er));
			break;
		case 2:
			dlmk.real(norm * (long double)dqag(avgD220,0.0,M_PI,0.0,1e-10,6,&ae,&nv,&er));
			break;
	}
	dlmk.imag(0.0);
	return dlmk;
}

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
 Name        : stvec.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class to evaluate representation of starting vector on basis
 ============================================================================
 */
#include "stvec.h"

int stvec::L1  = 0;
int stvec::mK1 = 0;
int stvec::mK  = 0;
int stvec::L2  = 0;
int stvec::K1  = 0;
int stvec::K2  = 0;
int stvec::nMult = 0;
int stvec::nMult2 = 0;
int stvec::tors = 0;
int stvec::nTors = 0;
double stvec::delta;
long double stvec::c20 = 0.0;
long double stvec::c22 = 0.0;
long double stvec::c40 = 0.0;
long double stvec::c42 = 0.0;
long double stvec::c44 = 0.0;
dcvector stvec::coef_int_pot;
potential_fb2 stvec::ufb2;
wigner stvec::w(60);

#undef JEMDICUB

/****************
 *  Constructor *
 ****************/

stvec::stvec()
{
	L1constrain = -1000000;
	isL1constrained = false;
	M1constrain = -1000000;
	isM1constrained = false;
	jjconstrain = -1000000;
	isjjconstrained = false;
	Mmconstrain = -1000000;
	isMmconstrained = false;
	initiated = false;
}

/************** 
 * Destructor *
 **************/

stvec::~stvec()
{
#ifdef WRITE_DESTROY_MESSAGE
	std::cout << "TASK " << rank << ": Cleared stvec object" << std::endl;
#endif
}

/***************
 *  Initiators *
 ***************/

void stvec::init(physics *p)
{
	rank = 0;
	phy = p;
	dynModel = phy->getDynamicsModel();
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
		potential u = phy->getPotentialCoefficients();
		c20 = u.c20.at(0); 
		c22 = u.c22.at(0);
		c40 = u.c40.at(0);
		c42 = u.c42.at(0);
		c44 = u.c44.at(0);
	}
	else if (!(dynModel.compare("fb1")) || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb")))
	{
		nTors = 1;
                coef_int_pot = phy->getFB1Coeff();
	}
        else if (!(dynModel.compare("fb2")))
        {
                nTors = 2;
                ufb2 = phy->getFB2Coeff();
        }

	sd = phy->getSpectralDensities();
	projections = new VECTOR_double[sd.size()<<1];
	npmkk1 = dvector(sd.size()<<1,0.0);
	normpm = dvector(sd.size()<<1,0.0);
	for (unsigned int i = 0; i < sd.size()<<1; i+=2)
	{
		KK =  sd[i>>1] / 5 - 2;
		KK1 = sd[i>>1] % 5 - 2;
		npmkk1[i] = 2.0 + (KK == KK1 ? 2.0 : 0.0);
		npmkk1[i] += (KK == 0 ?  1.0 : 0.0);
		npmkk1[i] += (KK1 == 0 ? 1.0 : 0.0);
		npmkk1[i] += (KK == -KK1 ? (!(KK%2) ? 2.0 : -2.0) : 0.0);
		npmkk1[i] *= 0.5;
		npmkk1[i+1] = 2.0 + (KK == KK1 ? 2.0 : 0.0);
		npmkk1[i+1] -= (KK == 0 ?  1.0 : 0.0);
		npmkk1[i+1] -= (KK1 == 0 ? 1.0 : 0.0);
		npmkk1[i+1] -= (KK == -KK1 ? (!(KK%2) ? 2.0 : -2.0) : 0.0);
		npmkk1[i+1] *= 0.5;
		if (npmkk1[i] > ZERO) normpm[i] = 1.0/sqrt(npmkk1[i]);
		else normpm[i] = npmkk1[i] = 0.0;
		if (npmkk1[i+1] > ZERO) normpm[i+1] = 1.0/sqrt(npmkk1[i+1]);
		else normpm[i+1] = npmkk1[i+1] = 0.0;
	}
	std::cout << "TASK " << rank << ": Created stvec object" << std::endl;
	initiated = true;
	return;
}

void stvec::init(physics *p, int r)
{
	init(p);
	rank = r;
	return;
}

/*********************************************
 * Method for object update during a fitting *
 *********************************************/

void stvec::update(void)
{
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
		potential u = phy->getPotentialCoefficients();
	        c20 = u.c20.at(0);
	        c22 = u.c22.at(0);
	        c40 = u.c40.at(0);
	        c42 = u.c42.at(0);
	        c44 = u.c44.at(0);
	}
	else if (!(dynModel.compare("fb1")) || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb")))
	{
		coef_int_pot = phy->getFB1Coeff();
	}
	else if (!(dynModel.compare("fb2")))
	{
		ufb2 = phy->getFB2Coeff();
	}
	
	npmkk1 = dvector(sd.size()<<1,0.0);
	normpm = dvector(sd.size()<<1,0.0);
	for (unsigned int i = 0; i < sd.size()<<1; i+=2)
	{
		KK =  sd[i>>1] / 5 - 2;
		KK1 = sd[i>>1] % 5 - 2;
		npmkk1[i] = 2.0 + (KK == KK1 ? 2.0 : 0.0);
		npmkk1[i] += (KK == 0 ?  1.0 : 0.0);
		npmkk1[i] += (KK1 == 0 ? 1.0 : 0.0);
		npmkk1[i] += (KK == -KK1 ? (!(KK%2) ? 2.0 : -2.0) : 0.0);
		npmkk1[i] *= 0.5;
		npmkk1[i+1] = 2.0 + (KK == KK1 ? 2.0 : 0.0);
		npmkk1[i+1] -= (KK == 0 ?  1.0 : 0.0);
		npmkk1[i+1] -= (KK1 == 0 ? 1.0 : 0.0);
		npmkk1[i+1] -= (KK == -KK1 ? (!(KK%2) ? 2.0 : -2.0) : 0.0);
		npmkk1[i+1] *= 0.5;
		if (npmkk1[i] > ZERO) normpm[i] = 1.0/sqrt(npmkk1[i]);
		else normpm[i] = npmkk1[i] = 0.0;
		if (npmkk1[i+1] > ZERO) normpm[i+1] = 1.0/sqrt(npmkk1[i+1]);
		else normpm[i+1] = npmkk1[i+1] = 0.0;
	}

	return;
}

void stvec::updatePotential(void)
{
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
		potential u = phy->getPotentialCoefficients();
	        c20 = u.c20.at(0);
	        c22 = u.c22.at(0);
	        c40 = u.c40.at(0);
	        c42 = u.c42.at(0);
	        c44 = u.c44.at(0);
	}
	return;
}

/********************************
 * Return initialization status *
 ********************************/

bool stvec::hasBeenInit(void)
{
	return initiated;
}

/*************************************
 *  Project starting vector on basis *
 *************************************/

void stvec::projectOnBasis(basis *bas, int nr, int gr)
{
#ifdef WRITE_ALL	
#ifdef WRITE_STVEC
	std::cout << "TASK " << rank << ": Porjecting starting vector on basis..." << std::endl;
#endif
#endif
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
		projectOnBasis_SRLS(bas,nr,gr);
	else if (!(dynModel.compare("fb1")) || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb")))
	{
		projectOnBasis_FB1(bas,nr,gr);
	}
	else if (!(dynModel.compare("fb2")))
		projectOnBasis_FB2(bas,nr,gr);
#ifdef WRITE_ALL
        std::cout << "TASK " << rank << ": ... starting vector projections calculated successfully." << std::endl;
#endif
}

/********************
 * FB1 projections *
 ********************/

void stvec::projectOnBasis_FB1(basis *bas, int nr, int gr)
{	
	nrows = nr;
	globalRow = gr;
	ivector idx;
	nbf = bas[0].getNumberOfBasisFunctions();
	int LL,MM,KK,NN,jj;
	int TK, TK1;
	int nv = 0, er = 0;
	double ae = 0.0;
	double normalization, multKK, multNN, deltaFactorK, deltaFactorK1;
    	dvector internalIntegralReal(phy[0].getNmax()+1,0.0);
    	dvector internalIntegralImag(phy[0].getNmax()+1,0.0);
	
	for (int i = (sd.size()<<1)-2; i >=0; i-=2)
	{		
		projections[i] = VECTOR_double(nbf,0.0);
		projections[i+1] = VECTOR_double(nbf,0.0);
	}
	
	unsigned int s = sd.size();
	
	/**********************************
	 * Precalculate useful quantities *
	 **********************************/

	double peqNorm;	
	int nmax = phy[0].getNmax();
 	peqNorm=1.0;
    	for (int n_value = 0; n_value <= nmax; n_value++){
		nMult = n_value;
		internalIntegralReal[n_value]=(double)dqag(internal_integrand_real,0.0,2.0*M_PI,0.0,1e-10,6,&ae,&nv,&er)*peqNorm;
		if (fabs(internalIntegralReal[n_value])<1.0e-13) internalIntegralReal[n_value]=0.0;
      		internalIntegralImag[n_value]=(double)dqag(internal_integrand_imag,0.0,2.0*M_PI,0.0,1e-10,6,&ae,&nv,&er)*peqNorm;
      		if (fabs(internalIntegralImag[n_value])<1.0e-13) internalIntegralImag[n_value]=0.0;
	}

	/*****************************************
	 * Calculate the projection on the basis *
	 *****************************************/
	
	for (int i = 0; i < nrows; i++){
		
		/**********************
		 * Load basis indexes *
		 **********************/
		
		idx = bas[0].getBasisFunction(i+globalRow);
		NN = idx[0]; LL = idx[1]; MM = idx[2]; KK = idx[3]; jj = idx[4];
		
		/*********************************************
		 * Evaluate the projections for Tpm(2,0,K,K')*
		 *********************************************/
		
		normalization = (!KK && !NN ? 0.5 :  SQRT_ONE_OVER_TWO);
		normalization *= SQRT_FIVE; // (2L+1)^(1/2), but here L = 2 only

		/*****************
		 * Evaluate sign *
		 *****************/
			
		multKK =  (!(KK%2)  ? 1.0 : -1.0);
		multNN =  ( NN>=0  ? 1.0 : -1.0);
			
		for (unsigned int v = 0; v < s<<1; v+=2)
		{
				
			TK =  (int)((unsigned int)sd[v>>1] / 5 - 2);
			TK1 = (int)((unsigned int)sd[v>>1] % 5 - 2);

			/********************************** 
			 * Calculate and store projection *
			 **********************************/

			deltaFactorK =  (KK == TK ?  (KK == -TK  ? 1.0 + multKK : 1.0) : (KK == -TK  ? multKK : 0.0));
			deltaFactorK1 = (KK == TK1 ? (KK == -TK1 ? 1.0 + multKK : 1.0) : (KK == -TK1 ? multKK : 0.0));

			projections[v][i+globalRow] = (double)normalization*normpm[v]*(deltaFactorK + deltaFactorK1)* \
							( jj == 1 ? internalIntegralReal[abs(NN)] : multNN*internalIntegralImag[abs(NN)]);
			if (fabs(projections[v][i+globalRow]) <= ZERO)  projections[v][i+globalRow] = 0.0;

			deltaFactorK =  (KK == TK ?  (KK == -TK  ? 1.0 - multKK : 1.0) : (KK == -TK  ? -multKK : 0.0));
			deltaFactorK1 = (KK == TK1 ? (KK == -TK1 ? 1.0 - multKK : 1.0) : (KK == -TK1 ? -multKK : 0.0));

			projections[v+1][i+globalRow] = (double)normalization*normpm[v+1]*(deltaFactorK + deltaFactorK1)* \
							( jj == 1 ? multNN*internalIntegralImag[abs(NN)] : -internalIntegralReal[abs(NN)]);
			if (fabs(projections[v+1][i+globalRow]) <= ZERO)  projections[v+1][i+globalRow] = 0.0;
				
		} /* close cycle over spectral densities */
			
	} /* Close cycle over i */

	return;
}

/********************
 * FB2 projections *
 ********************/

#include "jemdi.h"
double dn1, dn2;
void stvec::projectOnBasis_FB2(basis *bas, int nr, int gr)
{	
	nrows = nr;
	globalRow = gr;
	ivector idx;
	nbf = bas->getNumberOfBasisFunctions();
	int LL,MM,KK,NN1,NN2,jj;
	int TK, TK1;
	int nv = 0, er = 0;
	double ae = 0.0;
	double normalization, multKK, multNN1, multNN2, deltaFactorK, deltaFactorK1;
	
	for (int i = (sd.size()<<1)-2; i >=0; i-=2)
	{		
		projections[i]   = VECTOR_double(nbf,0.0);
		projections[i+1] = VECTOR_double(nbf,0.0);
	}
	
	unsigned int s = sd.size();

#ifdef JEMDICUB
	/*****************************************/
        /* SETUP INTEGRAL CALCULATION WITH JEMDI */
        /*****************************************/

        int nX = 2; // Number of variables
        int M = 5; // Number of repetitions
        int Ntrj = 50; // Number of trajectories per repetition
        int Nstep = 50000; // Number of steps per trajectory
        double val, err; // These two numbers will contain the result of the integral (val) and the error estimation (err)
        double *lb = new double[nX]; // pointer to double array containing lower bounds
        double *ub = new double[nX]; // pointer to double array containing upper bounds
#else
	int maxEval;
	double *lb = new double[2];
	double *ub = new double[2];
	double *mm = new double[2];
	double val, err, relErr;

	phy->getCubatureParams(&maxEval, &relErr);
#endif
	/**********************************
	 * Precalculate useful quantities *
	 **********************************/

	int n1max = phy->getNmax(1);
	int n2max = phy->getNmax(2);
	int ndim1 = 2 * n1max + 1;
	int ndim2 = 2 * n2max + 1;
	int n2UpLim;

    	dvector internalIntegralReal(ndim1 * ndim2, 0.0);
    	dvector internalIntegralImag(ndim1 * ndim2, 0.0);
	dcomplex torsionalIntegral;

	lb[0] = lb[1] = -M_PI; //0.0;
	ub[0] = ub[1] =  M_PI; //2.0 * M_PI;
	
	for (int n1_value = -n1max; n1_value <= 0; n1_value++)
	{
#ifndef JEMDICUB
		mm[0] = (double)n1_value;
#endif
		dn1 = (double)n1_value;
		n2UpLim = n1_value < 0 ? n2max : 0;
		for (int n2_value = -n2max; n2_value <= n2UpLim; n2_value++)
		{
			dn2 = (double)n2_value;
#ifndef JEMDICUB
			mm[1] = (double)n2_value;
			hcubature(1, fb2_integrand_real, (void *)mm, 2, lb, ub, maxEval, 1.0e-15, relErr, ERROR_L2, &val, &err);
			internalIntegralReal[(n1_value + n1max) * ndim2 + (n2_value + n2max)] = internalIntegralReal[(-n1_value + n1max) * ndim2 + (-n2_value + n2max)] = val;
			std::cout << "REAL; " << n1_value << ", " << n2_value << ", " << val << ", " << err << std::endl;

			hcubature(1, fb2_integrand_imag, (void *)mm, 2, lb, ub, maxEval, 1.0e-15, relErr, ERROR_L2, &val, &err);
			std::cout << "IMAG; " << n1_value << ", " << n2_value << ", " << val << ", " << err << std::endl;
			internalIntegralImag[(-n1_value + n1max) * ndim2 + (-n2_value + n2max)] = -(internalIntegralImag[(n1_value + n1max) * ndim2 + (n2_value + n2max)] = val);
#else
			jemg(fb2_integrand_real_j, nX, lb, ub, M, Ntrj, Nstep, &val, &err);
			internalIntegralReal[(n1_value + n1max) * ndim2 + (n2_value + n2max)] = internalIntegralReal[(-n1_value + n1max) * ndim2 + (-n2_value + n2max)] = val;
			jemg(fb2_integrand_imag_j, nX, lb, ub, M, Ntrj, Nstep, &val, &err);
			internalIntegralImag[(-n1_value + n1max) * ndim2 + (-n2_value + n2max)] = -(internalIntegralImag[(n1_value + n1max) * ndim2 + (n2_value + n2max)] = val);
#endif
		}
	}

	/*****************************************
	 * Calculate the projection on the basis *
	 *****************************************/
	
	for (int i = 0; i < nrows; i++){
		
		/**********************
		 * Load basis indexes *
		 **********************/
		
		idx = bas[0].getBasisFunction(i+globalRow);
		NN1 = idx[0]; NN2 = idx[1]; LL = idx[2]; MM = idx[3]; KK = idx[4]; jj = idx[5];
		
		/*********************************************
		 * Evaluate the projections for Tpm(2,0,K,K')*
		 *********************************************/
		
		normalization = (!KK && !NN1 && !NN2 ? 0.5 :  SQRT_ONE_OVER_TWO);
		normalization *= SQRT_FIVE; // (2L+1)^(1/2), but here L = 2 only

		/*****************
		 * Evaluate sign *
		 *****************/
			
		multKK =  (!(KK % 2) ? 1.0 : -1.0);
		multNN1 = (NN1 >= 0  ? 1.0 : -1.0);
		multNN2 = (NN2 >= 0  ? 1.0 : -1.0);

		torsionalIntegral = dcomplex(internalIntegralReal[(NN1 + n1max) * ndim2 + (NN2 + n2max)], internalIntegralImag[(NN1 + n1max) * ndim2 + (NN2 + n2max)]);
			
		for (unsigned int v = 0; v < s<<1; v+=2)
		{
				
			TK =  (int)((unsigned int)sd[v>>1] / 5 - 2);
			TK1 = (int)((unsigned int)sd[v>>1] % 5 - 2);

			/********************************** 
			 * Calculate and store projection *
			 **********************************/

			deltaFactorK =  (KK == TK  ? (KK == -TK  ? 1.0 + multKK : 1.0) : (KK == -TK  ? multKK : 0.0));
			deltaFactorK1 = (KK == TK1 ? (KK == -TK1 ? 1.0 + multKK : 1.0) : (KK == -TK1 ? multKK : 0.0));

			projections[v][i+globalRow] = (double)normalization*normpm[v]*(deltaFactorK + deltaFactorK1)* \
							( jj == 1 ? torsionalIntegral.real() : torsionalIntegral.imag());

			if (fabs(projections[v][i+globalRow]) <= ZERO)  projections[v][i+globalRow] = 0.0;

			deltaFactorK =  (KK == TK ?  (KK == -TK  ? 1.0 - multKK : 1.0) : (KK == -TK  ? -multKK : 0.0));
			deltaFactorK1 = (KK == TK1 ? (KK == -TK1 ? 1.0 - multKK : 1.0) : (KK == -TK1 ? -multKK : 0.0));

			projections[v+1][i+globalRow] = (double)normalization*normpm[v+1]*(deltaFactorK + deltaFactorK1)* \
							( jj == 1 ? torsionalIntegral.imag() : -torsionalIntegral.real());

			if (fabs(projections[v+1][i+globalRow]) <= ZERO)  projections[v+1][i+globalRow] = 0.0;
				
		} /* close cycle over spectral densities */
			
	} /* Close cycle over i */

	return;
}

// USEFUL ROUTINES FOR FBn MODELS

double stvec::fb1_peq(double x) // NB: this calculates exp[-U(t)/2]  !!!!
{
	int i;
	int n_coef_int_pot = coef_int_pot.size();
  	double y = (double)coef_int_pot[0].real();
	for(i = 1 ; i < n_coef_int_pot; i++)
		y += 2.0 * (double)coef_int_pot[i].real() * cos((double)i*x) + 2.0 * (double)coef_int_pot[i].imag() * sin((double)i*x);
	y = exp(0.5*y);
	return (y);
}

double stvec::fb2_peq(const double *x, int nx) // NB: this calculates exp[-U(x1,x2)/2]  !!!!
{
	int n2UpLim;
	double angle, peq, mult, ca, sa;
	dcomplex y = dcomplex(0.0, 0.0);

	for (int n1 = -ufb2.n1Max; n1 <= 0; n1++)
	{
		n2UpLim = n1 < 0 ? ufb2.n2Max : 0;
		for (int n2 = -ufb2.n2Max; n2 <= n2UpLim; n2++)
		{
			mult = n1 == 0 && n2 == 0 ? 0.5 : 1.0;
			angle = (double)n1 * x[0] + (double)n2 * x[1];
			ca = cos(angle); sa = sin(angle);
			y += mult * ( ufb2.c[( n1 + ufb2.n1Max) * ufb2.dim2 + ( n2 + ufb2.n2Max)] * dcomplex(ca,-sa) + \
				      ufb2.c[(-n1 + ufb2.n1Max) * ufb2.dim2 + (-n2 + ufb2.n2Max)] * dcomplex(ca, sa) );
		}
	}
	peq = exp(0.5*y.real());
	return (peq);
}


double stvec::peq_ratio(double t1, double t2, double t1p, double t2p)  // Calculates exp[-(U(x')-U(x))/2]
{
	int n2UpLim;
	double x, xp, mtw, mult, ca, sa;
	dcomplex y = dcomplex(0.0, 0.0);

	for (int n1 = -ufb2.n1Max; n1 <= 0; n1++)
	{
		n2UpLim = n1 < 0 ? ufb2.n2Max : 0;
		for (int n2 = -ufb2.n2Max; n2 <= n2UpLim; n2++)
		{
			mult = n1 == 0 && n2 == 0 ? 0.5 : 1.0;
			x  = (double)n1 * t1  + (double)n2 * t2;
			xp = (double)n1 * t1p + (double)n2 * t2p;
			ca  = cos(xp) - cos(x);
			sa  = sin(xp) - sin(x);
			y += mult * ( ufb2.c[( n1 + ufb2.n1Max) * ufb2.dim2 + ( n2 + ufb2.n2Max)] * dcomplex(ca,-sa) + \
				      ufb2.c[(-n1 + ufb2.n1Max) * ufb2.dim2 + (-n2 + ufb2.n2Max)] * dcomplex(ca, sa) );
		}
	}
	mtw = exp(0.5 * y.real());
	return (mtw);
}

double stvec::internal_integrand_real(double x)
{
	double ff;
	ff = cos((double)nMult * x);
	ff *= fb1_peq(x);
	return (ff);
}

double stvec::internal_integrand_imag(double x)
{
	double ff;
	ff = sin((double)nMult * x);
	ff *= fb1_peq(x);
	return (ff);
}


int stvec::fb2_integrand_real(unsigned int ndim, const double *x, void *params, unsigned int fdim, double *ff)
{
	double *mm = (double *)params;
	ff[0] = cos(mm[0] * x[0] + mm[1] * x[1]);
	ff[0] *= fb2_peq(x,2);
	return 0;
}

double stvec::fb2_integrand_real_j(double *x, int *nX)
{
	double ff = cos(dn1 * x[0] + dn2 * x[1]);
	ff *= fb2_peq(x,2);
	return ff;
}

int stvec::fb2_integrand_imag(unsigned int ndim, const double *x, void *params, unsigned int fdim, double *ff)
{
	double *mm = (double *)params;
	ff[0] = sin(mm[0] * x[0] + mm[1] * x[1]);
	ff[0] *= fb2_peq(x,2);
	return 0;
}

double stvec::fb2_integrand_imag_j(double *x, int *nX)
{
	double ff = sin(dn1 * x[0] + dn2 * x[1]);
	ff *= fb2_peq(x,2);
	return ff;
}

/********************
 * SRLS projections *
 ********************/

void stvec::projectOnBasis_SRLS(basis *bas, int nr, int gr)
{
	nrows = nr;
	globalRow = gr;
	ivector idx;
	nbf = bas->getNumberOfBasisFunctions();
	int L1,M1,K1,L2,M2,K2,jj,Mm,Mp;
	bool L1check, M1check, jjcheck, Mmcheck;
	double normalization, multK1, multKK, multKK1;
	dvector betaIntegral(5,0.0);
	
	for (int i = (sd.size()<<1)-2; i >=0; i-=2)
	{		
		projections[i] = VECTOR_double(nbf,0.0);
		projections[i+1] = VECTOR_double(nbf,0.0);
	}
	
	/*************************
	 * Determine constraints *
	 *************************/
	
	ivector bounds(2);
	bounds = bas->getL1bounds();
	if (bounds[0] == bounds[1]) {L1constrain = bounds[0]; isL1constrained = true;}
	bounds = bas->getM1bounds();
	if (bounds[0] == bounds[1]) {M1constrain = bounds[0]; isM1constrained = true;}
	bounds = bas->getjjbounds();
	if (bounds[0] == bounds[1]) {jjconstrain = bounds[0]; isjjconstrained = true;}
	bounds = bas->getMmbounds();
	if (bounds[0] == bounds[1]) {Mmconstrain = bounds[0]; isMmconstrained = true;}

	unsigned int s = sd.size();
	
	/**********************************
	 * Precalculate useful quantities *
	 **********************************/

	dvector sqrtL2;
	for (int i = 0; i <= phy->getProbeLMax(); i++) sqrtL2.push_back(sqrt((double)(2*i+1)));

	int integrandFunction;
	bool is_c44_zero = (fabs(c44) <= ZERO);
	bool are_c22_and_c42_zero = (fabs(c22) <= ZERO && fabs(c42) <= ZERO);
	bool are_c20_and_c40_zero = (fabs(c20) <= ZERO && fabs(c40) <= ZERO);

	switch (is_c44_zero)
	{
		case true:
			switch (are_c22_and_c42_zero)
			{
				case true:
					switch (are_c20_and_c40_zero)
					{
						case true:
							/* No coupling potential - the integral is analytic */
							integrandFunction = 0;
							break;
						default:
							/* Integral with [c22=0 & c42=0 & c44=0] - integral in gamma is a constant */
							integrandFunction = 1;
							break;
					}
					break;
				default:
					/* Integral with c44=0 - integral in gamma is a modified Bessel function of first kind */
					integrandFunction = 2;
					break;
			}
			break;
		default:
			switch (are_c22_and_c42_zero)
			{
				case true:
					/* Integral with [c22=0 & c42=0] - integral in gamma is a modified Bessel function of first kind */
					integrandFunction = 3;
					break;
				default:
					/* Integral with [(c22!=0 | c42!=0) & c44!=0] - integral in gamma is not analytic */
					integrandFunction = 4;
					break;
			}
			break;		
	}

	/*****************************************
	 * Calculate the projection on the basis *
	 *****************************************/
	
	for (int i = 0; i < nrows; i++){
		
		/**********************
		 * Load basis indexes *
		 **********************/
		
		idx = bas->getBasisFunction(i+globalRow);
		L1 = idx[0]; M1 = idx[1]; K1 = idx[2];
		L2 = idx[3]; M2 = idx[4]; K2 = idx[5];
		jj = idx[6]; Mm = idx[7]; Mp = idx[8];
		
		/************************************************
		 * Check for constraints in the quantum numbers *
		 ************************************************/
		
		L1check = L1 == 2;  ////((isL1constrained && L1 == L1constrain) || !isL1constrained);
		M1check = M1 == 0;  ////((isM1constrained && M1 == M1constrain) || !isM1constrained);
		jjcheck = ((isjjconstrained && jj == jjconstrain) || !isjjconstrained);
		Mmcheck = K1==M2;
		
		/*********************************************
		 * Evaluate the projections for Tpm(2,0,K,K')*
		 *********************************************/
		
		if (L1check && M1check && Mmcheck){
			
			normalization = (!K1 && !K2 && !M2 ? 0.5 :  ONE_OVER_SQRT_TWO);
			normalization *= sqrtL2[L2];

			multK1 =  (!(K1%2)  ? 1.0 : -1.0);
			
			betaIntegral[0] = (!(abs(-2-K2)%2) ? integrate(L1,-K1, 2,L2,K1,K2,integrandFunction) : 0.0);
			betaIntegral[1] = (!(abs(-1-K2)%2) ? integrate(L1,-K1, 1,L2,K1,K2,integrandFunction) : 0.0);
			betaIntegral[2] = (!(abs( 0-K2)%2) ? integrate(L1,-K1, 0,L2,K1,K2,integrandFunction) : 0.0);
			betaIntegral[3] = (!(abs( 1-K2)%2) ? integrate(L1,-K1,-1,L2,K1,K2,integrandFunction) : 0.0);
			betaIntegral[4] = (!(abs( 2-K2)%2) ? integrate(L1,-K1,-2,L2,K1,K2,integrandFunction) : 0.0);

			for (unsigned int v = 0; v < s<<1; v+=2)
			{
				
				KK =  (int)((unsigned int)sd[v>>1] / 5 - 2);
				KK1 = (int)((unsigned int)sd[v>>1] % 5 - 2);

				/*****************
				 * Evaluate sign *
				 *****************/

				multKK =  (!(KK%2)  ? 1.0 : -1.0);
				multKK1 = (!(KK1%2) ? 1.0 : -1.0);
			
				/********************************** 
				 * Calculate and store projection *
				 **********************************/

				projections[v][i+globalRow] = ( jj == 1 ? (double)(normalization*normpm[v]*multK1*(multKK*betaIntegral[KK+2]+betaIntegral[-KK+2]+multKK1*betaIntegral[KK1+2]+betaIntegral[-KK1+2])) : 0.0);
				if (fabs(projections[v][i+globalRow]) <= ZERO)  projections[v][i+globalRow] = 0.0;
				projections[v+1][i+globalRow] = ( jj == -1 ? (double)(normalization*normpm[v+1]*multK1*(multKK*betaIntegral[KK+2]-betaIntegral[-KK+2]+multKK1*betaIntegral[KK1+2]-betaIntegral[-KK1+2])) : 0.0);
				if (fabs(projections[v+1][i+globalRow]) <= ZERO)  projections[v+1][i+globalRow] = 0.0;
				
			} /* close cycle over spectral densities */
			
	 	} /* Close check of constraints */
	
	} /* Close cycle over i */
	
#ifdef WRITE_ALL
#ifdef WRITE_STVEC
	std::cout << "TASK " << rank << ": ... starting vector projections calculated successfully." << std::endl;
#endif
#endif
	
	return;
}

/*********************************
 *  Calculation of beta integral *
 *********************************/

extern "C" {
 #include "cquadpak.h"
}

inline double stvec::integrate(int L1, int mK1, int mK, int L2, int K1, int K2, int func)
{
	
	int nv = 0, er = 0;
	double ae = 0.0;
	
	double integral = 0.0;
	
	/* 
	 * The code distinguish among different cases
	 * to calculate fastly the integral in special cases
	 */ 
	this->L1 = L1;
	this->mK1 = mK1;
	this->mK = mK;
	this->L2 = L2;
	this->K1 = K1;
	this->K2 = K2;
	
	switch (func)
	{
		case 0:
		{
			/* No coupling potential - the integral is analytic */
			integral = ((L1==L2 && mK+K2==0) ? (!((K1-mK)%2) ? 1.0/((double)(2*L1+1)) : -1.0/((double)(2*L1+1))) : 0.0);
			break;
		}
		case 1:
		{
			/* Integral with [c22=0 && c42=0 && c44=0] - integral in gamma is a constant */
			if (mK+K2==0)
				integral = dqag(integrand_no_bessel,0.0,M_PI,1.0e-10,1.0e-10,6,&ae,&nv,&er);
			break;
		}
		case 2:
		{
			/* Integral with c44=0 - integral in gamma is a modified Bessel function of first kind */
			if (!(abs(K2+mK)%2))
				integral = dqag(integrand_one_bessel_u2,0.0,M_PI,1.0-10,1.0e-10,6,&ae,&nv,&er);
			break;
		}
		case 3:
		{
			/* Integral with [c22=0 && c42=0] - integral in gamma is a modified Bessel function of first kind */
			if (!(abs(K2+mK)%4))
				integral = dqag(integrand_one_bessel_u4,0.0,M_PI,1.0e-10,1e-10,6,&ae,&nv,&er);
			break;
		}
		case 4:
		{
			/* Integral with [(c22!=0 | c42!=0) && c44!=0] - integral in gamma is not analytic */
			integral = dqag(integrand_inf_bessel,0.0,M_PI,1.0e-10,1.0e-10,6,&ae,&nv,&er);
			break;
		}
		case 5:
		{
			/* Integral for calculation of order parameters */
			integral = 2.0 * M_PI * dqag(integrandOrderParams,0.0,M_PI,1.0e-10,1.0e-10,6,&ae,&nv,&er);
			break;
		}
	}

	return integral;
	
}

/******************************************************
 *  Integrand function - gamma integral is a constant *
 ******************************************************/

double stvec::integrand_no_bessel(double x)
{
	double u0 = (double)(-0.5*(c20*w.getReducedWignerMatrix(2,0,0,x)+c40*w.getReducedWignerMatrix(4,0,0,x)));
	double f = sin(x);
	f  *= w.getReducedWignerMatrix(L1,mK1,mK,x);
	f *= w.getReducedWignerMatrix(L2,K1,K2,x);
	f *= exp(u0);
	
	return f;
}

/************************************************************************************
 *  Integrand function - gamma integral is a modified bessel function of first kind *
 ************************************************************************************/

double stvec::integrand_one_bessel_u2(double x)
{
	int order = abs(K2 + mK) / 2;
	double f = 0.0;
	double u0 = -0.5*(double)(c20*w.getReducedWignerMatrix(2,0,0,x)+c40*w.getReducedWignerMatrix(4,0,0,x));
	double u2 = -(double)(c22*w.getReducedWignerMatrix(2,0,2,x)+c42*w.getReducedWignerMatrix(4,0,2,x));

	f = sin(x);
	f *= (double)w.getReducedWignerMatrix(L1,mK1,mK,x);
	f *= (double)w.getReducedWignerMatrix(L2,K1,K2,x);
	f *= exp(u0);
	f *= bessi(order, fabs(u2));
	f *= (u2 < 0.0 ? (order%2 == 0 ? 1.0 : -1.0) : 1.0);
		
	return f;
}

/************************************************************************************
 *  Integrand function - gamma integral is a modified bessel function of first kind *
 ************************************************************************************/

double stvec::integrand_one_bessel_u4(double x)
{	
	int order = abs(K2 + mK) / 4; 
	double f = 0.0;
	double u0 = -0.5*(double)(c20*w.getReducedWignerMatrix(2,0,0,x)+c40*w.getReducedWignerMatrix(4,0,0,x));
	double u4 = -c44*(double)w.getReducedWignerMatrix(4,0,4,x);
		
	
	f = sin(x);
	f *= (double)w.getReducedWignerMatrix(L1,mK1,mK,x);
	f *= (double)w.getReducedWignerMatrix(L2,K1,K2,x);
	f *= exp(u0);
	f *= bessi(order, fabs(u4));
	f *= (u4 < 0.0 ? (order%2 == 0 ? 1.0 : -1.0) : 1.0);
		
	return f;	
}
/********************************************************************************************
 *  Integrand function - gamma integral is a sum of modified bessel functions of first kind *
 ********************************************************************************************/

double stvec::integrand_inf_bessel(double x)
{
	long double f = 0.0;
	std::cout << std::endl << std::endl << "ERROR: [(c22!=0 | c42 !=0) & c44 !=0] logical expression evaluated to TRUE. This combination of c22, c42 and c44 parameters is still not implemented in C++OPPS" << std::endl << std::endl;
	exit(1);
	return f;
}

/****************************************************************************************
 * mpiFlag = TRUE  scatter all local projections to global projections of all processes *
 * mpiFlag = FALSE make global projections point to local projections                   *
 ****************************************************************************************/

void stvec::scatterProjections(bool mpiFlag, int nproc)
{
	switch (mpiFlag)
	{
		case (true):
		{
#ifdef _MPI_
			int pos;
			int procRows[nproc];
			int sProcRows[nproc];
			for (int i=0; i<nproc; i++) sProcRows[i] = nrows;
			MPI_Alltoall(sProcRows,1,MPI_INT,procRows,1,MPI_INT,MPI_COMM_WORLD);
			MPI_Request srequest, rrequest;
			MPI_Status status;
			double* rb;
			double* sb;
			sb = (double *)calloc(nrows, sizeof(double));
			for (unsigned int i = 0; i < (sd.size()<<1); i++)
			{
				pos = 0;
				for (int j = 0; j < nproc; j++)
				{
						/* send local buffer to j */
						for (int u = 0; u < nrows; u++) sb[u] = projections[i].p_[globalRow + u];
						MPI_Isend(sb,nrows,MPI_DOUBLE,j,99,MPI_COMM_WORLD,&srequest);
						/* recive to local buffer from j */
						rb = &projections[i].p_[pos];
						MPI_Irecv(rb,procRows[j],MPI_DOUBLE,j,99,MPI_COMM_WORLD,&rrequest);
						MPI_Wait(&srequest,&status);
						MPI_Wait(&rrequest,&status);
						pos += procRows[j];
				}
			}
#endif
			break;
		}
		case (false):
		{
			break;
		}
	}
	
       	
	/********************
	 * Normalize vector *
	 ********************/
	
	double vecnorm;
	for (unsigned int i = 0; i < (sd.size()<<1); i++)
	{
		vecnorm = sqrt(dot(projections[i],projections[i]));
		if (vecnorm > ZERO)
		{
			vecnorm = 1.0/vecnorm;
			for (int j = 0; j < nbf; j++) projections[i][j] *= vecnorm;
		}
		else
			npmkk1[i] = 0.0;
		/* DEBUGGING */
		/*
		std::cout << "NORM: " << npmkk1[i] << std::endl;
		*/
	}
	return;	
}




/******************************
 * Calculate Order Parameters *
 ******************************/
double stvec::integrandOrderParams(double x)
{
	//using  L1 to pass j, L2 to pass k. Note - k is the order.
	int order = L2 / 2;

	double u0beta = - (double)(c20 * w.getReducedWignerMatrix(2,0,0,x)+ 
			c40 * w.getReducedWignerMatrix(4,0,0,x));
	double u2beta = - 2.0*(double)(c22 * w.getReducedWignerMatrix(2,0,2,x)+
			c42 * w.getReducedWignerMatrix(4,0,2,x));
	double u4beta = - 2.0*(double)(c44 * w.getReducedWignerMatrix(4,0,4,x));

	// Calc returned function
	double f  = sin(x);
	f *= (double)w.getReducedWignerMatrix(L1,0,L2,x);
	f *= exp(u0beta);
	f *= bessi(order, fabs(u2beta));
	f *= (u2beta < 0.0 ? (order%2 == 0 ? 1.0 : -1.0) : 1.0);
	f *= bessi0(fabs(u4beta));
	return f;
}	


dvector stvec::calculateOrdersParams()
  {
	// vector for result reporting. [S20, S22, Sxx, Syy, Szz]
	dvector ordersParams;
	
	// ORDPARM points at the right integrand function in the integrate method
	int ORDPARM = 5;
	
	// Calculate S2/0
	double S20 = integrate(2,0,0,0,0,0,ORDPARM) / integrate(0,0,0,0,0,0,ORDPARM);
	ordersParams.push_back(S20);
	
	// Calculate S2/2
	double S22 = 2.0*integrate(2,0,0,2,0,0,ORDPARM) / integrate(0,0,0,0,0,0,ORDPARM);
	ordersParams.push_back(S22);
	 
	// Calculate Sxx
	double Sxx = 0.5 * (sqrt(1.5) * S22 - S20);
	ordersParams.push_back(Sxx);
	
	// Calculate Syy
	double Syy = -0.5 * (sqrt(1.5) * S22 + S20);
	ordersParams.push_back(Syy);
	
	// Szz equals S20
	ordersParams.push_back(S20);
	
	return ordersParams;

  }


/***************************
 *  Return all projections *
 ***************************/

VECTOR_double* stvec::getProjections(void)
{
	return projections;
}

/**************************************************************
 * Return projection of observable with index i < 2*sd.size() *
 **************************************************************/

VECTOR_double stvec::getProjections(int i)
{
		return projections[i];
}

VECTOR_double stvec::getProjections(unsigned int i)
{
		return projections[i];
}

/********************************************
 *  Return all small n constants n +/- K,K' *
 ********************************************/

dvector stvec::getNpmKK1(void)
{
	return npmkk1;
}

/***********************************************************
 *  Return the small n constant with index i < 2*sd.size() *
 ***********************************************************/

double stvec::getNpmKK1(int i)
{
		return npmkk1[i];
}

double stvec::getNpmKK1(unsigned int i)
{
		return npmkk1[i];
}

/***************************
 *  Writes the projections *
 ***************************/

std::string stvec::toString(basis *bas)
{
	ivector idx;
	int N,L,M,K,j,N1,N2;
	int L1,M1,K1,L2,M2,K2,Mm,Mp;
	std::ostringstream ostr;
	
	/****************************************
	 * Non zero projection of J(2,0,KK,KK1) *
	 ****************************************/

	for (unsigned s = 0; s < 2*sd.size(); s+=2)
	{
		
		unsigned int y = s/2;
		
		KK  = sd[y] / 5 - 2;
		KK1 = sd[y] % 5 - 2;
		
		ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
		ostr << "TASK " << rank << ":           Porjections for T+ L = 2, M = 0, K = " << KK << ", K' = " << KK1 << std::endl;
		ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
		if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
		{
			ostr << "TASK " << rank << ": norm = " << npmkk1[s] << std::endl;
			ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << ": \tL1\tM1\tK1\tL2\tM2\tK2\tprj" << std::endl;
		}
		else if (!(dynModel.compare("fb1")))
		{
			ostr << "TASK " << rank << ": norm = " << npmkk1[s] << std::endl;
			ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << ": \tN\tL\tM\tK\tj\tprj" << std::endl;
		}
		else if (!(dynModel.compare("fb2")))
		{
			ostr << "TASK " << rank << ": norm = " << npmkk1[s] << std::endl;
			ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << ": \tN1\tN2\tL\tM\tK\tj\tprj" << std::endl;
		}
		for (int i = 0; i < bas[0].getNumberOfBasisFunctions(); i++)
		{
			if (fabs(projections[s][i]) > ZERO)
			{
				idx = bas[0].getBasisFunction(i);
				if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
				{
					L1 = idx[0]; M1 = idx[1]; K1 = idx[2];
					L2 = idx[3]; M2 = idx[4]; K2 = idx[5];
					Mm = idx[6]; Mp = idx[7];
					ostr << "TASK " << rank << ": " << i+1 << "\t" << L1 << "\t" << M1 << "\t" << K1 << "\t" << L2
						 << "\t" << M2 << "\t" << K2 << "\t" << setiosflags(ios::fixed) << setprecision(4) 
						<< showpos << projections[s][i] << noshowpos << std::endl;
				}
				else if (!(dynModel.compare("fb1")))
				{
					N = idx[0]; L = idx[1]; M = idx[2]; K = idx[3]; j = idx[4];
					ostr << "TASK " << rank << ": " << i+1 << "\t" << N << "\t" << L << "\t" << M << "\t" << K
						 << "\t" << j << "\t" << setiosflags(ios::fixed) << setprecision(4) 
						<< showpos << projections[s][i] << noshowpos << std::endl;
				}
				else if (!(dynModel.compare("fb2")))
				{
					N1 = idx[0]; N2= idx[1]; L = idx[2]; M = idx[3]; K = idx[4]; j = idx[5];
					ostr << "TASK " << rank << ": " << i+1 << "\t" << N1 << "\t" << N2  << "\t" << L << "\t" << M << "\t" << K
						 << "\t" << j << "\t" << setiosflags(ios::fixed) << setprecision(4) 
						<< showpos << projections[s][i] << noshowpos << std::endl;
				}
			}
		}
		ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
		ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
		ostr << "TASK " << rank << ":           Porjections for T- L = 2, M = 0, K = " << KK << ", K' = " << KK1 << std::endl;
		ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
		if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
		{
			ostr << "TASK " << rank << ": norm = " << npmkk1[s] << std::endl;
			ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << ": \tL1\tM1\tK1\tL2\tM2\tK2\tprj" << std::endl;
		}
		else if (!(dynModel.compare("fb1")))
		{
			ostr << "TASK " << rank << ": norm = " << npmkk1[s] << std::endl;
			ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << ": \tN\tL\tM\tK\tj\tprj" << std::endl;
		}
		else if (!(dynModel.compare("fb1")))
		{
			ostr << "TASK " << rank << ": norm = " << npmkk1[s] << std::endl;
			ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;
			ostr << "TASK " << rank << ": \tN1\tN2\tL\tM\tK\tj\tprj" << std::endl;
		}
		for (int i = 0; i < bas[0].getNumberOfBasisFunctions(); i++)
		{
			if (fabs(projections[s+1][i]) > ZERO)
			{
				idx = bas[0].getBasisFunction(i);
				if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
				{
					L1 = idx[0]; M1 = idx[1]; K1 = idx[2];
					L2 = idx[3]; M2 = idx[4]; K2 = idx[5];
					Mm = idx[6]; Mp = idx[7];
					ostr << "TASK " << rank << ": " << i+1 << "\t" << L1 << "\t" << M1 << "\t" << K1 << "\t" << L2 
						<< "\t" << M2 << "\t" << K2 << "\t" << setiosflags(ios::fixed) << setprecision(4) 
						<< showpos << projections[s+1][i] << noshowpos << std::endl;
				}
				else if (!(dynModel.compare("fb1")))
				{
					N = idx[0]; L = idx[1]; M = idx[2]; K = idx[3]; j = idx[4];
					ostr << "TASK " << rank << ": " << i+1 << "\t" << N << "\t" << L << "\t" << M << "\t" << K
						 << "\t" << j << "\t" << setiosflags(ios::fixed) << setprecision(4) 
						<< showpos << projections[s+1][i] << noshowpos << std::endl;
				}
				else if (!(dynModel.compare("fb2")))
				{
					N1 = idx[0]; N2= idx[1]; L = idx[2]; M = idx[3]; K = idx[4]; j = idx[5];
					ostr << "TASK " << rank << ": " << i+1 << "\t" << N1 << "\t" << N2  << "\t" << L << "\t" << M << "\t" << K
						 << "\t" << j << "\t" << setiosflags(ios::fixed) << setprecision(4) 
						<< showpos << projections[s+1][i] << noshowpos << std::endl;
				}
			}
		}
		ostr << "TASK " << rank << ": ----------------------------------------------------------------" << std::endl;

	}
	return ostr.str();
}

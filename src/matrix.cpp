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
 Name        : matrix.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class to store the matrix associated to the diffusive operator
 ============================================================================
 */

#include "matrix.h"

/***************
 * Constructor *
 ***************/

matrix::matrix()
{
	isInit =  false;
}

/************** 
 * Destructor *
 **************/

matrix::~matrix()
{
	if (A.val_.p_ != NULL) A.val_.ref_ = 1;
	if (A.colind_.p_ != NULL) A.colind_.ref_ = 1;
#ifdef WRITE_DESTROY_MESSAGE
	std::cout << "TASK " << rank << ": Cleared matrix object" << std::endl;
#endif
}

/**********************************************************
 * Initiators:                                            *
 * basis *b = I - pointer to a basis object;              *
 * int r    = I - rank of task;                           *
 * int nr   = I - number of rows to be handled by object; *
 * int g    = I - starting row in the full matrix.        *
 **********************************************************/

void matrix::init(basis *b, physics *p, s3j *s)
{
	isInit = true;
	bas = b;
	phy = p;
	trj = s;
	rank = 0;
	nrows = ncols = bas->getNumberOfBasisFunctions();
	global_row = 0;
	
	/*************************
	 * Initialize the matrix *
	 *************************/
	
	A = CompRow_Mat_double();
	A.dim_[0] = A.dim_[1] = nrows;
	A.nz_ = 0;
	A.rowptr_.newsize(A.dim_[0]+1);
	A.rowptr_(0) = 0;
	A.rowptr_(A.dim_[0]) = 0;
	A.val_.dim_ = 0;
	A.val_.p_ = NULL;
	A.colind_.dim_ = 0;
	A.colind_.p_ = NULL;
	
	diag = VECTOR_double(nrows,0.0);
	
	/****************************/
	/* Setup matrix calculation */
	/****************************/

	setupCalculation();

	/***************************/
	/* Initialization complete */
	/***************************/
	
	std::cout << "TASK " << rank << ": Created matrix object" << std::endl;
	
	return;
}
void matrix::init(basis *b, physics *p, s3j *s, int r, int nr, int gr, int nc, int gc)
{
	isInit = true;
	bas = b;
	phy = p;
	trj = s;
	rank = r;
	nrows = nr;
	global_row = gr;
	ncols = nc;
	global_col = gc;
	
	/*************************
	 * Initialize the matrix *
	 *************************/
	
	A = CompRow_Mat_double();
	A.dim_[0] = nrows; A.dim_[1] = ncols;
	A.nz_ = 0;
	A.rowptr_.newsize(A.dim_[0]+1);
	A.rowptr_(0) = 0;
	A.rowptr_(A.dim_[0]) = 0;
	A.val_.dim_ = 0;
	A.val_.p_ = NULL;
	A.colind_.dim_ = 0;
	A.colind_.p_ = NULL;
	
	diag = VECTOR_double(nrows,0.0);
	
	/****************************/
	/* Setup matrix calculation */
	/****************************/

	setupCalculation();

	/***************************/
	/* Initialization complete */
	/***************************/
	
	std::cout << "TASK " << rank << ": Created matrix object" << std::endl;

	return;
}

void matrix::setupCalculation()
{
	/*****************************/
	/* Obtain model for dynamics */
	/*****************************/

	dynModel = phy->getDynamicsModel();

	/***************************************
	 * Precalculate some useful quantities *
	 ***************************************/
	
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
		/* sqrt[(2L-3)!/24*(2L+2)!] */
		
		L2Mult.clear();
		L2Mult.push_back(0.0);
		for (int i = 1; i <= phy->getProbeLMax(); i++)
			L2Mult.push_back(sqrt(exp(trj->getLnFac(2*i+3)-trj->getLnFac(2*i-2))*ONE_OVER_TWENTYFOUR));
		
		/* c+ and c- coefficients */
		
		int Lmax = max(2,phy->getProbeLMax());
		cp = dvector((Lmax+1)*(Lmax+2)+1,0.0);
		cm = dvector((Lmax+1)*(Lmax+2)+1,0.0);
		for (int i = 0; i <= Lmax; i++)
		{
			for (int j = -i; j <= i; j++)
			{
				cp[(i+1)*(i+1)+j] = cplm(i,j);
				cm[(i+1)*(i+1)+j] = cmlm(i,j);
			}
		}

		/* fa and fb coefficients */
		
		fa = new dvector[6];
		fb = new dvector[6];
		for (int i = 0; i < 6; i++) { fa[i] = dvector(9*9*9,0.0); fb[i] = dvector(9*9*9,0.0); }
		storeFaFb();
	}

	else if (!(dynModel.compare("fb1")) || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb")))
	{
		if (f11.size() > 0) f11.clear();
		storeFIIcoeff(1);
	}

	else if (!(dynModel.compare("fb2")))
	{
		// NB: coefficients entering matrix elements are calculated "on-the-fly"
		//if (f22.size() > 0) f22.clear();
		storeFIIcoeff(2);
	}

	/********************/
	/* Diffusion tensor */
	/********************/

	setDiften();

	/****************/
	/* Finish setup */
	/****************/

	return;
}

/************************
 * Set diffusion tensor *
 ************************/

void matrix::setDiften(void)
{
	/******************
	 * Model booleans *
	 ******************/

	bool useSRLS = !(dynModel.compare("srls")) || !(dynModel.compare("ts-srls"));
	bool useFB1 =  !(dynModel.compare("fb1"))  || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb"));
	bool useFB2 =  !(dynModel.compare("fb2"));

	/************************
	 * Obtain physical data *
	 ************************/

	if (useSRLS)
	{
		probeD[0] = (dcomplex)phy->getProbeDlm(0, 0);
		probeD[1] = (dcomplex)phy->getProbeDlm(2,-2);
		probeD[2] = (dcomplex)phy->getProbeDlm(2,-1);
		probeD[3] = (dcomplex)phy->getProbeDlm(2, 0);
		probeD[4] = (dcomplex)phy->getProbeDlm(2, 1);
		probeD[5] = (dcomplex)phy->getProbeDlm(2, 2);
	
		proteinD[0] = (dcomplex)phy->getProteinDlm(0, 0);
		proteinD[1] = (dcomplex)phy->getProteinDlm(2,-2);
		proteinD[2] = (dcomplex)phy->getProteinDlm(2,-1);
		proteinD[3] = (dcomplex)phy->getProteinDlm(2, 0);
		proteinD[4] = (dcomplex)phy->getProteinDlm(2, 1);
		proteinD[5] = (dcomplex)phy->getProteinDlm(2, 2);
	}

	else if (useFB1)
	{
		dvector fullD = phy->getFB1Diften();
		dxx = fullD[0]; dyy = fullD[1]; dzz = fullD[2];
		d11 = fullD[3];
		dx1 = fullD[4]; dy1 = fullD[5]; dz1 = fullD[6];
		dpr = dxx + dyy;
		dmr = dxx - dyy;
		dp1.real(0.5*dx1); dp1.imag(-0.5*dy1);
		dm1.real(0.5*dx1); dm1.imag( 0.5*dy1);
	}

	else if (useFB2)
	{
		dvector fullD = phy->getFB2Diften();
		dxx = fullD[0]; dyy = fullD[1];  dzz = fullD[2];
		d11 = fullD[3]; d22 = fullD[4];  d12 = fullD[5];
		dx1 = fullD[6]; dy1 = fullD[7];  dz1 = fullD[8];
		dx2 = fullD[9]; dy2 = fullD[10]; dz2 = fullD[11];

		dpr = dxx + dyy;
		dmr = dxx - dyy;

		dp1.real(0.5*dx1); dp1.imag(-0.5*dy1);
		dm1.real(0.5*dx1); dm1.imag( 0.5*dy1);
                                                        
		dp2.real(0.5*dx2); dp2.imag(-0.5*dy2);
		dm2.real(0.5*dx2); dm2.imag( 0.5*dy2);
	}

	return;
}

/*********************************************
 * Method for object update during a fitting *
 *********************************************/

void matrix::update(void)
{
	//A = CompRow_Mat_double();
	A.dim_[0] = nrows; A.dim_[1] = ncols;
	A.nz_ = 0;
	A.rowptr_.newsize(A.dim_[0]+1);
	A.rowptr_(0) = 0;
	A.rowptr_(A.dim_[0]) = 0;
	A.val_.dim_ = 0;
	A.val_.p_ = NULL;
	A.colind_.dim_ = 0;
	A.colind_.p_ = NULL;
	
	diag = VECTOR_double(nrows,0.0);
        
	setupCalculation();

	return;
}

/*****************************************************
 * Calculate and store fa and fb coefficients - SRLS *
 *****************************************************/

void matrix::storeFaFb(void)
{
	
	/***************************************
	 * Chek if object has been initialized *
	 ***************************************/
	
	if (!isInit)
	{
		std::cout << std::endl << "ERROR in matrix::storeFaFb(). Matrix object created but not initialized." << std::endl << std::endl;
		exit(0);
	}
	int mu, mup, j, nu, nup, pos;
	double mult,sign,c1,c2;
	double ca00, ca20, ca2p1, ca2m1, ca2p2, ca2m2;
	double cb00, cb20, cb2p1, cb2m1, cb2p2, cb2m2;
	ldvector coeff = phy->getCoeff();
	
	/*********************************************************** 
	 * Determine the limits of variability of indexes mu and j *
	 ***********************************************************/

	if (fabs(coeff[14]) <= ZERO)
	{
		if (fabs(coeff[12]) <= ZERO && fabs(coeff[5]) <= ZERO)
		{
			if (fabs(coeff[3]) <= ZERO && fabs(coeff[10]) <= ZERO)
			{
				mulim = -1;
				jlim = -1;
			}
			else
			{
				mulim = 0;
				jlim = (fabs(coeff[10]) <= ZERO ? 4 : 8);
			}
		}
		else
		{
			mulim = 2;
			jlim = ((fabs(coeff[10]) <= ZERO && fabs(coeff[12]) <= ZERO) ? 4 : 8);  
		}
	}
	else
	{
		mulim = 4;
		jlim = 8;
	}

	int nuidx, nupidx;

	for (mu = -mulim; mu <= mulim; mu++){
		for (mup = -mulim; mup <= mulim; mup++){
			for (j = 0; j <= jlim; j++){
				
				ca00 = ca20 = 0.0; ca2p1 = ca2m1 = 0.0; ca2p2 = ca2m2 = 0.0;
				cb00 = cb20 = 0.0; cb2p1 = cb2m1 = 0.0; cb2p2 = cb2m2 = 0.0;
				
				for (nu = abs(mu); nu <= 4; nu++){
					
					nuidx = (nu+1)*(nu+1);
					
					c1 = (double)coeff[(nu*(nu+1)>>1)+abs(mu)];
					
					if (j == nu && mup == 0)
					{
						ca00  -= c1*(cplm(nu,mu)*cmlm(nu,mu+1)+cmlm(nu,mu)*cplm(nu,mu-1)+2.0*(double)(mu*mu));
						ca20  += c1*(cplm(nu,mu)*cmlm(nu,mu+1)+cmlm(nu,mu)*cplm(nu,mu-1)-4.0*(double)(mu*mu)); 
						ca2p1 += c1*(double)(2*mu+1)*cplm(nu,mu);
						ca2m1 += c1*(double)(2*mu-1)*cmlm(nu,mu);
						ca2p2 += -2.0*c1*cplm(nu,mu)*cplm(nu,mu+1);
						ca2m2 += -2.0*c1*cmlm(nu,mu)*cmlm(nu,mu-1);
						
						cb00  += -2.0*c1*(double)(nu*(nu+1));
						cb2p1 += c1*cmlm(nu,0);
						cb2m1 += c1*cplm(nu,0);
						cb2p2 += -2.0*c1*cmlm(nu,0)*cmlm(nu,-1);
						cb2m2 += -2.0*c1*cplm(nu,0)*cplm(nu,1);
					}
					
					for (nup = max(abs(mup),abs(nu-j)); nup <= min(4,(nu+j)); nup++)
					{
						nupidx = (nup+1)*(nup+1);
						
						c2 = (double)coeff[(nup*(nup+1)>>1)+abs(mup)];
						sign = (!((mu+mup)%2) ? 1.0 : -1.0);
						
						mult = sign*c1*c2*(double)trj->getTrj(nu,nup,j,0,0,0)*(double)(2*j+1);
											
						ca00 += mult*((double)trj->getTrj(nu,nup,j,mu+1,mup-1,-(mu+mup))*cplm(nu,mu)*cmlm(nup,mup)+
								(double)trj->getTrj(nu,nup,j,mu,mup,-(mu+mup))*(double)(mu*mup));
						ca20 += mult*(-(double)trj->getTrj(nu,nup,j,mu+1,mup-1,-(mu+mup))*cplm(nu,mu)*cmlm(nup,mup)+
								(double)trj->getTrj(nu,nup,j,mu,mup,-(mu+mup))*2.0*(double)(mu*mup));
						ca2p1 += mult*cplm(nup,mup)*(double)trj->getTrj(nu,nup,j,mu,mup+1,-(mu+mup+1))*(double)mu; 
						ca2m1 += mult*cmlm(nup,mup)*(double)trj->getTrj(nu,nup,j,mu,mup-1,-(mu+mup-1))*(double)mu;
						ca2p2 += mult*cplm(nu,mu)*cplm(nup,mup)*(double)trj->getTrj(nu,nup,j,mu+1,mup+1,-(mu+mup+2));
						ca2m2 += mult*cmlm(nu,mu)*cmlm(nup,mup)*(double)trj->getTrj(nu,nup,j,mu-1,mup-1,-(mu+mup-2));
						
						mult = sign*c1*c2*(double)(2*j+1);
						
						cb00  += mult*(double)trj->getTrj(nu,nup,j,-1, 1, 0)*(double)trj->getTrj(nu,nup,j,mu,mup,-(mu+mup))*cplm(nu,0)*cmlm(nup,0);
						cb2p2 += mult*(double)trj->getTrj(nu,nup,j,-1,-1, 2)*cmlm(nu,0)*cmlm(nup,0)*(double)trj->getTrj(nu,nup,j,mu,mup,-(mu+mup));
						cb2m2 += mult*(double)trj->getTrj(nu,nup,j, 1, 1,-2)*cplm(nu,0)*cplm(nup,0)*(double)trj->getTrj(nu,nup,j,mu,mup,-(mu+mup));
					}
				}
				
				ca00  *=  0.25*ONE_OVER_SQRT_THREE;
				ca20  *= -0.25*ONE_OVER_SQRT_SIX;
				ca2p1 *= -0.25;
				ca2m1 *=  0.25;
				ca2p2 *= -0.125;
				ca2m2 *= -0.125;
					
				cb00  *=  0.25*ONE_OVER_SQRT_THREE;
				cb20   =  ONE_OVER_SQRT_TWO*cb00;
				cb2p1 *=  0.25;
				cb2m1 *=  0.25;
				cb2p2 *= -0.125;
				cb2m2 *= -0.125;

				pos = (mu+4)*9*9+(mup+4)*9+j;
				
				fa[0][pos] = ca00;
				fa[1][pos] = ca2m2;
				fa[2][pos] = ca2m1;
				fa[3][pos] = ca20;
				fa[4][pos] = ca2p1;				
				fa[5][pos] = ca2p2;
				
				fb[0][pos] = cb00;
				fb[1][pos] = cb2m2;
				fb[2][pos] = cb2m1;
				fb[3][pos] = cb20;
				fb[4][pos] = cb2p1;				
				fb[5][pos] = cb2p2;
				
			}
		}
	}

	return;
}

/**********************************************
 * Calculate and store fii coefficients - FBn *
 **********************************************/

void matrix::storeFIIcoeff(int nTors)
{
	int nu1,nu2,nu1p,nu2p;
	dcomplex c1, c2;
	
	dcomplex cval;
	dcomplex czero(0.0,0.0);
	
	switch (nTors)
	{

		/******************************************
		 * FB1: U(t1) = Sum_n1 c_n1 exp(-i n1 t1) *
		 ******************************************/
		case 1:
		{
			dcvector coeff = phy->getFB1Coeff();
	
			int ncipmo = 2 * coeff.size() - 1;
			int ncippo = coeff.size() - 1;
			int dim = ncipmo;
	
			f11 = dcvector(dim*dim,czero);
	
			for (nu1 = -ncippo; nu1 <= ncippo; nu1++){
				c1 = (nu1 >= 0 ? coeff[abs(nu1)] : conj(coeff[abs(nu1)]));
				for (nu2 = -ncippo; nu2 <= ncippo; nu2++){
					c2 = (nu2 >= 0 ? coeff[abs(nu2)] : conj(coeff[abs(nu2)]));
					cval = (double)(nu1*nu2) * c1 * c2;
					cval = cval + (!nu2 ? 2.0*(double)(nu1*nu1) * c1 : czero);
					f11[(nu1+ncippo)*dim+(nu2+ncippo)] = cval;
				}
			}
			break;
		}
	
		/*****************************************************************
		 * FB2: U(t1,t2) = Sum_n1 Sum_n2 c_n1,n2 exp[-i (n1 t1 + n2 t2)] *
		 *****************************************************************/
		case 2:
		{
			ufb2 = phy->getFB2Coeff();
			break;
		}
	}
	
	return;
}


/*********************************
 * Calculate the matrix elements *
 *********************************/

void matrix::buildMatrix(void)
{
	
	/***************************************
	 * Chek if object has been initialized *
	 ***************************************/
	
	if (!isInit)
	{
		std::cout << std::endl << "ERROR in matrix::buildMatrix(). Matrix object created but not initialized." << std::endl << std::endl;
		exit(0);
	}

	/******************
	 * Model booleans *
	 ******************/
	
	bool useSRLS = !(dynModel.compare("srls")) || !(dynModel.compare("ts-srls"));
	bool useFB1 =  !(dynModel.compare("fb1"))  || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb"));
	bool useFB2 =  !(dynModel.compare("fb2"));

	/**************************************************************************
	 * Every task will build its slice of matrix following this simple rules: *
	 * - if i+j is even calculate Aij                                         *
	 * - if i+j is odd, calculate Aji                                         *
	 * which better permit to distribute the computational load.              *
	 **************************************************************************/
	
	int deltaj, jj;
	double norm,element,mult;
	ivector colIndexes;
	dcomplex elementp, elementm;
	dvector rowElements;
	oldi = -1;

	for (int i = global_row; i < global_row+nrows; i++){

		colIndexes.clear();
		rowElements.clear();
		
		/********************
		 * Diagonal element *
		 ********************/
		
		if (useSRLS)
		{
			elementp = buildElement_SRLS(i,i,1,&norm,&mult,&deltaj,&jj);
			elementm = buildElement_SRLS(i,i,-1,&norm,&mult,&deltaj,&jj);
		}
		else if (useFB1)
		{
			elementp = buildElement_FB1(i,i,1,&norm,&mult,&deltaj,&jj);
			elementm = buildElement_FB1(i,i,-1,&norm,&mult,&deltaj,&jj);
		}
		else if (useFB2)
		{
			elementp = buildElement_FB2(i,i,1,&norm,&mult,&deltaj,&jj);
			elementm = buildElement_FB2(i,i,-1,&norm,&mult,&deltaj,&jj);
		}

		diag[i-global_row] = norm * (!deltaj ? elementp.real()+mult*elementm.real() : -(double)jj*(elementp.imag()+mult*elementm.imag()));

		/****************************
		 * Elements before diagonal *
		 ****************************/
		
		for (int j = i%2; j < i; j+=2)
		{
			if (useSRLS)
			{
				elementp = buildElement_SRLS(i,j,1,&norm,&mult,&deltaj,&jj);
				elementm = buildElement_SRLS(i,j,-1,&norm,&mult,&deltaj,&jj);
			}
			else if (useFB1)
			{
				elementp = buildElement_FB1(i,j,1,&norm,&mult,&deltaj,&jj);
				elementm = buildElement_FB1(i,j,-1,&norm,&mult,&deltaj,&jj);
			}                                  
			else if (useFB2)                   
			{                                  
				elementp = buildElement_FB2(i,j,1,&norm,&mult,&deltaj,&jj);
				elementm = buildElement_FB2(i,j,-1,&norm,&mult,&deltaj,&jj);
			}

			element = norm * ( !deltaj ? elementp.real()+mult*elementm.real() : -(double)jj*(elementp.imag()+mult*elementm.imag()));
			if (fabs(element) > ZERO){
				rowElements.push_back(element);
				colIndexes.push_back(j);
			}	
		}

		
		/***************************
		 * Elements after diagonal *
		 ***************************/
		
		for (int j = i+1; j < ncols; j+=2)
		{
			if (useSRLS)
			{
				elementp = buildElement_SRLS(i,j,1,&norm,&mult,&deltaj,&jj);
				elementm = buildElement_SRLS(i,j,-1,&norm,&mult,&deltaj,&jj);
			}
			else if (useFB1)
			{
				elementp = buildElement_FB1(i,j,1,&norm,&mult,&deltaj,&jj);
				elementm = buildElement_FB1(i,j,-1,&norm,&mult,&deltaj,&jj);
			}                                  
			else if (useFB2)                   
			{                                  
				elementp = buildElement_FB2(i,j,1,&norm,&mult,&deltaj,&jj);
				elementm = buildElement_FB2(i,j,-1,&norm,&mult,&deltaj,&jj);
			}

			element = norm * ( !deltaj ? elementp.real()+mult*elementm.real() : -(double)jj*(elementp.imag()+mult*elementm.imag()));
			if (fabs(element) > ZERO){
				rowElements.push_back(element);
				colIndexes.push_back(j);
			}
		}
		
		/*****************************
		 * Add the row to the matrix *
		 *****************************/
		
		grow(&A,i-global_row,rowElements,colIndexes);
		// FOR DEBUGGING
		//std::cout << "MATRIX ROW: " << i << " / " << nrows << std::std::endl;
	}

	return;
	
}

/*****************************************
 * Build a matrix element for SRLS model *
 *****************************************/

dcomplex matrix::buildElement_SRLS(int i, int j, int s, double* nnorm, double* nmult, int* deltaj, int* jidx)
{
	
	/***************************************
	 * Chek if object has been initialized *
	 ***************************************/
	
	if (!isInit)
	{
		std::cout << std::endl << "ERROR in matrix::buildElement_SRLS(). Matrix object created but not initialized." << std::endl << std::endl;
		exit(0);
	}
		
	/*******************************
	 * Initialization of variables *
	 *******************************/

	int pos;
	int L1p, M1p, K1p, L2p, M2p, K2p, jjp;
	bool DL1, DM1, DK1, DL2, DM2, DK2, D1, D2, K2mK2pLE2, M2mM2pLE2, K1mK1pLE2;
	ivector idx;
	double mult;
	dcomplex el, el0, el2, factor;
	
	el0 = dcomplex(0.0,0.0);
	el2 = dcomplex(0.0,0.0);
	factor = dcomplex(0.0,0.0);
	
	/*************************
	 * Extract basis indexes *
	 *************************/
	
	switch (i-oldi)
	{
		case 0:
		{
			break;
		}
		default:
		{
			idx = bas->getBasisFunction(i);
			iL1 = idx[0]; iM1 = idx[1]; iK1 = idx[2];
			iL2 = idx[3]; iM2 = idx[4]; iK2 = idx[5];
			ijj = idx[6];
			oldi = i;
		}
	}
	
	idx = bas->getBasisFunction(j);
	L1p = idx[0]; M1p = idx[1]; K1p = s*idx[2];
	L2p = idx[3]; M2p = s*idx[4]; K2p = s*idx[5];
	jjp = idx[6];
	
	/*************************
	 * Normalization factors *
	 *************************/
	
	*nnorm = 2.0;
	*nnorm *= ((iK1 == 0 && iK2 == 0 && iM2 == 0) ? 0.5 : ONE_OVER_SQRT_TWO);
	*nnorm *= ((K1p == 0 && K2p == 0 && M2p == 0) ? 0.5 : ONE_OVER_SQRT_TWO);
	*nmult = (!((K1p+K2p+M2p)%2) ? (double)jjp : -(double)jjp);
	*deltaj = ijj - jjp;
	*jidx = ijj;
	
	/***************************************
	 * Evaluate useful logical expressions *
	 ***************************************/ 
	
	DL1 = !(iL1-L1p); DM1 = !(iM1-M1p); DK1 = !(iK1-K1p);
	DL2 = !(iL2-L2p); DM2 = !(iM2-M2p); DK2 = !(iK2-K2p);
	D1 = DL1 && DM1 && DK1;
	D2 = DL2 && DM2 && DK2;
	K2mK2pLE2 = abs(iK2-K2p) <= 2;
	M2mM2pLE2 = abs(iM2-M2p) <= 2;
	K1mK1pLE2 = abs(iK1-K1p) <= 2;
	
	/********************************************************
	 * Potential independent part of the diffusive operator *
	 ********************************************************/
	
	if (D1)
	{
		
		/**********************************
		 * 1) Ja(0,0) Jb(0,0) and Jc(0,0) *
		 **********************************/

		if (D2)
		{
#ifdef MAKEJA
			el0 += dcomplex(-probeD[0].real() * ONE_OVER_SQRT_THREE * (double)(iL2*(iL2+1)),0.0);
#endif
#ifdef MAKEJB
			el0 += dcomplex(-proteinD[0].real() * ONE_OVER_SQRT_THREE * (double)(iL2*(iL2+1)),0.0);
#endif
#ifdef MAKEJC
			el0 += dcomplex(-proteinD[0].real() * ONE_OVER_SQRT_THREE * (double)(iL1*(iL1+1)),0.0);
#endif
		}
	}
		
	/**********************************
	 * 2) Ja(2,m) Jb(2,m) and Jc(2,m) *
	 **********************************/

	if (D1 && DL2)
	{
#ifdef MAKEJA
		if (DM2 && K2mK2pLE2)
		{
			mult = (!((iL2-iK2)%2) ? 1.0 : -1.0);
			factor = dcomplex(mult * L2Mult[iL2] * trj->getTrj(iL2,2,iL2,-iK2,iK2-K2p,K2p),factor.imag());
			el2 += probeD[3+iK2-K2p] * factor;
		}
#endif
#ifdef MAKEJB
		if (DK2 && M2mM2pLE2)
		{
			mult = (!((iL2-iM2)%2) ? 1.0 : -1.0);
			factor = dcomplex(mult * L2Mult[iL2] * trj->getTrj(iL2,2,iL2,-iM2,iM2-M2p,M2p),factor.imag());
			el2 += proteinD[3+M2p-iM2] * factor;
		}
#endif
	}
#ifdef MAKEJC			
	if (D2 && DL1 && DM1 && K1mK1pLE2)
	{
		mult = (!((iL1-iK1)%2) ? 1.0 : -1.0);
		factor = dcomplex(mult * L2Mult[iL1] * trj->getTrj(iL1,2,iL1,-iK1,iK1-K1p,K1p),factor.imag());
		el2 += proteinD[3+iK1-K1p] * factor;
	}
#endif
	
	/*********
	 * 3) Jd *
	 *********/
	
#ifdef MAKEJD
	int L1idx = (iL1+1)*(iL1+1);
	int L2idx = (iL2+1)*(iL2+1);
	
	if (DL1 && DM1 && DL2 && DK2)
	{
		if (DK1 && DM2){
			el0 += dcomplex(2.0 * ONE_OVER_SQRT_THREE * proteinD[0].real() * (double)(iK1*iM2),0.0);
			factor = dcomplex(-4.0 * ONE_OVER_SQRT_SIX * (double)(iK1*iM2),factor.imag());
			el2 += proteinD[3] * factor;
		}
		else if ((iK1-K1p) == -1 && (iM2-M2p) == -1)
		{
			el0 += dcomplex(ONE_OVER_SQRT_THREE * proteinD[0].real() * cm[L1idx+iK1+1] * cm[L2idx+iM2+1],0.0);
			factor = dcomplex(ONE_OVER_SQRT_SIX * cm[L1idx+iK1+1] * cm[L2idx+iM2+1],factor.imag());
			el2 += proteinD[3] * factor;
		}
		else if ((iK1-K1p) == 1 && (iM2-M2p) == 1)
		{
			el0 += dcomplex(ONE_OVER_SQRT_THREE * proteinD[0].real() * cp[L1idx+iK1-1] * cp[L2idx+iM2-1],0.0);
			factor = dcomplex(ONE_OVER_SQRT_SIX * cp[L1idx+iK1-1] * cp[L2idx+iM2-1],factor.imag());
			el2 += proteinD[3] * factor;
		}
		else if ((iK1-K1p) == 1 && DM2)
		{
			factor = dcomplex(cp[L1idx+iK1-1] * (double)iM2,factor.imag());
			el2 += proteinD[4] * factor;
		}
		else if ((iK1-K1p) == -1 && DM2)
		{
			factor = dcomplex(-cm[L1idx+iK1+1] * (double)iM2,factor.imag());
			el2 += proteinD[2] * factor;
		}
		else if (DK1 && (iM2-M2p) == -1)
		{
			factor = dcomplex((double)iK1 * cm[L2idx+iM2+1]);
			el2 += proteinD[4] * factor;
		}
		else if (DK1 && (iM2-M2p) == 1)
		{
			factor = dcomplex(-(double)iK1 * cp[L2idx+iM2-1]);
			el2 += proteinD[2] * factor;
		}
		else if ((iK1-K1p) == 1 && (iM2-M2p) == -1)
		{
			factor = dcomplex(-cp[L1idx+iK1-1] * cm[L2idx+iM2+1],factor.imag());
			el2 += proteinD[5] * factor;
		}
		else if ((iK1-K1p) == -1 && (iM2-M2p) == 1)
		{
			factor = dcomplex(-cm[L1idx+iK1+1] * cp[L2idx+iM2-1],factor.imag());
			el2 += proteinD[1] * factor;
		}
    }
#endif
	
    el = el0 + el2;
    dcomplex tmpel = el;

	/******************************************************
	 * Potential dependent part of the diffusive operator *
	 ******************************************************/

#ifdef MAKEFA_OR_MAKEFB
	if (D1)
	{

		int jlow = max(max(max(0,abs(iK2-K2p)),abs(iM2-M2p)),abs(iL2-L2p));
		int jhigh = min(jlim,(iL2+L2p));
		
		dvector trjL2M2(5,0.0);
		dvector trjL2K2(jlim+1,0.0);
		for (int j = 0; j <= jlim; j++) trjL2K2[j] = trj->getTrj(iL2,j,L2p,-iK2,iK2-K2p,K2p);
		
		mult = (!((iM2-iK2)%2) ? 1.0 : -1.0)*sqrt((double)((2*iL2+1)*(2*L2p+1)));
			
		for (int mu = -mulim; mu <= mulim; mu+=2)
		{
			for (int mup = -mulim; mup <= mulim; mup+=2)
			{
				for (int j = jlow; j <= jhigh; j++)
				{
					
					pos = (mu+4)*9*9+(mup+4)*9+j;

					/*******
					 *  Fa *
					 *******/
#ifdef MAKEFA
					if (DM2)
					{
						trjL2M2[2] = trj->getTrj(iL2,j,L2p,-iM2,0,iM2);
						switch (iK2-K2p-mu-mup)
						{
							case (-2):
							{
								factor = dcomplex(mult * fa[1].at(pos) * trjL2M2[2] * trjL2K2[j],factor.imag());
								el += probeD[1] * factor;
								break;
							}
							case (-1):
							{
								factor = dcomplex(mult * fa[2].at(pos) * trjL2M2[2] * trjL2K2[j],factor.imag());
								el += probeD[2] * factor;
								break;
							}
							case (0):
							{
								factor = dcomplex(mult * fa[0].at(pos) * trjL2M2[2] * trjL2K2[j],factor.imag());
								el += probeD[0] * factor;
								factor = dcomplex(mult * fa[3].at(pos) * trjL2M2[2] * trjL2K2[j],factor.imag());
								el += probeD[3] * factor;
								break;
							}
							case (1):
							{
								factor = dcomplex(mult * fa[4].at(pos) * trjL2M2[2] * trjL2K2[j],factor.imag());
								el += probeD[4] * factor;
								break;
							}
							case (2):
							{
								factor = dcomplex(mult * fa[5].at(pos) * trjL2M2[2] * trjL2K2[j],factor.imag());
								el += probeD[5] * factor;
								break;
							}
						}
					}
#endif			
#ifdef MAKEFB
					/*******
					 *  Fb *
					 *******/
					
					if ((mu+mup == iK2-K2p))
					{
						trjL2M2[0] = trj->getTrj(iL2,j,L2p,-iM2,-2,M2p);
						trjL2M2[1] = trj->getTrj(iL2,j,L2p,-iM2,-1,M2p);
						trjL2M2[2] = trj->getTrj(iL2,j,L2p,-iM2, 0,M2p);
						trjL2M2[3] = trj->getTrj(iL2,j,L2p,-iM2, 1,M2p);
						trjL2M2[4] = trj->getTrj(iL2,j,L2p,-iM2, 2,M2p);
						
						switch(iM2-M2p)
						{
							case (-2):
							{
								factor = dcomplex(mult * fb[5].at(pos) * trjL2M2[0] * trjL2K2[j],factor.imag()); 
								el += proteinD[5] * factor;
								break;
							}
							case (-1):
							{
								factor = dcomplex(mult * fb[4].at(pos) * trjL2M2[1] * trjL2K2[j],factor.imag()); 
								el += proteinD[4] * factor;
								break;
							}
							case (0):
							{
								factor = dcomplex(mult * fb[0].at(pos) * trjL2M2[2] * trjL2K2[j],factor.imag()); 
								el += proteinD[0] * factor;
								factor = dcomplex(mult * fb[3].at(pos) * trjL2M2[2] * trjL2K2[j],factor.imag()); 
								el += proteinD[3] * factor;
								break;
							}
							case (1):
							{
								factor = dcomplex(mult * fb[2].at(pos) * trjL2M2[3] * trjL2K2[j],factor.imag()); 
								el += proteinD[2] * factor;
								break;
							}
							case (2):
							{
								factor = dcomplex(mult * fb[1].at(pos) * trjL2M2[4] * trjL2K2[j],factor.imag()); 
								el += proteinD[1] * factor;
								break;
							}
						}
					}
#endif
					
				}
			}
		}
		
	}
#endif
	//if (fabs(el.real()-tmpel.real())>0.0 | fabs(el.imag()-tmpel.imag())>0.0) std::cout << "EL " << iK2 << " " << K2p << std::endl;
		
	return el;
}

/*******************************
 * Return Sqrt[L*(L+1)-M(M+1)] *
 *******************************/

inline double matrix::cplm(int cL, int cM)
{
	if ( (cM < -cL) | (cM > cL) )
		return 0.0;
	return sqrt((double)(cL*(cL+1)-cM*(cM+1)));
}

/*******************************
 * Return Sqrt[L*(L+1)-M(M-1)] *
 *******************************/

inline double matrix::cmlm(int cL, int cM)
{
	if ( (cM < -cL) | (cM > cL) )
		return 0.0;
	return sqrt((double)(cL*(cL+1)-cM*(cM-1)));
}

/****************************************
 * Build a matrix element for FB1 model *
 ****************************************/

dcomplex matrix::buildElement_FB1(int i, int j, int s, double* nnorm, double* nmult, int* deltaj, int* jidx)
{
	/***************************************
	 * Chek if object has been initialized *
	 ***************************************/
	
	if (!isInit)
	{
		std::cout << std::endl << "ERROR in matrix::buildElement_FB1(int i, int j, int s). Matrix object created but not initialized." << std::endl << std::endl;
		exit(0);
	}
	
	/*******************************
	 * Initialization of variables *
	 *******************************/

	int pos;
	int N12, L2, M2, K2, j2;
	ivector idx;
	double mult;
	
	/*************************
	 * Extract basis indexes *
	 *************************/
	
	switch (i - oldi)
	{
		case 0:
		{
			break;
		}
		default:
		{
			idx = bas->getBasisFunction(i);
			iN11 = idx[0]; iL1 = idx[1]; iM1 = idx[2]; iK1 = idx[3]; ijj = idx[4];
			oldi = i;
		}
	}
	
	idx = bas->getBasisFunction(j);
	N12 = s*idx[0]; L2 = idx[1]; M2 = idx[2]; K2 = s*idx[3]; j2 = idx[4];

	/****************
	 * Various info *
	 ****************/
	
	*nnorm = 2.0;
	*nnorm *= ((iK1 == 0 && iN11 == 0) ? 0.5 : ONE_OVER_SQRT_TWO);
	*nnorm *= ((K2 == 0 && N12 == 0)   ? 0.5 : ONE_OVER_SQRT_TWO);
	*nmult = (!((L2+K2)%2) ? (double)j2 : -(double)j2);
	*deltaj = ijj - j2;
	*jidx = ijj;

	/*******************************************
	 * [1] Pure rotational part (no potential) *
	 *******************************************/	

	double grot = 0.0;

 	if ((iL1 == L2) && (iM1 == M2) /*&& (ijj == j2)*/ && (iN11 == N12))
	{

		switch(iK1-K2)
		{
			case 0:
        			grot = dzz*(double)(iK1*iK1) + .25*dpr*(CMENO(L2,K2)*CPIU(L2,K2-1)+CPIU(L2,K2)*CMENO(L2,K2+1));
        			break;
      			case 2:
 				grot = .25*dmr*CPIU(L2,K2)*CPIU(L2,K2+1);
				break;
			case -2:
				grot = .25*dmr*CMENO(L2,K2)*CMENO(L2,K2-1);
				break;
		}
	}

	/**************************
	 * [2] Pure internal part *
	 **************************/

  	int ncipmo = phy[0].getFB1NCoeff() - 1;
  	int dim = 2 * ncipmo + 1;
	dcomplex gint(0.0,0.0);

 	if ((iL1 == L2) && (iM1 == M2) && (iK1 == K2))
	{
		int dN1N2=iN11-N12;

		// Potential independent part

		if (dN1N2 == 0)
			gint.real((double)(iN11*iN11));

		// Potential dependent part

		for (int nu1 = -ncipmo; nu1 <= ncipmo; nu1++){
			for (int nu2 = -ncipmo; nu2 <= ncipmo; nu2++){
				if ((nu1+nu2)==dN1N2)
				{
					gint = gint - 0.25 * f11[(nu1+ncipmo)*dim+(nu2+ncipmo)];
				}
      			}
    		}
		gint = d11 * gint;
	}

	/**********************************
	 * [3] Rotational - internal part *
	 **********************************/

	dcomplex groin(0.0,0.0);

	if ((iL1==L2) && (iM1==M2) && (iN11==N12))
	{
		switch(iK1-K2)
		{
			case 0:
				groin.real(2.0*(double)(iN11*iK1)*dz1);
				break;
			case 1:
				groin = 2.0*(double)iN11*CPIU(iL1,K2)*dp1;
				break;
			case -1:
				groin = 2.0*(double)iN11*CMENO(iL1,K2)*dm1;
				break;
		}
	}

	/*********************/
	/* Sum up everything */
	/*********************/

	dcomplex el = grot + gint + groin;

	return el;

}

/****************************************
 * Build a matrix element for FB2 model *
 ****************************************/

dcomplex matrix::buildElement_FB2(int i, int j, int s, double* nnorm, double* nmult, int* deltaj, int* jidx)
{
	/***************************************
	 * Chek if object has been initialized *
	 ***************************************/
	
	if (!isInit)
	{
		std::cout << std::endl << "ERROR in matrix::buildElement_FB2(int i, int j, int s). Matrix object created but not initialized." << std::endl << std::endl;
		exit(0);
	}
	
	/*******************************
	 * Initialization of variables *
	 *******************************/

	int pos;
	int N12, N22, L2, M2, K2, j2;
	ivector idx;
	double mult;
	dcomplex czero(0.0,0.0);
	
	/*************************
	 * Extract basis indexes *
	 *************************/
	
	switch (i - oldi)
	{
		case 0:
		{
			break;
		}
		default:
		{
			idx = bas->getBasisFunction(i);
			iN11 = idx[0]; iN21 = idx[1]; iL1 = idx[2]; iM1 = idx[3]; iK1 = idx[4]; ijj = idx[5];
			oldi = i;
		}
	}
	
	idx = bas->getBasisFunction(j);
	N12 = s*idx[0]; N22 = s*idx[1]; L2 = idx[2]; M2 = idx[3]; K2 = s*idx[4]; j2 = idx[5];

	/****************
	 * Various info *
	 ****************/
	
	*nnorm = 2.0;
	*nnorm *= ((!iK1 && !iN11 && !iN21) ? 0.5 : ONE_OVER_SQRT_TWO);
	*nnorm *= ((!K2  && !N12  && !N22) ?  0.5 : ONE_OVER_SQRT_TWO);
	*nmult = (!((L2+K2)%2) ? (double)j2 : -(double)j2);
	*deltaj = ijj - j2;
	*jidx = ijj;

	/*******************************************
	 * [1] Pure rotational part (no potential) *
	 *******************************************/	

	dcomplex grot(0.0,0.0);

 	if ((iL1 == L2) && (iM1 == M2) && (iN11 == N12) && (iN21 == N22))
	{
		switch(iK1-K2)
		{
			case 0:
        			grot.real(dzz*(double)(iK1*iK1) + .25*dpr*(CMENO(L2,K2)*CPIU(L2,K2-1)+CPIU(L2,K2)*CMENO(L2,K2+1)));
        			break;
      			case 2:
 				grot.real(.25*dmr*CPIU(L2,K2)*CPIU(L2,K2+1));
				break;
			case -2:
				grot.real(.25*dmr*CMENO(L2,K2)*CMENO(L2,K2-1));
				break;
		}
		if (fabs(grot.real()) < ZERO*1.0e2) grot = 0.0;
	}

	/**************************
	 * [2] Pure internal part *
	 **************************/

	dcomplex g11(0.0,0.0), g22(0.0,0.0), g12(0.0,0.0), gint(0.0,0.0), ctmp1, ctmp2;

 	if ((iL1 == L2) && (iM1 == M2) && (iK1 == K2))
	{
		int delta1 = iN11 - N12;
		int delta2 = iN21 - N22;

		// Potential independent part

		if (delta1 == 0 && delta2 == 0)
		{
			g11.real((double)(iN11 * iN11));
			g22.real((double)(iN21 * iN21));
			g12.real(2.0 * (double)(iN11 * iN21));
		}

		// Potential dependent part

		for (int it1 = -ufb2.n1Max; it1 <= ufb2.n1Max; it1++)
		{
			if ( abs(delta1 - it1) <= ufb2.n1Max)
			{
				for (int it2 = -ufb2.n2Max; it2 <= ufb2.n2Max; it2++)
				{
					if ( abs(delta2 - it2) <= ufb2.n1Max)
					{
						ctmp1 = ufb2.c[(it1 + ufb2.n1Max) * ufb2.dim2 + (it2 + ufb2.n2Max)];
						ctmp2 = ufb2.c[(delta1 - it1 + ufb2.n1Max) * ufb2.dim2 + (delta2 - it2 + ufb2.n2Max)];

						/* 11 */
						g11 = g11 - 0.25 * (double)(it1 * (delta1 - it1)) * ctmp1 * ctmp2 + \
						      ( (delta1 == it1 && delta2 == it2) ? -0.5 * (double)(it1 * it1) * ctmp1 : czero );

						/* 22 */
						g22 = g22 - 0.25 * (double)(it2 * (delta2 - it2)) * ctmp1 * ctmp2 + \
							      ( (delta1 == it1 && delta2 == it2) ? -0.5 * (double)(it2 * it2) * ctmp1 : czero );

						/* 12 */
						g12 = g12 - 0.5 * (double)(it1 * (delta2 - it2)) * ctmp1 * ctmp2 + \
						      ( (delta1 == it1 && delta2 == it2) ? -1.0 * (double)(it1 * it2) * ctmp1 : czero );
					}
				}
			}
		}

		// Final update

		g11 = d11 * g11;
		if (fabs(g11.real()) < ZERO*1.0e2) g11.real(0.0);
		if (fabs(g11.imag()) < ZERO*1.0e2) g11.imag(0.0);
		g22 = d22 * g22;
		if (fabs(g22.real()) < ZERO*1.0e2) g22.real(0.0);
		if (fabs(g22.imag()) < ZERO*1.0e2) g22.imag(0.0);
		g12 = d12 * g12;
		if (fabs(g12.real()) < ZERO*1.0e2) g12.real(0.0);
		if (fabs(g12.imag()) < ZERO*1.0e2) g12.imag(0.0);

		gint = g11 + g22 + g12;
		if (fabs(gint.real()) < ZERO*1.0e2) gint.real(0.0);
		if (fabs(gint.imag()) < ZERO*1.0e2) gint.imag(0.0);
	}

	/**********************************
	 * [3] Rotational - internal part *
	 **********************************/

	dcomplex gr1(0.0,0.0), gr2(0.0,0.0), groin(0.0,0.0);

	if ((iL1 == L2) && (iM1 == M2) && (iN11 == N12) && (iN21 == N22))
	{
		switch(iK1-K2)
		{
			case 0:
			{
				gr1 = dcomplex(2.0 * (double)(iN11 * iK1) * dz1, 0.0);
				gr2 = dcomplex(2.0 * (double)(iN21 * iK1) * dz2, 0.0);
				break;
			}
			case 1:
			{
				gr1 = 2.0 * (double)iN11 * CPIU(L2, K2) * dp1;
				gr2 = 2.0 * (double)iN21 * CPIU(L2, K2) * dp2;
				break;
			}
			case -1:
			{
				gr1 = 2.0 * (double)iN11 * CMENO(L2, K2) * dm1;
				gr2 = 2.0 * (double)iN21 * CMENO(L2, K2) * dm2;
				break;
			}
		}
	}

	if (fabs(gr1.real()) < ZERO*1.0e2) gr1.real(0.0);
	if (fabs(gr1.imag()) < ZERO*1.0e2) gr1.imag(0.0);

	if (fabs(gr2.real()) < ZERO*1.0e2) gr2.real(0.0);
	if (fabs(gr2.imag()) < ZERO*1.0e2) gr2.imag(0.0);

	groin = gr1 + gr2;
	if (fabs(groin.real()) < ZERO*1.0e2) groin.real(0.0);
	if (fabs(groin.imag()) < ZERO*1.0e2) groin.imag(0.0);

	/***********************************/
	/* UPDATE THE TOTAL MATRIX ELEMENT */
	/***********************************/

	dcomplex el = grot + gint + groin;
	if (fabs(el.real()) < ZERO*1.0e2) el.real(0.0);
	if (fabs(el.imag()) < ZERO*1.0e2) el.imag(0.0);

	return el;

}

/******************************************************************************
 * Grow matrix by one row:                                                    *
 * cmatrix *mat      : I - matrix to be grown / O - matrix with the row added *
 * int row           : I - index of row to be added                           *
 * cdvector elements : I - vector of real and imaginary elements of the row   *
 * civector columns  : I - vector of real and immaginary column indexs        * 
 ******************************************************************************/

void matrix::grow(CompRow_Mat_double *mat, int row, dvector elements, ivector columns)
{
	if (elements.size() > 0){
		mat->rowptr_[row] = mat->nz_;
		mat->nz_ += elements.size();
		mat->val_.dim_ = mat->nz_;
		mat->colind_.dim_ = mat->nz_;
		if (row == 0)
		{
			mat->val_.p_ = (double *)calloc(mat->nz_, sizeof(double));
			mat->colind_.p_ = (int *)calloc(mat->nz_, sizeof(int));
		}
		else
		{
			mat->val_.p_ = (double *)realloc(mat->val_.p_, sizeof(double)*mat->nz_);
			mat->colind_.p_ = (int *)realloc(mat->colind_.p_, sizeof(int)*mat->nz_);
		}		
		for (unsigned int i = 0; i < elements.size(); i++) {
			mat->val_(i+mat->rowptr_(mat->dim_[0])) = elements[i];
			mat->colind_(i+mat->rowptr_(mat->dim_[0])) = columns[i];
		}
		mat->rowptr_(mat->dim_[0]) += elements.size();
	}
	else
		mat->rowptr_[row] = mat->nz_;	

	return;
}

/********************************************
 * Return the number of rows stored by task *
 ********************************************/

int matrix::getNrows(void)
{
	return nrows;
}

/****************************************************
 * Set the number of rows that the task must handle *
 ****************************************************/

void matrix::setNrows(int nr)
{
	nrows = nr;
	return;
}

/***************************************************** 
 * Return the starting global row in the full matrix *
 *****************************************************/ 

int matrix::getGlobalRow(void)
{
	return global_row;
}

/**************************************************
 * Set the starting global row in the full matrix *
 **************************************************/

void matrix::setGlobalRow(int g)
{
	global_row = g;
	return;
}

/***********************************************************
 * Return the number of non zero elements of the real part *
 ***********************************************************/

int matrix::getRealNZ(void)
{
	return realNZ;
}

/****************************************************************
 * Return the number of non zero elements of the imaginary part *
 ****************************************************************/

int matrix::getImagNZ(void)
{
	return imagNZ;
}

/**********************************************************
 * Return a pointer to the out-of-diagonal part of matrix *
 **********************************************************/

CompRow_Mat_double* matrix::getAptr(void)
{
	CompRow_Mat_double *Aptr;
	Aptr = &A;
	return Aptr;
}

/**************************************************
 * Return a pointer to the diagonal of the matrix *
 **************************************************/

VECTOR_double* matrix::getDiagPtr(void)
{
	VECTOR_double *Dptr;
	Dptr = &diag;
	return Dptr;
}

/********************************
 * Print informations on matrix *
 ********************************/

std::string matrix::toString(int realFlag, int preFlag)
{
	std::ostringstream ostr;
	std::ostringstream tmp;
	
	if (realFlag)
	{
		ostr << "TASK " << rank << ": Diagonal Matrix Elements" << std::endl;
		ostr << "--------------------------------------------------------" << std::endl;
		ostr << setw(14) << "row" << setw(18) << "column" << setw(24) << "element" << std::endl;
		for (int i = 0; i < nrows; i++)
			ostr << setw(14) << i+global_row+1 << setw(18) << i+global_row+1 << setw(24) << setiosflags(ios::scientific) << setprecision(12) << diag[i] << std::endl;
		ostr << "--------------------------------------------------------" << std::endl;
		
		ostr << "TASK " << rank << ": Out of Diagonal Matrix Elements" << std::endl;
		ostr << "--------------------------------------------------------" << std::endl;
		ostr << setw(14) << "row" << setw(18) << "column" << setw(24) << "element" << std::endl;
		ostr << "--------------------------------------------------------" << std::endl;
		ostr << A;
		ostr << "--------------------------------------------------------" << std::endl;
	}
	
	if ( preFlag && ( !(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")) ) ) // SRLS only
	{
		int w = 20;
		int p = 4;
		int size = 0;
		int oldSize = 0;
		
		ostr << "TASK " << rank << ": Coefficients fa" << std::endl;
		tmp << "(mu,mup,j)";
		oldSize = tmp.str().size();
		ostr << "TASK " << rank << ": " << setfill('-') << setw(w*7-10) << "" << setfill(' ') << std::endl;
		ostr << "TASK " << rank << ": " << tmp.str() << setw(w-3) << "fa(0,0)" << setw(w) << "fa(2,-2)" << setw(w) << "fa(2,-1)" << setw(w) << "fa(2,0)" << setw(w) << "fa(2,1)" << setw(w) << "fa(2,2)" << std::endl;
		ostr << "TASK " << rank << ": " << setfill('-') << setw(w*7-10) << "" << setfill(' ') << std::endl;
		int pos = 0;
		tmp.str("");
		for (int mu = -4; mu <= 4; mu++)
		{
			for (int mup = -4; mup <= 4; mup++)
			{
				for (int j = 0; j <= 8; j++)
				{
					tmp << "(" << mu << "," << mup << "," << j << ")";
					size = tmp.str().size();
					ostr << "TASK " << rank << ": " << setfill(' ') << tmp.str() << setw(w+(oldSize-size)) << showpos << setiosflags(ios::scientific) << setprecision(p) << fa[0].at(pos) << setw(w) << setprecision(p) << fa[1].at(pos) << setw(w) << setprecision(p) << fa[2].at(pos) << setw(w) << setprecision(p) << fa[3].at(pos) << setw(w) << setprecision(p) << fa[4].at(pos) << setw(w) << setprecision(p) << fa[5].at(pos) << noshowpos << std::endl;
					pos ++;
					tmp.str("");
				}
			}
		}
		ostr << "TASK " << rank << ": " << setfill('-') << setw(w*7-10) << "" << setfill(' ') << std::endl;
	
		ostr << "TASK " << rank << ": Coefficients fb" << std::endl;
		tmp << "(mu,mup,j)";
		oldSize = tmp.str().size();
		ostr << "TASK " << rank << ": " << setfill('-') << setw(w*7-10) << "" << setfill(' ') << std::endl;
		ostr << "TASK " << rank << ": " << tmp.str() << setw(w-3) << "fb(0,0)" << setw(w) << "fb(2,-2)" << setw(w) << "fb(2,-1)" << setw(w) << "fb(2,0)" << setw(w) << "fb(2,1)" << setw(w) << "fb(2,2)" << std::endl;
		ostr << "TASK " << rank << ": " << setfill('-') << setw(w*7-10) << "" << setfill(' ') << std::endl;
		pos = 0;
		tmp.str("");
		for (int mu = -4; mu <= 4; mu++)
		{
			for (int mup = -4; mup <= 4; mup++)
			{
				for (int j = 0; j <= 8; j++)
				{
					tmp << "(" << mu << "," << mup << "," << j << ")";
					size = tmp.str().size();
					ostr << "TASK " << rank << ": " << setfill(' ') << tmp.str() << setw(w+(oldSize-size)) << showpos << setiosflags(ios::scientific) << setprecision(p) << fb[0].at(pos) << setw(w) << setprecision(p) << fb[1].at(pos) << setw(w) << setprecision(p) << fb[2].at(pos) << setw(w) << setprecision(p) << fb[3].at(pos) << setw(w) << setprecision(p) << fb[4].at(pos) << setw(w) << setprecision(p) << fb[5].at(pos) << noshowpos << std::endl;
					pos ++;
					tmp.str("");
				}
			}
		}
		ostr << "TASK " << rank << ": " << setfill('-') << setw(w*7-10) << "" << setfill(' ') << std::endl;

	}
	
	return ostr.str();

}

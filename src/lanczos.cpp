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
 Name        : lanczos.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Class to perform Lanczos tridiagonalization and calculation 
               of spectral density
 ============================================================================
 */
#include "lanczos.h"

extern std::string path, project;

/***************
 * Constructor *
 ***************/

lanczos::lanczos()
{
}

/**************
 * Destructor *
 **************/

lanczos::~lanczos()
{
#ifdef WRITE_DESTORY_MESSAGE
	std::cout << "TASK " << rank << ": Cleared lanczos object" << std::endl;
#endif
}

/**************
 * Initiators *
 **************/

void lanczos::init(int n, stvec* sv, matrix* m, ivector js)
{
	rank = 0;
	nstep = nstep0 = n;
	v = sv;
	mat = m;
	sd = js;
	jnumber = sd.size();
	alpha = dvvector(2*jnumber);
	beta  = dvvector(2*jnumber);
	lastRun = oldRun = false;
	return;
}

void lanczos::init(int n, stvec* sv, matrix *m, ivector js, int r, int nr)
{
	rank = r;
	nproc = nr;
	nstep = nstep0 = n;
	v = sv;
	mat = m;
	sd = js;
	jnumber = sd.size();
	alpha = dvvector(jnumber<<1);
	beta  = dvvector(jnumber<<1);
	lastRun = oldRun = false;
	return;
}

/************************************************** 
 * Method to update object in a fitting procedure *
 **************************************************/

void lanczos::update(void)
{
	alpha = dvvector(jnumber<<1);
	beta  = dvvector(jnumber<<1);
	return;
}

/*********************
 * Run Lanczos cycle *
 *********************/

void lanczos::runLaczos(void)
{
	
	VECTOR_double *x = new VECTOR_double[jnumber<<1];
	VECTOR_double *y = new VECTOR_double[jnumber<<1];
	VECTOR_double *diag;
	CompRow_Mat_double* A;
	
	A = mat->getAptr();
	diag = mat->getDiagPtr();
#ifdef _MPI_
	nbf = 0;
	int local_nbf = diag[0].size();
	MPI_Allreduce(&local_nbf,&nbf,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#else
	nbf = diag[0].size();
#endif
	VECTOR_double u = VECTOR_double(nbf,0.0);
	VECTOR_double ug = VECTOR_double(nbf,0.0);
	
	/***********************/
	/* Send - Receive info */
	/***********************/

#ifdef _MPI_
//	int *grArray;
//	grArray = (int *)calloc(nproc, sizeof(int));
//	grArray[rank] = mat->getGlobalRow();
//	MPI_Allreduce(&grArray[0], &grArray[0], nproc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//	int *nrArray;
//	nrArray = (int *)calloc(nproc, sizeof(int));
//	nrArray[rank] = mat->getNrows();
//	MPI_Allreduce(&nrArray[0], &nrArray[0], nproc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

	/******************************************
	 * Determine which j have to be evaluated *
	 ******************************************/
	
	for (int i = 0; i < (jnumber<<1); i+=2)
	{
		if (v->hasBeenInit() && v->getNpmKK1(i) > ZERO)
		{
			x[i] = v->getProjections(i);
			y[i] = VECTOR_double(nbf,0.0);
		}
		
		if (v->hasBeenInit() && v->getNpmKK1(i+1) > ZERO)
		{
			x[i+1] = v->getProjections(i+1);
			y[i+1] = VECTOR_double(nbf,0.0);
		}
		 
	}
	
	/*****************
	 * Lanczos cycle *
	 *****************/

	for (int i = 0; i < nstep; i++)
	{
		
		/*********************************************
		 * Cycle over symmetrized spectral densities *
		 *********************************************/
		
		for (int j = 0; j < (jnumber<<1); j+=2)
		{

			/******
			 * j+ *
			 ******/
			
			if (v->hasBeenInit() && v->getNpmKK1(j) > ZERO && (!i || (i && fabs(beta[j].back())>ZERO)))
			{

					/*******************************************
					 * Perform spzerse matrix - vector product *
					 ********************************************/
					
					smvm(A,diag,&x[j],&u);

#ifdef _MPI_
					MPI_Allreduce(&u[0],&ug[0],nbf,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
					ug = u;
#endif
					
					/************
					 * Update y *
					 ************/
				
					y[j] += ug;
				
					/*******************
					 * Calculate alpha *
					 *******************/

					alpha[j].push_back(ddotu(&x[j],&y[j]));
					
					/*************************
					 * Update Y = Y - alphaX *
					 *************************/

					daxpy(&x[j],&y[j],-alpha[j].back());
					
					/******************
					 * Calculate beta *
					 ******************/

					beta[j].push_back(sqrt(ddotu(&y[j],&y[j])));
	
					/*********************
					 * Rearrange vectors *
					 *********************/
			
					dscsw(&y[j],&x[j],beta[j].back());
			
			} // j+ != 0 condition
			
			/******
			 * j- *
			 ******/

			if (v->hasBeenInit() && v->getNpmKK1(j+1) > ZERO  && (!i || (i && fabs(beta[j+1].back())>ZERO)))
			{

					/*******************************************
					 * Perform spzerse matrix - vector product *
					 ********************************************/
				
					smvm(A,diag,&x[j+1],&u);
	
#ifdef _MPI_
					MPI_Allreduce(&u[0],&ug[0],nbf,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
					ug = u;
#endif
			
					/************
					 * Update y *
					 ************/
					
					y[j+1] += ug;
					
					/*******************
					 * Calculate alpha *
					 *******************/
					
					alpha[j+1].push_back(ddotu(&x[j+1],&y[j+1]));
					
					/************
					 * Update y *
					 ************/
					
					daxpy(&x[j+1],&y[j+1],-alpha[j+1].back());
					
					/******************
					 * Calculate beta *
					 ******************/
					
					beta[j+1].push_back(sqrt(ddotu(&y[j+1],&y[j+1])));
					
					/*********************
					 * Rearrange vectors *
					 *********************/
					
					dscsw(&y[j+1],&x[j+1],beta[j+1].back());
				
			} // j- != 0 condition

		} // End cycle over j
	
	// DEBUGGING
	//std::cout << "LANCZOS STEP " << i+1 << " / " << nstep << std::std::endl;
	} // End of Lanczos cycle

	lastRun = true;
	oldRun = false;

	return;
}

/***************************************
 * Sparse Matrix Vector Multiplication *
 ***************************************/

void lanczos::smvm(CompRow_Mat_double* A, VECTOR_double* diag, VECTOR_double* b, VECTOR_double* u)
{
	int global_row = mat->getGlobalRow();

	VECTOR_double z;
	VECTOR_double w(nbf,0.0);
	
	/********************* 
	 * Zero out u vector *
	 *********************/
	
	for (int i = 0; i < nbf; i++) u[0](i) = 0.0;
	
	/************************
	 * Out of diagonal part *
	 ************************/
	
	z.dim_ = mat->getNrows();
	z.p_ = new double(z.dim_);
	z.p_ = &u[0].p_[global_row];
	z.ref_ = 1;
	// MATRIX IS STORED IN 'CHESS-LIKE' PATTERN THUS THE MATRIX
	// VECTOR MULTIPLICATION IS SPLIT INTO 2 STEPS
        if (A[0].nz_ > 0)
        {
	  // STEP 1: COLUMNS CONTRIBUTION
          z = A[0]*b[0]; // THIS UPDATES u
	  // STEP 2: ROWS CONTRIBUTION
	  z.p_ = &b[0].p_[global_row];
	  w = A[0].trans_mult(z);
	  u[0] += w;
        }
	
	/*****************
	 * Diagonal part *
	 *****************/
	
	for (int i = 0; i < mat->getNrows(); i++)
		u[0](i+global_row) += diag[0](i) * b[0](i+global_row);

	return;
}

/********************
 * Internal product *
 ********************/

double lanczos::ddotu(VECTOR_double* x, VECTOR_double* y)
{
	double dot = 0.0;
	for (int i = 0; i < nbf; i++)
		dot += x[0](i)*y[0](i);
	return dot;
}

/******************************
 * Vector scaling and summing *
 * Y = aX + Y                 *
 ******************************/

void lanczos::daxpy(VECTOR_double* x, VECTOR_double* y, double a)
{
	for (int i = 0; i < nbf; i++)
		y[0](i) = a*x[0](i) + y[0](i);
	return;
}

/***********************
 * Swap vectors        *
 * Y = X/a and X = -aY *
 ***********************/

void lanczos::dscsw(VECTOR_double* x, VECTOR_double* y, double a)
{
	double scaley = 1.0/a;
	double scalex = -a;
	double tmp;
	for (int i = 0; i < nbf; i++)
	{
		tmp = y[0](i);
		y[0](i) = scaley * x[0](i);
		x[0](i) = scalex * tmp;
	}
	return;
}

/******************************************************************
 * Calculate a spectral density given alpha, beta and frequencies *
 ******************************************************************/

dcvector lanczos::calculateSpectralDensity(dvector a, dvector b, dvector freq, double intrinsic, bool returnFirstDerivative)
{
	int nab = a.size();
	int nfreq = freq.size();
	dcvector jw;
	int iz,k;
	dcomplex tv,td,s,x;
	dcomplex rone,ione;
	rone = dcomplex(1.0, 0.0);
	ione = dcomplex(0.0,1.0);
	
	/******************************************************************
	 * Compute 0th derivative of continued fraction with respect to z *
	 ******************************************************************/
	
	if (!returnFirstDerivative)
	{
		for (iz = 0; iz < nfreq; iz++)
		{
			x = dcomplex(intrinsic, freq[iz]);
			s = rone / (a[nab-1] + x);
			tv = s * b[nab-2] * b[nab-2];
			for (k = nab-2; k >= 1; k--)
			{
				s = rone / (a[k] + x - tv);
				tv = s * b[k-1] * b[k-1];
			}
			s = rone / (a[0] + x - tv);
			tv = s;
			jw.push_back(tv);
		}	
	}
	/*******************************************************************
	 *  Compute 1st derivative of continued fraction with respect to z *
	 *******************************************************************/
	else
	{
		for (iz = 0; iz < nfreq; iz++)
		{
			x = dcomplex(intrinsic, freq[iz]);
			s = rone / (a[nab-1] + x);
			tv = s * b[nab-2] * b[nab-2];
			td = -tv * s;
			for (k = nab-2; k >=1; k--)
			{
				s = rone / (a[k] + x - tv);
				tv = s * b[k-1] * b[k-1];
				td = -tv * (rone - td) * s;
			}
			s = rone / (a[0] + x - tv);
			tv = s;
			jw.push_back(ione * tv * (rone - td) * s);
		}
	}
	
	return jw;
}

/*****************************************
 * Calculates all the spectral densities *
 *****************************************/

dcvector* lanczos::calculateSpectralDensities(dvector freq, double intrinsic, bool returnFirstDerivative, bool eigensystem, ldvector omegaDip, ldvector omegaDip2, ldvector omegaCsa, int nH)
{

	if (!lastRun && !oldRun)
	{
		std::cout << "TASK " << rank << ": ERROR : in lanczos::calculateSpectralDensities. Lanczos has not been run yet" << std::endl;
		exit(1);
	}
#ifdef WRITE_ALL
#ifdef WRITE_LANCZOS_WARNINGS
	if (!lastRun && oldRun)
		std::cout << "TASK " << rank << ": WARNING: in lanczos::calculateSpectralDensities. Alpha and Beta correspond to an old run Lanczos" << std::endl;
#endif
#endif
	dcvector *jw;
	dcvector czero(5,complex<double>(0.0,0.0));
	
	jw = new dcvector[2*jnumber];

	dvvector ew;
	dvvector cf;
	dvvector sf;
	dvvector spd(2*jnumber);
	dvector ww;
	
	std::string separator;
#ifdef _LINUX_
	separator = "/";
#else
	separator = "\\";
#endif

	fstream outFile;
	std::string fileName;
	fstream cfFile;
	std::string cfname;
	int NT;
	double minT, maxT;
	int NW;
	double minW, maxW;
	bool wout = !rank && eigensystem;
	if (wout)
	{
		fileName = path + separator + project + "_tw.par";
		outFile.open(fileName.c_str(),ios::in);
		if (outFile.is_open())
		{
			outFile >> NT;
			outFile >> minT;
			outFile >> maxT;
			outFile >> NW;
			outFile >> minW;
			outFile >> maxW;
			outFile.close();
			fileName = path + separator + project + ".eig";
			outFile.open(fileName.c_str(),ios::out);
		}
		else
		{
			std::cout << std::endl << std::endl  << "ERROR: cannot access file " << fileName << std::endl << std::endl;
			exit(1);
		}
	}
	
	int k,k1;
	char kch[4], k1ch[4];
	for (int i = 0; i < 2*jnumber; i+=2)
	{
		if (v->hasBeenInit() && v->getNpmKK1(i) > ZERO)
		{
			jw[i] = calculateSpectralDensity(alpha[i],beta[i],freq,intrinsic,false);
			if (wout)
			{
				k =  sd[i>>1] / 5 - 2;
				k1 = sd[i>>1] % 5 - 2;
				
				/**************************************************
				 * Calculate eigenvalues and weights for C+(k,k1) *
				 **************************************************/
				
				ew = tqli(alpha[i].size(),alpha[i],beta[i]);
				outFile << "*** C+, " << k << ", " << k1 << " ***" << std::endl << std::endl;
				for (unsigned int h = 0; h < alpha[i].size(); h++)
					outFile << ew[0].at(h) << "\t" << ew[1].at(h) << std::endl;
				outFile << "**********************************************" << std::endl << std::endl;
				
				/*******************************************************
				 * Calculate the time correlation function C+(k,k1)(t) *
				 *******************************************************/
				
				sprintf(kch,"%d",k);
				sprintf(k1ch,"%d",k1);
				cf = corrFunc(ew[0],ew[1],NT,minT,maxT);
				cfname = path + separator + project + "_C+_" + kch + "_" + k1ch + ".dat";
				cfFile.open(cfname.c_str(),ios::out);
				for (unsigned int h = 0; h < cf[0].size(); h++)
					cfFile << cf[0].at(h) << " " << cf[1].at(h) << std::endl;
				cfFile.close();
				
				/**********************************************
				 * Calculate the spectral density j+(k,k1)(w) *
				 **********************************************/
				
				sf = specDens(ew[0],ew[1],NW,minW,maxW);
				cfname = path + separator + project + "_j+_" + kch + "_" + k1ch + ".dat";
				cfFile.open(cfname.c_str(),ios::out);
				for (unsigned int h = 0; h < sf[0].size(); h++)
					cfFile << sf[0].at(h) << " " << sf[1].at(h) << std::endl;
				cfFile.close();
				spd[i] = sf[1];
				ww = sf[0];
			}
		}
		else
		{
			jw[i] = czero;
			if (wout)
				spd[i] = dvector(NW,0.0);
		}
		if (v->hasBeenInit() && v->getNpmKK1(i+1) > ZERO)
		{
			jw[i+1] = calculateSpectralDensity(alpha[i+1],beta[i+1],freq,intrinsic,false);
			if (wout)
			{
				k =  sd[i>>1] / 5 - 2;
				k1 = sd[i>>1] % 5 - 2;
				
				/**************************************************
				 * Calculate eigenvalues and weights for C-(k,k1) *
				 **************************************************/

				ew = tqli(alpha[i+1].size(),alpha[i+1],beta[i+1]);
				outFile << "*** C-, " << k << ", " << k1 << " ***" << std::endl << std::endl;
				for (unsigned int h = 0; h < alpha[i+1].size(); h++)
					outFile << ew[0].at(h) << "\t" << ew[1].at(h) << std::endl;
				outFile << "**********************************************" << std::endl << std::endl;
				
				/*******************************************************
				 * Calculate the time correlation function C-(k,k1)(t) *
				 *******************************************************/
				
				sprintf(kch,"%d",k);
				sprintf(k1ch,"%d",k1);
				cf = corrFunc(ew[0],ew[1],NT,minT,maxT);
				cfname = path + separator + project + "_C-_" + kch + "_" + k1ch + ".dat";
				cfFile.open(cfname.c_str(),ios::out);
				for (unsigned int h = 0; h < cf[0].size(); h++)
					cfFile << cf[0].at(h) << " " << cf[1].at(h) << std::endl;
				cfFile.close();
				
				/**********************************************
				 * Calculate the spectral density j-(k,k1)(w) *
				 **********************************************/
				
				sf = specDens(ew[0],ew[1],NW,minW,maxW);
				cfname = path + separator + project + "_j-_" + kch + "_" + k1ch + ".dat";
				cfFile.open(cfname.c_str(),ios::out);
				for (unsigned int h = 0; h < sf[0].size(); h++)
					cfFile << sf[0].at(h) << " " << sf[1].at(h) << std::endl;
				cfFile.close();
				spd[i+1] = sf[1];
				ww = sf[0];
			}
		}
		else
		{
			jw[i+1] = czero;
			if (wout)
				spd[i+1] = dvector(NW,0.0);
		}
	}
	
	/********************************
	 * Calculate Jxx(w) if required *
	 ********************************/
	
	if (wout)
	{

		int ikk1, k, k1, ikk, ik1k1, z, z1;
		unsigned int y, h;
		dvector* smallj = new dvector[jnumber];

		char smjn[2];
		std::string fileName;
		fstream file;
		
		/************************
		 * Cycle over small j's *
		 ************************/
		
		for (int i = 0; i < 2*jnumber; i+=2) {
			y = i>>1;				
			smallj[y].clear();
			ikk1 = sd[y];
			k = ikk1/5;  ikk   = (abs(k-2)+2)*5+(abs(k-2)+2);
			k1 = ikk1%5; ik1k1 = (abs(k1-2)+2)*5+(abs(k1-2)+2);
			for (int j = 0; j < jnumber; j++) {
				if (sd[j] == ikk)   z  = 2*j;
				if (sd[j] == ik1k1) z1 = 2*j;
			}

			sprintf(smjn,"%d",y);
			fileName = path + separator + project + "_smallj_" + smjn + ".dat";
			file.open(fileName.c_str(),ios::out);

			for (h = 0; h < ww.size(); h++) {
				smallj[y].push_back(4.0*(v->getNpmKK1(i)*spd[i].at(h)	+ v->getNpmKK1(i+1)*spd[i+1].at(h)));
				smallj[y].back() -= (v->getNpmKK1(z)*spd[z].at(h)	+ v->getNpmKK1(z+1)*spd[z+1].at(h));
				smallj[y].back() -= (v->getNpmKK1(z1)*spd[z1].at(h) + v->getNpmKK1(z1+1)*spd[z1+1].at(h));
				smallj[y].back() *= 0.125;
				file << ww.at(h) << " " << smallj[y].back() << std::endl;
			}
			file.close();
		}
					
		/*****************************************
		 * Calculate weights for CSA interaction *
		 *****************************************/
		
		dcomplex *dk0csa, *dk0dip, *dk0dip2;
		ldcomplex dlmk;
		wigner w(10);
		
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
				dk0dip2[0] += (dcomplex)(w.getWignerMatrix(2,-2,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
				dk0dip2[1] += (dcomplex)(w.getWignerMatrix(2,-1,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
				dk0dip2[2] += (dcomplex)(w.getWignerMatrix(2, 0,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
				dk0dip2[3] += (dcomplex)(w.getWignerMatrix(2, 1,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
				dk0dip2[4] += (dcomplex)(w.getWignerMatrix(2, 2,i,omegaDip[0],omegaDip[1],omegaDip[2]) * dlmk);
			}
		}
		
		std::string name = path + separator + project + "_Jdip.dat";
		fstream jdipFile(name.c_str(),ios::out);
		name = path + separator + project + "_Jcsa.dat";
		fstream jcsaFile(name.c_str(),ios::out);
		name = path + separator + project + "_Jdip2.dat";
		fstream jdip2File(name.c_str(),ios::out);
		name = path + separator + project + "_Jdip1-dip2.dat";
		fstream jdip1dip2File(name.c_str(),ios::out);

		double bigJdip, bigJcsa, bigJdip2, bigKdip;

		double f_0_0_dip   = (conj(dk0dip[2])*dk0dip[2]).real();
		double f_0_p1_dip  = (conj(dk0dip[2])*dk0dip[3]).real();
		double f_0_p2_dip  = (conj(dk0dip[2])*dk0dip[4]).real();
		double f_m1_0_dip  = (conj(dk0dip[1])*dk0dip[2]).real();
		double f_m1_m1_dip = (conj(dk0dip[1])*dk0dip[1]).real();
		double f_m1_p2_dip = (conj(dk0dip[1])*dk0dip[4]).real();
		double f_p1_p1_dip = (conj(dk0dip[3])*dk0dip[3]).real();
		double f_p1_p2_dip = (conj(dk0dip[3])*dk0dip[4]).real();
		double f_m1_p1_dip = (conj(dk0dip[1])*dk0dip[3]).real();
		double f_m2_m1_dip = (conj(dk0dip[0])*dk0dip[1]).real();
		double f_m2_m2_dip = (conj(dk0dip[0])*dk0dip[0]).real();
		double f_p2_p2_dip = (conj(dk0dip[4])*dk0dip[4]).real();
		double f_m2_0_dip  = (conj(dk0dip[0])*dk0dip[2]).real();
		double f_m2_p1_dip = (conj(dk0dip[0])*dk0dip[3]).real();
		double f_m2_p2_dip = (conj(dk0dip[0])*dk0dip[4]).real();

		double f_0_0_csa   = (conj(dk0csa[2])*dk0csa[2]).real();
		double f_0_p1_csa  = (conj(dk0csa[2])*dk0csa[3]).real();
		double f_0_p2_csa  = (conj(dk0csa[2])*dk0csa[4]).real();
		double f_m1_0_csa  = (conj(dk0csa[1])*dk0csa[2]).real();
		double f_m1_m1_csa = (conj(dk0csa[1])*dk0csa[1]).real();
		double f_m1_p2_csa = (conj(dk0csa[1])*dk0csa[4]).real();
		double f_p1_p1_csa = (conj(dk0csa[3])*dk0csa[3]).real();
		double f_p1_p2_csa = (conj(dk0csa[3])*dk0csa[4]).real();
		double f_m1_p1_csa = (conj(dk0csa[1])*dk0csa[3]).real();
		double f_m2_m1_csa = (conj(dk0csa[0])*dk0csa[1]).real();
		double f_m2_m2_csa = (conj(dk0csa[0])*dk0csa[0]).real();
		double f_p2_p2_csa = (conj(dk0csa[4])*dk0csa[4]).real();
		double f_m2_0_csa  = (conj(dk0csa[0])*dk0csa[2]).real();
		double f_m2_p1_csa = (conj(dk0csa[0])*dk0csa[3]).real();
		double f_m2_p2_csa = (conj(dk0csa[0])*dk0csa[4]).real();


		double f_0_0_dip2   = 0.0;
		double f_0_p1_dip2  = 0.0;
		double f_0_p2_dip2  = 0.0;
		double f_m1_0_dip2  = 0.0;
		double f_m1_m1_dip2 = 0.0;
		double f_m1_p2_dip2 = 0.0;
		double f_p1_p1_dip2 = 0.0;
		double f_p1_p2_dip2 = 0.0;
		double f_m1_p1_dip2 = 0.0;
		double f_m2_m1_dip2 = 0.0;
		double f_m2_m2_dip2 = 0.0;
		double f_p2_p2_dip2 = 0.0;
		double f_m2_0_dip2  = 0.0;
		double f_m2_p1_dip2 = 0.0;
		double f_m2_p2_dip2 = 0.0;

		double f_0_0_d1d2   = 0.0;
		double f_0_p1_d1d2  = 0.0;
		double f_0_p2_d1d2  = 0.0;
		double f_m1_0_d1d2  = 0.0;
		double f_m1_m1_d1d2 = 0.0;
		double f_m1_p2_d1d2 = 0.0;
		double f_p1_p1_d1d2 = 0.0;
		double f_p1_p2_d1d2 = 0.0;
		double f_m1_p1_d1d2 = 0.0;
		double f_m2_m1_d1d2 = 0.0;
		double f_m2_m2_d1d2 = 0.0;
		double f_p2_p2_d1d2 = 0.0;
		double f_m2_0_d1d2  = 0.0;
		double f_m2_p1_d1d2 = 0.0;
		double f_m2_p2_d1d2 = 0.0;

		if (nH == 2)
		{
			f_0_0_dip2   = (conj(dk0dip2[2])*dk0dip2[2]).real();
			f_0_p1_dip2  = (conj(dk0dip2[2])*dk0dip2[3]).real();
			f_0_p2_dip2  = (conj(dk0dip2[2])*dk0dip2[4]).real();
			f_m1_0_dip2  = (conj(dk0dip2[1])*dk0dip2[2]).real();
			f_m1_m1_dip2 = (conj(dk0dip2[1])*dk0dip2[1]).real();
			f_m1_p2_dip2 = (conj(dk0dip2[1])*dk0dip2[4]).real();
			f_p1_p1_dip2 = (conj(dk0dip2[3])*dk0dip2[3]).real();
			f_p1_p2_dip2 = (conj(dk0dip2[3])*dk0dip2[4]).real();
			f_m1_p1_dip2 = (conj(dk0dip2[1])*dk0dip2[3]).real();
			f_m2_m1_dip2 = (conj(dk0dip2[0])*dk0dip2[1]).real();
			f_m2_m2_dip2 = (conj(dk0dip2[0])*dk0dip2[0]).real();
			f_p2_p2_dip2 = (conj(dk0dip2[4])*dk0dip2[4]).real();
			f_m2_0_dip2  = (conj(dk0dip2[0])*dk0dip2[2]).real();
			f_m2_p1_dip2 = (conj(dk0dip2[0])*dk0dip2[3]).real();
			f_m2_p2_dip2 = (conj(dk0dip2[0])*dk0dip2[4]).real();

			f_0_0_d1d2   = (conj(dk0dip[2])*dk0dip2[2]).real();
			f_0_p1_d1d2  = (conj(dk0dip[2])*dk0dip2[3]).real();
			f_0_p2_d1d2  = (conj(dk0dip[2])*dk0dip2[4]).real();
			f_m1_0_d1d2  = (conj(dk0dip[1])*dk0dip2[2]).real();
			f_m1_m1_d1d2 = (conj(dk0dip[1])*dk0dip2[1]).real();
			f_m1_p2_d1d2 = (conj(dk0dip[1])*dk0dip2[4]).real();
			f_p1_p1_d1d2 = (conj(dk0dip[3])*dk0dip2[3]).real();
			f_p1_p2_d1d2 = (conj(dk0dip[3])*dk0dip2[4]).real();
			f_m1_p1_d1d2 = (conj(dk0dip[1])*dk0dip2[3]).real();
			f_m2_m1_d1d2 = (conj(dk0dip[0])*dk0dip2[1]).real();
			f_m2_m2_d1d2 = (conj(dk0dip[0])*dk0dip2[0]).real();
			f_p2_p2_d1d2 = (conj(dk0dip[4])*dk0dip2[4]).real();
			f_m2_0_d1d2  = (conj(dk0dip[0])*dk0dip2[2]).real();
			f_m2_p1_d1d2 = (conj(dk0dip[0])*dk0dip2[3]).real();
			f_m2_p2_d1d2 = (conj(dk0dip[0])*dk0dip2[4]).real();

		}

		for (h = 0; h < ww.size(); h++)
		{
			bigJdip = bigJcsa = 0.0;
			bigJdip2 = bigKdip = 0.0;
			
			/***************
			 * K = K' part *
			 ***************/
			
			/* (0, 0) */
			bigJdip  += smallj[0].at(h) * f_0_0_dip;
			bigJcsa  += smallj[0].at(h) * f_0_0_csa;
			bigJdip2 += smallj[0].at(h) * f_0_0_dip2;
			bigKdip  += smallj[0].at(h) * f_0_0_d1d2;
			/* (-1,-1) + (1,1) */
			bigJdip  += smallj[1].at(h) * ( f_m1_m1_dip  + f_p1_p1_dip  );
			bigJcsa  += smallj[1].at(h) * ( f_m1_m1_csa  + f_p1_p1_csa  );
			bigJdip2 += smallj[1].at(h) * ( f_m1_m1_dip2 + f_p1_p1_dip2 );
			bigKdip  += smallj[1].at(h) * ( f_m1_m1_d1d2 + f_p1_p1_d1d2 );
			/* (-2,-2) + (2,2) */
			bigJdip  += smallj[2].at(h) * ( f_m2_m2_dip  + f_p2_p2_dip  );
			bigJcsa  += smallj[2].at(h) * ( f_m2_m2_csa  + f_p2_p2_csa  );
			bigJdip2 += smallj[2].at(h) * ( f_m2_m2_dip2 + f_p2_p2_dip2 );
			bigKdip  += smallj[2].at(h) * ( f_m2_m2_d1d2 + f_p2_p2_d1d2 );

			/****************
			 * K != K' part *
			 ****************/
					
			if (sd.size() > 3)
			{
				/* (-2, 2) */
				bigJdip  += 2.0 * smallj[3].at(h) * f_m2_p2_dip;
				bigJcsa  += 2.0 * smallj[3].at(h) * f_m2_p2_csa;
				bigJdip2 += 2.0 * smallj[3].at(h) * f_m2_p2_dip2;
				bigKdip  += 2.0 * smallj[3].at(h) * f_m2_p2_d1d2;
				/* (-1, 1) */
				bigJdip  += 2.0 * smallj[4].at(h) * f_m1_p1_dip;
				bigJcsa  += 2.0 * smallj[4].at(h) * f_m1_p1_csa;
				bigJdip2 += 2.0 * smallj[4].at(h) * f_m1_p1_dip2;
				bigKdip  += 2.0 * smallj[4].at(h) * f_m1_p1_d1d2;
				/* (-2, 0) + (0, 2) */
				bigJdip  += 2.0 * smallj[5].at(h) * ( f_m2_0_dip  + f_0_p2_dip  );
				bigJcsa  += 2.0 * smallj[5].at(h) * ( f_m2_0_csa  + f_0_p2_csa  );
				bigJdip2 += 2.0 * smallj[5].at(h) * ( f_m2_0_dip2 + f_0_p2_dip2 );
				bigKdip  += 2.0 * smallj[5].at(h) * ( f_m2_0_d1d2 + f_0_p2_d1d2 );
				if (sd.size() > 6)
				{
					/* (-1, 2) + (-2, 1) */
					bigJdip  += 2.0 * smallj[6].at(h) * ( f_m1_p2_dip  + f_m2_p1_dip  );
					bigJcsa  += 2.0 * smallj[6].at(h) * ( f_m1_p2_csa  + f_m2_p1_csa  );
					bigJdip2 += 2.0 * smallj[6].at(h) * ( f_m1_p2_dip2 + f_m2_p1_dip2 );
					bigKdip  += 2.0 * smallj[6].at(h) * ( f_m1_p2_d1d2 + f_m2_p1_d1d2 );
					/* (0, 1) + (-1, 0) */
					bigJdip  += 2.0 * smallj[7].at(h) * ( f_0_p1_dip   + f_m1_0_dip   );
					bigJcsa  += 2.0 * smallj[7].at(h) * ( f_0_p1_csa   + f_m1_0_csa   );
					bigJdip2 += 2.0 * smallj[7].at(h) * ( f_0_p1_dip2  + f_m1_0_dip2  );
					bigKdip  += 2.0 * smallj[7].at(h) * ( f_0_p1_d1d2  + f_m1_0_d1d2  );
					/* (1, 2) + (-2, -1) */
					bigJdip  += 2.0 * smallj[8].at(h) * ( f_p1_p2_dip  + f_m2_m1_dip  );
					bigJcsa  += 2.0 * smallj[8].at(h) * ( f_p1_p2_csa  + f_m2_m1_csa  );
					bigJdip2 += 2.0 * smallj[8].at(h) * ( f_p1_p2_dip2 + f_m2_m1_dip2 );
					bigKdip  += 2.0 * smallj[8].at(h) * ( f_p1_p2_d1d2 + f_m2_m1_d1d2 );
				}
			}
			
			jdipFile << ww.at(h) << " " << bigJdip << std::endl;
			jcsaFile << ww.at(h) << " " << bigJcsa << std::endl;
			if (nH == 2)
			{
				jdip2File     << ww.at(h) << " " << bigJdip2 << std::endl;
				jdip1dip2File << ww.at(h) << " " << bigKdip  << std::endl;
			}
			if (h == 0 && !rank)
				std::cout << std::endl << std::endl << "CORRELATION TIME = " << 1.0e12 * bigJdip / scale << " ps" << std::endl << std::endl;
		}
		
		jdipFile.close();
		jcsaFile.close();
		jdip2File.close();
		jdip1dip2File.close();
	}
	
	lastRun = false;
	oldRun = true;
	
	if (wout)
		outFile.close();
	
	return jw;
}

/**********************************************************
 * Return eigenvalues and weigths of a tridiagonal matrix *
 **********************************************************/

#define SIGN(a,b) ((b)>=0.0 ? fabs((a)) : -fabs((a)))
#define pythag(a,b) sqrt((a)*(a)+(b)*(b))

dvvector lanczos::tqli(unsigned int n, dvector d, dvector e)
{
	
	if (d.size() != n || e.size() != n)
	{
		std::cout << std::endl << std::endl << "ERROR in lancoz:tqli(int n, dvector d, dvector e) : d and e must be of size n. Size(d) = " << d.size() << "; Size(e) = " << e.size()<<"." <<std::endl << std::endl;
		exit(1);
	}
	dvvector zw(2);
	
	int i, iter, l, m, maxIter = 100;
	double temp;
	double b, c, f, g, p, r, s;
	dvector z(n,0.0);
	
	/****************************
	 * Initialize output vector *
	 ****************************/
	
	z.at(0) = 1.0;
	
	/*************************
	 * Special case of n = 1 *
	 *************************/
	
	if (n == 1)
	{
		zw[0] = d;
		zw[1] = z;
		return zw;
	}
	
	/*************************
	 * Loop over eigenvalues *
	 *************************/

	e.at(n-1) = 0.0;

	for (l = 0; l < (int)n; l++)
	{
		iter = 0;
		
		do{
			for (m = l; m < (int)(n-1); m++)
			{
				temp = fabs(d.at(m)) + fabs(d.at(m+1));
				if((fabs(e.at(m)) + temp) == temp) break;
			}
			if (m != l)
			{
				if (iter++ == maxIter)
				{
					std::cout <<std::endl<<std::endl<< "ERROR in lanczos::tqli - max iterations reached"<<std::endl<<std::endl;
					exit(1);
				}
				g = 0.5 * (d.at(l+1)-d.at(l)) / e.at(l);
				r = pythag(g, 1.0);
				g = d.at(m) - d.at(l) + e.at(l)/(g+SIGN(r,g));
				s = c = 1.0;
				p = 0.0;
				for (i = m-1; i >= l; i--)
				{
					f = s*e.at(i);
					b = c*e.at(i);
					e.at(i+1) = (r = pythag(f,g));
					if (fabs(r) <= ZERO)
					{
						d.at(i+1) -= p;
						e.at(m) = 0.0;
						break;
					}
					s = f/r;
					c = g/r;
					g = d.at(i+1)-p;
					r = (d.at(i)-g)*s + 2.0*c*b;
					d.at(i+1) = g + (p=s*r);
					g = c*r - b;
					f = z.at(i+1);
					z.at(i+1) = s*z.at(i) + c*f;
					z.at(i) = c*z.at(i) - s*f;
				}
				if (fabs(r) <= ZERO && i >= l) continue;
				d.at(l) -= p;
				e.at(l) = g;
				e.at(m) = 0.0;
			}
		} while (m != l);
	}
	
	/*****************************************
	 * Insertion sort (decreasing by weight) *
	 *****************************************/
	
	for (l = 0; l < (int)n; l++) z.at(l) = z.at(l)*z.at(l);

	for (l = 0; l < (int)(n-1); l++)
	{
		s = z.at(l);
		for (m = l+1; m < (int)n; m++)
		{
			p = z.at(m);
			if (p > s)
			{
				s = p;
				f = d.at(m);
				d.at(m) = d.at(l);
				d.at(l) = f;
				f = z.at(m);
				z.at(m) = z.at(l);
				z.at(l) = f;
			}
		}
	}
	
	zw[0] = d;
	zw[1] = z;
	
	return zw;
}

/**************************************
 * Calculate the correlation function *
 **************************************/

dvvector lanczos::corrFunc(dvector eigenvalue, dvector weight, int nt, double mint, double maxt)
{
	if (eigenvalue.size() != weight.size())
	{
		std::cout << std::endl << std::endl << "ERROR in lanczos::corrFunc(dvector eigenvalue, dvector weight, int timeLimit) - Size(eigenvalue) = " << eigenvalue.size() << " is different from Size(weight) = " << weight.size() << std::endl << std::endl;
		exit(1);
	}
	dvvector cf(2);
	double dt = (maxt-mint)/(double)(nt-1);
	double t;
	for (int i = 0; i < nt; i++)
	{
		t = mint+dt*(double)i;
		cf[0].push_back(t);
		cf[1].push_back(0.0);
		for (unsigned int j = 0; j < eigenvalue.size(); j++)
			cf[1].back() += weight.at(j)*exp(-eigenvalue.at(j)*t);
	}
	return cf;
}

/**********************************
 * Calculate the spectrla density *
 **********************************/

dvvector lanczos::specDens(dvector eigenvalue, dvector weight, int nw, double minw, double maxw)
{
	if (eigenvalue.size() != weight.size())
	{
		std::cout << std::endl << std::endl << "ERROR in lanczos::corrFunc(dvector eigenvalue, dvector weight, int freqLimit) - Size(eigenvalue) = " << eigenvalue.size() << " is different from Size(weight) = " << weight.size() << std::endl << std::endl;
		exit(1);
	}
	dvvector cf(2);
	double dw = (maxw - minw) / (double)(nw - 1);
	double w;
	double tau;
	for (int i = 0; i < nw; i++)
	{
		w = minw + dw*(double)i;
		cf[0].push_back(w);
		cf[1].push_back(0.0);
		for (unsigned int j = 0; j < eigenvalue.size(); j++)
		{
			tau = 1.0 / eigenvalue.at(j);
			cf[1].back() += weight.at(j)*tau/(1.0+w*w*tau*tau);
		}
	}
	return cf;
}

/*****************************
 * Return a pointer to alpha *
 *****************************/

dvvector* lanczos::getAlphaPtr(void)
{
	dvvector *aptr;
	aptr = &alpha;
	return aptr;
}

/****************************
 * Return a pointer to beta *
 ****************************/

dvvector* lanczos::getBetaPtr(void)
{
	dvvector *bptr;
	bptr = &beta;
	return bptr;
}

/************************
 * Write alpha and beta *
 ************************/

std::string lanczos::toString(bool writeAlpha, bool writeBeta, int ikk1)
{
	std::ostringstream ostr;
	if (writeAlpha && ikk1<50 && ikk1 < (int)alpha.size())
	{
		ostr << "TASK " << rank << ": " << setfill('-') << setw(12) << "" << setfill(' ') << std::endl;
		ostr << "TASK " << rank << ": " << " alpha(" << ikk1 << ")" << std::endl;
		ostr << "TASK " << rank << ": " << setfill('-') << setw(12) << "" << setfill(' ') << std::endl;
		for (unsigned int i = 0; i < alpha.at(ikk1).size(); i++)
			ostr << /*"TASK " << rank << ": " << */scientific << setprecision(4) << showpos << alpha[ikk1][i] << std::endl;
		ostr << "TASK " << rank << ": " << setfill('-') << setw(12) << "" << setfill(' ') << std::endl;
	}
	
	if (writeBeta && ikk1<50 && ikk1 < (int)beta.size())
	{
		ostr << "TASK " << rank << ": " << setfill('-') << setw(12) << "" << setfill(' ') << std::endl;
		ostr << "TASK " << rank << ": " << " beta(" << ikk1 << ")" << std::endl;
		ostr << "TASK " << rank << ": " << setfill('-') << setw(12) << "" << setfill(' ') << std::endl;
		for (unsigned int i = 0; i < beta.at(ikk1).size(); i++)
			ostr << /*"TASK " << rank << ": " << */scientific << setprecision(4) << showpos << beta[ikk1][i] << std::endl;
		ostr << "TASK " << rank << ": " << setfill('-') << setw(12) << "" << setfill(' ') << std::endl;
	}
	return ostr.str();
}

/********************************************/
/* Set scaling factor for correlation times */
/********************************************/

void lanczos::setScale(double s)
{
	scale = s;
	return;
}

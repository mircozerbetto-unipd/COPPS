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

#include "copps.h"

int S_OUT_OF_RANGE;
dvvector T1T2NOE;

int observables_minpack(void *p, int m, int n, const double *par, double *fvec, int iflag)
{

	double *e = (double *)p;
	
	clock_t start, stop;

	potential u, v;
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
		if (expdata.compHasPotential(expdata.getActiveComponent()))
			u = expdata.getPotentialOfComp(expdata.getActiveComponent());
		else
			u = data.getPotentialCoefficients();

		v.npop = 1;
		v.c20 = ldvector(1,0.0);
		v.c22 = ldvector(1,0.0);
		v.c40 = ldvector(1,0.0);
		v.c42 = ldvector(1,0.0);
		v.c44 = ldvector(1,0.0);
	}
	else
		u.npop = 1;

        dvector B0;
        dvvector tmp_T1T2NOE;

	for (int npop = 0; npop < u.npop; npop++)
        {

		if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
		{
	                v.c20.at(0) = u.c20.at(npop);
	                v.c22.at(0) = u.c22.at(npop);
	                v.c40.at(0) = u.c40.at(npop);
	                v.c42.at(0) = u.c42.at(npop);
	                v.c44.at(0) = u.c44.at(npop);
		}

		/************************
		 * Update physical data *
		 ************************/

		updatePhysicalData(par,m,n,p,&v);

		/******************
		 * Update objects *
		 ******************/
		
		updateObjects(v);
		
		/***************************** 
		 * Calculate starting vector *
		 *****************************/

		if (!fitStep || (fitStep && data.isPotentialFit()))
		{
			start = clock();
#ifdef WRITE_ALL
			if (!mpi_rank) std::cout << "TASK 0: Projecting starting vector on basis functions..." << std::endl;
#endif
			Tkk1.projectOnBasis(&basisFunctions,nrows,global_row);
			Tkk1.scatterProjections(true,mpi_ntasks);
#ifdef WRITE_ALL
#ifdef WRITE_STVEC
			if (!mpi_rank)
				std::cout << Tkk1.toString(&basisFunctions);
#endif
#endif
			stop = clock();
		}
			
		/*****************************
		 * Calculate matrix elements *
		 *****************************/
		
		start = clock();
#ifdef WRITE_ALL
			if (!mpi_rank) std::cout << "TASK 0: Building matrix elements..." << std::endl;
#endif
		mat.buildMatrix();
#ifdef WRITE_ALL
#ifdef WRITE_MATRIX
		if (!mpi_rank)
			std::cout << mat.toString(1,0);
#endif
#endif
		stop = clock();
		
		/***************
		 * Run Lanczos *
		 ***************/
		
		start = clock();
		lcz.runLaczos();
#ifdef WRITE_ALL
#ifdef WRITE_LANCZOS
			std::cout << lcz.toString(true,true,0);
			std::cout << lcz.toString(true,true,1);
			std::cout << lcz.toString(true,true,2);
#endif
#endif
		stop = clock();
		
		/******************************
		 * Calculate relaxation times *
		 ******************************/

		start = clock();
		tmp_T1T2NOE = rel.getT1T2NOE(false);
		B0 = rel.getField();
		stop = clock();

		for (unsigned int i = 0; i < B0.size(); i++)
		{
                        /* store 1/T1 */
			tmp_T1T2NOE[i][0] = 1.0 / tmp_T1T2NOE[i][0];
			/* store 1/T2 */
			tmp_T1T2NOE[i][1] = 1.0 / tmp_T1T2NOE[i][1];
			/* store NOE/T1 */
			tmp_T1T2NOE[i][2] *= tmp_T1T2NOE[i][0];
		}

		if (!npop)
			T1T2NOE = tmp_T1T2NOE;
		else
		{
			for (int i = 0; i < tmp_T1T2NOE.size(); i++)
			{
				
				for (int j = 0; j < tmp_T1T2NOE[i].size(); j++)
					T1T2NOE[i][j] += tmp_T1T2NOE[i][j];
			}
		}

	} // cycle over populations

	/**************************************/
	/* Calculate average over populations */
	/**************************************/
	
	double invNpop = 1.0 / (double)u.npop;
	for (unsigned int i = 0; i < B0.size(); i++)
	{
		/* T1 */
		T1T2NOE[i][0] *= invNpop;
		T1T2NOE[i][0] = 1.0 / T1T2NOE[i][0];
		/* T2 */
		T1T2NOE[i][1] *= invNpop;
		T1T2NOE[i][1] = 1.0 / T1T2NOE[i][1];
		/* NOE */
		T1T2NOE[i][2] *= invNpop;
		T1T2NOE[i][2] *= T1T2NOE[i][0];
		if (data.getNHydrogens() == 2)
		{
			/* CCRT1 */
			T1T2NOE[i][3] *= invNpop;
			/* CCRT2 */
			T1T2NOE[i][4] *= invNpop;
		}
	}

	/************************
	 * Calculate chi square *
	 ************************/
	
	double chisq = calculateChiSquare(T1T2NOE,e,fvec);

	/**************************/
	/* Update fitStep counter */
	/**************************/

	fitStep++;

	/****************************/
	/* Ouput fitting parameters */
	/****************************/

	bool errorsConverged = outputData(T1T2NOE,e,chisq,true);

	/*************************/
	/* Check exit conditions */
	/*************************/

	if (chisq < 1.0e-3)
	{
		fitout << std::endl << "IMPORTANT : FITTING HAS BEEN STOPPED BECAUSE REDUCED CHI-SQUARE < 1.0e-3." << std::endl << std::endl;
		iflag = -1;
		return iflag;
	}
	if (nOrderParametersStopChecks == n+1)
	{
		fitout << std::endl << "IMPORTANT : FITTING HAS BEEN STOPPED BECAUSE VARIANCE OF ORDER PARAMETERS IS < " << data.getOrderParametersTolerance() << " AFTER " << (n+1) << " FIT STEPS.\nPLEASE, LOOK AT RESULTS CAREFULLY. IN CASE OF UNACCEPTABLE RESULTS CONSIDER RESTARTING\nTHE SIMULATION WITH A STRICKTER TOLERANCE FOR THE VARIANCE OF ORDER PARAMETERS." << std::endl << std::endl;
		iflag = -1;
		return iflag;
	}
	if (errorsConverged)
	{
		fitout << std::endl << "FITTING HAS BEEN STOPPED BECAUSE ALL THE THEORETICAL ERRORS ARE BELOW THE EXPERIMENTAL ONES." << std::endl << std::endl;
		iflag = -1;
		return iflag;
	}

	/***************/
	/* MPI Barrier */
	/***************/

#ifdef _MPI_	
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	return 0;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void observables_levmar (double *p, double *sim, int n, int m, void *otherData)
{

        double *e = (double *)otherData;
        double *fvec = new double[m];

        clock_t start, stop;

        potential u, v;
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
	        if (expdata.compHasPotential(expdata.getActiveComponent()))
	                u = expdata.getPotentialOfComp(expdata.getActiveComponent());
	        else
	                u = data.getPotentialCoefficients();
	
	        v.npop = 1;
	        v.c20 = ldvector(1,0.0);
	        v.c22 = ldvector(1,0.0);
	        v.c40 = ldvector(1,0.0);
	        v.c42 = ldvector(1,0.0);
	        v.c44 = ldvector(1,0.0);
	}
	else
		u.npop = 1;

        dvector B0;
        dvvector tmp_T1T2NOE;

        for (int npop = 0; npop < u.npop; npop++)
        {
	       
		if(!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
		{
			v.c20.at(0) = u.c20.at(npop);
			v.c22.at(0) = u.c22.at(npop);
			v.c40.at(0) = u.c40.at(npop);
			v.c42.at(0) = u.c42.at(npop);
			v.c44.at(0) = u.c44.at(npop);
		}

		/************************
		 * Update physical data *
		 ************************/

		updatePhysicalData(p,m,n,e,&v);

		/******************
		 * Update objects *
		 ******************/

		updateObjects(v);

	       /***************************** 
		 * Calculate starting vector *
		 *****************************/

		if (!fitStep || (fitStep && data.isPotentialFit()))
		{
			start = clock();
#ifdef WRITE_ALL
			if (!mpi_rank) std::cout << "TASK 0: Projecting starting vector on basis functions..." << std::endl;
#endif
			Tkk1.projectOnBasis(&basisFunctions,nrows,global_row);    
			Tkk1.scatterProjections(true,mpi_ntasks);
#ifdef WRITE_ALL
#ifdef WRITE_STVEC
			if (!mpi_rank)
				std::cout << Tkk1.toString(&basisFunctions);
#endif
#endif
			stop = clock();
		}
			
		/*****************************
		 * Calculate matrix elements *
		 *****************************/
		
		start = clock();
#ifdef WRITE_ALL
			if (!mpi_rank) std::cout << "TASK 0: Building matrix elements..." << std::endl;
#endif
		mat.buildMatrix();
#ifdef WRITE_ALL
#ifdef WRITE_MATRIX
		if (!mpi_rank)
			std::cout << mat.toString(1,0);
#endif
#endif
		stop = clock();
		
		/***************
		 * Run Lanczos *
		 ***************/
		
		start = clock();
		lcz.runLaczos();
#ifdef WRITE_ALL
#ifdef WRITE_LANCZOS
		std::cout << lcz.toString(true,true,0);
		std::cout << lcz.toString(true,true,1);
		std::cout << lcz.toString(true,true,2);
#endif
#endif
		stop = clock();
		
		/******************************
		 * Calculate relaxation times *
		 ******************************/

                start = clock();
                tmp_T1T2NOE = rel.getT1T2NOE(false);
                B0 = rel.getField();
                stop = clock();

		for (unsigned int i = 0; i < B0.size(); i++)
                {
                        /* store 1/T1 */
                        tmp_T1T2NOE[i][0] = 1.0 / tmp_T1T2NOE[i][0];
                        /* store 1/T2 */
                        tmp_T1T2NOE[i][1] = 1.0 / tmp_T1T2NOE[i][1];
                        /* store NOE/T1 */
                        tmp_T1T2NOE[i][2] *= tmp_T1T2NOE[i][0];
                }

                if (!npop)
                        T1T2NOE = tmp_T1T2NOE;
                 else
                 {
                        for (int i = 0; i < tmp_T1T2NOE.size(); i++)
                        {
                                for (int j = 0; j < tmp_T1T2NOE[i].size(); j++)
                                        T1T2NOE[i][j] += tmp_T1T2NOE[i][j];
                        }
                 }

        } // cycle over populations

        /**************************************/
        /* Calculate average over populations */
        /**************************************/

        double invNpop = 1.0 / (double)u.npop;
	for (unsigned int i = 0; i < B0.size(); i++)
        {
                /* T1 */
                T1T2NOE[i][0] *= invNpop;
                T1T2NOE[i][0] = 1.0 / T1T2NOE[i][0];
                /* T2 */
                T1T2NOE[i][1] *= invNpop;
                T1T2NOE[i][1] = 1.0 / T1T2NOE[i][1];
                /* NOE */
                T1T2NOE[i][2] *= invNpop;
                T1T2NOE[i][2] *= T1T2NOE[i][0];
		if (data.getNHydrogens() == 2)
		{
			/* CCRT1 */
			T1T2NOE[i][3] *= invNpop;
			/* CCRT2 */
			T1T2NOE[i][4] *= invNpop;
		}
        }

        /************************
         * Calculate chi square *
         ************************/
	
	double chisq = calculateChiSquare(T1T2NOE,e,fvec);
 
	/********************/
	/* Update sim array */
	/********************/

        for (unsigned int i = 0; i < B0.size(); i++)
	{
		sim[i*3+0]  = T1T2NOE[i][0];
		sim[i*3+1]  = T1T2NOE[i][1];
		sim[i*3+2]  = T1T2NOE[i][2];
	}

	/**************************/
	/* Update fitStep counter */
	/**************************/
               
        fitStep++;

	/*******************************/
	/* Output fitting pararameters */
	/*******************************/
        
	outputData(T1T2NOE,e,chisq,true);

	/***************/
	/* MPI Barrier */
	/***************/

#ifdef _MPI_    
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        for (int ii = 0; ii < n; ii++)
		std::cout << p[ii] << std::endl;
	std::cout << "===================" << std::endl;

        return;
}

/*************************************/
/* Function called by POWELL routine */
/*************************************/

extern "C"{
	double func_(double *p)
	{
		int powell_iflag = 0;
	        double *powell_fvec = new double[powell_m];
		int om = observables_minpack(powell_otherData, powell_m, powell_n, p, powell_fvec, powell_iflag);
		double fret = 0.0;
		if (S_OUT_OF_RANGE) fret = 1.0e12;
		else if (om == -1) fret = 0.0;
		else fret = calculateChiSquare(T1T2NOE, powell_otherData, powell_fvec);
#ifdef _MPI_	
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		return fret;
	}
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

int updatePhysicalData(const double *par, int m, int n, void *otherData, potential *v)
{
        int nfit     = (!(fn % 2) ? data.getNFit()   : data.getNFit2());
        bool *fitpar = (!(fn % 2) ? data.getFitPar() : data.getFitPar2());
        int i = 0;
        int iDZ1=0, iDZ2=0;

	S_OUT_OF_RANGE = 0;

	int nanFound = 0;
	for (i = 0; i < nfit; i++)
	{
		if ( numeric_limits<double>::has_quiet_NaN && par[i] == numeric_limits<double>::quiet_NaN() )
			nanFound ++;
	}
	if (nanFound)
	{
		S_OUT_OF_RANGE = 1;
		return 1;
	}

	// If potential coefficients have been specified in the .exp file fitting
	// if automatically disabled
	if (expdata.compHasPotential(expdata.getActiveComponent()))
	{
		fitpar[18] = false;
		fitpar[19] = false;
		fitpar[20] = false;
		fitpar[21] = false;
		fitpar[22] = false;
	}
	
	// Process fitpar array to understand which parameters are going to be fit
	i = 0;
        // Protein diffusion
        if (fitpar[0])  {if(!cdz1) data.setProteinDyy(exp((long double)par[i])); else iDZ1=i; ++i;}
        if (fitpar[1])  {data.setProteinDxx((long double)(par[i])*data.getProteinDyy()); ++i;}
        if (fitpar[2])  {data.setProteinDzz((long double)(par[i])*data.getProteinDyy()); ++i;}
        // Omega V
        if (fitpar[6])  {data.setValueOf("protein_alpha", (long double)par[i]); ++i;}
        if (fitpar[7])  {data.setValueOf("protein_beta",  (long double)par[i]); ++i;}
        if (fitpar[8])  {data.setValueOf("protein_gamma", (long double)par[i]); ++i;}
        // Omega O
        if (fitpar[9])  {data.setValueOf("probe_alpha", (long double)par[i]);   ++i;}
        if (fitpar[10]) {data.setValueOf("probe_beta",  (long double)par[i]);   ++i;}
        if (fitpar[11]) {data.setValueOf("probe_gamma", (long double)par[i]);   ++i;}
        // Omega D
        if (fitpar[12]) {data.setValueOf("dipolar_alpha", (long double)par[i]); ++i;}
        if (fitpar[13]) {data.setValueOf("dipolar_beta",  (long double)par[i]); ++i;}
        if (fitpar[14]) {data.setValueOf("dipolar_gamma", (long double)par[i]); ++i;}
        // Omega D2
        if (fitpar[26]) {data.setValueOf("dipolar_alpha2", (long double)par[i]); ++i;}
        if (fitpar[27]) {data.setValueOf("dipolar_beta2",  (long double)par[i]); ++i;}
/////////////////////////////////////////
// ONLY TEMPORARY
        if (fitpar[27] &&  !(dynModel.compare("ts-srls"))) data.setDB2((long double)par[i-1]);
/////////////////////////////////////////
        if (fitpar[28]) {data.setValueOf("dipolar_gamma2", (long double)par[i]); ++i;}
        // Probe diffusion
        if (fitpar[3])  {if(!cdz2) data.setProbeDyy(exp((long double)par[i])); else iDZ2=i; ++i;}
        if (fitpar[4])  {data.setProbeDxx((long double)(par[i])*data.getProbeDyy()); ++i;}
        if (fitpar[5])  {data.setProbeDzz((long double)(par[i])*data.getProbeDyy()); ++i;}
        // Omega C
        if (fitpar[15]) {data.setValueOf("csa_alpha", (long double)par[i]);     ++i;}
        if (fitpar[16]) {data.setValueOf("csa_beta",  (long double)par[i]);     ++i;}
        if (fitpar[17]) {data.setValueOf("csa_gamma", (long double)par[i]);     ++i;}
        // Poteintial coefficients
////	if (fitpar[18] && fitpar[19]) // For rhombic potential S20 and S22 are fitted in place of c20 and c22 
////	{
////		double *Spar = new double[2];
////		Spar[0] = par[i]; ++i;
////		Spar[1] = par[i]; ++i;
////		if ( fabs(Spar[0]) < STOL && fabs(Spar[1]) < STOL)
////		{
////			data.setValueOf("c20",0.0);
////			data.setValueOf("c22",0.0);
////		}
////		else if (fabs(Spar[1])>=SQRT_THREE_OVER_TWO || Spar[0] <= -0.5 || Spar[0] >= 1.0)
////			S_OUT_OF_RANGE = 1;
////		else
////		{
////			int nS = 2;
////			int infoS;
////			int lwaS = 20;
////			double tolS = 0.01;
////			double *xS = new double[2];
////			xS[0] = 0.0;;
////			xS[1] = 0.0;;
////			double foutS[2];
////			double waS[lwaS];
////			int hS = hybrd1 (StoC, (void *)Spar, nS, xS, foutS, tolS, waS, lwaS);
////			if ( (std::numeric_limits<double>::has_quiet_NaN && xS[0] == std::numeric_limits<double>::quiet_NaN()) ||
////			     (std::numeric_limits<double>::has_quiet_NaN && xS[1] == std::numeric_limits<double>::quiet_NaN())    )
////			{
////				data.setValueOf("c20",31.0);
////				data.setValueOf("c22",31.0);
////				if (!mpi_rank) std::cout << "NaN detected in S -> c conversion" << std::endl;
////				S_OUT_OF_RANGE = 1;
////			}
////			else if (fabs(xS[0]) <= 30 && fabs(xS[1]) <= 30)
////			{
////				data.setValueOf("c20",(long double)xS[0]);
////				data.setValueOf("c22",(long double)xS[1]);
////			}
////			else
////				S_OUT_OF_RANGE = 1;
////		}
////	}
////	else
////	{
		if (fitpar[18]) {data.setValueOf("c20", (long double)par[i]); ++i;}
		if (fitpar[19])
		{
			if (ratio22)
				data.setValueOf("c22", (long double)(par[i]*(double)data.getValueOf("c20")));
			else
				data.setValueOf("c22", (long double)(par[i]));
			++i;
		}
////	}
        if (fitpar[20]) {data.setValueOf("c40", (long double)par[i]);           ++i;}
        if (fitpar[21])
	{
		if (ratio42)
			data.setValueOf("c42", (long double)(par[i]*(double)data.getValueOf("c40")));
		else
			data.setValueOf("c42", (long double)(par[i]));
		++i;
	}
        if (fitpar[22])
	{
		if (ratio44)
			data.setValueOf("c44", (long double)(par[i]*(double)data.getValueOf("c40")));
		else
			data.setValueOf("c44", (long double)(par[i]));
		++i;
	}
        // Conformational exchange rate
        if (fitpar[23]) {data.setValueOf("Rexchange", (long double)exp(par[i])+1.0e-3); ++i;}
	if (fitpar[25])
	{
		if (!(dynModel.compare("fb1")))
			data.scaleFB1diften(data.getScale()/(long double)(par[i]*par[i]));
		else if (!(dynModel.compare("fb2")))
			data.scaleFB2diften(data.getScale()/(long double)(par[i]*par[i]));
		++i;
	}
	// Population in ts-srls model
	if (fitpar[29]) {data.setPopulation(cos(par[i]) * cos(par[i])); ++i;} // population = cos^2(par)
	// Jump frequency in ts-srls model
	if (fitpar[30]) {data.setJumpFrequency(exp(par[i])); ++i;} // jump freq = exp(par)
	// hch sigma
	if (fitpar[31]) {data.setHchSigma(exp(par[i])); ++i;}

        /*************************************
         * Check constrains among parameters *
         *************************************/

	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
		// Protein's D constraints 
		
		if (cdz1)
		{
			long double DX_OVER_DY = data.getProteinDxx()/data.getProteinDyy();
			long double DZ_OVER_DY = data.getProteinDzz()/data.getProteinDyy();
			long double DY = (long double)exp(par[iDZ1]);
			data.setProteinDxx(DX_OVER_DY*DY);
			data.setProteinDzz(DZ_OVER_DY*DY);
			data.setProteinDyy(DY);
		}
		else if (cdy1)
			data.setProteinDxx(data.getProteinDyy());
		
		// Probe's D constraints 
		
		if (cdz2)
		{
			long double DX_OVER_DY = data.getProbeDxx()/data.getProbeDyy();
			long double DZ_OVER_DY = data.getProbeDzz()/data.getProbeDyy();
			long double DY = (long double)exp(par[iDZ2]);
			data.setProbeDxx(DX_OVER_DY*DY);
			data.setProbeDzz(DZ_OVER_DY*DY);
			data.setProbeDyy(DY);
		}
		else if (cdy2)
			data.setProbeDxx(data.getProbeDyy());

		// Update OmegaD

		if ((fitpar[6] || fitpar[7] || fitpar[8]) & data.getConstrainDV())
		{
			dvector o1 = dvector(3,0.0);
			o1[0] = (double)data.getValueOf("protein_alpha");
			o1[1] = (double)data.getValueOf("protein_beta");
			o1[2] = (double)data.getValueOf("protein_gamma");
			eul.setOmega1(o1);
			eul.omega1ToE1();
			eul.calculateE12();
			eul.E12ToOmega12();
			o1 = eul.getOmega12();
			data.setValueOf("dipolar_alpha",(long double)o1[0]);
			data.setValueOf("dipolar_beta", (long double)o1[1]);
			data.setValueOf("dipolar_gamma",(long double)o1[2]);
		}

		// Update OmegaV

		if ((fitpar[12] || fitpar[13] || fitpar[14]) & data.getConstrainVD())
		{
			dvector o12 = dvector(3,0.0);
			o12[0] = (double)data.getValueOf("dipolar_alpha");
			o12[1] = (double)data.getValueOf("dipolar_beta");
			o12[2] = (double)data.getValueOf("dipolar_gamma");
			eul.setOmega12(o12);
			eul.omega12ToE12();
			eul.calculateE1();
			eul.E1ToOmega1();
			dvector o1 = eul.getOmega1();
			data.setValueOf("protein_alpha",(long double)o1[0]);
			data.setValueOf("protein_beta", (long double)o1[1]);
			data.setValueOf("protein_gamma",(long double)o1[2]);
		}
		
		// Update potential vector
		if (!expdata.compHasPotential(expdata.getActiveComponent()))
		{
			v->c20.at(0) = data.getValueOf("c20");
			v->c22.at(0) = data.getValueOf("c22");
			v->c40.at(0) = data.getValueOf("c40");
			v->c42.at(0) = data.getValueOf("c42");
			v->c44.at(0) = data.getValueOf("c44");
		}
	}

	return 0;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void updateObjects(potential v)
{
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
	        data.transformProteinD();
	        data.transformProbeD();
		data.setPotentialCoefficients(v);
	}
        
	if (!fitStep || (fitStep && data.isPotentialFit()))
	        Tkk1.update();
        mat.update();
        lcz.update();
	if (!(dynModel.compare("ts-srls")))
	{
		dvector od_1  = data.getOmegaDipolarState(1);
		dvector od_2  = data.getOmegaDipolarState(2);
		dvector oc_1  = data.getOmegaCsaState(1);
		dvector oc_2  = data.getOmegaCsaState(2);
		dvector od2_1 = data.getOmegaDipolar2State(1);
		dvector od2_2 = data.getOmegaDipolar2State(2);
		rel.update_tssrls(od_1,oc_1,od2_1,od_2,oc_2,od2_2,data.getPopulation(),data.getJumpFrequency(),rel.getField());
		rel.setHchSigma(data.getHchSigma());
	}
	else
        	rel.update(data.getOmegaDipolar(), data.getOmegaDipolar2(), data.getOmegaCSA(), rel.getField(), data.getPopulation(), data.getJumpFrequency());

        rel.setRex(data.getValueOf("Rexchange"));

	return; 
}

/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

double calculateChiSquare(dvvector T1T2NOE, double *e, double *fvec)
{
        double chisq = 0.0;
	dvector B0 = rel.getField();

	/********************************************/
	/* Rectification based on Diffusion tensors */
	/********************************************/

	double Drect = 0.0;

	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{

		double D2;
		double D2x = (double)data.getProbeDxx();
		double D2y = (double)data.getProbeDyy();
		double D2z = (double)data.getProbeDzz();
		double minD = data.getDmin();
		double K= data.getKrect();
		if (D2x < minD || D2y < minD || D2z < minD)
		{
			D2 = (D2x < D2y ? (D2x < D2z ? D2x : D2z) : (D2y < D2z ? D2y : D2z));
			Drect = K*(D2-minD)*(D2-minD);
			std::cout << "Chisq has been rectified with K = " << K << std::endl;
		}

//		if (D2z < D2y)
//		{
//			Drect += K*(D2z-D2y)*(D2z-D2y);
//			std::cout << "Chisq has been rectified" << std::endl;
//		}

	}

	if (!(dynModel.compare("ts-srls")))
	{
		if (fabs(data.getPopulation()) < ZERO || fabs(data.getPopulation()-1.0) < ZERO)
			Drect += 1.0e6;
	}

	int nd = expdata.getNExpData();
	int nd2 = 2 * nd;

        for (unsigned int i = 0; i < B0.size(); i++)
        {
		for (int j = 0; j < nd; j++)
		{
	                fvec[i*nd+j] = (T1T2NOE[i][j] - e[i*nd2+2*j]) / e[i*nd2+2*j+1] + Drect;
		        chisq += fvec[i*nd+j]*fvec[i*nd+j];
		}
	}
	
	/******************************************************************************/
	chisq /= dof;
	/******************************************************************************/

	return chisq;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

bool outputData(dvvector T1T2NOE, double *e, double chisq, bool isFitting)
{
	dvector B0 = rel.getField();
        //int nfit = data.getNFit();
        bool *fitpar = (!(fn % 2) ? data.getFitPar() : data.getFitPar2());
	bool errorsConverged = true;

        fitout <<"#######################################################" << std::endl;
        fitout << "               Component " << comp << " - Fit Step " << fitStep << std::endl;
        fitout <<"#######################################################" << std::endl<<std::endl;
        fitout << "* Parameters:" << std::endl<<std::endl;
        int i = 0;
        if (fitpar[1] || cdy1 || cdz1)  fitout << scientific << "  Protein Dxx = " << data.getValueOf("protein_dxx") * data.getScale() << " Hz" << std::endl;
        if (fitpar[0] || cdy1 || cdz1)  fitout << scientific << "  Protein Dyy = " << data.getValueOf("protein_dyy") * data.getScale() << " Hz" << std::endl;
        if (fitpar[2] || cdz1)  fitout << scientific << "  Protein Dzz = " << data.getValueOf("protein_dzz") * data.getScale() << " Hz" << std::endl;
        if (fitpar[4] || cdy2 || cdz2)  fitout << scientific << "  Probe Dxx = "   << data.getValueOf("probe_dxx")   * data.getScale() << " Hz" << std::endl;
        if (fitpar[3] || cdy2 || cdz2)  fitout << scientific << "  Probe Dyy = "   << data.getValueOf("probe_dyy")   * data.getScale() << " Hz" << std::endl;
        if (fitpar[5] || cdz2)  fitout << scientific << "  Probe Dzz = "   << data.getValueOf("probe_dzz")   * data.getScale() << " Hz" << std::endl;
	if (fitpar[6])  fitout << fixed      << "  Alpha V = " << data.getValueOf("protein_alpha") * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[7])  fitout << fixed      << "  Beta  V = " << data.getValueOf("protein_beta")  * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[8])  fitout << fixed      << "  Gamma V = " << data.getValueOf("protein_gamma") * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[9])  fitout << fixed      << "  Alpha O = " << data.getValueOf("probe_alpha")   * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[10]) fitout << fixed      << "  Beta  O = " << data.getValueOf("probe_beta")    * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[11]) fitout << fixed      << "  Gamma O = " << data.getValueOf("probe_gamma")   * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[12]) fitout << fixed      << "  Alpha D = " << data.getValueOf("dipolar_alpha") * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[13]) fitout << fixed      << "  Beta  D = " << data.getValueOf("dipolar_beta")  * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[14]) fitout << fixed      << "  Gamma D = " << data.getValueOf("dipolar_gamma") * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[15]) fitout << fixed      << "  Alpha C = " << data.getValueOf("csa_alpha")     * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[16]) fitout << fixed      << "  Beta  C = " << data.getValueOf("csa_beta")      * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[17]) fitout << fixed      << "  Gamma C = " << data.getValueOf("csa_gamma")     * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[26]) fitout << fixed      << "  Alpha D2 = " << data.getValueOf("dipolar_alpha2") * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[27]) fitout << fixed      << "  Beta  D2 = " << data.getValueOf("dipolar_beta2")  * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[28]) fitout << fixed      << "  Gamma D2 = " << data.getValueOf("dipolar_gamma2") * RAD_TO_DEG << " deg" << std::endl;
	if (fitpar[18]) fitout << fixed      << "  c20 = " << -data.getValueOf("c20") << " kT" << std::endl;
	if (fitpar[19]) fitout << fixed      << "  c22 = " << -data.getValueOf("c22") << " kT" << std::endl;
	if (fitpar[20]) fitout << fixed      << "  c40 = " << -data.getValueOf("c40") << " kT" << std::endl;
	if (fitpar[21]) fitout << fixed      << "  c42 = " << -data.getValueOf("c42") << " kT" << std::endl;
	if (fitpar[22]) fitout << fixed      << "  c44 = " << -data.getValueOf("c44") << " kT" << std::endl;
	if (fitpar[23]) fitout << fixed      << "  R exchange = " << data.getValueOf("Rexchange") << " Hz" << std::endl;
	if (fitpar[25])
	{
		if (!(dynModel.compare("fb1")))
		{
			dvector d = data.getFB1Diften();
			fitout << fixed << scientific << "dxx = " << d[0]*data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dyy = " << d[1]*data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dzz = " << d[2]*data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "d11 = " << d[3]*data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dx1 = " << d[4]*data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dy1 = " << d[5]*data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dz1 = " << d[6]*data.getScale() << " Hz" << std::endl;
		}
		else if (!(dynModel.compare("fb2")))
		{
			dvector d = data.getFB2Diften();
			fitout << fixed << scientific << "dxx = " << d[0] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dyy = " << d[1] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dzz = " << d[2] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "d11 = " << d[3] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "d22 = " << d[4] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "d12 = " << d[5] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dx1 = " << d[6] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dy1 = " << d[7] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dz1 = " << d[8] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dx2 = " << d[9] *data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dy2 = " << d[10]*data.getScale() << " Hz" << std::endl;
			fitout << fixed << scientific << "dz2 = " << d[11]*data.getScale() << " Hz" << std::endl;
		}
	}
	if (fitpar[29]) fitout << fixed << "  State 1 population = " << data.getPopulation() << std::endl;
        if (fitpar[30]) fitout << fixed << scientific << "  Jump frequency = " << data.getJumpFrequency() * data.getScale() << " Hz" << std::endl;
        if (fitpar[31]) fitout << fixed << scientific << "  HCH sigma = " << data.getHchSigma() << " deg" << std::endl;


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
	if (data.getConstrainDV()) //(fitpar[6] || fitpar[7] || fitpar[8])
	{
		fitout << fixed      << "  Alpha D = " << data.getValueOf("dipolar_alpha") * RAD_TO_DEG << " deg" << std::endl;
		fitout << fixed      << "  Beta  D = " << data.getValueOf("dipolar_beta")  * RAD_TO_DEG << " deg" << std::endl;
		fitout << fixed      << "  Gamma D = " << data.getValueOf("dipolar_gamma") * RAD_TO_DEG << " deg" << std::endl;
	}
	if (data.getConstrainVD()) //(fitpar[12] || fitpar[13] || fitpar[14])
	{
		fitout << fixed      << "  Alpha V = " << data.getValueOf("protein_alpha") * RAD_TO_DEG << " deg" << std::endl;
		fitout << fixed      << "  Beta  D = " << data.getValueOf("protein_beta")  * RAD_TO_DEG << " deg" << std::endl;
		fitout << fixed      << "  Gamma D = " << data.getValueOf("protein_gamma") * RAD_TO_DEG << " deg" << std::endl;
	}
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
        fitout << std::endl<< "* Relaxation data:" << std::endl;
        ////int w1 = 13;
        int w2 = 11;
        int w3 = 2;

        fitout << std::endl;
        fitout << "               field/MHz      theo        exp     diff(%)"<< std::endl;
        fitout << "  __________________________________________________________" << std::endl<<std::endl;

	int nd = 2 * expdata.getNExpData();
        double diff = 0.0;
	errorsConverged = true;

        for (i = 0; i < (int)B0.size(); i++)
        {
		/* T1 */
		if (expdata.hasExpValueOf("T1"))
		{
		        diff = (T1T2NOE[i][0] - e[i*nd+0])*100.0 / e[i*nd+0];
			errorsConverged = errorsConverged & (fabs(diff) < 100.0*e[i*nd+1]/e[i*nd+0]);
			fitout << "    T1/ms:      " << fixed << setprecision(3) << setw(w2-2*w3) << B0[i]*1.0e-6 \
				<< fixed << setprecision(3) << setw(w2) << T1T2NOE[i][0]*1.0e3 << fixed << setprecision(3) << setw(w2) << e[i*nd+0]*1.0e3 << fixed << setprecision(3) << setw(w2) << diff << std::endl;
		}
		/* T2 */
		if (expdata.hasExpValueOf("T2"))
		{
		        diff = (T1T2NOE[i][1] - e[i*nd+2])*100.0 / e[i*nd+2];
			errorsConverged = errorsConverged & (fabs(diff) < 100.0*e[i*nd+3]/e[i*nd+2]); 
			fitout << "    T2/ms:      " << fixed << setprecision(3) << setw(w2-2*w3) << B0[i]*1.0e-6 \
				<< fixed << setprecision(3) << setw(w2) << T1T2NOE[i][1]*1.0e3 << fixed << setprecision(3) << setw(w2) << e[i*nd+2]*1.0e3 << fixed << setprecision(3) << setw(w2) << diff << std::endl;
		}
		/* NOE */
		if (expdata.hasExpValueOf("NOE"))
		{
		        diff = (T1T2NOE[i][2] - e[i*nd+4])*100.0 / e[i*nd+4];
			errorsConverged = errorsConverged & (fabs(diff) < 100.0*e[i*nd+5]/e[i*nd+4]); 
			fitout << "    NOE:        " << fixed << setprecision(3) << setw(w2-2*w3) << B0[i]*1.0e-6 \
				<< fixed << setprecision(3) << setw(w2) << T1T2NOE[i][2] << fixed << setprecision(3) << setw(w2) << e[i*nd+4] << fixed << setprecision(3) << setw(w2) << diff << std::endl;
		}
		/* CCRRT1 */
		if (data.getNHydrogens() == 2 && expdata.hasExpValueOf("CCRRT1"))
		{
		        diff = (T1T2NOE[i][3] - e[i*nd+6])*100.0 / e[i*nd+6];
			errorsConverged = errorsConverged & (fabs(diff) < 100.0*e[i*nd+7]/e[i*nd+6]); 
			fitout << "    CCRRT1/Hz:  " << fixed << setprecision(3) << setw(w2-2*w3) << B0[i]*1.0e-6 \
				<< fixed << setprecision(3) << setw(w2) << T1T2NOE[i][3] << fixed << setprecision(3) << setw(w2) << e[i*nd+6] << fixed << setprecision(3) << setw(w2) << diff << std::endl;
		}
		/* CCRRT2 */
		if (expdata.hasExpValueOf("CCRRT2"))
		{
		        diff = (T1T2NOE[i][4] - e[i*nd+8])*100.0 / e[i*nd+8];
			errorsConverged = errorsConverged & (fabs(diff) < 100.0*e[i*nd+9]/e[i*nd+8]); 
			fitout << "    CCRRT2/Hz:  " << fixed << setprecision(3) << setw(w2-2*w3) << B0[i]*1.0e-6 \
				<< fixed << setprecision(3) << setw(w2) << T1T2NOE[i][4] << fixed << setprecision(3) << setw(w2) << e[i*nd+8] << fixed << setprecision(3) << setw(w2) << diff << std::endl;
		}
		
		fitout << "  __________________________________________________________" << std::endl<<std::endl;
	}
	if (isFitting)
        	fitout << std::endl<< "* Reduced Chi-square = " << chisq << std::endl;
	else
        	fitout << std::endl<< "* Chi-square (dof = 1) = " << chisq << std::endl;

	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
	{
		dvector orderParameters = Tkk1.calculateOrdersParams();
		if (fitpar[18] || fitpar[19] || fitpar[20] || fitpar[21] || fitpar[22])
			checkOrderParametersConvergence(orderParameters[0],orderParameters[1],data.getOrderParametersTolerance());

		// Order parameters are not echoed if potential is not being fit
		if (!fitpar[18] && !fitpar[19] && !fitpar[20] && !fitpar[21] && !fitpar[22])
		{
			fitout << std::endl;
			fitout << "Order parameters:"<< std::endl <<std::endl;
			fitout << fixed << setprecision(4) << "S(2,0) = " << orderParameters[0] << std::endl;
			fitout << fixed << setprecision(4) << "S(2,2) = " << orderParameters[1] << std::endl;
			fitout << fixed << setprecision(4) << "Sxx    = " << orderParameters[2] << std::endl;
			fitout << fixed << setprecision(4) << "Syy    = " << orderParameters[3] << std::endl;
			fitout << fixed << setprecision(4) << "Szz    = " << orderParameters[4] << std::endl;
			fitout << std::endl;
		}
	}
        fitout << std::endl<< "* Elapsed time = " << getcputime() << "s" << std::endl<<std::endl;
        
        if (!mpi_rank)
        {
                fstream file_op((path+"/"+project+"_copps.output").c_str(),ios::out);
                file_op << fitout.str();
                file_op.close();
		if (!isFitting)
			std::cout << fitout.str();
        }

	return errorsConverged;
}

void checkOrderParametersConvergence(double s20, double s22, double tol)
{
	double check20;
	double check22 = 0.0;
	
	if (nOrderParametersStopChecks == 0)
	{
		oldS20 = s20;
		oldS22 = s22;
	}

	check20 = 100.0 * fabs( (s20 - oldS20) / oldS20 );
	if (fabs(data.getValueOf("c22")) > 0.0)
		check22 = 100.0 * fabs( (s22 - oldS22) / oldS22 );
	double maxErr = MAX(check20, check22);

	if (maxErr < tol)
	{
		nOrderParametersStopChecks ++;
		double nc = (double)nOrderParametersStopChecks;
		oldS20 = (nc*oldS20 + s20) / (nc+1.0);
		oldS22 = (nc*oldS22 + s22) / (nc+1.0);
	}
	else
	{
		nOrderParametersStopChecks = 0;
		oldS20 = s20;
		oldS22 = s22;
	}

	return;
}

int StoC(void *p, int n, const double *x, double *fout, int iflag)
{
	double *S = (double *)p;
	data.setValueOf("c20",(long double)x[0]);
	data.setValueOf("c22",(long double)x[1]);
	Tkk1.updatePotential();
	dvector op = Tkk1.calculateOrdersParams();
	fout[0] = S[0]-op[0];
	fout[1] = S[1]-op[1];
	return 0;
}


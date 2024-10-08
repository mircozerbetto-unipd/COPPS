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
 Name        : physics.h
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Header of physics class
 ============================================================================
 */
#ifndef PHYSICS_H_
#define PHYSICS_H_

#include <cstdlib>
#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <math.h>

#include "constants.h"
#include "tensor.h"
#include "tensor4.h"
#include "tensor5.h"
#include "euler.h"
#include "prep.h"

#ifdef _MPI_
#include <mpi.h>
#endif

#include "types.h"

// TOTAL NUMBER OF POSSIBLE FITTING PARAMETERS
#define NPARS 32

extern int OUT_CONTROL;

class physics
{
	public:
		
		physics();
		virtual ~physics();
		void init(void);
		void init(int);
		
		// Methods for input
		void readInputFile(std::string);
		void setValueOf(std::string,long double);
		void setPotentialCoefficients(potential);
		
		// Methods for output
		std::string toString(void);
		long double getValueOf(std::string);
		potential getPotentialCoefficients(void);
		bool getPotentialFromExpFile(void);

		
		// MPI mehods
		void setRank(int);
		int getRank(void);
		
		// Specific methods		
		// Diffusion tensor methods
		long double getProteinDxx(void);
		void setProteinDxx(long double);
		long double getProteinDyy(void);
		void setProteinDyy(long double);
		long double getProteinDzz(void);
		void setProteinDzz(long double);
		void scaleProteinD(long double);
		void transformProteinD(void);
		ldcomplex getProteinDlm(int, int);
		long double getProbeDxx(void);
		void setProbeDxx(long double);
		long double getProbeDyy(void);
		void setProbeDyy(long double);
		long double getProbeDzz(void);
		void setProbeDzz(long double);
		void scaleProbeD(long double);
		void transformProbeD(void);
		ldcomplex getProbeDlm(int, int);
		cldvector getProteinDSph(void);
		cldvector getProbeDSph(void);
		void scaleFB1diften(long double);
		void scaleFB2diften(long double);
		int getNCoeff(void);
		ldvector getCoeff(void);
		ldvector getOmegaDipolar(void);
		ldvector getOmegaDipolar2(void);
		ldvector getOmegaCSA(void);
		long double getScale(void);
		dvector getField(void);
		int getLanczosNStep(void);

		// Write out calculation parameters
		int getProbeLMax(void);
		bool getFitting(void);
		int getNFit(void);
		bool* getFitPar(void);
		int getNFit2(void);
		bool* getFitPar2(void);
		std::string getFittingMethod(void);
		
		ivector getSpectralDensities(void);

		dvector getNuclearData(void);
		
		std::string getDynamicsModel(void);

		// TS-SRLS related methods
		dvector getOmegaDipolarState(int);
		dvector getOmegaCsaState(int);
		dvector getOmegaDipolar2State(int);
		double getPopulation(void);
		void setPopulation(double);
		double getJumpFrequency(void);
		void setJumpFrequency(double);
//////////////////////////////////
// ONLY TEMPORARY
		void setDB2(double);
		double getHchSigma(void);
		void setHchSigma(double);
//////////////////////////////////
	
		// FB1 related methods
		int getNmax(void);
		int getNmax(int);
		int getFB1NCoeff(void);
		dcvector getFB1Coeff(void);
		dvector getFB1Diften(void);

		// FB2 related methods
		potential_fb2 getFB2Coeff();
		int getFB2CoeffDim(int);
		dcomplex getFB2CoeffAt(int, int);
		dvector getFB2Diften(void);
		void getCubatureParams(int *, double *);

		// 3S-FB related methods
		dvector _3SFB_getOmegaDipolar1Site(int);
		dvector _3SFB_getOmegaDipolar2Site(int);
		dvector _3SFB_getOmegaCSASite(int);
		double  _3SFB_getPopulation(int);
		double  _3SFB_getSqPopulation(int);
		dvector _3SFB_getSqPopulations(void);
		void    _3SFB_setPopulation(int, double);
		double  _3SFB_getJumpFrequency(int);
		dvector _3SFB_getJumpFrequencies(void);
		void    _3SFB_setJumpFrequency(int, double);

		// Other methods
		bool fitProteinDiffWithConstraints;
		bool fitProbeDiffWithConstraints;

		bool isPotentialFit();

		long double getBkpDipolarAngle(std::string);
	
		bool getConstrainVD();
		bool getConstrainDV();

		double getDmin();
		double getKrect();
		
		double getOrderParametersTolerance();
		void setOrderParametersTolerance(double);

		int getNHydrogens(void);

		bool readSmallJs;
		int getNw(void);

	private:
		bool isFB1;		

		bool potentialFromExpFile;

		int rank; // MPI rank

		int nstep;
		int nfit, nfit2;
		bool fitting;
		bool potentialFit;
		bool fitpar[NPARS], fitpar2[NPARS];
		bool constrain_OmegaV_OmegaD; // when fitting OmegaD
		bool constrain_OmegaD_OmegaV; // when fitting OmegaV
		double fitTol1, fitTol2; // Tolerance. Introduced as external param.
		double order_parameters_tolerance;
		double Dmin, Krect; // D rectification parameters
		std::string fitMethod;

		std::string dynModel; // string indicating the dynamicl model

		int probeLmax;    // truncation parameter for probe rotational space in SRLS
		int Nmax;         // truncation parameter for internal rotation in FB1
		int N1max, N2max; // truncation parameter for internal rotation in FB2
		int t1Max, t2Max; // truncation parameters of the internal potential in the FB2 model

		ivector spectralDensities; // Indexes numbering spectral densities jkk' to calculate

		dvector field; // Contains the frequencies at which caltulating NMR parameters

		long double temperature, viscosity, effective_radius;
		tensor proteinD, probeD; // Diffusion tensors in SRLS model
		tensor4 fb1_diften; // Diffusion tensor in FB1 model [dxx, dyy, dzz, d11, dx1, dy1, dz1]
		tensor5 fb2_diften; // Diffusion tensor in FB2 model [dxx, dyy, dzz, d11, d22, d12, dx1, dy1, dz1, dx2, dy2, dz2]

		int nHydrogens; // Number of hydrogen atoms attached to nucleus
		double gyromag, bondLength, deltaCSA; // Nuclear values
                long double Rexchange; // Rate of exchange

		long double scale; // Scaling factor for frequencies

		int n_coeff;   // number of potential coefficients (in all models)
		long double *coeff; // SRLS coefficients
		potential Ucoeff;   // SRLS coefficients
		dcvector fb1_coeff; // vector of (n>=0) coefficients of the internal potential in FB1 model
		potential_fb2 fb2_coeff; // FB2 model potential coefficients
		int maxEval; // Max number of evalutions in cubature integration routine
		double relErr; // Target relative error in cubature integration routine

		long double omega_D[3];     // MF -> DF
		long double omega_D_bkp[3]; // a backup of omega_D input
		long double omega_D2[3];    // DF -> D2F (for CH2 probes)
		long double omega_CSA[3];   // DF -> CF

		/* TS-SRLS */
		double omegaDip_1[3], omegaDip_2[3];
		double omegaCsa_1[3], omegaCsa_2[3];
		double omegaDip2_1[3], omegaDip2_2[3];
		double population, jumpFrequency, hchSigma;

		/* 3S-FB */
		double _3SFB_OmegaDip1[3][3];
		double _3SFB_OmegaDip2[3][3];
		double _3SFB_OmegaCSA[3][3];
		double _3SFB_population[3], _3SFB_sqPopulation[3];
		double _3SFB_jumpFrequency[3];

		/* Other */
		bool acfFound;
		int Nt, Nw;
		double ti, tf, wi, wf;

		std::ostringstream outputString;
};

#endif /*PHYSICS_H_*/

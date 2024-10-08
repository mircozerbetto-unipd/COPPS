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
 Name        : copps.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Main of C++OPPS
 ============================================================================
 */

#include "copps.h"

int mpi_rank, mpi_ntasks;
int nrows, global_row;
int fitStep, comp;
int OUT_CONTROL;
double dof;
double oldD01, oldD02;
std::string dynModel;
ostringstream fitout;
physics data;
s3j symbols3j;
basis basisFunctions;
stvec Tkk1;
matrix mat;
lanczos lcz;
relax rel;
experimental expdata;
string path, project;
euler eul;
int nOrderParametersStopChecks;
double oldS20, oldS22;
int fn;

/* Powell related quantities */
int powell_n;
int powell_m;
int powell_np;
double *powell_otherData;
double *powell_fvec;
double *powell_sim;
/*****************************/

/* LAPACK ROUTINES */
extern "C"
{
	void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
	void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}
/*******************/

/***************************************************************/
int main(int argc, char *argv[])
/***************************************************************/
{
#ifdef _MPI_
	
	int mpi_err;
	
	/*******************
	 *  Initialize MPI *
	 *******************/
	
	if ((mpi_err = MPI_Init(&argc,&argv)) != MPI_SUCCESS){
		std::cout << "\n\n COOPS ERROR : Cannot initialize MPI environment.\n\n";
		return 1;
	}
	
	/****************************
	 *  Retrive number of tasks *
	 ****************************/
	
	if ((mpi_err = MPI_Comm_size(MPI_COMM_WORLD,&mpi_ntasks)) != MPI_SUCCESS){
		std::cout << "\n\nCOOPS ERROR : Cannot retrive number of tasks.\n\n";
		return 1;
	}
	
	/*****************
	 *  Retrive rank *
	 *****************/
	
	if ((mpi_err = MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank)) != MPI_SUCCESS){
		std::cout << "\n\nCOOPS ERROR : Cannot retrive rank.\n\n";
		return 1;
	}

#else
	
	mpi_rank = 0;
	mpi_ntasks = 1;
	
#endif
	
	/***********************/
	/* Version and license */
	/***********************/

	char VERSION[8];
	sprintf(VERSION,"2.1");
	if(!mpi_rank)
	{
		std::cout << std::endl; system("date"); std::cout << std::endl;
		char gpl_info[2];

		printf("\n\n");
		printf("/*****************************************************************\n");
		printf(" * C++OPPS %s, Copyright (C) 2008 Mirco Zerbetto                *\n",VERSION);
		printf(" * C++OPPS %s comes with ABSOLUTELY NO WARRANTY; for details    *\n",VERSION);
		printf(" * type `copps w'.  This is free software, and you are welcome   *\n");
		printf(" * to redistribute it under certain conditions; type `copps c'   *\n");
		printf(" * for details. Type 'copps v' for version info.                 *\n");
		printf(" *****************************************************************/\n");

		if(argc==2){
			sprintf(gpl_info,"%s",argv[1]);
			if (!strcmp(gpl_info,"w")){
				printf("\n\nWARRANTY NOTES:\n");
				printf("This program is distributed in the hope that it will be useful,\n");
				printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
				printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
				printf("GNU General Public License for more details.\n\n");
				exit(0);
			}
			else if (!strcmp(gpl_info,"c")){
				printf("\n\nREDISTRIBUTION NOTES:\n");
				printf("This program is free software; you can redistribute it and/or\n");
				printf("modify it under the terms of the GNU General Public License\n");
				printf("as published by the Free Software Foundation; either version 2\n");
				printf("of the License, or any later version.\n\n");
				exit(0);
			}
			else if (!strcmp(gpl_info,"v")){
				printf("\n\nC++OPPS VERSION IS %s\n\n",VERSION);
				exit(0);
			}
		}
	}

	/********
	 * Code *
	 ********/
	
	/*********************************************
	 * Check how many arguments have been passed *
	 *********************************************/
	
	if (argc < 3)
	{
		std::cout << std::endl << "ERROR. Project ABSOLUTE PATH and NAME have to be passed in calling COPPS" << std::endl << std::endl;
		exit(0);
	}
	
	/********************************************
	 * Determine absolute path and project name *
	 ********************************************/
	
	path.append(argv[1]);
	project.append(argv[2]);
	if (!mpi_rank)
	{
		std::cout << std::endl << std::endl <<"Echo passed arguments:\n* C++OPPS job '" << project << "'\n* job path '" << path << std::endl;
		if(argc == 4)
			std::cout << "* C / J flag " << argv[3];
		std::cout << std::endl ;
	}

	/*********************************************
	 * Remove the .finished file in project path *
	 *********************************************/
	
	if (!mpi_rank){
		fstream finish;
		finish.open((path+"/.finished").c_str(),ios::in);
		if( finish.is_open() )
		{
			system(("rm -rf "+path+"/.finished").c_str());
			std::cout << ".finished file found and succesfully removed..."<<std::endl<<std::endl;
		}
		else
			std::cout << ".finished file not found..."<<std::endl<<std::endl;
		finish.close();
	}

	/**************************************
	 * Determine which data must be given *
	 **************************************/
	
	// NB: OUT_CONTROL may be changed inside physics.cpp

	if (argc >= 4)
		sscanf(argv[3],"%d",&OUT_CONTROL);
	else
		OUT_CONTROL = 0;
	if (OUT_CONTROL != 0 && OUT_CONTROL != 1)
	{
		std::cout << std::endl << std::endl << "Error : ouput code " << OUT_CONTROL << " not recognized. Please, choose between: " << std::endl << " 0 - calculation on relaxation times only;" << std::endl << " 1 - calculation of relaxatin times, eigenvalues and weights, correlation functions and spectral densities;" << std::endl << std::endl;
		exit(1);
	}

	/**********
	 * Timing *
	 **********/
	
	clock_t start, stop;
	////double tosec = (1.0/(double)CLOCKS_PER_SEC);
	
	/**************************************************
	 * Initiate physics object and read physical data *
	 **************************************************/
	
	start = clock();
	data.init(mpi_rank);
#ifdef _WINDWOS_
	string fileName = path + "\\" + project + "_copps.input";
#else
	string fileName = path + "/" + project + "_copps.input";
#endif
	data.readInputFile(fileName);
#ifdef WRITE_ALL
#ifdef WRITE_PHYSICS
	if (!mpi_rank)
	{
		std::cout << "TASK 0: Echoing physical and calculation parameters:" << std::endl;
		std::cout << data.toString();
	}
#endif
#endif
	stop = clock();
	double physics_time = getcputime();

    	/****************************
    	 * Obtain experimental data *
    	 ****************************/

#ifdef _WINDWOS_
	fileName = path + "\\" + project + "_copps.exp";
#else
	fileName = path + "/" + project + "_copps.exp";
#endif
    	expdata.init(fileName,mpi_rank);
    	bool haveExpFile = expdata.readData(data.getFitting());
#ifdef WRITE_ALL
	if (data.getFitting())
		std::cout << expdata.toString();
#endif

	/***************************************
	 Obtain model of dynamics to be used *
	***************************************/

	dynModel = data.getDynamicsModel();

	/*********************** 
	 * Generate 3j symbols *
	 ***********************/

	start = clock();
	if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
		symbols3j.init(data.getProbeLMax(),mpi_rank);
	else
		symbols3j.init(20,mpi_rank);
	symbols3j.storeSymbols();
	stop = clock();
	double s3jstore_time = getcputime();
	
	/****************** 
	 * Generate basis *
	 ******************/
	
	start = clock();
	if (!data.readSmallJs)
	{
		if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
		{
			basisFunctions.init(2,2,0,0,-2,2,0,data.getProbeLMax(),-data.getProbeLMax(),data.getProbeLMax(),-data.getProbeLMax(),data.getProbeLMax(),mpi_rank);
			basisFunctions.isOmegaVFromGeometry = expdata.isOmegaVFromGeometry;
			ldvector ldv(6);
			ldv[0] = data.getValueOf("protein_alpha");
			ldv[1] = data.getValueOf("protein_beta");
			ldv[2] = data.getValueOf("protein_gamma");
			ldv[3] = data.getProteinDxx() - data.getProteinDyy();
			ldv[4] = abs(data.getProbeDlm(2,1))+abs(data.getProbeDlm(2,2))+abs(data.getProteinDlm(2,1))+abs(data.getProteinDlm(2,2));
			if (!data.getPotentialFromExpFile())
				ldv[5] = fabs(data.getValueOf("c22")) + fabs(data.getValueOf("c42")) + fabs(data.getValueOf("c44"));
			else
				ldv[5] = 1.0;
			basisFunctions.checkSymmetry(ldv,data.getFitPar());
		}
		else if (!(dynModel.compare("fb1")) || !(dynModel.compare("ts-fb1")) || !(dynModel.compare("3s-fb")))
		{
			basisFunctions.init(2,2,0,0,-2,2,0,0,0,0,-data.getNmax(),data.getNmax(),mpi_rank);
		}
		else if (!(dynModel.compare("fb2")))
		{
			basisFunctions.init(2,2,0,0,-2,2,0,0,0,0,data.getNmax(1),data.getNmax(2),mpi_rank);
		}
		basisFunctions.buildSpace(dynModel);
#ifdef WRITE_ALL
		std::cout << "TASK " << mpi_rank << ": Number of basis functions = " << basisFunctions.getNumberOfBasisFunctions() << std::endl;
#ifdef WRITE_BASIS
		std::cout << basisFunctions.toString(dynModel);
#endif
#endif
	}
	stop = clock();
	double basis_time = getcputime();

	/*********************************************************
	 * Cut the basis in slices basing on the number of tasks *
	 *********************************************************/
	
	int nbf, ncols;
	if (!data.readSmallJs)
	{
		nbf = basisFunctions.getNumberOfBasisFunctions();
		ncols = nbf;
		int quotient = nbf / mpi_ntasks;
		int reminder = nbf % mpi_ntasks;
		nrows = quotient + (mpi_rank < reminder ? 1 : 0);
		global_row = 0;
		for (int i = 0; i < mpi_rank; i++){
			global_row += quotient + (i < reminder ? 1 : 0);
		}
#ifdef WRITE_ALL
		std::cout << "TASK " << mpi_rank << ": Number of rows = " << nrows << "\t" << "Global matrix row = " << global_row << std::endl;
#endif
	}
	
	/*******************************
	 * Init starting vector object *
	 *******************************/
	
	Tkk1.init(&data,mpi_rank);
		
	/**********************
	 * Init matrix object *
	 **********************/
	
	if (!data.readSmallJs)
		mat.init(&basisFunctions,&data,&symbols3j,mpi_rank,nrows,global_row,ncols,0);
	
	/***********************
	 * Init lanczos object *
	 ***********************/
	
	if (!data.readSmallJs)
	{
		int nstep = data.getLanczosNStep();
		if (nstep <= 0)
			nstep = 15*basisFunctions.getNumberOfBasisFunctions() / 100;
		lcz.init(nstep,&Tkk1,&mat,data.getSpectralDensities(),mpi_rank,mpi_ntasks);
		lcz.setScale((double)data.getScale());
	}

	/*********************
	 * Init relax object *
	 *********************/
	
	dvector B0;
	if (!haveExpFile)
	    	B0 = data.getField();
	else
		B0 = expdata.getFieldOfComp(1);
	
	if (!(dynModel.compare("ts-srls")) || !(dynModel.compare("ts-fb1")))
	{
		dvector od_1  = data.getOmegaDipolarState(1);
		dvector od_2  = data.getOmegaDipolarState(2);
		dvector oc_1  = data.getOmegaCsaState(1);
		dvector oc_2  = data.getOmegaCsaState(2);
		dvector od2_1 = data.getOmegaDipolar2State(1);
		dvector od2_2 = data.getOmegaDipolar2State(2);
		rel.init(&Tkk1,&lcz,od_1,oc_1,od2_1,od_2,oc_2,od2_2,data.getPopulation(),data.getJumpFrequency(),B0,data.getScale(),data.getSpectralDensities(),data.getNHydrogens(),mpi_rank);
	}
	else if (!dynModel.compare("3s-fb"))
	{
		dvector od1_1 = data._3SFB_getOmegaDipolar1Site(1);
		dvector od2_1 = data._3SFB_getOmegaDipolar2Site(1);
		dvector oc_1  = data._3SFB_getOmegaCSASite(1);

		dvector od1_2 = data._3SFB_getOmegaDipolar1Site(2);
		dvector od2_2 = data._3SFB_getOmegaDipolar2Site(2);
		dvector oc_2  = data._3SFB_getOmegaCSASite(2);

		dvector od1_3 = data._3SFB_getOmegaDipolar1Site(3);
		dvector od2_3 = data._3SFB_getOmegaDipolar2Site(3);
		dvector oc_3  = data._3SFB_getOmegaCSASite(3);

		rel.init(&Tkk1,&lcz,od1_1,oc_1,od2_1,od1_2,oc_2,od2_2,od1_3,oc_3,od2_3,data._3SFB_getSqPopulations(),data._3SFB_getJumpFrequencies(),B0,data.getScale(),data.getSpectralDensities(),data.getNHydrogens(),mpi_rank);
	}
	else
	{
		ldvector omegaDip = data.getOmegaDipolar();
		ldvector omegaDip2 = data.getOmegaDipolar2();
		ldvector omegaCsa = data.getOmegaCSA();
		rel.init(&Tkk1,&lcz,omegaDip,omegaDip2,omegaCsa,B0,data.getScale(),data.getSpectralDensities(),data.getNHydrogens(),mpi_rank);
	}
	rel.setNuclearData(data.getNuclearData());
	rel.setHchSigma(data.getHchSigma());

	/**************************************************
	 * Select between one-shot calculation or fitting *
	 **************************************************/
	
	if (!data.getFitting() || data.readSmallJs)
	{
		double stvec_time   = 0.0;
		double matrix_time  = 0.0;
		double lanczos_time = 0.0;
		double relax_time   = 0.0;

		dvvector T1T2NOE;
		dvector orderParameters;
	
		if (data.readSmallJs)
		{
			rel.readSmallJs(data.getNw());
			T1T2NOE = rel.getT1T2NOE(true);

		}
		else
		{
			potential u, v;
			if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
			{
				if (haveExpFile && expdata.compHasPotential(1))
				u = expdata.getPotentialOfComp(1);
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
			
					data.setPotentialCoefficients(v);
					Tkk1.updatePotential();
					mat.update();
				}
		
				/*******************************
				 * Calculate Orders Parameters *
				 *******************************/
			
				if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
				{

					orderParameters = Tkk1.calculateOrdersParams();

#ifdef WRITE_ALL
#ifdef WRITE_ORDERPAR
					if (!mpi_rank)
					{
						std::cout << "TASK " << mpi_rank << ": " << std::endl;
						std::cout << "TASK " << mpi_rank << ": " << "Order parameters for population " << npop+1 << ":"<< std::endl;
						std::cout << "TASK " << mpi_rank << ": " << std::endl;
						std::cout << "TASK " << mpi_rank << ": " << fixed << setprecision(4) << "S(2,0) = " << orderParameters[0] << std::endl;
						std::cout << "TASK " << mpi_rank << ": " << fixed << setprecision(4) << "S(2,2) = " << orderParameters[1] << std::endl;
						std::cout << "TASK " << mpi_rank << ": " << fixed << setprecision(4) << "Sxx    = " << orderParameters[2] << std::endl;
						std::cout << "TASK " << mpi_rank << ": " << fixed << setprecision(4) << "Syy    = " << orderParameters[3] << std::endl;
						std::cout << "TASK " << mpi_rank << ": " << fixed << setprecision(4) << "Szz    = " << orderParameters[4] << std::endl;
						std::cout << "TASK " << mpi_rank << ": " << std::endl;
					}
#endif
#endif
				}

				/***************************** 
				 * Calculate starting vector *
				 *****************************/
			
				start = clock();
#ifdef WRITE_ALL
				if (!mpi_rank) std::cout << "TASK 0: Projecting starting vector on basis functions..." << std::endl;
#endif
				Tkk1.projectOnBasis(&basisFunctions,nrows,global_row);
#ifdef _MPI_
				Tkk1.scatterProjections(true,mpi_ntasks);
#else
				Tkk1.scatterProjections(false,1);
#endif
#ifdef WRITE_ALL
#ifdef WRITE_STVEC
				if (!mpi_rank)
					std::cout << Tkk1.toString(&basisFunctions);
#endif
#endif
				stop = clock();
				stvec_time += getcputime();
			
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
				matrix_time += getcputime();
			
				/***************
				 * Run Lanczos *
				 ***************/
				
				start = clock();
				lcz.update();
				lcz.runLaczos();
#ifdef WRITE_ALL
#ifdef WRITE_LANCZOS
				std::cout << lcz.toString(true,true,0);
				std::cout << lcz.toString(true,true,1);
				std::cout << lcz.toString(true,true,2);
				std::cout << lcz.toString(true,true,3);
				std::cout << lcz.toString(true,true,4);
				std::cout << lcz.toString(true,true,5);
#endif
#endif
				stop = clock();
				lanczos_time += getcputime();
			
				/******************************
				 * Calculate relaxation times *
				 ******************************/
			
				start = clock();
				tmp_T1T2NOE = rel.getT1T2NOE(false);
				stop = clock();
				relax_time += getcputime();

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

			} // std::endl of cycle over populations

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
		}

#ifdef WRITE_ALL
#ifdef WRITE_RELAX
		if (!mpi_rank)
		{
			if (!haveExpFile)
			{
				ostringstream toOut;
        			toOut << std::endl << "* Relaxation data:" << std::endl;
			        int w2 = 11;
			        int w3 = 2;

			        toOut << std::endl;
			        toOut << "             field/MHz      theo"<< std::endl;
			        toOut << "  ______________________________________" << std::endl<<std::endl;

        			for (unsigned int ib = 0; ib < B0.size(); ib++)
			        {
					/* T1 */
					toOut << "    T1/ms:      " << fixed << setprecision(3) << setw(w2-2*w3) << B0[ib]*1.0e-6 \
						<< fixed << setprecision(3) << setw(w2) << T1T2NOE[ib][0]*1.0e3 << std::endl;
					/* T2 */
					toOut << "    T2/ms:      " << fixed << setprecision(3) << setw(w2-2*w3) << B0[ib]*1.0e-6 \
						<< fixed << setprecision(3) << setw(w2) << T1T2NOE[ib][1]*1.0e3 << std::endl;
					/* NOE */
					toOut << "    NOE:        " << fixed << setprecision(3) << setw(w2-2*w3) << B0[ib]*1.0e-6 \
						<< fixed << setprecision(3) << setw(w2) << T1T2NOE[ib][2] << std::endl;
					if (data.getNHydrogens() == 2)
					{
						/* CCRRT1 */
						toOut << "    CCRRT1/Hz:  " << fixed << setprecision(3) << setw(w2-2*w3) << B0[ib]*1.0e-6 \
							<< fixed << setprecision(3) << setw(w2) << T1T2NOE[ib][3] << std::endl;
						/* CCRRTT2 */
						toOut << "    CCRRT2/Hz:  " << fixed << setprecision(3) << setw(w2-2*w3) << B0[ib]*1.0e-6 \
							<< fixed << setprecision(3) << setw(w2) << T1T2NOE[ib][4] << std::endl;
					}
		
					toOut << "  ______________________________________" << std::endl<<std::endl;
				}
				std::cout << toOut.str();

				/* Write the project_copps.output file */ 

				fstream file_op((path+"/"+project+"_copps.output").c_str(),ios::out);
				file_op << "C++OPPS VERSION " << VERSION << std::endl << std::endl;
				file_op << toOut.str();
				file_op << std::endl;
				if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
				{
					file_op << "Order parameters:"<< std::endl <<std::endl;
					file_op << fixed << setprecision(4) << "S(2,0) = " << orderParameters[0] << std::endl;
					file_op << fixed << setprecision(4) << "S(2,2) = " << orderParameters[1] << std::endl;
					file_op << fixed << setprecision(4) << "Sxx    = " << orderParameters[2] << std::endl;
					file_op << fixed << setprecision(4) << "Syy    = " << orderParameters[3] << std::endl;
					file_op << fixed << setprecision(4) << "Szz    = " << orderParameters[4] << std::endl;
				}
				file_op.close();
			}
			else
			{
				dof = 1;
				int m = 1;
				double *vexp = expdata.getExpDataOfComp(&m);
				m >>= 1;
				double* fvec = new double[m];
				double chisq = calculateChiSquare(T1T2NOE,vexp,fvec);
				outputData(T1T2NOE,vexp,chisq,false);
			}
		}
#endif
#endif

		/*********************
		 * Writes time table *
		 *********************/
		
		int ntimes = 8;
		double local_times[ntimes];
		double global_times[ntimes];
		local_times[0] = global_times[0] = physics_time;
		local_times[1] = global_times[1] = s3jstore_time - physics_time;
		local_times[2] = global_times[2] = basis_time - s3jstore_time;
		local_times[3] = global_times[3] = stvec_time - basis_time;
		local_times[4] = global_times[4] = matrix_time - stvec_time;
		local_times[5] = global_times[5] = lanczos_time - matrix_time;
		local_times[6] = global_times[6] = relax_time - lanczos_time;
		local_times[7] = global_times[7] = relax_time;
		
#ifdef _MPI_
		MPI_Allreduce(local_times,global_times,ntimes,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
		std::cout << "TASK " << mpi_rank << ": Task CPU time:" << std::endl;
		std::cout << "TASK " << mpi_rank << ": * physical data capture       : " << scientific << setprecision(2) << local_times[0]  << " s (" << fixed << setprecision(2) << 100.0*local_times[0]/local_times[7]  << " %)" << std::endl;
		std::cout << "TASK " << mpi_rank << ": * 3j symbols storage          : " << scientific << setprecision(2) << local_times[1] << " s (" << fixed << setprecision(2) << 100.0*local_times[1]/local_times[7] << " %)" << std::endl;
		std::cout << "TASK " << mpi_rank << ": * basis indexes construction  : " << scientific << setprecision(2) << local_times[2]    << " s (" << fixed << setprecision(2) << 100.0*local_times[2]/local_times[7]    << " %)" << std::endl;
		std::cout << "TASK " << mpi_rank << ": * starting vector projections : " << scientific << setprecision(2) << local_times[3]    << " s (" << fixed << setprecision(2) << 100.0*local_times[3]/local_times[7]    << " %)" << std::endl;
		std::cout << "TASK " << mpi_rank << ": * matrix construction         : " << scientific << setprecision(2) << local_times[4]   << " s (" << fixed << setprecision(2) << 100.0*local_times[4]/local_times[7]   << " %)" << std::endl;
		std::cout << "TASK " << mpi_rank << ": * Lanczos tridiagonalization  : " << scientific << setprecision(2) << local_times[5]  << " s (" << fixed << setprecision(2) << 100.0*local_times[5]/local_times[7]  << " %)" << std::endl;
		std::cout << "TASK " << mpi_rank << ": * T1, T2 and NOE calculation  : " << scientific << setprecision(2) << local_times[6]    << " s (" << fixed << setprecision(2) << 100.0*local_times[6]/local_times[7]  << " %)" << std::endl;
		std::cout << "TASK " << mpi_rank << ": * total running time          : " << scientific << setprecision(2) << local_times[7]    << " s (" << fixed << setprecision(2) << 100.0  << " %)" << std::endl;
		if (!mpi_rank & mpi_ntasks > 1)
		{
#ifdef _MPI_    	
			std::cout << "TASK " << mpi_rank << ": Max CPU time over tasks:" << std::endl;
			std::cout << "TASK " << mpi_rank << ": * physical data capture       : " << scientific << setprecision(2) << global_times[0] << " s " << std::endl;
			std::cout << "TASK " << mpi_rank << ": * 3j symbols storage          : " << scientific << setprecision(2) << global_times[1] << " s " << std::endl;
			std::cout << "TASK " << mpi_rank << ": * basis indexes construction  : " << scientific << setprecision(2) << global_times[2] << " s " << std::endl;
			std::cout << "TASK " << mpi_rank << ": * starting vector projections : " << scientific << setprecision(2) << global_times[3] << " s " << std::endl;
			std::cout << "TASK " << mpi_rank << ": * matrix construction         : " << scientific << setprecision(2) << global_times[4] << " s " << std::endl;
			std::cout << "TASK " << mpi_rank << ": * Lanczos tridiagonalization  : " << scientific << setprecision(2) << global_times[5] << " s " << std::endl;
			std::cout << "TASK " << mpi_rank << ": * T1, T2 and NOE calculation  : " << scientific << setprecision(2) << global_times[6] << " s " << std::endl;
			std::cout << "TASK " << mpi_rank << ": * total running time          : " << scientific << setprecision(2) << global_times[7] << " s " << std::endl;
#endif
		}
	}
	else // FIT
	{
		/*******************************
		 * Start cycle over components *
		 *******************************/

		fitout.str("");
		
		for (int c = 0; c < expdata.getNComp(); c++)
		{
		
			expdata.setActiveComponent(c+1);

			/***************************************
			 * Initial guess of fitting parameters *
			 ***************************************/
			
#ifdef _WINDWOS_
			string fileName = path + "\\" + project + "_copps.input";
#else
			string fileName = path + "/" + project + "_copps.input";
#endif
			data.readInputFile(fileName);

			/***************************/
			/* Start cycle over blocks */
			/***************************/

			int fnMax = data.getNFit2() ? MAX_BLOCK_FIT_CYCLES : 1;
			for (fn = 0; fn < fnMax; ++fn)
			{
			
				int fitBlock = fn % 2;
				int i;
				string fitRoutine = data.getFittingMethod();
				
				int n = !fitBlock ? data.getNFit() : data.getNFit2();
				double par[n];
				bool *fitpar = !fitBlock ? data.getFitPar() : data.getFitPar2();

				/****************************
				 * Update orientation of VF *
				 ****************************/
				if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
				{
					// Here the M1F -> VF angles given in the input are overwritten in order to fit local geometry
					data.setValueOf("protein_alpha",(long double)expdata.getAlphaV(c));	
					data.setValueOf("protein_beta",(long double)expdata.getBetaV(c));	
					data.setValueOf("protein_gamma",(long double)expdata.getGammaV(c));
					data.transformProteinD();
				}

				/********************************************/
				/* Store the rotation matrix from M1F to DF */
				/********************************************/
				if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
				{
					dvector o1 = dvector(3,0.0);
					o1[0] = (double)expdata.getAlphaV(c);
					o1[1] = (double)expdata.getBetaV(c);
					o1[2] = (double)expdata.getGammaV(c);
					dvector o12 = dvector(3,0.0);
					o12[0] = (double)data.getBkpDipolarAngle("alpha");
					o12[1] = (double)data.getBkpDipolarAngle("beta");
					o12[2] = (double)data.getBkpDipolarAngle("gamma");
					eul.setOmega1(o1);
					eul.setOmega12(o12);
					eul.omega1ToE1();
					eul.omega12ToE12();
					eul.calculateE2();
				}

				// Mod for Levmar library
				if (!fitRoutine.compare("levmar"))
				{
					data.setValueOf("protein_dyy",data.getValueOf("protein_dyy")+0.1);
					data.setValueOf("probe_dyy",data.getValueOf("probe_dyy")+0.1);
				}
				
				// If potential coefficients have been specified in the .exp file fitting
				// if automatically disabled
				if (expdata.compHasPotential(c+1))
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
				if (fitpar[0])  {par[i] = log((double)data.getValueOf("protein_dyy")); ++i;}
				if (fitpar[1])  {par[i] = (double)(data.getValueOf("protein_dxx")/data.getValueOf("protein_dyy")); ++i;}
				if (fitpar[2])  {par[i] = (double)(data.getValueOf("protein_dzz")/data.getValueOf("protein_dyy")); ++i;}
				// Omega V
				if (fitpar[6])  {par[i] = (double)data.getValueOf("protein_alpha");    ++i;}
				if (fitpar[7])  {par[i] = (double)data.getValueOf("protein_beta");     ++i;}
				if (fitpar[8])  {par[i] = (double)data.getValueOf("protein_gamma");    ++i;}
				// Omega O
				if (fitpar[9])  {par[i] = (double)data.getValueOf("probe_alpha");      ++i;}
				if (fitpar[10]) {par[i] = (double)data.getValueOf("probe_beta");       ++i;}
				if (fitpar[11]) {par[i] = (double)data.getValueOf("probe_gamma");      ++i;}
				// Omega D
				if (fitpar[12]) {par[i] = (double)data.getValueOf("dipolar_alpha");    ++i;}
				if (fitpar[13]) {par[i] = (double)(data.getValueOf("dipolar_beta"));   ++i;}
				if (fitpar[14]) {par[i] = (double)data.getValueOf("dipolar_gamma");    ++i;}
				// Omega D2
				if (fitpar[26]) {par[i] = (double)data.getValueOf("dipolar_alpha2");   ++i;}
				if (fitpar[27]) {par[i] = (double)(data.getValueOf("dipolar_beta2"));  ++i;}
				if (fitpar[28]) {par[i] = (double)data.getValueOf("dipolar_gamma2");   ++i;}
				// probe diffusion
				if (fitpar[3])  {par[i] = log((double)data.getValueOf("probe_dyy")); ++i;}
				if (fitpar[4])  {par[i] = (double)(data.getValueOf("probe_dxx")/data.getValueOf("probe_dyy"));   ++i;}
				if (fitpar[5])  {par[i] = (double)(data.getValueOf("probe_dzz")/data.getValueOf("probe_dyy"));   ++i;}
				// Omega C
				if (fitpar[15]) {par[i] = (double)data.getValueOf("csa_alpha");        ++i;}
				if (fitpar[16]) {par[i] = (double)data.getValueOf("csa_beta");         ++i;}
				if (fitpar[17]) {par[i] = (double)data.getValueOf("csa_gamma");        ++i;}
				// Potential coefficients
		////		if (fitpar[18] && fitpar[19]) // If fitting both c20 and c22 the program fits the order parameters
		////		{ 
		////			Tkk1.updatePotential();
		////			dvector op = Tkk1.calculateOrdersParams();
		////			par[i] = op[0]; ++i;
		////			par[i] = op[1]; ++i;
		////		}
		////		else
		////		{
					if (fitpar[18]) {par[i] = (double)data.getValueOf("c20");       ++i;}
					if (fitpar[19])
					{
						if (ratio22)
							par[i] = (double)(data.getValueOf("c22")/data.getValueOf("c20"));
						else
							par[i] = (double)data.getValueOf("c22");
						++i;
					}
		////		}
				if (fitpar[20]) {par[i] = (double)data.getValueOf("c40");              ++i;}
				if (fitpar[21])
				{
					if (ratio42)
						par[i] = (double)(data.getValueOf("c42")/data.getValueOf("c40")); 
					else
						par[i] = (double)data.getValueOf("c42");
					++i;
				}
				if (fitpar[22])
				{
					if (ratio44)
						par[i] = (double)(data.getValueOf("c44")/data.getValueOf("c40"));
					else
						par[i] = (double)data.getValueOf("c44");
					++i;
				}
				// Conformational exchange rate
				if (fitpar[23]) {par[i] = (double)log(data.getValueOf("Rexchange"));   ++i;}
				// Correction to FBn diffusion tensor
				if (fitpar[25]) {par[i] = sqrt(1.0);         ++i;}
				// Population in ts-srls model
				if (fitpar[29]) {par[i] = acos(sqrt(data.getPopulation())); ++i;} // population = cos^2(par)
				// Jump frequency in ts-srls model
				if (fitpar[30]) {par[i] = log(data.getJumpFrequency()); ++i;} // jump freq = exp(par)
				// hch sigma
				if (fitpar[31]) {par[i] = log(data.getHchSigma()); ++i;}
			
				/*********************
				 * Init relax object *
				 *********************/
				
				dvector B0(1,1.0e8);
				ldvector omegaDip = data.getOmegaDipolar();
				ldvector omegaDip2 = data.getOmegaDipolar2();
				ldvector omegaCsa = data.getOmegaCSA();
				rel.init(&Tkk1,&lcz,omegaDip,omegaDip2,omegaCsa,B0,data.getScale(),data.getSpectralDensities(),data.getNHydrogens(),mpi_rank);
				rel.setNuclearData(data.getNuclearData());
			
				fitStep = 0;
				comp = c+1;
				
				/***********************************************
				 * Obatain experimental data for component c+1 *
				 ***********************************************/
				
				int m = c + 1;
				rel.update(data.getOmegaDipolar(), data.getOmegaDipolar2(), data.getOmegaCSA(), expdata.getFieldOfComp(m), data.getPopulation(), data.getJumpFrequency());
				rel.setHchSigma(data.getHchSigma());
				double *vexp = expdata.getExpDataOfComp(&m);
				m >>= 1;

				/**********************/
				/* Degrees of freedom */
				/**********************/
			
				dof = (double)(m-n);
				if (dof <= 0) dof = 1;

				/*************************/
				/* Other initializations */
				/*************************/

				nOrderParametersStopChecks = 0;
				oldS20 = oldS22 = 100.0;

				if (!fitRoutine.compare("minpack"))
				{

		//			LMDIF CODE
		//
					/******************************************
					 * Set parameters for the fitting routine *
					 ******************************************/
		//    	
					int maxfev = 200*(n+1);
					int mode = 1;
					int nprint = 0;
					int nfev = 0;
					int ldfjac = m;
					int* ipvt = new int[n];
					double factor = 0.01;
					double xtol = data.getValueOf("fit_tolerance_1");
					double ftol = xtol;
					double gtol = xtol;
					double epsfcn = data.getValueOf("fit_tolerance_2");
					double* diag = new double[n];
					//for (int ii = 0; ii < n; ii++) diag[ii] = 1.0;
					//diag[0] = 0.01;
					double* fjac = new double[m*n];
					double* jacob = new double[m*n];
					double* qtf = new double[n];
					double* wa1 = new double[n];
					double* wa2 = new double[n];
					double* wa3 = new double[n];
					double* wa4 = new double[m];
					double* fvec = new double[m];
					if (!mpi_rank) fvec[0] = -1.0; else fvec[0] = 0.0;
		//    		
					/******************
					 * Launch fitting *
					 ******************/
	#ifdef _MPI_	
					MPI_Barrier(MPI_COMM_WORLD);
	#endif
		//
					std::cout << std::endl << std::endl << "NUMBER OF FITTING PARAMS = " << n << std::endl << std::endl;
		//
					int info = lmdif(observables_minpack, vexp, m, n, par, fvec, 
						ftol, xtol, gtol, maxfev, epsfcn, 
						diag, mode, factor, nprint, &nfev, 
						fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, jacob);
					double fnorm = enorm(m, fvec);
					printf("      final l2 norm of the residuals%15.7g\n\n", fnorm);
					printf("      number of function evaluations%10i\n\n", nfev);
					printf("      exit parameter                %10i\n\n", info);
		//
					/************************
					 * Write final analysis *
					 ************************/
					
					
					if (!mpi_rank)
					{
						fitout << "Fitting procedure termintated for component " << c+1 << std::endl;
						fitout << "Final l2 norm of the residual = " << fnorm << std::endl;
						fitout << "Exit parameter = " << info << std::endl;
						/* Output raw parameters */
						fitout << std::endl << "Raw fitting parameters: " <<  std::endl;
						for (int i = 0; i < n; ++i)
							fitout << "p(" <<i+1 << ") = " << par[i] << std::endl;
						/* Output Jacobian, J, and calculate J = W * J */
						// W is diagonal matrix containing inverses of weights
						fitout << std::endl << "Jacobian: " << std::endl;
						for (int i = 0; i < m; ++i)
						{
							for (int j = 0; j < n; ++j)
							{
								fitout << scientific << setprecision(2) << jacob[i*n+j] << "\t";
								jacob[i*n+j] /= vexp[i*2+1];
							}
							fitout << std::endl;
						}
						/* Calculate varaince-covariance matrix  V = inv(J^t W W J) */
						fitout << std::endl << "Variance-Covariance matrix: " << std::endl;
						double *V = new double[n*n];
						for (int i = 0; i < n; ++i)
						{
							for (int j = 0; j < n; ++j)
							{
								V[i*n+j] = 0.0;
								for (int k = 0; k < m; ++k)
									V[i*n+j] += jacob[k*n+i] * jacob[k*n+j];
							}
						}
						int *ipiv2 = new int[n+1];
						int lwork2 = n*n;
						double *work2 = new double[lwork2];
						int info2;
						dgetrf_(&n,&n,V,&n,ipiv2,&info2);
						dgetri_(&n,V,&n,ipiv2,work2,&lwork2,&info2);
						delete ipiv2;
						delete work2;
						for (int i = 0; i < n; ++i)
						{
							for (int j = 0; j < n; ++j)
								fitout << scientific << setprecision(2) << V[i*n+j] << "\t";
							fitout << std::endl;
						}
						/* Print the absolute and relative errors */
						fitout << std::endl << "Absolute and relative errors:" << std::endl;
						for (int i = 0; i < n; ++i)
							fitout << scientific << setprecision(2) << "s(" << i+1 << ") " << sqrt(V[i*n+i]) << "\t; " << 100.0*sqrt(V[i*n+i])/fabs(par[i]) << "%" << std::endl;

						/*************************************************************************
						 * Calculate Orders Parameters at the std::endl of fitting for every component *
						 *************************************************************************/

						fitout << std::endl;
						if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
						{
							dvector orderParameters = Tkk1.calculateOrdersParams();
							fitout << "Order parameters:"<< std::endl <<std::endl;
							fitout << fixed << setprecision(4) << "S(2,0) = " << orderParameters[0] << std::endl;
							fitout << fixed << setprecision(4) << "S(2,2) = " << orderParameters[1] << std::endl;
							fitout << fixed << setprecision(4) << "Sxx    = " << orderParameters[2] << std::endl;
							fitout << fixed << setprecision(4) << "Syy    = " << orderParameters[3] << std::endl;
							fitout << fixed << setprecision(4) << "Szz    = " << orderParameters[4] << std::endl;
							fitout << std::endl;
						}
	#ifdef WRITE_ALL
						std::cout << fitout.str();
	#endif
						fstream file_op((path+"/"+project+"_copps.output").c_str(),ios::out);
						file_op << "C++OPPS VERSION " << VERSION << std::endl << std::endl;
						file_op<< fitout.str();
						file_op.close();
					}
				}
				else if (!fitRoutine.compare("levmar"))
				{

				       /****************/
				       /* Setup LEVMAR */
				       /****************/
			
				       int maxIter = 5000;
				       double info[8];
				       double opt[5];
				       double *cov = (double*)calloc(n*n,sizeof(double));
				       double *measurements = (double *)calloc(m,sizeof(double));
				       for (int ii = 0; ii < 2*m; ii+=2)
					 measurements[ii/2] = vexp[ii];
					
			
				       /******************
					* Launch fitting *
					******************/
	#ifdef _MPI_    
				  MPI_Barrier(MPI_COMM_WORLD);
	#endif
		//
					std::cout << std::endl << std::endl << "NUMBER OF FITTING PARAMS = " << n << std::endl << std::endl;

					std::cout <<                 "======================================================" << std::endl;
					std::cout << std::endl << std::endl << "WARNING : USAGE OF LEVMAR IN C++OPPS IS STILL UNSTABLE" ;
					std::cout <<                 "======================================================" << std::endl << std::endl;
		//
					dlevmar_dif(observables_levmar,par,measurements,n,m,maxIter,opt,info,NULL,cov,vexp);
		//
				       /************************
					* Write final analysis *
					************************/
			
				       if (!mpi_rank)
				       {
						fitout <<                 "======================================================" << std::endl;
						fitout << std::endl << std::endl << "WARNING : USAGE OF LEVMAR IN C++OPPS IS STILL UNSTABLE" ;
						fitout <<                 "======================================================" << std::endl << std::endl;
					       fitout << std::endl << std::endl;;
					       fitout << "====================================" << std::endl;
					       fitout << "==         Fitting Results        ==" << std::endl;
					       fitout << "====================================";
					       fitout << std::endl << std::endl;
					       fitout << "* Convergence Parameters : ";
					       for (int z = 0; z < n; z++) 
						       fitout << "   " << par[z] ;
					       fitout << std::endl << std::endl << "* ||e||_2 at initial p : " << info[0];
					       fitout << std::endl << "* ||e||_2 at std::endling p  : " << info[1];
					       fitout << std::endl << "* Fitting steps        : " << info[5];
					       fitout << std::endl << "* Stop condition       : " << info[6];
					       fitout << std::endl << std::endl << "* Legstd::endl :";
					       fitout << std::endl << "   1 - stopped by small gradient J^T e" << std::endl << "   2 - stopped by small Dp (small change in parameters)" << std::endl << "   3 - stopped by itmax" << std::endl << "   4 - singular matrix. Restart from current p with increased mu" << std::endl << "   5 - no further error reduction is possible. Restart with increased mu" << std::endl << "   6 - stopped by small ||e||_2";
					       fitout << std::endl << std::endl << "* Covariance matrix : " << std::endl;
					       for (int z = 0; z< n; z++)
					       {
						       for (int z1 = 0; z1 < n; z1++)
						       {
							       fitout << "\t" << cov[z*n+z1];
						       }
						       fitout << std::endl;
					       }
					       fitout << std::endl << std::endl << "====================================" << std::endl;
			
			
					       /*************************************************************************
					       * Calculate Orders Parameters at the std::endl of fitting for every component *
					       *************************************************************************/
			
					       fitout << std::endl;
						if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
						{
							dvector orderParameters = Tkk1.calculateOrdersParams();
							fitout << "Order parameters:"<< std::endl <<std::endl;
							fitout << fixed << setprecision(4) << "S(2,0) = " << orderParameters[0] << std::endl;
							fitout << fixed << setprecision(4) << "S(2,2) = " << orderParameters[1] << std::endl;
							fitout << fixed << setprecision(4) << "Sxx    = " << orderParameters[2] << std::endl;
							fitout << fixed << setprecision(4) << "Syy    = " << orderParameters[3] << std::endl;
							fitout << fixed << setprecision(4) << "Szz    = " << orderParameters[4] << std::endl;
							fitout << std::endl;
						}
						fitout <<                 "======================================================" << std::endl;
						fitout << std::endl << std::endl << "WARNING : USAGE OF LEVMAR IN C++OPPS IS STILL UNSTABLE" ;
						fitout <<                 "======================================================" << std::endl << std::endl;
	#ifdef WRITE_ALL        
						std::cout << fitout.str();
	#endif  
					       fstream file_op((path+"/"+project+"_copps.output").c_str(),ios::out);
					       file_op << "C++OPPS VERSION " << VERSION << std::endl << std::endl;
					       file_op<< fitout.str();
					       file_op.close();
				       }
				}

				/**********/
				/* POWELL */
				/**********/

				else if (!fitRoutine.compare("powell"))
				{
					powell_np = NPARS;
					powell_n  = n;
					powell_m  = m;
					int iter = 0;
					int iflag = 0;
					double ftol = 1.0e-6;
					double atol = 1.0e-6;
					double fret = 0.0;
					powell_otherData = vexp;
					powell_fvec = new double[m];
					double *xi = new double[powell_np*powell_np];
					for (int in1 = 0; in1 < powell_np; in1++)
					{
						for (int in2 = 0; in2 < powell_np; in2++)
						{
							xi[in1*powell_np+in2] = 0.0;
						}
						xi[in1*powell_np+in1] = 1.0;
					}

	#ifdef _MPI_    
					MPI_Barrier(MPI_COMM_WORLD);
	#endif

					powell_(par, xi, &powell_n, &powell_np, &ftol, &atol, &iter, &fret, &iflag);

					if (!mpi_rank)
					{
						fitout << "Fitting procedure termintated for component " << c+1 << std::endl;
						fitout << "Final l2 norm of the residual = " << fret << std::endl;
						fitout << "Exit parameter = " << iflag << std::endl;
						
						/*************************************************************************
						 * Calculate Orders Parameters at the std::endl of fitting for every component *
						 *************************************************************************/
						
						fitout << std::endl;
					
						if (!(dynModel.compare("srls")) || !(dynModel.compare("ts-srls")))
						{
							dvector orderParameters = Tkk1.calculateOrdersParams();
							
							fitout << "Order parameters:"<< std::endl <<std::endl;
							fitout << fixed << setprecision(4) << "S(2,0) = " << orderParameters[0] << std::endl;
							fitout << fixed << setprecision(4) << "S(2,2) = " << orderParameters[1] << std::endl;
							fitout << fixed << setprecision(4) << "Sxx    = " << orderParameters[2] << std::endl;
							fitout << fixed << setprecision(4) << "Syy    = " << orderParameters[3] << std::endl;
							fitout << fixed << setprecision(4) << "Szz    = " << orderParameters[4] << std::endl;
							fitout << std::endl;
						}
	#ifdef WRITE_ALL
						std::cout << fitout.str();
	#endif
						fstream file_op((path+"/"+project+"_copps.output").c_str(),ios::out);
						file_op << "C++OPPS VERSION " << VERSION << std::endl << std::endl;
						file_op<< fitout.str();
						file_op.close();
					}
				}
				else
				{
					if (!mpi_rank)
						std::cout << std::endl << std::endl << "ERROR: " << fitRoutine << " fitting routine not implemented. Plaese, choose between minpack and levmar." << std::endl << std::endl;
					exit(1);
			 
				}

			}
		}	
   	}
			
	/*********************************************
	 * Create the .finished file in project path *
	 *********************************************/
	
	if (!mpi_rank)
		system(("touch "+path+"/.finished").c_str());

	/*********/
	/* Quote */
	/*********/

#ifdef _MPI_
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	if (!mpi_rank) quote();

#ifdef _MPI_
	
	/**************** 
	 * Finalize MPI *
	 ****************/
	
	if ((mpi_err = MPI_Finalize()) != MPI_SUCCESS){
		if (!mpi_rank) std::cout << "\n\nCOOPS ERROR : Cannot finalize MPI environment.\n\n";
		return 1;
	}
	
#endif

if(!mpi_rank) system("date");
    return 0; 
}


/* TIMING ROUTINE */

double getcputime(void)
{
	struct timeval tim;
	struct rusage ru;        
	getrusage(RUSAGE_SELF, &ru);        
	tim=ru.ru_utime;        
	double t=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;        
	tim=ru.ru_stime;        
	t+=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;        
	return t;
}


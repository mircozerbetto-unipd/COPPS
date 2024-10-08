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
 Name        : copps.h
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : Main header of C++OPPS
 ============================================================================
 */

#ifndef COPPS_H_
#define COPPS_H_

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <math.h>
#include <limits>

#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>

#include "prep.h"

#ifdef _MPI_
#include "mpi.h"
#endif

#include "constants.h"
#include "s3j.h"
#include "tensor.h"
#include "physics.h"
#include "wigner.h"
#include "basis.h"
#include "stvec.h"
#include "matrix.h"
#include "lanczos.h"
#include "relax.h"
#include "experimental.h"
#include "euler.h"

#include "cminpack_jac.h"
#include "lm.h"

#include "types.h"

#define MAX_BLOCK_FIT_CYCLES 4 // Runs 2 fits per block (a total of 4: block1 - block2 - block1 - block2)

extern bool cdy1, cdz1, cdy2, cdz2;
extern bool ratio22, ratio42, ratio44;
extern int mpi_rank, mpi_ntasks;
extern int nrows, global_row;
extern int fitStep, comp;
extern int OUT_CONTROL;
extern double dof;
extern double oldD01,oldD02;
extern std::string dynModel;
extern std::ostringstream fitout;
extern physics data;
extern s3j symbols3j;
extern basis basisFunctions;
extern stvec Tkk1;
extern matrix mat;
extern lanczos lcz;
extern relax rel;
extern experimental expdata;
extern std::string path, project;
extern euler eul;
extern int nOrderParametersStopChecks;
extern double oldS20, oldS22;
extern int fn;

/* POWELL F77 interface */
extern int powell_n;
extern int powell_np;
extern int powell_m;
extern double *powell_fvec;
extern double *powell_otherData;
extern double *powell_sim;

extern "C"{
	void powell_(double p[], double xi[], int *n, int *np, double *ftol, double *atol, int *iter, double *fret, int *iflag);
}

/* Observables.cpp routines */
int observables_minpack(void *p, int m, int n, const double *par, double *fvec, int iflag);
void observables_levmar (double *, double *, int, int, void *);
extern "C"{
	double func_(double *);
}
int updatePhysicalData(const double *, int, int, void *, potential *);
void updateObjects(potential);
double calculateChiSquare(dvvector, double *, double *);
bool outputData(dvvector, double *, double, bool);
void checkOrderParametersConvergence(double, double, double);
int StoC(void *p, int n, const double *x, double *fout, int iflag);

double getcputime(void);
void quote(void);

#endif /*COPPS_H_*/

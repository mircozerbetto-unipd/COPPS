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
 Name        : s3j.h
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto
 Description : header of s3j class
 ============================================================================
 */
#ifndef S3J_H_
#define S3J_H_

#include <iostream>
#include <math.h>
#include <vector>

typedef std::vector<long double> ldvector;

#define S3J_EQUAL(a,b)		(fabs((a)-(b))<S3J_0)
#define S3J_MAX(a,b,c,ris)	(((a)>(b)?(ris=(a)):(ris=(b)))>(c)?ris:(ris=(c)))
#define S3J_MIN(a,b,c,ris)	(((a)<(b)?(ris=(a)):(ris=(b)))<(c)?ris:(ris=(c)))

class s3j {
	
	public:
		
		s3j();
		virtual ~s3j();
		void init(int);
		void init(int,int);
		
		int getNSymbolsStored(void);
		int getReggePar(void);
		void setReggePar(int);
		void storeSymbols(void);
		long double getTrj(int,int,int,int,int,int);
		long double getLnFac(int);
		
		//MPI
		int rank;
		int getRank(void);
		void setRank(int);
		
	private:
		
		int ReggePar;
		unsigned int n_s3j_stored;
		ldvector lnf;
		long double* s3j_vector;
		struct rocol {int r,c;};
		
		long double trj(long double,long double,long double,long double,long double,long double);
		void ebs (int *v, rocol *w, int n, char order);
		
};

#endif /*S3J_H_*/

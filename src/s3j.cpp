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
 Name        : s3j.cpp
 Author      : Mirco Zerbetto
 Version     : 2.2
 Copyright   : 2008 Mirco Zerbetto 
 Description : Class for calculation and storage of 3j symbols
 ============================================================================
 */
#include "s3j.h"

#include <cstdlib>

// Constructor
s3j::s3j()
{
	s3j_vector = NULL;
	n_s3j_stored = 0;
}

// Deconstructor
s3j::~s3j()
{
#ifdef WRITE_DESTROY_MESSAGE
	std::cout << "TASK " << rank << ": Cleared s3j object" << std::endl;
#endif
}

// Initiators
void s3j::init(int LRegge)
{
	rank = 0;
	std::cout << "TASK " << rank << ": Strarted 3js object" << std::endl;
	std::cout << "TASK " << rank << ": Allocating ln(x!) vector..." << std::endl;
	ReggePar = LRegge;
	if (lnf.size() > 0) lnf.clear();
	lnf.push_back(0.0);
	long double j1=1.0;
	for (int i=1; i<=(6*ReggePar+1); i++) {
	    lnf.push_back(lnf[i-1] + log(j1));
	    j1 += 1.0;
	}
	std::cout << "TASK " << rank << ": ... ln(x!) vector allocated successfully" << std::endl;
	return;
}

void s3j::init(int LRegge, int r)
{
	rank = r;
	std::cout << "TASK " << rank << ": Strarted 3js object" << std::endl;
	std::cout << "TASK " << rank << ": Allocating ln(x!) vector..." << std::endl;
	ReggePar = LRegge;
	if (lnf.size() > 0) lnf.clear();
	lnf.push_back(0.0);
	long double j1=1.0;
	for (int i=1; i<=(6*ReggePar+1); i++) {
	    lnf.push_back(lnf[i-1] + log(j1));
	    j1 += 1.0;
	}
	std::cout << "TASK " << rank << ": ... ln(x!) vector allocated successfully" << std::endl;
	return;
}

// Return the dimensions of the vector
int s3j::getNSymbolsStored(void){
	return n_s3j_stored;
}

// Return the Regge parameter
int s3j::getReggePar(void){
	return ReggePar;
}

// Set the Regge parameter
void s3j::setReggePar(int rp){
	ReggePar = rp;
}

// Store 3j symbols
void s3j::storeSymbols(void){
	
	if (s3j_vector != NULL) delete s3j_vector;
	
	std::cout << "TASK " << rank << ": Storing 3j symbols compatibily with memory." << std::endl;
	
	long double j1,j2,j3,m1,m2,m3;
	int L,X,T,B,S;
	unsigned long int cmax,mem;
	
	mem = 1024000000/16;
	
	/* Determines the length of the global array */
	cmax=(long int)(((double)(ReggePar)/120.0)*(double)(274+ReggePar*(225+ReggePar*(85+ReggePar*(15+ReggePar)))))+1;
		
	/* Allocates the array */
	/* 3j symbols will be allocated to a maximum of about 1 GB of RAM */
	
	if (cmax > mem){
		do {
			ReggePar = ((ReggePar/3)-1)/2 - 1;
			ReggePar = 3*(2*ReggePar+1);
			cmax=(long int)(((double)(ReggePar)/120.0)*(double)(274+ReggePar*(225+ReggePar*(85+ReggePar*(15+ReggePar)))))+1;
		}while (ReggePar > 1 && cmax > mem);
	}
	
	std::cout << "TASK " << rank << ": Starting 3j storage for Regge parameter = " << ReggePar << "..." << std::endl;

	s3j_vector = new long double[cmax];

	/* Stores in memory the 3j symbols */
	
	for (L=0;L<=ReggePar;L++){
		for (X=0;X<=L;X++){
			for (T=0;T<=X;T++){
				for (B=0;B<=T;B++){
					for (S=0;S<=B;S++){
						j1=(long double)(X+L+B-T)*0.5;
						m1=(long double)(L+B-T-X)*0.5;
						j2=(long double)(S+X-T+B)*0.5;
						m2=(long double)(S+X-T-B)*0.5;
						j3=(long double)(S+L)*0.5;
						m3=(long double)(2*T-S-L)*0.5;
						s3j_vector[n_s3j_stored]=trj(j1,j2,j3,m1,m2,m3);
						n_s3j_stored++;
					}
				}
			}
		}
	}
	
	std::cout << "TASK " << rank << ": ... storage completed successfully" << std::endl;
	
	return;
}

// Calculate the 3j symbol
long double s3j::trj(long double j1, long double j2, long double j3, long double m1, long double m2, long double m3){
	
	/*       Evaluates 3j symbol
	**
	** Author: Riccardo Gusmeroli (web.address@libero.it)
	**
	** Notes: 
	**     - defining S3J_TEST enables the compilation of a very small test suite.
	**     - the maximum allowed factorial is S3J_MAX_FACT (currently 25!). 
	**
	**
	** This program is free software; you can redistribute it and/or  
	** modify it under the terms of the GNU General Public License 
	** as published by the Free Software Foundation; either version 2  
	** of the License, or (at your option) any later version.
	** This program is distributed in the hope that it will be useful, 
	** but WITHOUT ANY WARRANTY; without even the implied warranty of 
	** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
	** GNU General Public License for more details.
	** You should have received a copy of the GNU General Public License 
	** along with this program; if not, write to the Free Software 
	** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA. */
	
	int k, kmin, kmax;
	int jpm1, jmm1, jpm2, jmm2, jpm3, jmm3;
	int j1pj2mj3, j3mj2pm1, j3mj1mm2;
	long double ris, mult;

	jpm1 = (int)(j1+m1);
	jpm2 = (int)(j2+m2);
	jpm3 = (int)(j3+m3);
	jmm1 = (int)(j1-m1);
	jmm2 = (int)(j2-m2);
	jmm3 = (int)(j3-m3);

	/* delta(m1+m2+m3,0) */
	if ((jpm1 - jmm1 + jpm2 - jmm2 + jpm3 - jmm3) != 0) return 0.0;

	if (m1 == 0.0 && m2 == 0.0 && m3 == 0.0) {
	  k = (int)(j1 + j2 + j3);
	  if (k%2 != 0) return 0.0;
	}   	 	

	/* j1+j2-j3 = (j1+j2-j3) + (m1+m2+m3) = jpm1+jpm2-jmm3 */
	j1pj2mj3 = jpm1 + jpm2 - jmm3;
	
	/* j3-j2+m1 = (j3-j2+m1) - (m1+m2+m3) = jmm3-jpm2 */
	j3mj2pm1 = jmm3 - jpm2;
	  
	/* j3-j1-m2 = (j3-j1-m2) + (m1+m2+m3) = jpm3-jmm1 */
	j3mj1mm2 = jpm3 - jmm1;
	  
	S3J_MAX(-j3mj2pm1, -j3mj1mm2,   0,         kmin);
	S3J_MIN(j1pj2mj3,  jmm1,       jpm2,       kmax);
	if (kmin > kmax) return 0.0;	
	  
	ris=0.0;
	mult=(kmin%2 == 0 ? 1.0 : -1.0);
	
	for (k=kmin; k<=kmax; k++) {
	  ris+=mult*exp(-(lnf[k]+lnf[j1pj2mj3-k]+lnf[jmm1-k]+lnf[jpm2-k]+lnf[j3mj2pm1+k]+lnf[j3mj1mm2+k]));
	  mult=-mult;
	}
	
	/* (-1)^(j1-j2-m3)=(-1)^(j1-j2-m3+m1+m2+m3)=(-1)^(jpm1-jmm2) */
	if ((jpm1-jmm2)%2 !=0 ) ris = -ris;
	  
	ris*=exp(0.5*(lnf[j1pj2mj3]+lnf[jpm1-jmm2+jpm3]+lnf[-jmm1+jpm2+jpm3]+lnf[jpm1]+lnf[jpm2]+lnf[jpm3]+lnf[jmm1]+lnf[jmm2]+lnf[jmm3]-lnf[jpm1+jpm2+jpm3+1]));
		
	return ris;

}

// Gets the 3js symbol corresponding to j1,j2,j3,m1,m2,m3
long double s3j::getTrj(int j1, int j2, int j3, int m1, int m2, int m3){
	  int R[9],MR[3][3];
	  int i1,i2;
	  long int L,X,T,B,S;
	  unsigned long int c;
	  rocol rcr[9],Sc,Lc;
	  long double tjs=0.0,segno=1.0;

	  if (abs(j1-j2)>j3 || abs(j2-j3)>j1 || abs(j1-j3)>j2) return (tjs);
	  if ((m1+m2+m3)!=0) return(tjs);
	  if (j1<abs(m1) || j2<abs(m2) || j3<abs(m3)) return (tjs);

	  MR[0][0]=R[0]=-j1+j2+j3;   rcr[0].r=0; rcr[0].c=0;
	  MR[0][1]= R[1]=j1-j2+j3;   rcr[1].r=0; rcr[1].c=1;
	  MR[0][2]= R[2]=j1+j2-j3;   rcr[2].r=0; rcr[2].c=2;
	  MR[1][0]= R[3]=j1-m1;      rcr[3].r=1; rcr[3].c=0;
	  MR[1][1]= R[4]=j2-m2;      rcr[4].r=1; rcr[4].c=1;
	  MR[1][2]= R[5]=j3-m3;      rcr[5].r=1; rcr[5].c=2;
	  MR[2][0]= R[6]=j1+m1;      rcr[6].r=2; rcr[6].c=0;
	  MR[2][1]= R[7]=j2+m2;      rcr[7].r=2; rcr[7].c=1;
	  MR[2][2]= R[8]=j3+m3;      rcr[8].r=2; rcr[8].c=2;

	  ebs(R,rcr,9,'d');
	  
	  L=(long int)R[0];
	  S=(long int)R[8];
	  
	  for (i1=8;i1>=0;i1--){
	    for (i2=0;i2<=8;i2++){
	      if (R[i1]==S && R[i2]==L && rcr[i1].r==rcr[i2].r){
		S=(long int)R[i1]; L=(long int)R[i2];
		Sc.r=rcr[i1].r; Sc.c=rcr[i1].c;
		Lc.r=rcr[i2].r; Lc.c=rcr[i2].c;
		goto _10;
	      }
	    }
	  }

	 for (i1=8;i1>=0;i1--){
	    for (i2=0;i2<=8;i2++){
	      if (R[i1]==S && R[i2]==L && rcr[i1].c==rcr[i2].c){
		S=(long int)R[i1]; L=(long int)R[i2];
		Sc.r=rcr[i1].r; Sc.c=rcr[i1].c;
		Lc.r=rcr[i2].r; Lc.c=rcr[i2].c;
		goto _20;
	      }
	    }
	  }


	 _10:

	  X=(Sc.r==1 ? MR[0][Sc.c] : MR[1][Sc.c]);
	  B=(Lc.r==1 ? MR[0][Lc.c] : MR[1][Lc.c]);
	  T=X+B+L+S-(long int)(j1+j2+j3);
	  if ((j1+j2+j3)%2!=0){
	    if (Sc.c<Lc.c) i1=Sc.c+Lc.c-1;
	    else i1=Sc.c-Lc.c;
	    i1+=(Sc.r==0 ? 0 : 1);
	    segno=(i1%2==0 ? 1.0 : -1.0);
	  }
	  goto _end;

	_20:

	  X=(Sc.c==1 ? (long int)MR[Sc.r][0] : (long int)MR[Sc.r][1]);
	  B=(Lc.c==1 ? (long int)MR[Lc.r][0] : (long int)MR[Lc.r][1]);
	  T=X+B+L+S-(long int)(j1+j2+j3);
	  if ((j1+j2+j3)%2!=0){
	    if (Sc.r<Lc.r) i1=Sc.r+Lc.r-1;
	    else i1=Sc.r-Lc.r;
	    i1+=(Sc.c==0 ? 0 : 1);
	    segno=(i1%2==0 ? 1.0 : -1.0);
	  }

	 _end:
	 
	  c=L*(24+L*(50+L*(35+L*(10+L))))/120+X*(6+X*(11+X*(6+X)))/24+T*(2+T*(3+T))/6+B*(B+1)/2+S;
	  if (c < n_s3j_stored)
	    return ((segno*s3j_vector[c]));
	  return (trj(j1,j2,j3,m1,m2,m3));

}

// Efficient bubble sorting algorithm
void s3j::ebs (int *v, rocol *w, int n, char order){
	
	int i=n-1, j, temp1, swapped=1, wcheck=1;
	rocol temp2;
	
	if (v==NULL) return;
	if (w==NULL) wcheck=0;
	
	switch (order){
	
	case 'a':{
		while ((swapped==1) && (i >= 1)){
			swapped = 0;
			for (j = 0; j <= (i-1); j++)
				if (v[j] > v[j+1]){
					temp1=v[j];
					v[j]=v[j+1];
					v[j+1]=temp1;
					if (wcheck==1){
						temp2=w[j];
						w[j]=w[j+1];
						w[j+1]=temp2;
					}
					swapped = 1; // swapping happened
				}
			i--;
			}
		}
	
	  case 'd':{ 
	    while ((swapped==1) && (i >= 1)){
	    	swapped = 0;
	    	for (j = 0; j <= (i-1); j++)
	    		if (v[j] < v[j+1]){
	    			temp1=v[j];
	    			v[j]=v[j+1];
	    			v[j+1]=temp1;
	    			if (wcheck==1){
	    				temp2=w[j];
	    				w[j]=w[j+1];
	    				w[j+1]=temp2;
	    			}
	    			swapped = 1; // swapping happened
	    		}
	    	i--;
	    	}
	    }
	  }

	  return;

}

// Get rank of task in MPI calculation
int s3j::getRank(void){
	return rank;
}

// Set rank of task in MPI calculation
void s3j::setRank(int r){
	rank = r;
	return;
}

/********************************************
 * Return the natural logatihm of integer n *
 ********************************************/

long double s3j::getLnFac(int n)
{
	long double f;
	if (n <= 6*ReggePar)
		f = lnf.at(n);
	else
	{
		f = lnf.at(6*ReggePar);
		long double x = (long double)(6*ReggePar + 1); 
		for (int i = 6*ReggePar+1; i <= n; i++) {
			f += log(x);
			x += 1.0;
		}
	}
	return f;	
}

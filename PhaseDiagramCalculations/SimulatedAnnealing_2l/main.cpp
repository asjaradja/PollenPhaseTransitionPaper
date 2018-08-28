#include <iostream>
#include <time.h>
#include <math.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cmath>
#include <random>
#include <fstream>
#include <complex>
#include <iomanip>
#include "nrutil.h"
#include "nrtypes.h"
#include "nr.h"
#include "wigxjpf.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace std;
const double pi = 3.14159265358979;
random_device rd;
const double tau = -1.0;
const double lambda4 = 1.0;

///Set lambda3 to be any value between -1.0 and 1.0//////////
const double lambda3 = 0.0;

////Set l_not to be any value between el_not and el_not +1//////////
double l_not = 3.2;

double* Make6DArray(int arraySizeA, int arraySizeB, int arraySizeC, int arraySizeD, int arraySizeE, int arraySizeF)
{

    int size;
    size=arraySizeA*arraySizeB*arraySizeC*arraySizeD*arraySizeE*arraySizeF;
    double *Array = new double[size];
    return Array;
}



double* Make4DArray(int arraySizeA, int arraySizeB, int arraySizeC, int arraySizeD)
{

    int size;
    size=arraySizeA*arraySizeB*arraySizeC*arraySizeD;
    double *Array = new double[size];
    return Array;
}


//Define size of Gaunt matrices
int Na = 2*l0+3;
int Nb = 3*l0+5;
int Nc = 3*l0+5;
int Nd = 6*l0+9;

int NaNb = (2*l0+3)*(3*l0+5);
int NaNbNc = (2*l0+3)*(3*l0+5)*(3*l0+5);
 

double* Gaunt_matrix_3rd_order = Make6DArray(2,2,2,Na,Na,Na);
double* Gaunt_matrix_4th_order_1 = Make6DArray(2,2,Na,Nb,Nc,Nd);
double* Gaunt_matrix_4th_order_2 = Make6DArray(2,2,Na,Nb,Nc,Nd);


int hs = 4*l0+4;
int z = 4*l0+4;
int zz = (4*l0+4)*(4*l0+4);
int zzz = (4*l0+4)*(4*l0+4)*(4*l0+4);
int bb = 2*l0+2;

double* Alpha_matrix = Make4DArray(z,z,z,z);

namespace NR
{
   void amebsa(Mat_IO_DP &p, Vec_IO_DP &y, Vec_O_DP &pb, DP &yb, const DP ftol, DP funk(Vec_I_DP &x, double *Gaunt_matrix, double *Alpha_matrix), int &iter, const DP temptr);
    DP amotsa(Mat_IO_DP &p, Vec_O_DP &y, Vec_IO_DP &psum, Vec_IO_DP &pb, DP &yb, DP funk(Vec_I_DP &x, double *Gaunt_matrix, double *Alpha_matrix), const int ihi, DP &yhi, const DP fac);
}

//Intialize tables for wigner3j coefficients
void wig_table_init(int max_two_j, int wigner_type);
void wig_temp_init(int max_two_j);        /* Single-threaded */
void wig_thread_temp_init(int max_two_j); /* Multi-threaded. */
double wig3jj(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3);
void wig_temp_free();  /* Per-thread when multi-threaded. */
void wig_table_free();

namespace
{
    inline void get_psum(Mat_I_DP &p, Vec_O_DP &psum)
    {
        int n,m;
        DP sum;
        
        int mpts=p.nrows();
        int ndim=p.ncols();
        for(n=0;n<ndim;n++)
        {
            for(sum=0.0,m=0;m<mpts;m++) sum+=p[m][n];
            psum[n]=sum;
        }
    }
}

double Gaunt(int l, int lp, int lpp, int m, int mp, int mpp)
{
    
    double GauntCoefficient = sqrt(((2*l+1)*(2*lp+1)*(2*lpp+1))/(4*pi))*wig3jj(2*l, 2*lp, 2*lpp, 2*m, 2*mp, 2*mpp)*wig3jj(2*l,2*lp,2*lpp,0,0,0);
    return GauntCoefficient;
}

complex<double>** Make2DArray(int arraySizeR, int arraySizeC)
{
    complex<double> **Array = new complex<double>*[arraySizeR];
    for(int i=0;i<arraySizeR;i++)
	Array[i] = new complex<double>[arraySizeC];
    return Array;
}


double** Make2DRArray(int arraySizeR, int arraySizeC)
{
    double **Array = new double*[arraySizeR];
    for(int i=0;i<arraySizeR;i++)
	Array[i] = new double[arraySizeC];
    return Array;
}

//DERIVATIVE CALCULATION - NUMERICAL APPROXIMATION
double ngrad(Vec_I_DP &x, double *g,  DP funk(Vec_I_DP &, double *Gauntm, double *Alpha))
{
	double gg=0;
	double del=1e-7;

DP ptempl[(4*l0+4)];
DP ptempr[(4*l0+4)];
DP ptempc[(4*l0+4)];
for (int i=0; i<4*l0+4; i++) { ptempl[i]=x[i]; ptempr[i]=x[i]; ptempc[i]=x[i]; }

for(int k=0; k<4*l0+4; k++)
{
        ptempr[k]=ptempr[k]+del;
        ptempl[k]=ptempl[k]-del;

        Vec_I_DP tempr(ptempr,4*l0+4);
        Vec_I_DP templ(ptempl,4*l0+4);

        g[k]=(funk(tempr,Gaunt_matrix_3rd_order, Alpha_matrix)-funk(templ,Gaunt_matrix_3rd_order, Alpha_matrix))/(2*del);
        ptempr[k]=x[k];
        ptempl[k]=x[k];

}


	return gg;
}

//GRADIENT CALCULATION
double grad(Vec_I_DP &x, double *g, double *Gaunt_matrix, double* alp)
{

	double gg=0;
	double t;

    double** re = Make2DRArray(2, int(l0+2));
    double** im = Make2DRArray(2, int(l0+2));

    for(int q=0;q<2;q++) {
	re[q][0] = x[q];
	im[q][0] = 1000000;		}	

    int count=2;
    
    for(int i=0;i<=1;i++)
    {
        for(int j=1;j<=l0+i;j++)
        {
           re[i][j] = x[count];
           im[i][j] = x[(2*l0+1+count)];
           count++;
        }
    }
	re[0][l0+1]=100000;
	im[0][l0+1]=100000;

    complex<double>** c = Make2DArray(2, l0+2);

    for(int q=0;q<=1;q++)
        c[q][0] = complex<double> (x[q],0);

    count=2;

    for(int i=0;i<=1;i++)
    {
        for(int j=1;j<=l0+i;j++)
        {
           c[i][j]=complex<double> (x[count],x[(2*l0+1+count)]);
           count++;
        }
    }
    c[0][(l0+1)]=complex<double> (100000,10000);

for (int d=0; d<2; d++) {

	g[d]=(tau+pow((l_not-(double(l0+d))),2))*re[d][0];

for(int d1=0;d1<2;d1++)
	for(int d2=0;d2<2;d2++)	{ g[d]+=lambda3*Gaunt_matrix[d+d1*2+4*d2+8*(l0+d)+8*Na*(l0+d1)+8*Na*Na*(l0+d2)]*re[d1][0]*re[d2][0]/2.0;
	for(int m=1;m<=l0+d1;m++) {
	t=real(c[d1][m]*c[d2][m]);
	g[d]+=lambda3*double(1-2*(m%2))*Gaunt_matrix[d+d1*2+d2*4+8*(l0+d+m)+8*Na*(l0+d1-m)+8*Na*Na*(l0+d2)]*t;
						}
				}
for(int d1=0;d1<2;d1++)
	for(int d2=0;d2<2;d2++)
		for(int d3=0;d3<2;d3++)  {
 g[d]+=lambda4*alp[(l0+d1*bb)+(l0+d2*bb)*z+(l0+d3*bb)*zz+(l0+d*bb)*zzz]*re[d1][0]*re[d2][0]*re[d3][0]/6.0;

	for(int m=1;m<=l0+min(d1,d2);m++) {
	t=lambda4*real(c[d1][m]*conj(c[d2][m]))/6.;
	g[d]+=(alp[d*bb+l0+(d3*bb+l0)*z+(d1*bb+l0+m)*zz+(d2*bb+l0-m)*zzz]*double(1-2*(m%2))+2*alp[d*bb+l0+(d1*bb+l0+m)*z+(d3*bb+l0)*zz+(d2*bb+l0-m)*zzz]+alp[d3*bb+l0+(d*bb+l0)*z+(d1*bb+l0+m)*zz+(d2*bb+l0-m)*zzz]*double(1-2*(m%2))+2*alp[d3*bb+l0+(d1*bb+l0+m)*z+(d*bb+l0)*zz+(d2*bb+l0-m)*zzz])*re[d3][0]*t;

				}

	for(int m=1;m<=l0+d1;m++) {
	for(int mb=1;mb<=l0+d3-m;mb++) {
		t=lambda4*real(c[d1][m]*c[d2][mb]*conj(c[d3][m+mb]))/6.;
		            g[d]+=alp[l0+bb*d+(l0+m+bb*d1)*z+(l0+mb+bb*d2)*zz+(l0-m-mb+bb*d3)*zzz]*double(1-2*(mb%2))*t;
	
	if(mb>=m) {

	g[d]+=t*double(2-int(m==mb))*(alp[d1*bb+l0+m+(d2*bb+l0+mb)*z+(d*bb+l0)*zz+(d3*bb+l0-m-mb)*zzz]+alp[d1*bb+l0+m+(d3*bb+l0-m-mb)*z+(d2*bb+l0+mb)*zz+(d*bb+l0)*zzz]*double(1-2*(m%2)));
			}
	}


	for(int mb=1;mb<=min(m-1,l0+d2);mb++) {
		t=lambda4*real(conj(c[d1][m])*c[d2][mb]*c[d3][m-mb])/6.;
		g[d]+=alp[d*bb+l0+(d1*bb+l0+m)*z+(d2*bb+l0-mb)*zz+(d3*bb+l0+mb-m)*zzz]*t;
	}


	for(int mb=m+1;mb<=l0+d2;mb++) {
		t=lambda4*real(conj(c[d1][m])*c[d2][mb]*conj(c[d3][mb-m]))/6.;
		g[d]+=(alp[d*bb+l0+(d1*bb+l0+m)*z+(d2*bb+l0-mb)*zz+(d3*bb+l0+mb-m)*zzz]*double(1-2*((m+mb)%2))+alp[d1*bb+l0+m+(d2*bb+l0-mb)*z+(d*bb+l0)*zz+(d3*bb+l0+mb-m)*zzz]*double(1-2*(m%2)))*t;
		}

	}
			}

}



for (int d=0;d<2;d++)  {

//Begin calculating dE/dRk
	for(int k=1;k<=l0+d;k++) {


        g[k+1+d*l0]=2*(tau+pow(l_not-(double(l0+d)),2))*re[d][k];


for(int d1=0;d1<2;d1++) 
for(int d2=0;d2<2;d2++) {
// Do cobic term stuff
if(k<=l0+d2)	g[k+1+d*l0]+=2*lambda3*Gaunt_matrix[d+2*d2+4*d1+8*(l0+d+k)+8*Na*(l0-k+d2)+8*Na*Na*(l0+d1)]*double(1-2*(k%2))*re[d1][0]*re[d2][k];
  
      for(int m=1;m<=l0+min(d1,d2)-k;m++) {
                t=re[d1][m]*re[d2][m+k]+im[d1][m]*im[d2][k+m];
        g[k+1+d*l0]+=2*lambda3*(Gaunt_matrix[d1+2*d+4*d2+8*(l0+d1+m)+(l0+d+k)*8*Na+(l0+d2-m-k)*8*Na*Na]+Gaunt_matrix[d2+2*d1+d*4+8*(l0+d2+k+m)+(l0+d1-m)*8*Na+(l0+d-k)*8*Na*Na])*double(1-2*((m+k)%2))*t/3.;
        }
// The m,m-k single loop
        for(int m=k+1;m<=l0+d1;m++) {
                t=re[d1][m]*re[d2][m-k]+im[d1][m]*im[d2][m-k];
        g[k+1+d*l0]+=2*lambda3*Gaunt_matrix[d1+2*d+4*d2+8*(l0+d1+m)+(l0+d-k)*8*Na+(l0+d2+k-m)*8*Na*Na]*double(1-2*(m%2))*t/3.;
        }
// The m,k-m single loop
        for(int m=1;m<=k-1;m++) {
                t=re[d1][m]*re[d2][k-m]+im[d1][m]*im[d2][k-m];
        g[k+1+d*l0]+=2*lambda3*Gaunt_matrix[d+2*d1+4*d2+8*(l0+d+k)+(l0+d1-m)*8*Na+(l0+d2+m-k)*8*Na*Na]*double(1-2*(k%2))*t/3.0;
        if(m<=((int) floor(double(k)/2.))) {
        if(k==2*m) g[k+1+d*l0]+=lambda3*Gaunt_matrix[d1+2*d2+4*d+8*(l0+d1+m)+(l0+d2+k-m)*8*Na+(l0+d-k)*8*Na*Na]*double(1-2*(k%2))*t/3.0;
        else g[k+1+d*l0]+=2*lambda3*Gaunt_matrix[d1+2*d2+d*4+8*(l0+d1+m)+(l0+d2+k-m)*8*Na+(l0+d-k)*8*Na*Na]*double(1-2*(k%2))*t/3.0;
                                        }
                                }

	} 	


	for(int d1=0;d1<2;d1++) 
	for(int d2=0;d2<2;d2++)
	for(int d3=0;d3<2;d3++)  {

// Do quartic term stuff
// The m,m+k single loop
   if(k<=l0+d3)     g[k+1+d*l0]+=lambda4*(2*alp[l0+d1*bb+(l0+k+d*bb)*z+(l0+d2*bb)*zz+(l0+d3*bb-k)*zzz]+alp[l0+d*bb+k+(l0+d3*bb-k)*z+(l0+d1*bb)*zz+(l0+d2*bb)*zzz]*double(1-2*(k%2)))*re[d1][0]*re[d2][0]*re[d3][k]/3.;



        for(int m=1;m<=min(l0+d2-k,l0+d1);m++) {
	if(k+m<=l0+d2) {
                t=re[d1][m]*re[d2][m+k]+im[d1][m]*im[d2][k+m];
        g[k+1+d*l0]+=lambda4*re[d3][0]*(alp[l0+d3*bb+(l0+d*bb+k)*z+(l0+d1*bb+m)*zz+(l0+d2*bb-k-m)*zzz]*double(1-2*(m%2))+alp[l0+d*bb+k+(l0+d1*bb+m)*z+(l0+d2*bb-k-m)*zz+(l0+d3*bb)*zzz]+alp[l0+d3*bb+(l0+d1*bb+m)*z+(l0+d*bb+k)*zz+(l0+d2*bb-k-m)*zzz]*double(1-2*(k%2)))*t/3.;
        if(m>k) g[k+1+d*l0]+=lambda4*re[d3][0]*double(1-2*(k%2))*alp[l0+d2*bb+k+m+(l0+d*bb-k)*z+(l0+d1*bb-m)*zz+(l0+d3*bb)*zzz]*t/3.;
        if(m<k) g[k+1+d*l0]+=lambda4*re[d3][0]*double(1-2*(m%2))*alp[l0+d2*bb+k+m+(l0+d1*bb-m)*z+(l0+d*bb-k)*zz+(l0+d3*bb)*zzz]*t/3.;
       if(m==k) g[k+1+d*l0]+=lambda4*re[d3][0]*double(1-2*(k%2))*alp[l0+d2*bb+k+k+(l0+d*bb-k)*z+(l0+d1*bb-k)*zz+(l0+d3*bb)*zzz]*t/6.;
       if(m==k) g[k+1+d*l0]+=lambda4*re[d3][0]*double(1-2*(k%2))*alp[l0+d2*bb+k+k+(l0+d1*bb-k)*z+(l0+d*bb-k)*zz+(l0+d3*bb)*zzz]*t/6.;
 }
        }
// The m,m-k single loop
        for(int m=k+1;m<=min(l0+d1,l0+d1+k);m++) {
                t=re[d1][m]*re[d2][m-k]+im[d1][m]*im[d2][m-k];
        g[k+1+d*l0]+=lambda4*re[d3][0]*(alp[l0+d3*bb+(l0+d*bb+k)*z+(l0+d1*bb-m)*zz+(l0+d2*bb-k+m)*zzz]*double(1-2*((k+m)%2))+2*alp[l0+d3*bb+(l0+d1*bb+m)*z+(l0+d*bb-k)*zz+(l0+d2*bb-m+k)*zzz]+alp[l0+d*bb+k+(l0+d1*bb-m)*z+(l0+d3*bb)*zz+(l0+d2*bb+m-k)*zzz]*double(1-2*(k%2)))*t/6.;
        }
// The m,k-m single loop
        for(int m=1;m<=min(k-1,l0+d1);m++) {
                t=re[d1][m]*re[d2][k-m]-im[d1][m]*im[d2][k-m];
        g[k+1+d*l0]+=lambda4*re[d3][0]*(alp[l0+d3*bb+(l0+d*bb-k)*z+(l0+d1*bb+m)*zz+(l0+d2*bb+k-m)*zzz]+2*alp[l0+d3*bb+(l0+d1*bb+m)*z+(l0+d2*bb+k-m)*zz+(l0+d*bb-k)*zzz]*double(1-2*((m+k)%2))+alp[l0+d1*bb+m+(l0+d*bb-k)*z+(l0+d3*bb)*zz+(l0+d2*bb+k-m)*zzz]*double(1-2*(m%2)))*t/6.0;
        if(m<=((int) floor(double(k)/2.))) {
        if(k==2*m) {
        g[k+1+d*l0]+=lambda4*re[d3][0]*(alp[l0+d1*bb+m+(l0+d2*bb+k-m)*z+(l0+d*bb-k)*zz+(l0+d3*bb)*zzz]+alp[l0+d*bb+k+(l0+d1*bb-m)*z+(l0+d2*bb+m-k)*zz+(l0+d3*bb)*zzz]*double(1-2*(m%2)))*t/6.;
                } else
                {
      g[k+1+d*l0]+=2*lambda4*re[d3][0]*(alp[l0+d1*bb+m+(l0+d2*bb+k-m)*z+(l0+d*bb-k)*zz+(l0+d3*bb)*zzz]+alp[l0+d*bb+k+(l0+d1*bb-m)*z+(l0+d2*bb+m-k)*zz+(l0+d3*bb)*zzz]*double(1-2*(m%2)))*t/6.;
                }
                                        }
                                }
// Start the double loops
// The double loops from the new gradient
// m+k+mb
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=1;mb<=l0+d2;mb++) {
	if(m+k+mb<=l0+d3) {


		t=lambda4*((re[d1][m]*re[d2][mb]-im[d1][m]*im[d2][mb])*re[d3][m+k+mb]+(re[d1][m]*im[d2][mb]+im[d1][m]*re[d2][mb])*im[d3][m+k+mb])/6.;
	g[k+1+d*l0]+=alp[l0+d1*bb+m+(l0+d*bb+k)*z+(l0+d2*bb+mb)*zz+(l0+d3*bb-m-mb-k)*zzz]*double(1-2*(mb%2))*t;
	if(mb>m) g[k+1+d*l0]+=2*(alp[l0+d1*bb-m+(l0+d3*bb+m+mb+k)*z+(l0+d2*bb-mb)*zz+(l0+d*bb-k)*zzz]*double(1-2*(m%2))+alp[l0+d1*bb+m+(l0+d2*bb+mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb-m-mb-k)*zzz]*double(1-2*(k%2)))*t;
	if(mb==m)  g[k+1+d*l0]+=(alp[l0+d1*bb-m+(l0+d3*bb+m+mb+k)*z+(l0+d2*bb-mb)*zz+(l0+d*bb-k)*zzz]*double(1-2*(m%2))+alp[l0+d1*bb+m+(l0+d2*bb+mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb-m-mb-k)*zzz]*double(1-2*(k%2)))*t/2.;
	if(mb==m)  g[k+1+d*l0]+=(alp[l0+d2*bb-m+(l0+d3*bb+m+mb+k)*z+(l0+d1*bb-mb)*zz+(l0+d*bb-k)*zzz]*double(1-2*(m%2))+alp[l0+d2*bb+m+(l0+d1*bb+mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb-m-mb-k)*zzz]*double(1-2*(k%2)))*t/2.;
	}	}	}

// k-m-mb
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=m;mb<=l0+d2;mb++) {
	if((k-m-mb)>=1) {
		          t=lambda4*((re[d1][m]*re[d2][mb]-im[d1][m]*im[d2][mb])*re[d3][k-m-mb]-(re[d1][m]*im[d2][mb]+im[d1][m]*re[d2][mb])*im[d3][k-m-mb])/6.0;
		if(m!=mb) t=lambda4*((re[d1][m]*re[d2][mb]-im[d1][m]*im[d2][mb])*re[d3][k-m-mb]-(re[d1][m]*im[d2][mb]+im[d1][m]*re[d2][mb])*im[d3][k-m-mb])/3.0;
		g[k+1+d*l0]+=(alp[l0+d1*bb+m+(l0+d*bb-k)*z+(l0+d2*bb+mb)*zz+(l0+d3*bb+k-m-mb)*zzz]*double(1-2*(m%2))+alp[l0+d1*bb+m+(l0+d2*bb+mb)*z+(l0+d*bb-k)*zz+(l0+d3*bb+k-m-mb)*zzz]*double(1-2*((k+m+mb)%2)))*t;
			}   }	}	


// k+m-mb 
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=1;mb<=min(k+m-1,l0+d2);mb++) {
	if((k+m-mb)<=l0+d3 && (k+m-mb)>=1) {

	    t=lambda4*((re[d1][m]*re[d2][mb]+im[d1][m]*im[d2][mb])*re[d3][k+m-mb]+(im[d1][m]*re[d2][mb]-im[d2][mb]*re[d1][m])*im[d3][k+m-mb])/6.0;
	if(m>k) g[k+1+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb+m)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-m+mb)*zzz]+2*alp[l0+d*bb+k+(l0+d2*bb-mb)*z+(l0+d1*bb+m)*zz+(l0+d3*bb+mb-k-m)*zzz]*double(1-2*((m+mb)%2)))*t;
	if(m<k) g[k+1+d*l0]+=(alp[l0+d1*bb+m+(l0+d*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-m+mb)*zzz]+2*alp[l0+d1*bb+m+(l0+d2*bb-mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb+mb-k-m)*zzz]*double(1-2*((k+mb)%2)))*t;
      if(m==k)  g[k+1+d*l0]+=(alp[l0+d1*bb+k+(l0+d*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-k+mb)*zzz]+2*alp[l0+d1*bb+k+(l0+d2*bb-mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb+mb-k-k)*zzz]*double(1-2*((k+mb)%2)))*t/2.;
      if(m==k)  g[k+1+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-k+mb)*zzz]+2*alp[l0+d*bb+k+(l0+d2*bb-mb)*z+(l0+d1*bb+k)*zz+(l0+d3*bb+mb-k-k)*zzz]*double(1-2*((k+mb)%2)))*t/2.;

			}
	}
	}

// m+mb-k
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=m;mb<=l0+d2;mb++) {
	if((m+mb-k)<=l0+d3 && (m+mb-k)>=1) {

   	   t=lambda4*((re[d1][m]*re[d2][mb]-im[d1][m]*im[d2][mb])*re[d3][m+mb-k]+(im[d1][m]*re[d2][mb]+re[d1][m]*im[d2][mb])*im[d3][m+mb-k])/6.0;
 if(m!=mb) t=lambda4*((re[d1][m]*re[d2][mb]-im[d1][m]*im[d2][mb])*re[d3][m+mb-k]+(im[d1][m]*re[d2][mb]+re[d1][m]*im[d2][mb])*im[d3][m+mb-k])/3.0;

	g[k+1+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb-m)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k+m+mb)*zzz]*double(1-2*((k+mb)%2))+alp[l0+d1*bb+m+(l0+d3*bb+k-m-mb)*z+(l0+d2*bb+mb)*zz+(l0+d*bb-k)*zzz]*double(1-2*((k+m)%2))+alp[l0+d1*bb+m+(l0+d2*bb+mb)*z+(l0+d*bb-k)*zz+(l0+d3*bb+k-m-mb)*zzz])*t;
			}
	}
	}


// mb-m-k
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=k+m+1;mb<=l0+d2;mb++) {
		if((mb-m-k)>=1) {
         t=lambda4*((re[d1][m]*re[d2][mb]+im[d1][m]*im[d2][mb])*re[d3][mb-m-k]+(re[d1][m]*im[d2][mb]-im[d1][m]*re[d2][mb])*im[d3][mb-m-k])/6.;
 
	if(m>k) g[k+1+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb+m)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-m+mb)*zzz]*double(1-2*((k+m+mb)%2))+2*alp[l0+d*bb+k+(l0+d2*bb-mb)*z+(l0+d1*bb+m)*zz+(l0+d3*bb+mb-m-k)*zzz]*double(1-2*(k%2)))*t;
	if(m<k) g[k+1+d*l0]+=(alp[l0+d1*bb+m+(l0+d*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-m+mb)*zzz]*double(1-2*((k+m+mb)%2))+2*alp[l0+d1*bb+m+(l0+d2*bb-mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb+mb-m-k)*zzz]*double(1-2*(m%2)))*t;
       if(m==k) g[k+1+d*l0]+=(alp[l0+d1*bb+k+(l0+d*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-k+mb)*zzz]*double(1-2*((k+k+mb)%2))+2*alp[l0+d1*bb+k+(l0+d2*bb-mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb+mb-k-k)*zzz]*double(1-2*(k%2)))*t/2.;
       if(m==k) g[k+1+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-k+mb)*zzz]*double(1-2*((k+k+mb)%2))+2*alp[l0+d*bb+k+(l0+d2*bb-mb)*z+(l0+d1*bb+k)*zz+(l0+d3*bb+mb-k-k)*zzz]*double(1-2*(k%2)))*t/2.;

				}
	}
	}


		}

		} // end of the Rk loop

}


int bbb=2*l0+1;

for (int d=0;d<2;d++)  {

//Begin calculating dE/dIk
	for(int k=1;k<=l0+d;k++) {


        g[k+bb+d*l0]=2*(tau+pow(l_not-(double(l0+d)),2))*im[d][k];


for(int d1=0;d1<2;d1++) 
for(int d2=0;d2<2;d2++) {
// Do cobic term stuff
if(k<=l0+d2)	g[k+bb+d*l0]+=2*lambda3*Gaunt_matrix[d+2*d2+4*d1+8*(l0+d+k)+8*Na*(l0-k+d2)+8*Na*Na*(l0+d1)]*double(1-2*(k%2))*re[d1][0]*im[d2][k];
  
      for(int m=1;m<=l0+min(d1,d2)-k;m++) {
                t=re[d1][m]*im[d2][m+k]-im[d1][m]*re[d2][k+m];
        g[k+bb+d*l0]+=2*lambda3*(Gaunt_matrix[d1+2*d+4*d2+8*(l0+d1+m)+(l0+d+k)*8*Na+(l0+d2-m-k)*8*Na*Na]+Gaunt_matrix[d2+2*d1+d*4+8*(l0+d2+k+m)+(l0+d1-m)*8*Na+(l0+d-k)*8*Na*Na])*double(1-2*((m+k)%2))*t/3.;
        }
// The m,m-k single loop
        for(int m=k+1;m<=l0+d1;m++) {
                t=im[d1][m]*re[d2][m-k]-re[d1][m]*im[d2][m-k];
        g[k+bb+d*l0]+=2*lambda3*Gaunt_matrix[d1+2*d+4*d2+8*(l0+d1+m)+(l0+d-k)*8*Na+(l0+d2+k-m)*8*Na*Na]*double(1-2*(m%2))*t/3.;
        }
// The m,k-m single loop
        for(int m=1;m<=k-1;m++) {
                t=re[d1][m]*im[d2][k-m]+im[d1][m]*re[d2][k-m];
        g[k+bb+d*l0]+=2*lambda3*Gaunt_matrix[d+2*d1+4*d2+8*(l0+d+k)+(l0+d1-m)*8*Na+(l0+d2+m-k)*8*Na*Na]*double(1-2*(k%2))*t/3.0;
        if(m<=((int) floor(double(k)/2.))) {
        if(k==2*m) g[k+bb+d*l0]+=lambda3*Gaunt_matrix[d1+2*d2+4*d+8*(l0+d1+m)+(l0+d2+k-m)*8*Na+(l0+d-k)*8*Na*Na]*double(1-2*(k%2))*t/3.0;
        else g[k+bb+d*l0]+=2*lambda3*Gaunt_matrix[d1+2*d2+d*4+8*(l0+d1+m)+(l0+d2+k-m)*8*Na+(l0+d-k)*8*Na*Na]*double(1-2*(k%2))*t/3.0;
                                        }
                                }

	} 	


	for(int d1=0;d1<2;d1++) 
	for(int d2=0;d2<2;d2++)
	for(int d3=0;d3<2;d3++)  {

// Do quartic term stuff
// The m,m+k single loop
   if(k<=l0+d3)     g[k+bb+d*l0]+=lambda4*(2*alp[l0+d1*bb+(l0+k+d*bb)*z+(l0+d2*bb)*zz+(l0+d3*bb-k)*zzz]+alp[l0+d*bb+k+(l0+d3*bb-k)*z+(l0+d1*bb)*zz+(l0+d2*bb)*zzz]*double(1-2*(k%2)))*re[d1][0]*re[d2][0]*im[d3][k]/3.;



        for(int m=1;m<=min(l0+d2-k,l0+d1);m++) {
	if(k+m<=l0+d2) {
                t=re[d1][m]*im[d2][m+k]-im[d1][m]*re[d2][k+m];
        g[k+bb+d*l0]+=lambda4*re[d3][0]*(alp[l0+d3*bb+(l0+d*bb+k)*z+(l0+d1*bb+m)*zz+(l0+d2*bb-k-m)*zzz]*double(1-2*(m%2))+alp[l0+d*bb+k+(l0+d1*bb+m)*z+(l0+d2*bb-k-m)*zz+(l0+d3*bb)*zzz]+alp[l0+d3*bb+(l0+d1*bb+m)*z+(l0+d*bb+k)*zz+(l0+d2*bb-k-m)*zzz]*double(1-2*(k%2)))*t/3.;
        if(m>k) g[k+bb+d*l0]+=lambda4*re[d3][0]*double(1-2*(k%2))*alp[l0+d2*bb+k+m+(l0+d*bb-k)*z+(l0+d1*bb-m)*zz+(l0+d3*bb)*zzz]*t/3.;
        if(m<k) g[k+bb+d*l0]+=lambda4*re[d3][0]*double(1-2*(m%2))*alp[l0+d2*bb+k+m+(l0+d1*bb-m)*z+(l0+d*bb-k)*zz+(l0+d3*bb)*zzz]*t/3.;
       if(m==k) g[k+bb+d*l0]+=lambda4*re[d3][0]*double(1-2*(k%2))*alp[l0+d2*bb+k+k+(l0+d*bb-k)*z+(l0+d1*bb-k)*zz+(l0+d3*bb)*zzz]*t/6.;
       if(m==k) g[k+bb+d*l0]+=lambda4*re[d3][0]*double(1-2*(k%2))*alp[l0+d2*bb+k+k+(l0+d1*bb-k)*z+(l0+d*bb-k)*zz+(l0+d3*bb)*zzz]*t/6.; }
        }
// The m,m-k single loop
        for(int m=k+1;m<=min(l0+d1,l0+d1+k);m++) {
                t=im[d1][m]*re[d2][m-k]-re[d1][m]*im[d2][m-k];
        g[k+bb+d*l0]+=lambda4*re[d3][0]*(alp[l0+d3*bb+(l0+d*bb+k)*z+(l0+d1*bb-m)*zz+(l0+d2*bb-k+m)*zzz]*double(1-2*((k+m)%2))+2*alp[l0+d3*bb+(l0+d1*bb+m)*z+(l0+d*bb-k)*zz+(l0+d2*bb-m+k)*zzz]+alp[l0+d*bb+k+(l0+d1*bb-m)*z+(l0+d3*bb)*zz+(l0+d2*bb+m-k)*zzz]*double(1-2*(k%2)))*t/6.;
        }
// The m,k-m single loop
        for(int m=1;m<=min(k-1,l0+d1);m++) {
                t=re[d1][m]*im[d2][k-m]+im[d1][m]*re[d2][k-m];
        g[k+bb+d*l0]+=lambda4*re[d3][0]*(alp[l0+d3*bb+(l0+d*bb-k)*z+(l0+d1*bb+m)*zz+(l0+d2*bb+k-m)*zzz]+2*alp[l0+d3*bb+(l0+d1*bb+m)*z+(l0+d2*bb+k-m)*zz+(l0+d*bb-k)*zzz]*double(1-2*((m+k)%2))+alp[l0+d1*bb+m+(l0+d*bb-k)*z+(l0+d3*bb)*zz+(l0+d2*bb+k-m)*zzz]*double(1-2*(m%2)))*t/6.0;
        if(m<=((int) floor(double(k)/2.))) {
        if(k==2*m) {
        g[k+bb+d*l0]+=lambda4*re[d3][0]*(alp[l0+d1*bb+m+(l0+d2*bb+k-m)*z+(l0+d*bb-k)*zz+(l0+d3*bb)*zzz]+alp[l0+d*bb+k+(l0+d1*bb-m)*z+(l0+d2*bb+m-k)*zz+(l0+d3*bb)*zzz]*double(1-2*(m%2)))*t/6.;
                } else
                {
      g[k+bb+d*l0]+=2*lambda4*re[d3][0]*(alp[l0+d1*bb+m+(l0+d2*bb+k-m)*z+(l0+d*bb-k)*zz+(l0+d3*bb)*zzz]+alp[l0+d*bb+k+(l0+d1*bb-m)*z+(l0+d2*bb+m-k)*zz+(l0+d3*bb)*zzz]*double(1-2*(m%2)))*t/6.;
                }
                                        }
                                }
// Start the double loops
// The double loops from the new gradient
// m+k+mb
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=1;mb<=l0+d2;mb++) {
	if(m+k+mb<=l0+d3) {


		t=lambda4*((re[d1][m]*re[d2][mb]-im[d1][m]*im[d2][mb])*im[d3][m+k+mb]-(re[d1][m]*im[d2][mb]+im[d1][m]*re[d2][mb])*re[d3][m+k+mb])/6.;
	g[k+bb+d*l0]+=alp[l0+d1*bb+m+(l0+d*bb+k)*z+(l0+d2*bb+mb)*zz+(l0+d3*bb-m-mb-k)*zzz]*double(1-2*(mb%2))*t;
	if(mb>m) g[k+bb+d*l0]+=2*(alp[l0+d1*bb-m+(l0+d3*bb+m+mb+k)*z+(l0+d2*bb-mb)*zz+(l0+d*bb-k)*zzz]*double(1-2*(m%2))+alp[l0+d1*bb+m+(l0+d2*bb+mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb-m-mb-k)*zzz]*double(1-2*(k%2)))*t;
	if(mb==m)  g[k+bb+d*l0]+=(alp[l0+d1*bb-m+(l0+d3*bb+m+mb+k)*z+(l0+d2*bb-mb)*zz+(l0+d*bb-k)*zzz]*double(1-2*(m%2))+alp[l0+d1*bb+m+(l0+d2*bb+mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb-m-mb-k)*zzz]*double(1-2*(k%2)))*t/2.;
	if(mb==m)  g[k+bb+d*l0]+=(alp[l0+d2*bb-m+(l0+d3*bb+m+mb+k)*z+(l0+d1*bb-mb)*zz+(l0+d*bb-k)*zzz]*double(1-2*(m%2))+alp[l0+d2*bb+m+(l0+d1*bb+mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb-m-mb-k)*zzz]*double(1-2*(k%2)))*t/2.;
	}	}	}

// k-m-mb
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=m;mb<=l0+d2;mb++) {
	if((k-m-mb)>=1) {
		          t=lambda4*((re[d1][m]*re[d2][mb]-im[d1][m]*im[d2][mb])*im[d3][k-m-mb]+(re[d1][m]*im[d2][mb]+im[d1][m]*re[d2][mb])*re[d3][k-m-mb])/6.0;
		if(m!=mb) t=lambda4*((re[d1][m]*re[d2][mb]-im[d1][m]*im[d2][mb])*im[d3][k-m-mb]+(re[d1][m]*im[d2][mb]+im[d1][m]*re[d2][mb])*re[d3][k-m-mb])/3.0;
		g[k+bb+d*l0]+=(alp[l0+d1*bb+m+(l0+d*bb-k)*z+(l0+d2*bb+mb)*zz+(l0+d3*bb+k-m-mb)*zzz]*double(1-2*(m%2))+alp[l0+d1*bb+m+(l0+d2*bb+mb)*z+(l0+d*bb-k)*zz+(l0+d3*bb+k-m-mb)*zzz]*double(1-2*((k+m+mb)%2)))*t;
			}   }	}	


// k+m-mb 
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=1;mb<=min(k+m-1,l0+d2);mb++) {
	if((k+m-mb)<=l0+d3 && (k+m-mb)>=1) {

	    t=lambda4*((re[d1][m]*re[d2][mb]+im[d1][m]*im[d2][mb])*im[d3][k+m-mb]+(re[d1][m]*im[d2][mb]-re[d2][mb]*im[d1][m])*re[d3][k+m-mb])/6.0;
	if(m>k) g[k+bb+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb+m)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-m+mb)*zzz]+2*alp[l0+d*bb+k+(l0+d2*bb-mb)*z+(l0+d1*bb+m)*zz+(l0+d3*bb+mb-k-m)*zzz]*double(1-2*((m+mb)%2)))*t;
	if(m<k) g[k+bb+d*l0]+=(alp[l0+d1*bb+m+(l0+d*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-m+mb)*zzz]+2*alp[l0+d1*bb+m+(l0+d2*bb-mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb+mb-k-m)*zzz]*double(1-2*((k+mb)%2)))*t;
      if(m==k)  g[k+bb+d*l0]+=(alp[l0+d1*bb+k+(l0+d*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-k+mb)*zzz]+2*alp[l0+d1*bb+k+(l0+d2*bb-mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb+mb-k-k)*zzz]*double(1-2*((k+mb)%2)))*t/2.;
      if(m==k)  g[k+bb+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-k+mb)*zzz]+2*alp[l0+d*bb+k+(l0+d2*bb-mb)*z+(l0+d1*bb+k)*zz+(l0+d3*bb+mb-k-k)*zzz]*double(1-2*((k+mb)%2)))*t/2.;

			}
	}
	}

// m+mb-k
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=m;mb<=l0+d2;mb++) {
	if((m+mb-k)<=l0+d3 && (m+mb-k)>=1) {

   	   t=lambda4*((im[d1][m]*im[d2][mb]-re[d1][m]*re[d2][mb])*im[d3][m+mb-k]+(im[d1][m]*re[d2][mb]+re[d1][m]*im[d2][mb])*re[d3][m+mb-k])/6.0;
 if(m!=mb) t=lambda4*((im[d1][m]*im[d2][mb]-re[d1][m]*re[d2][mb])*im[d3][m+mb-k]+(im[d1][m]*re[d2][mb]+re[d1][m]*im[d2][mb])*re[d3][m+mb-k])/3.0;

	g[k+bb+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb-m)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k+m+mb)*zzz]*double(1-2*((k+mb)%2))+alp[l0+d1*bb+m+(l0+d3*bb+k-m-mb)*z+(l0+d2*bb+mb)*zz+(l0+d*bb-k)*zzz]*double(1-2*((k+m)%2))+alp[l0+d1*bb+m+(l0+d2*bb+mb)*z+(l0+d*bb-k)*zz+(l0+d3*bb+k-m-mb)*zzz])*t;
			}
	}
	}


// mb-m-k
	for(int m=1;m<=l0+d1;m++) {
	for(int mb=k+m+1;mb<=l0+d2;mb++) {
		if((mb-m-k)>=1) {
         t=lambda4*((re[d1][m]*im[d2][mb]-im[d1][m]*re[d2][mb])*re[d3][mb-m-k]-(re[d1][m]*re[d2][mb]+im[d1][m]*im[d2][mb])*im[d3][mb-m-k])/6.;
 
	if(m>k) g[k+bb+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb+m)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-m+mb)*zzz]*double(1-2*((k+m+mb)%2))+2*alp[l0+d*bb+k+(l0+d2*bb-mb)*z+(l0+d1*bb+m)*zz+(l0+d3*bb+mb-m-k)*zzz]*double(1-2*(k%2)))*t;
	if(m<k) g[k+bb+d*l0]+=(alp[l0+d1*bb+m+(l0+d*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-m+mb)*zzz]*double(1-2*((k+m+mb)%2))+2*alp[l0+d1*bb+m+(l0+d2*bb-mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb+mb-m-k)*zzz]*double(1-2*(m%2)))*t;
       if(m==k) g[k+bb+d*l0]+=(alp[l0+d1*bb+k+(l0+d*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-k+mb)*zzz]*double(1-2*((k+k+mb)%2))+2*alp[l0+d1*bb+k+(l0+d2*bb-mb)*z+(l0+d*bb+k)*zz+(l0+d3*bb+mb-k-k)*zzz]*double(1-2*(k%2)))*t/2.;
       if(m==k) g[k+bb+d*l0]+=(alp[l0+d*bb+k+(l0+d1*bb+k)*z+(l0+d2*bb-mb)*zz+(l0+d3*bb-k-k+mb)*zzz]*double(1-2*((k+k+mb)%2))+2*alp[l0+d*bb+k+(l0+d2*bb-mb)*z+(l0+d1*bb+k)*zz+(l0+d3*bb+mb-k-k)*zzz]*double(1-2*(k%2)))*t/2.;

				}
	}
	}


		}

		} // end of the Ik loop

}

	
	
	t=0;
	for (int i=0;i<4*l0+4;i++) {
		t+=g[i]*g[i];
	}
	gg=sqrt(t);
	return gg;
}

//CALCULATE HESSIAN MATRIX - NUMERICAL APPROXIMATION
double nhessian(Vec_I_DP &x, double *hess,  DP funk(Vec_I_DP &, double *Gauntm, double *Alpha))
{

        DP ptempl[(4*l0+4)];
        DP ptempr[(4*l0+4)];
        DP ptempc[(4*l0+4)];
        DP ptempc2[(4*l0+4)];
        for (int i=0; i<4*l0+4; i++)
        {
                ptempl[i]=x[i];
                ptempr[i]=x[i];
                ptempc[i]=x[i];
                ptempc2[i]=x[i];
        }
        double del=1e-5;

        for(int k=0; k<4*l0+4; k++)
        {
                ptempr[k]=ptempr[k]+del;
                ptempl[k]=ptempl[k]-del;

                Vec_I_DP tempr(ptempr,4*l0+4);
                Vec_I_DP templ(ptempl,4*l0+4);
                Vec_I_DP tempc(ptempc,4*l0+4);
                hess[k+k*hs]=(funk(tempr,Gaunt_matrix_3rd_order, Alpha_matrix)+funk(templ,Gaunt_matrix_3rd_order, Alpha_matrix)-2*funk(tempc,Gaunt_matrix_3rd_order, Alpha_matrix))/(del*del);
                ptempr[k]=x[k];
                ptempl[k]=x[k];

        }

        for(int k=0; k<4*l0+4; k++) {
         for(int g=k+1; g<4*l0+4; g++) {
                ptempr[k]=ptempr[k]+del;
                ptempr[g]=ptempr[g]+del;
                ptempl[k]=ptempl[k]-del;
                ptempl[g]=ptempl[g]-del;
                ptempc[k]=ptempc[k]+del;
                ptempc[g]=ptempc[g]-del;
                ptempc2[k]=ptempc2[k]-del;
                ptempc2[g]=ptempc2[g]+del;

                Vec_I_DP tempr(ptempr,4*l0+4);
                Vec_I_DP templ(ptempl,4*l0+4);
                Vec_I_DP tempc(ptempc,4*l0+4);
                Vec_I_DP tempc2(ptempc2,4*l0+4);
        hess[k+g*hs]=(funk(tempr,Gaunt_matrix_3rd_order, Alpha_matrix)+funk(templ,Gaunt_matrix_3rd_order, Alpha_matrix)-funk(tempc,Gaunt_matrix_3rd_order, Alpha_matrix)-funk(tempc2,Gaunt_matrix_3rd_order, Alpha_matrix))/(4*del*del);
        hess[g+k*hs]=hess[k+g*hs];


                ptempr[k]=x[k];
                ptempr[g]=x[g];
                ptempl[k]=x[k];
                ptempl[g]=x[g];
                ptempc[k]=x[k];
                ptempc[g]=x[g];
                ptempc2[k]=x[k];
                ptempc2[g]=x[g];


                }
                }


        return 0;
}

//FUNCTION EVALUATION
DP funk(Vec_I_DP &x, double *Gauntm, double *Alpha)
{
	double H_gradient=0;
	double H_thirdorder =0;
	double H_fourthorder=0;
	double EPSILON = 1e-12;
	double H=0;

    complex<double>** c = Make2DArray(2, l0+2);

    for(int q=0;q<=1;q++)
    	c[q][0] = complex<double> (x[q],0);

    int count=2;
    
    for(int i=0;i<=1;i++)
    {
        for(int j=1;j<=l0+i;j++)
        {
           c[i][j]=complex<double> (x[count],x[(2*l0+1+count)]);
           count++;
        }
    }
    c[0][(l0+1)]=complex<double> (0,0);


   
    //Gradient term
    for(int l=l0; l<=l0+1; l++)
    {
        H_gradient += (1/2.0)*((pow(l-l_not,2))+tau)*real(c[(l-l0)][0]*c[(l-l0)][0]);
        for(int m=1; m<=l; m++)
        {
            H_gradient += ((pow(l-l_not,2)+tau)*(norm(c[(l-l0)][m])));
        }
    } 


   //Third Order Term

if (lambda3>EPSILON) {
 ///3rd order term in Hamiltonian
    int FirstIndex,SecondIndex,ThirdIndex;
 
    for(int el_1=l0;el_1<=l0+1;el_1++)
    {
	FirstIndex = el_1-l0;
        for(int el_2=l0;el_2<=l0+1;el_2++)
        {
	    SecondIndex = el_2-l0;
            for(int el_3=l0;el_3<=min(l0+1,el_1+el_2);el_3++)
            {
		ThirdIndex = el_3-l0;
                if((el_3>=abs(el_1-el_2))&&((el_1+el_2+el_3)%2==0))
                {
                    H_thirdorder += (lambda3/6.0)*Gauntm[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*el_1 + 8*Na*el_2 + 8*Na*Na*el_3]*real(c[ThirdIndex][0]*c[SecondIndex][0]*c[FirstIndex][0]);
                    for(int m2=1;m2<=el_2;m2++)
                    {
                            H_thirdorder+=4.0*(lambda3/6.0)*Gauntm[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m2+el_1) + 8*Na*(el_2-m2) + 8*Na*Na*el_3]*double(1-2*(m2%2))*real(c[ThirdIndex][0]*c[FirstIndex][m2]*conj(c[SecondIndex][m2]));
			    for(int m1=1;m1<=el_1;m1++)
                            {
                                    if(abs(m1+m2)<=el_3)
                                            H_thirdorder+=(lambda3/6.0)*Gauntm[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m1+el_1) + 8*Na*(m2+el_2) + 8*Na*Na*(-m1-m2+el_3)]*double(1-2*((m1+m2)%2))*2*(real(c[FirstIndex][m1]*c[SecondIndex][m2]*conj(c[ThirdIndex][m1+m2])));
                            }
                            for(int m1=1;m1<=m2;m1++)
                            {
                                    if((m2-m1)<=el_3)
                                            H_thirdorder+=(lambda3/6.0)*Gauntm[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m1+el_1) + 8*Na*(-m2+el_2) + 8*Na*Na*(m2-m1+el_3)]*double(1-2*(m2%2))*2*(real(c[FirstIndex][m1]*conj(c[SecondIndex][m2])*c[ThirdIndex][m2-m1]));
                            }
                            for(int m1=m2+1;m1<=el_1;m1++)
                            {
                                    if((m2-m1)<=el_3)
                                            H_thirdorder+=(lambda3/6.0)*Gauntm[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m1+el_1) + 8*Na*(-m2+el_2) + 8*Na*Na*(m2-m1+el_3)]*double(1-2*(m1%2))*2*(real(c[FirstIndex][m1]*conj(c[SecondIndex][m2])*conj(c[ThirdIndex][m1-m2])));
                            }
                    }
                }
            }
        }
    }

		}


// Fourth order term


//H_fourthorder=0;
	
	H_fourthorder=real(Alpha[(l0)+(l0)*z+(l0)*zz+(l0)*zzz]*pow(c[0][0],4)+Alpha[(bb+l0)+(bb+l0)*z+(bb+l0)*zz+(bb+l0)*zzz]*pow(c[1][0],4)+4*Alpha[(l0)+(bb+l0)*z+(bb+l0)*zz+(bb+l0)*zzz]*c[0][0]*pow(c[1][0],3)+6*Alpha[(l0)+(l0)*z+(bb+l0)*zz+(bb+l0)*zzz]*c[0][0]*c[0][0]*c[1][0]*c[1][0]+4*Alpha[(l0)+(l0)*z+(l0)*zz+(bb+l0)*zzz]*c[1][0]*pow(c[0][0],3));
	H_fourthorder*=lambda4/24.;



for(int d1=0;d1<2;d1++) {
for(int d2=0;d2<2;d2++) {
for(int d3=0;d3<2;d3++) {
for(int d4=0;d4<2;d4++) {


for(int m=1; m<=l0+min(d3,d4);m++)  
	H_fourthorder+=(lambda4/6.)*(Alpha[d1*bb+l0+(d2*bb+l0)*z+(d3*bb+l0+m)*zz+(d4*bb+l0-m)*zzz]*double(1-2*(m%2))+2*Alpha[d1*bb+l0+(d3*bb+l0+m)*z+(d2*bb+l0)*zz+(d4*bb+l0-m)*zzz])*real(c[d1][0]*c[d2][0]*c[d3][m]*conj(c[d4][m]));


for(int m=1; m<=l0+d2;m++) {  
for(int mb=1; mb<=l0+d3;mb++) {
if(m+mb<=l0+d4)
	H_fourthorder+=(lambda4/6.)*Alpha[d1*bb+l0+(d2*bb+l0+m)*z+(d3*bb+l0+mb)*zz+(d4*bb+l0-m-mb)*zzz]*double(1-2*(mb%2))*real(c[d1][0]*c[d2][m]*c[d3][mb]*conj(c[d4][m+mb]));
}
for(int mb=m; mb<=l0+d3;mb++) {
if(m+mb<=l0+d4)
	H_fourthorder+=double(2-int(m==mb))*(lambda4/6.)*(Alpha[d2*bb+l0+m+(d3*bb+l0+mb)*z+(d1*bb+l0)*zz+(d4*bb+l0-m-mb)*zzz]+Alpha[d2*bb+l0+m+(d4*bb+l0-m-mb)*z+(d3*bb+l0+mb)*zz+(d1*bb+l0)*zzz]*double(1-2*(m%2)))*real(c[d1][0]*c[d2][m]*c[d3][mb]*conj(c[d4][m+mb]));
}
for(int mb=1;mb<=min(m-1,l0+d3);mb++) {
if(m-mb<=l0+d4)
	H_fourthorder+=(lambda4/6.)*(Alpha[d1*bb+l0+(d2*bb+l0+m)*z+(d3*bb+l0-mb)*zz+(d4*bb+l0+mb-m)*zzz])*real(c[d1][0]*conj(c[d2][m])*c[d3][mb]*c[d4][m-mb]);
}
for(int mb=m+1;mb<=l0+d3;mb++) {
if(mb-m<=l0+d4)
	H_fourthorder+=(lambda4/6.)*(Alpha[d1*bb+l0+(d2*bb+l0+m)*z+(d3*bb+l0-mb)*zz+(d4*bb+l0+mb-m)*zzz]*double(1-2*((m+mb)%2))+Alpha[d2*bb+l0+m+(d3*bb+l0-mb)*z+(d1*bb+l0)*zz+(d4*bb+l0+mb-m)*zzz]*double(1-2*(m%2)))*real(c[d1][0]*conj(c[d2][m])*c[d3][mb]*conj(c[d4][mb-m]));
}
}
for(int m=1; m<=l0+d1;m++) {  
for(int mh=m;mh<=l0+d2;mh++) {

	for(int mb=1;mb<=l0+d3;mb++) {
	if(m+mh+mb<=l0+d4)
H_fourthorder+=(lambda4/12.)*double(2-int(m==mh))*(Alpha[d1*bb+l0+m+(d2*bb+l0+mh)*z+(d3*bb+l0+mb)*zz+(d4*bb+l0-m-mb-mh)*zzz])*double(1-2*(mb%2))*real(c[d1][m]*c[d2][mh]*c[d3][mb]*conj(c[d4][m+mh+mb]));
		}
	for(int mb=1;mb<=min(m+mh-1,l0+d3);mb++) {
	if(mh+m-mb<=l0+d4 && mh+m-mb>=1)

 H_fourthorder+=(lambda4/12.)*double(2-int(m==mh))*(Alpha[d1*bb+l0+m+(d2*bb+l0+mh)*z+(d3*bb+l0-mb)*zz+(d4*bb+l0-m-mh+mb)*zzz]+2*Alpha[d1*bb+l0+m+(d3*bb+l0-mb)*z+(d2*bb+l0+mh)*zz+(d4*bb+l0+mb-m-mh)*zzz]*double(1-2*((mh+mb)%2)))*real((c[d1][m])*(c[d2][mh])*conj(c[d3][mb])*conj(c[d4][mh+m-mb]));
		}
	for(int mb=m+mh+1;mb<=l0+d3;mb++) {
	if(mb-m-mh<=l0+d4 && mb-m-mh>=1) 

 H_fourthorder+=(lambda4/12.)*double(2-int(m==mh))*(Alpha[d1*bb+l0+m+(d2*bb+l0+mh)*z+(d3*bb+l0-mb)*zz+(d4*bb+l0-m-mh+mb)*zzz]*double(1-2*((m+mb+mh)%2))+2*Alpha[d1*bb+l0+m+(d3*bb+l0-mb)*z+(d2*bb+l0+mh)*zz+(d4*bb+l0+mb-m-mh)*zzz]*double(1-2*(m%2)))*real((c[d1][m])*(c[d2][mh])*conj(c[d3][mb])*(c[d4][mb-m-mh]));
		}

			     }
			}


		} } } } 


    H = H_gradient+H_thirdorder+H_fourthorder;

    for(int b = 0; b<2; b++)
	delete [] c[b];
    delete [] c;


	return H;
}


int counter=1;

DP tt;
DP NR::amotsa(Mat_IO_DP &p, Vec_O_DP &y, Vec_IO_DP &psum, Vec_O_DP &pb, DP &yb, DP funk(Vec_I_DP &x, double *Gaunt_matrix, double *Alpha_matrix), const int ihi, DP &yhi, const DP fac)
{
    ofstream yfluc ("yfluc.txt", ios_base::app);
    ofstream y_try ("ytry.txt", ios_base::app);
    mt19937 gen(100);

   
    int j;
    DP fac1,fac2,yflu,ytry;
    
    int ndim=p.ncols();
    Vec_DP ptry(ndim);
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for(j=0;j<ndim;j++)
        ptry[j]=psum[j]*fac1-p[ihi][j]*fac2; 
    ytry=funk(ptry,Gaunt_matrix_3rd_order, Alpha_matrix);
    if(ytry<=yb)
    {
        for(j=0;j<ndim;j++) pb[j]=ptry[j];
        yb=ytry;
    }
    yflu=ytry-tt*log(double(gen())/double(gen.max()));
    if(yflu<yhi)
    {
        y[ihi]=ytry;
        yhi=yflu;
        for(j=0;j<ndim;j++)
        {
            psum[j]+=ptry[j]-p[ihi][j];
            p[ihi][j]=ptry[j];
        }
    }

    counter+=1;

    yfluc.close(); 
    y_try.close();
    return yflu;

}

void NR::amebsa(Mat_IO_DP &p, Vec_IO_DP &y, Vec_O_DP &pb, DP &yb, const DP ftol, DP funk(Vec_I_DP &x, double *Gaunt_matrix, double *Alpha_matrix), int &iter, const DP temptr)
{
    mt19937 gen(rd());
    
    int i,ihi,ilo,j,n;
    DP rtol,yhi,ylo,ynhi,ysave,yt,ytry;
    
    int mpts=p.nrows();
    int ndim=p.ncols();
    Vec_DP psum(ndim);
    tt=-temptr;
    get_psum(p,psum);
    for(;;)
    {
        ilo=0;
        ihi=1;
        ynhi=ylo=y[0]+tt*log(double(gen())/double(gen.max()));
        yhi=y[1]+tt*log(double(gen())/double(gen.max()));
        ynhi=ylo=y[0]+tt*log(double(gen())/double(gen.max()));
        yhi=y[1]+tt*log(double(gen())/double(gen.max()));
        if(ylo>yhi)
        {
            ihi=0;
            ilo=1;
            ynhi=yhi;
            yhi=ylo;
            ylo=ynhi;
            ilo=1;
            ynhi=yhi;
            yhi=ylo;
            ylo=ynhi;
        }
        for(i=3;i<=mpts;i++)
        {
            yt=y[i-1]+tt*log(double(gen())/double(gen.max()));
            if(yt<=ylo)
            {
                ilo=i-1;
                ylo=yt;
            }
            if(yt>yhi)
            {
                ynhi=yhi;
                ihi=i-1;
                yhi=yt;
            }
            else if(yt>ynhi)
            {
                ynhi=yt;
            }
        }
        rtol=2.0*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo));
        if(rtol<ftol||iter<0)
        {
            SWAP(y[0],y[ilo]);
            for(n=0;n<ndim;n++)
            {
                SWAP(p[0][n],p[ilo][n]);
            }
            break;
            
        }
        iter-=2;
        ytry=amotsa(p,y,psum,pb,yb,funk,ihi,yhi,-1.0);
        if(ytry<=ylo)
        {
            ytry=amotsa(p,y,psum,pb,yb,funk,ihi,yhi,2.0);
            
        }
        else if(ytry>=ynhi)
        {
            ysave=yhi;
            ytry=amotsa(p,y,psum,pb,yb,funk,ihi,yhi,0.5);
            if(ytry>=ysave)
            {
                for(i=0;i<mpts;i++)
                {
                    if(i != ilo)
                    {
                        for(j=0;j<ndim;j++)
                        {
                            psum[j]=0.5*(p[i][j]+p[ilo][j]);
                            p[i][j]=psum[j];
                        }
                        y[i]=funk(psum,Gaunt_matrix_3rd_order, Alpha_matrix);
                    }
                }
                iter-=ndim;
                get_psum(p,psum);
            }
            
        }
        else ++iter;
    }
}



int main()
{
    clock_t t1;
    t1 = clock();   
    ofstream T_H,endcms,cms,parameters,time; 
    T_H.open("T_H.txt");
    endcms.open("endcms.txt");
    cms.open("cms.txt");
    parameters.open("parameters.txt");
    time.open("time.txt");
    mt19937 gen(500); //creates same set of random numbers
    wig_table_init(2*100,9); //define using largest 2*3j value and largest wigner type (i.e. 3j, 6j, 9j) - however this particular 2*100, 9 is commonly used becuase it doesn't take a lot of memory and yet it works for most programs
    wig_temp_init(2*100);
    
    double magg;
    double hess_array[(int)(4*l0+4)*(4*l0+4)];

    //Gaunt matrix initialization for Hamiltonian 3rd order term
    int FirstIndex,SecondIndex,ThirdIndex;
    for(int el_1=l0;el_1<=l0+1;el_1++)
    {   
	FirstIndex = el_1-l0;
        for(int el_2=l0;el_2<=l0+1;el_2++)
        {   
	    SecondIndex = el_2-l0;
            for(int el_3=l0;el_3<=l0+1;el_3++)
            {   
		ThirdIndex = el_3-l0;
                if((el_3<=abs(el_1+el_2))&&(el_3>=abs(el_1-el_2))&&((el_1+el_2+el_3)%2==0))
                {
                    Gaunt_matrix_3rd_order[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*el_1 + 8*Na*el_2 + 8*Na*Na*el_3] = Gaunt(el_1,el_2,el_3,0,0,0);
                    for(int m2=1;m2<=el_2;m2++)
                    {   
                            Gaunt_matrix_3rd_order[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m2+el_1) + 8*Na*(el_2-m2) + 8*Na*Na*el_3] = Gaunt(el_1,el_2,el_3,m2,-m2,0);
                            for(int m1=1;m1<=el_1;m1++)
                            {   
                                    if(abs(m1+m2)<=el_3)
                                            Gaunt_matrix_3rd_order[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m1+el_1) + 8*Na*(m2+el_2) + 8*Na*Na*(-m1-m2+el_3)] = Gaunt(el_1,el_2,el_3,m1,m2,-m1-m2);
                            }   
                            for(int m1=1;m1<=m2;m1++)
                            {   
                                    if((m2-m1)<=el_3)
                                           Gaunt_matrix_3rd_order[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m1+el_1) + 8*Na*(-m2+el_2) + 8*Na*Na*(m2-m1+el_3)] = Gaunt(el_1,el_2,el_3,m1,-m2,m2-m1);
			    }   
                            for(int m1=m2+1;m1<=el_1;m1++)
                            {   
                                    if((m2-m1)<=el_3)
                                           Gaunt_matrix_3rd_order[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m1+el_1) + 8*Na*(-m2+el_2) + 8*Na*Na*(m2-m1+el_3)] = Gaunt(el_1,el_2,el_3,m1,-m2,m2-m1);
                            }   
                    }   
                }   
            }   
        }   
    }
        //Gaunt matrix initialization for Hamiltonian 4th order term
int m1, m2, m3, m4;
int l1, l2, l3, l4;
double tempg;
for (int i=0; i<z*z*z*z; i++) {

	m1=i%z;
	m2=(i/z)%z;
	m3=(i/zz)%z;
	m4=(i/zzz)%z;

	if(m1>2*l0) {l1=l0+1; m1-=3*l0+2;} else {l1=l0; m1-=l0;}
	if(m2>2*l0) {l2=l0+1; m2-=3*l0+2;} else {l2=l0; m2-=l0;}
	if(m3>2*l0) {l3=l0+1; m3-=3*l0+2;} else {l3=l0; m3-=l0;}
	if(m4>2*l0) {l4=l0+1; m4-=3*l0+2;} else {l4=l0; m4-=l0;}
	
	tempg=0;
  for (int ll=0; ll<=2*l0+3; ll++) {

	tempg+=Gaunt(l1,l2,ll,m1,m2,-m1-m2)*Gaunt(l3,l4,ll,m3,m4,-m3-m4);
	}
	Alpha_matrix[i]=tempg;
}

   
    double t;
   
    
    parameters<<"l: "<<l0<<" "<<l0+1<<endl;
    parameters<<"el_not: "<<l_not<<endl;
    parameters<<"tau: "<<tau<<endl;
    parameters<<"lambda3: "<<lambda3<<endl;
    parameters<<"lambda4: "<<lambda4<<endl;
    parameters<<endl;
//Define staring point values
    DP p_array[4*l0+5][4*l0+4];   
    parameters<<"Starting p matrix"<<endl;
    //random starting cm's
    for(int a =0;a<4*l0+5;a++)
    {
	for(int b=0;b<4*l0+4;b++)
    	{
        	p_array[a][b]=sqrt(abs(tau/lambda4))*(double(gen())/double(gen.max())-0.5)*2;
		parameters<<p_array[a][b]<<" ";
        }
	parameters<<endl;
    }
    parameters<<endl;
    Mat_IO_DP p(*p_array,4*l0+5,4*l0+4);
    
    
//Set corresponding initial y-values
    DP psi_array[4*l0+4];
    DP y_array[4*l0+5];
    DP y_array2[4*l0+5];
    parameters<<"Starting y value: "<<endl;
    for(int xi=0; xi<4*l0+5; xi++)
    {
        for(int xii=0; xii<4*l0+4;xii++)
        {
            psi_array[xii]=p_array[xi][xii];
	}


        Vec_I_DP psi(psi_array,4*l0+4);


	y_array[xi]=funk(psi, Gaunt_matrix_3rd_order, Alpha_matrix);
	parameters<<y_array[xi]<< ", ";
    }
    parameters<<endl;
    Vec_IO_DP y(y_array,4*l0+5);
    
    double pb_array[(int)(4*l0+4)];
    Vec_O_DP pb(pb_array,4*l0+4);
    const DP ftol = 1.0e-8;
    parameters<<"FTolerance: "<<ftol<<endl;
    
    DP yb = 10;
    parameters<<"yb starting: "<<yb<<endl; 
    int iter= 10000;
    parameters<<"Number of iterations: "<<iter<<endl;
    const DP temptr = 1;
    parameters<<"Starting temperature: "<<temptr<<endl;
    double Hamiltonian_values[100000]; 

    NR::amebsa(p,y,pb,yb,ftol,funk,iter,temptr);

     double twoelgrad[4*l0+4];
     double ntwoelgrad[4*l0+4];
     double mag2elgrad;

   
    double tinterval = 0.1;
    parameters<<"Temperature interval: "<<tinterval<<endl;
    

/////////////ANNEALING SCHEDULE/////////////////

	    for(t=temptr*0.5;t>temptr*0.1;t-=tinterval*0.05)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_H<<t<<" "<<yb<<endl;

		cms<<pb[0]<<" 0"<<endl;
		int endcount=2;
		for(int i=0;i<l0;i++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}

		cms<<pb[1]<<" 0"<<endl;
		for(int j=0;j<=l0;j++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}
		cms<<endl;
	    }
	    
	    for(t=temptr*0.1;t>=0.01;t-=tinterval*0.01)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_H<<t<<" "<<yb<<endl;
			
		cms<<pb[0]<<" 0"<<endl;
		int endcount=2;
		for(int i=0;i<l0;i++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}

		cms<<pb[1]<<" 0"<<endl;
		for(int j=0;j<=l0;j++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}
		cms<<endl;
	    }

	    for(t=temptr*0.01;t>=0.001;t-=tinterval*0.001)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_H<<t<<" "<<yb<<endl;
			
		cms<<pb[0]<<" 0"<<endl;
		int endcount=2;
		for(int i=0;i<l0;i++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}

		cms<<pb[1]<<" 0"<<endl;
		for(int j=0;j<=l0;j++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}
		cms<<endl;    
	    }
	  
	    for(t=temptr*0.001;t>=0.0001;t-=tinterval*0.0001)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_H<<t<<" "<<yb<<endl;

		cms<<pb[0]<<" 0"<<endl;
		int endcount=2;
		for(int i=0;i<l0;i++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}

		cms<<pb[1]<<" 0"<<endl;
		for(int j=0;j<=l0;j++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}
		cms<<endl;
	    }

	    for(t=temptr*0.0001;t>=0.00001;t-=tinterval*0.00001)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_H<<t<<" "<<yb<<endl;

		cms<<pb[0]<<" 0"<<endl;
		int endcount=2;
		for(int i=0;i<l0;i++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}

		cms<<pb[1]<<" 0"<<endl;
		for(int j=0;j<=l0;j++)
		{
		    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
		    endcount++;
		}
		cms<<endl;
	    }
  
	for(t=temptr*0.00001;t>=0.0;t-=tinterval*0.000001)
   	{
        	iter=10000;
        	NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
        	T_H<<t<<" "<<yb<<endl;
        
        	cms<<pb[0]<<" 0"<<endl;
        	int endcount=2;
        	for(int i=0;i<l0;i++)
        	{
       		     cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
        	    endcount++;
       		}

        	cms<<pb[1]<<" 0"<<endl;
        	for(int j=0;j<=l0;j++)
        	{
        	    cms<<pb[endcount]<<" "<<pb[2*l0+1+endcount]<<endl;
        	    endcount++;
        	}
        	cms<<endl;
	}



    endcms<<pb[0]<<" 0"<<endl;
    int endcount_2=2;
    for(int i=0;i<l0;i++)
    {
        endcms<<pb[endcount_2]<<" "<<pb[2*l0+1+endcount_2]<<endl;
        endcount_2++;
    }

    endcms<<pb[1]<<" 0"<<endl;
    for(int j=0;j<=l0;j++)
    {
        endcms<<pb[endcount_2]<<" "<<pb[2*l0+1+endcount_2]<<endl;
        endcount_2++;
    }
   magg=nhessian(pb,hess_array,funk);
   Eigen::MatrixXd hhh(hs,hs);
   for (int m=0; m<hs*hs; m++) {
            int x = m%(4*l0+4);
            int y = (m/hs)%hs;
    hhh(x,y)=hess_array[m];
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es(hhh);

   parameters << "Eigenvalues: " << endl << es.eigenvalues() << endl;
   parameters<<"Lowest point derivative value: "<<grad(pb,twoelgrad,Gaunt_matrix_3rd_order,Alpha_matrix)<<endl; 
   parameters<<"Lowest_yb: "<<yb<<endl;
   T_H.close();
   endcms.close();
   cms.close();
   parameters.close();

   //Free Gaunt Matricies
   delete[] Gaunt_matrix_3rd_order;
   delete[] Gaunt_matrix_4th_order_1;
   delete[] Gaunt_matrix_4th_order_2;

   t = clock() - t;
   time<<"Computation Time: "<<t/CLOCKS_PER_SEC/60.0<<" minutes";
   time.close();

   return 0; 
}

#include <iostream>
#include <math.h>
#include <time.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cmath>
#include <random>
#include <fstream>
#include <complex>
#include <chrono>
#include <iomanip>
#include "nrutil.h"
#include "nrtypes.h"
#include "nr.h"
#include "wigxjpf.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace std;
const double pi = 3.1415926535897;
random_device rd;

const double tau = -1.0;
const double lambda4 = 1.0;

///Set lambda3 to be any value between -1.0 and 1.0//////////
const double lambda3 = 0.0;

double* Make1DArray(int arraySizeA)
{	
    int size;
    size=arraySizeA;
    double *Array = new double[size];
    return Array;
}


double* Make3DArray(int arraySizeA, int arraySizeB, int arraySizeC)
{

    int size;
    size=arraySizeA*arraySizeB*arraySizeC;
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

double* Hamiltonian_values = Make1DArray(10000);
double* Gradient_values = Make1DArray(10000);
    
//these need to be changed
int Na = 2*el_not+1;
int Nb = 2*el_not+1;
int Nc = 2*el_not+1;

int na = 2*el_not+1;
int nb = 2*el_not+1;
int nc = 2*el_not+1;
int nd = 2*el_not+1; 

int NaNb = (2*el_not+1)*(2*el_not+1);
int NaNbNc = (2*el_not+1)*(2*el_not+1)*(2*el_not+1);


int hs = (2*el_not+1);
int nanb = (2*el_not+1)*(2*el_not+1);
int nanbnc = (2*el_not+1)*(2*el_not+1)*(2*el_not+1);

double* Gaunt_matrix_3rd_order = Make3DArray(Na,Nb,Nc);
double* Gaunt_matrix_4th_order = Make4DArray(na,nb,nc,nd);

namespace NR
{
    void amebsa(Mat_IO_DP &p, Vec_IO_DP &y, Vec_O_DP &pb, DP &yb, const DP ftol, DP funk(Vec_I_DP &, double* Gaunt_matrix, double* Gaunt_4th), int &iter, const DP temptr);
    DP amotsa(Mat_IO_DP &p, Vec_O_DP &y, Vec_IO_DP &psum, Vec_IO_DP &psp, DP &yb, DP funk(Vec_I_DP &, double* Gaunt_matrix, double* Gaunt_4th), const int ihi, DP &yhi, const DP fac);
}

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

//DERIVATIVE CALCULATION
double grad(Vec_I_DP &x, double *g, double *Gaunt_matrix, double* alp)
{

	double gg=0;
	double t;
	vector<double> re;
	re.resize(el_not+1);
	vector<double> im;
	im.resize(el_not+1);
	re[0]=-1000000;
	im[0]=-1000000;
	for (int i=1; i<=el_not;i++){
		re[i]=x[i];
		im[i]=x[i+el_not];
		}

// dE/dc[0] starts
	g[0]=tau*x[0]+lambda3*Gaunt_matrix[el_not+Na*el_not+NaNb*el_not]*x[0]*x[0]/2.0+lambda4*alp[el_not+el_not*na+el_not*nanb+el_not*nanbnc]*pow(x[0],3)/6.0;
	for(int m=1;m<=el_not;m++) {
	t=re[m]*re[m]+im[m]*im[m];
	g[0]+=lambda3*double(1-2*(m%2))*Gaunt_matrix[el_not+m+Na*(el_not-m)+NaNb*el_not]*t;
	g[0]+=lambda4*(2*alp[el_not+(el_not+m)*na+(el_not)*nanb+(el_not-m)*nanbnc]+double(1-2*(m%2))*alp[el_not+m+(el_not-m)*na+(el_not)*nanb+(el_not)*nanbnc])*x[0]*t/3.0;
	// m+mb
	for(int mb=1;mb<=el_not-m;mb++) {
		t=lambda4*((re[m]*re[mb]-im[m]*im[mb])*re[m+mb]+(re[m]*im[mb]+im[m]*re[mb])*im[m+mb])/6.;
		g[0]+=alp[el_not+(el_not+m)*na+(el_not+mb)*nanb+(el_not-m-mb)*nanbnc]*double(1-2*(mb%2))*t;
		if(mb>m) g[0]+=2*(alp[el_not+m+(el_not+mb)*na+(el_not-m-mb)*nanb+(el_not)*nanbnc]+alp[el_not+m+mb+(el_not-m)*na+(el_not-mb)*nanb+el_not*nanbnc]*double(1-2*(m%2)))*t;
		if(mb==m)  g[0]+=(alp[el_not+m+(el_not+mb)*na+(el_not-m-mb)*nanb+(el_not)*nanbnc]+alp[el_not+m+mb+(el_not-m)*na+(el_not-mb)*nanb+el_not*nanbnc]*double(1-2*(m%2)))*t;
		}
	// m-mb
	for(int mb=1;mb<=m-1;mb++) {
		t=lambda4*((re[m]*re[mb]+im[m]*im[mb])*re[m-mb]+(im[m]*re[mb]-re[m]*im[mb])*im[m-mb])/6.;
		g[0]+=alp[el_not+(el_not+m)*na+(el_not-mb)*nanb+(el_not-m+mb)*nanbnc]*t;
	}
	// mb-m
	for(int mb=m+1;mb<=el_not;mb++) {
		t=lambda4*((re[m]*re[mb]+im[m]*im[mb])*re[mb-m]+(im[mb]*re[m]-re[mb]*im[m])*im[mb-m])/6.;
		g[0]+=(alp[el_not+(el_not+m)*na+(el_not-mb)*nanb+(el_not-m+mb)*nanbnc]*double(1-2*(mb%2))+alp[el_not+m+(el_not-mb)*na+(el_not)*nanb+(el_not+mb-m)*nanbnc])*double(1-2*(m%2))*t;
		}
	}


// dE/dc[0] ends and begins dE/dR[m]

	for(int k=1; k<=el_not; k++) {

	g[k]=2*tau+2*lambda3*Gaunt_matrix[el_not+k+Na*(el_not-k)+NaNb*el_not]*double(1-2*(k%2))*x[0];
        g[k]+=lambda4*(2*alp[el_not+(el_not+k)*na+el_not*nanb+(el_not-k)*nanbnc]+alp[el_not+k+(el_not-k)*na+el_not*nanb+el_not*nanbnc]*double(1-2*(k%2)))*x[0]*x[0]/3.;
	g[k]*=re[k];

// The m,m+k single loop
	for(int m=1;m<=el_not-k;m++) {
		t=re[m]*re[m+k]+im[m]*im[k+m];
	g[k]+=2*lambda3*(Gaunt_matrix[el_not+m+(el_not+k)*Na+(el_not-m-k)*NaNb]+Gaunt_matrix[el_not+k+m+(el_not-m)*Na+(el_not-k)*NaNb])*double(1-2*((m+k)%2))*t/3.;
	g[k]+=lambda4*x[0]*(alp[el_not+(el_not+k)*na+(el_not+m)*nanb+(el_not-k-m)*nanbnc]*double(1-2*(m%2))+alp[el_not+k+(el_not+m)*na+(el_not-k-m)*nanb+el_not*nanbnc]+alp[el_not+(el_not+m)*na+(el_not+k)*nanb+(el_not-k-m)*nanbnc]*double(1-2*(k%2)))*t/3.;
	if(m>k) g[k]+=lambda4*x[0]*double(1-2*(k%2))*alp[el_not+k+m+(el_not-k)*na+(el_not-m)*nanb+el_not*nanbnc]*t/3.;
	if(m<k) g[k]+=lambda4*x[0]*double(1-2*(m%2))*alp[el_not+k+m+(el_not-m)*na+(el_not-k)*nanb+el_not*nanbnc]*t/3.;
       if(m==k) g[k]+=lambda4*x[0]*double(1-2*(k%2))*alp[el_not+k+k+(el_not-k)*na+(el_not-k)*nanb+el_not*nanbnc]*t/3.;
	}


// The m,m-k single loop
	for(int m=k+1;m<=el_not;m++) {
		t=re[m]*re[m-k]+im[m]*im[m-k];
	g[k]+=2*lambda3*Gaunt_matrix[el_not+m+(el_not-k)*Na+(el_not+k-m)*NaNb]*double(1-2*(m%2))*t/3.;
	g[k]+=lambda4*x[0]*(alp[el_not+(el_not+k)*na+(el_not-m)*nanb+(el_not-k+m)*nanbnc]*double(1-2*((k+m)%2))+2*alp[el_not+(el_not+m)*na+(el_not-k)*nanb+(el_not-m+k)*nanbnc]+alp[el_not+k+(el_not-m)*na+(el_not)*nanb+(el_not+m-k)*nanbnc]*double(1-2*(k%2)))*t/6.;
	}
// The m,k-m single loop
	for(int m=1;m<=k-1;m++) {
		t=re[m]*re[k-m]-im[m]*im[k-m];
	g[k]+=2*lambda3*Gaunt_matrix[el_not+k+(el_not-m)*Na+(el_not+m-k)*NaNb]*double(1-2*(k%2))*t/3.0;
	g[k]+=lambda4*x[0]*(alp[el_not+(el_not-k)*na+(el_not+m)*nanb+(el_not+k-m)*nanbnc]+2*alp[el_not+(el_not+m)*na+(el_not+k-m)*nanb+(el_not-k)*nanbnc]*double(1-2*((m+k)%2))+alp[el_not+m+(el_not-k)*na+el_not*nanb+(el_not+k-m)*nanbnc]*double(1-2*(m%2)))*t/6.0;	
	if(m<=((int) floor(double(k)/2.))) {
	if(k==2*m) {
	g[k]+=lambda3*Gaunt_matrix[el_not+m+(el_not+k-m)*Na+(el_not-k)*NaNb]*double(1-2*(k%2))*t/3.0;
	g[k]+=lambda4*x[0]*(alp[el_not+m+(el_not+k-m)*na+(el_not-k)*nanb+(el_not)*nanbnc]+alp[el_not+k+(el_not-m)*na+(el_not+m-k)*nanb+el_not*nanbnc]*double(1-2*(m%2)))*t/6.;
		} else
		{
      g[k]+=2*lambda3*Gaunt_matrix[el_not+m+(el_not+k-m)*Na+(el_not-k)*NaNb]*double(1-2*(k%2))*t/3.0;
      g[k]+=2*lambda4*x[0]*(alp[el_not+m+(el_not+k-m)*na+(el_not-k)*nanb+(el_not)*nanbnc]+alp[el_not+k+(el_not-m)*na+(el_not+m-k)*nanb+el_not*nanbnc]*double(1-2*(m%2)))*t/6.;
		}
					}
				}
// Start the double loops
// m+k+mb

	for(int m=1;m<=el_not;m++) {
	for(int mb=1;mb<=el_not;mb++) {
	if(m+k+mb<=el_not) {

		t=lambda4*((re[m]*re[mb]-im[m]*im[mb])*re[m+k+mb]+(re[m]*im[mb]+im[m]*re[mb])*im[m+k+mb])/6.;
	g[k]+=alp[el_not+m+(el_not+k)*na+(el_not+mb)*nanb+(el_not-m-mb-k)*nanbnc]*double(1-2*(mb%2))*t;
	if(mb>m) g[k]+=2*(alp[el_not-m+(el_not+m+mb+k)*na+(el_not-mb)*nanb+(el_not-k)*nanbnc]*double(1-2*(m%2))+alp[el_not+m+(el_not+mb)*na+(el_not+k)*nanb+(el_not-m-mb-k)*nanbnc]*double(1-2*(k%2)))*t;
	if(mb==m)  g[k]+=(alp[el_not-m+(el_not+m+mb+k)*na+(el_not-mb)*nanb+(el_not-k)*nanbnc]*double(1-2*(m%2))+alp[el_not+m+(el_not+mb)*na+(el_not+k)*nanb+(el_not-m-mb-k)*nanbnc]*double(1-2*(k%2)))*t;
	}	}	}

// k-m-mb
	for(int m=1;m<=el_not;m++) {
	for(int mb=m;mb<=el_not;mb++) {
	if((k-m-mb)>=1 && (m+mb)<=(k-1)) {
		t=lambda4*((re[m]*re[mb]-im[m]*im[mb])*re[k-m-mb]-(re[m]*im[mb]+im[m]*re[mb])*im[k-m-mb])/6.0;
		if(m!=mb) t=lambda4*((re[m]*re[mb]-im[m]*im[mb])*re[k-m-mb]-(re[m]*im[mb]+im[m]*re[mb])*im[k-m-mb])/3.0;
		g[k]+=(alp[el_not+m+(el_not-k)*na+(el_not+mb)*nanb+(el_not+k-m-mb)*nanbnc]*double(1-2*(m%2))+alp[el_not+m+(el_not+mb)*na+(el_not-k)*nanb+(el_not+k-m-mb)*nanbnc]*double(1-2*((k+m+mb)%2)))*t;
			}   }	}	

// k+m-mb 
	for(int m=1;m<=el_not;m++) {
	for(int mb=1;mb<=min(k+m-1,el_not);mb++) {
	if((k+m-mb)<=el_not && (k+m-mb)>=1) {

	    t=lambda4*((re[m]*re[mb]+im[m]*im[mb])*re[k+m-mb]+(im[m]*re[mb]-im[mb]*re[m])*im[k+m-mb])/12.;
	if(m>k) g[k]+=2*(alp[el_not+k+(el_not+m)*na+(el_not-mb)*nanb+(el_not-k-m+mb)*nanbnc]+2*alp[el_not+k+(el_not-mb)*na+(el_not+m)*nanb+(el_not+mb-k-m)*nanbnc]*double(1-2*((m+mb)%2)))*t;
	if(m<k) g[k]+=2*(alp[el_not+k+(el_not+m)*na+(el_not-mb)*nanb+(el_not-k-m+mb)*nanbnc]+2*alp[el_not+m+(el_not-mb)*na+(el_not+k)*nanb+(el_not+mb-k-m)*nanbnc]*double(1-2*((k+mb)%2)))*t;
      if(m==k)  g[k]+=2*(alp[el_not+k+(el_not+k)*na+(el_not-mb)*nanb+(el_not-k-k+mb)*nanbnc]+2*alp[el_not+k+(el_not-mb)*na+(el_not+k)*nanb+(el_not+mb-k-k)*nanbnc]*double(1-2*((k+mb)%2)))*t;

			}
	}
	}
// m+mb-k
	for(int m=1;m<=el_not;m++) {
	for(int mb=m;mb<=el_not;mb++) {
	if((m+mb-k)<=el_not && (m+mb-k)>=1) {

   	   t=lambda4*((re[m]*re[mb]-im[m]*im[mb])*re[m+mb-k]+(re[m]*im[mb]+im[m]*re[mb])*im[m+mb-k])/6.0;
 if(m!=mb) t=lambda4*((re[m]*re[mb]-im[m]*im[mb])*re[m+mb-k]+(re[m]*im[mb]+im[m]*re[mb])*im[m+mb-k])/3.0;

	g[k]+=(alp[el_not+k+(el_not-m)*na+(el_not-mb)*nanb+(el_not-k+m+mb)*nanbnc]*double(1-2*((k+mb)%2))+alp[el_not+m+(el_not+k-m-mb)*na+(el_not+mb)*nanb+(el_not-k)*nanbnc]*double(1-2*((k+m)%2))+alp[el_not+m+(el_not+mb)*na+(el_not-k)*nanb+(el_not+k-m-mb)*nanbnc])*t;
			}
	}
	}
	for(int m=1;m<=el_not;m++) {
	for(int mb=k+m+1;mb<=el_not;mb++) {
		if(mb-m-k>=1) {
         t=lambda4*((re[m]*re[mb]+im[m]*im[mb])*re[mb-m-k]+(re[m]*im[mb]-im[m]*re[mb])*im[mb-m-k])/6.;
 
	if(m>k) g[k]+=(alp[el_not+k+(el_not+m)*na+(el_not-mb)*nanb+(el_not-k-m+mb)*nanbnc]*double(1-2*((k+m+mb)%2))+2*alp[el_not+k+(el_not-mb)*na+(el_not+m)*nanb+(el_not+mb-m-k)*nanbnc]*double(1-2*(k%2)))*t;
	if(m<k) g[k]+=(alp[el_not+m+(el_not+k)*na+(el_not-mb)*nanb+(el_not-k-m+mb)*nanbnc]*double(1-2*((k+m+mb)%2))+2*alp[el_not+m+(el_not-mb)*na+(el_not+k)*nanb+(el_not+mb-m-k)*nanbnc]*double(1-2*(m%2)))*t;
       if(m==k) g[k]+=(alp[el_not+k+(el_not+k)*na+(el_not-mb)*nanb+(el_not-k-k+mb)*nanbnc]*double(1-2*((k+k+mb)%2))+2*alp[el_not+k+(el_not-mb)*na+(el_not+k)*nanb+(el_not+mb-k-k)*nanbnc]*double(1-2*(k%2)))*t;

				}
	}
	}
			} // end of the Rk loop


//Begin calculating dE/dIk

	for(int k=1;k<=el_not;k++) {

        g[k+el_not]=2*tau+2*lambda3*Gaunt_matrix[el_not+k+Na*(el_not-k)+NaNb*el_not]*double(1-2*(k%2))*x[0];
        g[k+el_not]+=lambda4*(2*alp[el_not+(el_not+k)*na+el_not*nanb+(el_not-k)*nanbnc]+alp[el_not+k+(el_not-k)*na+el_not*nanb+el_not*nanbnc]*double(1-2*(k%2)))*x[0]*x[0]/3.;
        g[k+el_not]*=im[k];

// The m,m+k single loop
        for(int m=1;m<=el_not-k;m++) {
                t=re[m]*im[m+k]-im[m]*re[k+m];
        g[k+el_not]+=2*lambda3*(Gaunt_matrix[el_not+m+(el_not+k)*Na+(el_not-m-k)*NaNb]+Gaunt_matrix[el_not+k+m+(el_not-m)*Na+(el_not-k)*NaNb])*double(1-2*((m+k)%2))*t/3.;
        g[k+el_not]+=lambda4*x[0]*(alp[el_not+(el_not+k)*na+(el_not+m)*nanb+(el_not-k-m)*nanbnc]*double(1-2*(m%2))+alp[el_not+k+(el_not+m)*na+(el_not-k-m)*nanb+el_not*nanbnc]+alp[el_not+(el_not+m)*na+(el_not+k)*nanb+(el_not-k-m)*nanbnc]*double(1-2*(k%2)))*t/3.;
        if(m>k) g[k+el_not]+=lambda4*x[0]*double(1-2*(k%2))*alp[el_not+k+m+(el_not-k)*na+(el_not-m)*nanb+el_not*nanbnc]*t/3.;
        if(m<k) g[k+el_not]+=lambda4*x[0]*double(1-2*(m%2))*alp[el_not+k+m+(el_not-m)*na+(el_not-k)*nanb+el_not*nanbnc]*t/3.;
       if(m==k) g[k+el_not]+=lambda4*x[0]*double(1-2*(k%2))*alp[el_not+k+k+(el_not-k)*na+(el_not-k)*nanb+el_not*nanbnc]*t/3.;
        }


// The m,m-k single loop
        for(int m=k+1;m<=el_not;m++) {
                t=im[m]*re[m-k]-re[m]*im[m-k];
        g[k+el_not]+=2*lambda3*Gaunt_matrix[el_not+m+(el_not-k)*Na+(el_not+k-m)*NaNb]*double(1-2*(m%2))*t/3.;
        g[k+el_not]+=lambda4*x[0]*(alp[el_not+(el_not+k)*na+(el_not-m)*nanb+(el_not-k+m)*nanbnc]*double(1-2*((k+m)%2))+2*alp[el_not+(el_not+m)*na+(el_not-k)*nanb+(el_not-m+k)*nanbnc]+alp[el_not+k+(el_not-m)*na+(el_not)*nanb+(el_not+m-k)*nanbnc]*double(1-2*(k%2)))*t/6.;
        }
// The m,k-m single loop
        for(int m=1;m<=k-1;m++) {
                t=re[m]*im[k-m]+im[m]*re[k-m];
        g[k+el_not]+=2*lambda3*Gaunt_matrix[el_not+k+(el_not-m)*Na+(el_not+m-k)*NaNb]*double(1-2*(k%2))*t/3.0;
        g[k+el_not]+=lambda4*x[0]*(alp[el_not+(el_not-k)*na+(el_not+m)*nanb+(el_not+k-m)*nanbnc]+2*alp[el_not+(el_not+m)*na+(el_not+k-m)*nanb+(el_not-k)*nanbnc]*double(1-2*((m+k)%2))+alp[el_not+m+(el_not-k)*na+el_not*nanb+(el_not+k-m)*nanbnc]*double(1-2*(m%2)))*t/6.0;
        if(m<=((int) floor(double(k)/2.))) {
        if(k==2*m) {
        g[k+el_not]+=lambda3*Gaunt_matrix[el_not+m+(el_not+k-m)*Na+(el_not-k)*NaNb]*double(1-2*(k%2))*t/3.0;
        g[k+el_not]+=lambda4*x[0]*(alp[el_not+m+(el_not+k-m)*na+(el_not-k)*nanb+(el_not)*nanbnc]+alp[el_not+k+(el_not-m)*na+(el_not+m-k)*nanb+el_not*nanbnc]*double(1-2*(m%2)))*t/6.;
                } else
                {
      g[k+el_not]+=2*lambda3*Gaunt_matrix[el_not+m+(el_not+k-m)*Na+(el_not-k)*NaNb]*double(1-2*(k%2))*t/3.0;
      g[k+el_not]+=2*lambda4*x[0]*(alp[el_not+m+(el_not+k-m)*na+(el_not-k)*nanb+(el_not)*nanbnc]+alp[el_not+k+(el_not-m)*na+(el_not+m-k)*nanb+el_not*nanbnc]*double(1-2*(m%2)))*t/6.;
                }
                                        }
                                }


// Start the double loops
// The double loops from the new gradient
// m+k+mb

	for(int m=1;m<=el_not;m++) {
	for(int mb=1;mb<=el_not;mb++) {
	if(m+k+mb<=el_not) {

		t=lambda4*((re[m]*re[mb]-im[m]*im[mb])*im[m+k+mb]-(re[m]*im[mb]+im[m]*re[mb])*re[m+k+mb])/6.;
	g[k+el_not]+=alp[el_not+m+(el_not+k)*na+(el_not+mb)*nanb+(el_not-m-mb-k)*nanbnc]*double(1-2*(mb%2))*t;
	if(mb>m) g[k+el_not]+=2*(alp[el_not-m+(el_not+m+mb+k)*na+(el_not-mb)*nanb+(el_not-k)*nanbnc]*double(1-2*(m%2))+alp[el_not+m+(el_not+mb)*na+(el_not+k)*nanb+(el_not-m-mb-k)*nanbnc]*double(1-2*(k%2)))*t;
	if(mb==m)  g[k+el_not]+=(alp[el_not-m+(el_not+m+mb+k)*na+(el_not-mb)*nanb+(el_not-k)*nanbnc]*double(1-2*(m%2))+alp[el_not+m+(el_not+mb)*na+(el_not+k)*nanb+(el_not-m-mb-k)*nanbnc]*double(1-2*(k%2)))*t;
	}	}	}

// k-m-mb
	for(int m=1;m<=el_not;m++) {
	for(int mb=m;mb<=el_not;mb++) {
	if((k-m-mb)>=1) {
		          t=lambda4*((re[m]*re[mb]-im[m]*im[mb])*im[k-m-mb]+(re[m]*im[mb]+im[m]*re[mb])*re[k-m-mb])/6.0;
		if(m!=mb) t=lambda4*((re[m]*re[mb]-im[m]*im[mb])*im[k-m-mb]+(re[m]*im[mb]+im[m]*re[mb])*re[k-m-mb])/3.0;
		g[k+el_not]+=(alp[el_not+m+(el_not-k)*na+(el_not+mb)*nanb+(el_not+k-m-mb)*nanbnc]*double(1-2*(m%2))+alp[el_not+m+(el_not+mb)*na+(el_not-k)*nanb+(el_not+k-m-mb)*nanbnc]*double(1-2*((k+m+mb)%2)))*t;
			}   }	}	

// k+m-mb 
	for(int m=1;m<=el_not;m++) {
	for(int mb=1;mb<=min(k+m-1,el_not);mb++) {
	if((k+m-mb)<=el_not && (k+m-mb)>=1) {

	    t=lambda4*((re[m]*re[mb]+im[m]*im[mb])*im[k+m-mb]+(re[m]*im[mb]-re[mb]*im[m])*re[k+m-mb])/6.0;
	if(m>k) g[k+el_not]+=(alp[el_not+k+(el_not+m)*na+(el_not-mb)*nanb+(el_not-k-m+mb)*nanbnc]+2*alp[el_not+k+(el_not-mb)*na+(el_not+m)*nanb+(el_not+mb-k-m)*nanbnc]*double(1-2*((m+mb)%2)))*t;
	if(m<k) g[k+el_not]+=(alp[el_not+m+(el_not+k)*na+(el_not-mb)*nanb+(el_not-k-m+mb)*nanbnc]+2*alp[el_not+m+(el_not-mb)*na+(el_not+k)*nanb+(el_not+mb-k-m)*nanbnc]*double(1-2*((k+mb)%2)))*t;
      if(m==k)  g[k+el_not]+=(alp[el_not+k+(el_not+k)*na+(el_not-mb)*nanb+(el_not-k-k+mb)*nanbnc]+2*alp[el_not+k+(el_not-mb)*na+(el_not+k)*nanb+(el_not+mb-k-k)*nanbnc]*double(1-2*((k+mb)%2)))*t;

			}
	}
	}
// m+mb-k
	for(int m=1;m<=el_not;m++) {
	for(int mb=m;mb<=el_not;mb++) {
	if((m+mb-k)<=el_not && (m+mb-k)>=1) {

   	   t=lambda4*((im[m]*im[mb]-re[m]*re[mb])*im[m+mb-k]+(im[m]*re[mb]+re[m]*im[mb])*re[m+mb-k])/6.0;
 if(m!=mb) t=lambda4*((im[m]*im[mb]-re[m]*re[mb])*im[m+mb-k]+(im[m]*re[mb]+re[m]*im[mb])*re[m+mb-k])/3.0;

	g[k+el_not]+=(alp[el_not+k+(el_not-m)*na+(el_not-mb)*nanb+(el_not-k+m+mb)*nanbnc]*double(1-2*((k+mb)%2))+alp[el_not+m+(el_not+k-m-mb)*na+(el_not+mb)*nanb+(el_not-k)*nanbnc]*double(1-2*((k+m)%2))+alp[el_not+m+(el_not+mb)*na+(el_not-k)*nanb+(el_not+k-m-mb)*nanbnc])*t;
			}
	}
	}

// mb-m-k
	for(int m=1;m<=el_not;m++) {
	for(int mb=k+m+1;mb<=el_not;mb++) {
		if((mb-m-k)>=1) {
         t=lambda4*((re[m]*im[mb]-im[m]*re[mb])*re[mb-m-k]-(re[m]*re[mb]+im[m]*im[mb])*im[mb-m-k])/6.;
 
	if(m>k) g[k+el_not]+=(alp[el_not+k+(el_not+m)*na+(el_not-mb)*nanb+(el_not-k-m+mb)*nanbnc]*double(1-2*((k+m+mb)%2))+2*alp[el_not+k+(el_not-mb)*na+(el_not+m)*nanb+(el_not+mb-m-k)*nanbnc]*double(1-2*(k%2)))*t;
	if(m<k) g[k+el_not]+=(alp[el_not+m+(el_not+k)*na+(el_not-mb)*nanb+(el_not-k-m+mb)*nanbnc]*double(1-2*((k+m+mb)%2))+2*alp[el_not+m+(el_not-mb)*na+(el_not+k)*nanb+(el_not+mb-m-k)*nanbnc]*double(1-2*(m%2)))*t;
       if(m==k) g[k+el_not]+=(alp[el_not+k+(el_not+k)*na+(el_not-mb)*nanb+(el_not-k-k+mb)*nanbnc]*double(1-2*((k+k+mb)%2))+2*alp[el_not+k+(el_not-mb)*na+(el_not+k)*nanb+(el_not+mb-k-k)*nanbnc]*double(1-2*(k%2)))*t;

				}
	}
	}

		} // end of the Ik loop


	t=0;
	for (int i=0;i<=2*el_not;i++) {
		t+=g[i]*g[i];
	}
	gg=sqrt(t);

	return gg;
}

//HESSIAN MATRIX CALCULATION NUMERICAL APPROXIMATION
double nhessian(Vec_I_DP &x, double *hess,  DP funk(Vec_I_DP &, double* Gaunt_matrix, double* Gaunt_4th))
{

        DP ptempl[(2*el_not+1)];
        DP ptempr[(2*el_not+1)];
        DP ptempc[(2*el_not+1)];
        DP ptempc2[(2*el_not+1)];
        for (int i=0; i<=2*el_not; i++)
        {
                ptempl[i]=x[i];
                ptempr[i]=x[i];
                ptempc[i]=x[i];
                ptempc2[i]=x[i];
        }
        double del=1e-5;

        for(int k=0; k<=2*el_not; k++)
        {
                ptempr[k]=ptempr[k]+del;
                ptempl[k]=ptempl[k]-del;

                Vec_I_DP tempr(ptempr,2*el_not+1);
                Vec_I_DP templ(ptempl,2*el_not+1);
                Vec_I_DP tempc(ptempc,2*el_not+1);
                hess[k+k*hs]=(funk(tempr,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order)+funk(templ,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order)-2*funk(tempc,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order))/(del*del);
                ptempr[k]=x[k];
                ptempl[k]=x[k];

        }

        for(int k=0; k<=2*el_not; k++) {
         for(int g=k+1; g<=2*el_not; g++) {
                ptempr[k]=ptempr[k]+del;
                ptempr[g]=ptempr[g]+del;
                ptempl[k]=ptempl[k]-del;
                ptempl[g]=ptempl[g]-del;
                ptempc[k]=ptempc[k]+del;
                ptempc[g]=ptempc[g]-del;
                ptempc2[k]=ptempc2[k]-del;
                ptempc2[g]=ptempc2[g]+del;
		
		Vec_I_DP tempr(ptempr,2*el_not+1);
                Vec_I_DP templ(ptempl,2*el_not+1);
                Vec_I_DP tempc(ptempc,2*el_not+1);
                Vec_I_DP tempc2(ptempc2,2*el_not+1);
        hess[k+g*hs]=(funk(tempr,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order)+funk(templ,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order)-funk(tempc,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order)-funk(tempc2,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order))/(4*del*del);
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

//HESSIAN MATRIX CALCULATION 
double hessian(Vec_I_DP &x, double *hess, double *g3, double* alp)
{

        double gg=0;
        double t;
        vector<double> re;
        re.resize(el_not+1);
        vector<double> im;
        im.resize(el_not+1);
        re[0]=-1000000;
        im[0]=-1000000;
        for (int i=1; i<=el_not;i++){
                re[i]=x[i];
                im[i]=x[i+el_not];
                }

        double ta;

// diagonal components
// 0 0 component
        hess[0+0*hs] = tau+lambda3*g3[el_not+el_not*Na+el_not*NaNb]*x[0]+(lambda4/2.)*alp[el_not+el_not*na+el_not*nanb+el_not*nanbnc]*x[0]*x[0];
        for (int m=1; m<=el_not; m++) hess[0+0*hs]+=(lambda4/3.)*(alp[el_not+el_not*na+(el_not+m)*nanb+(el_not-m)*nanbnc]*double(1-2*(m%2))+2*alp[el_not+(el_not+m)*na+(el_not)*nanb+(el_not-m)*nanbnc])*(re[m]*re[m]+im[m]*im[m]);
// k k component for k>=1
        for (int k=1; k<=el_not; k++) {
        ta=2*tau+2*lambda3*g3[el_not+k+(el_not-k)*Na+el_not*NaNb]*double(1-2*(k%2))*x[0]+(lambda4/3.)*(alp[el_not+el_not*na+(el_not+k)*nanb+(el_not-k)*nanbnc]*double(1-2*(k%2))+2*alp[el_not+(el_not+k)*na+(el_not)*nanb+(el_not-k)*nanbnc])*x[0]*x[0];
        hess[k+k*hs]=ta;
        hess[k+el_not+(k+el_not)*hs]=ta;

        if(2*k<=el_not) {
        hess[k+k*hs]+=(2*lambda3*g3[el_not+k+(el_not+k)*Na+(el_not-2*k)*NaNb]+(2*x[0]*lambda4/3.)*(alp[el_not+k+(el_not+k)*na+(el_not)*nanb+(el_not-2*k)*nanbnc]+2*alp[el_not+(el_not+k)*na+(el_not+k)*nanb+(el_not-2*k)*nanbnc]*double(1-2*(k%2))))*re[2*k];
hess[k+el_not+(k+el_not)*hs]-=(2*lambda3*g3[el_not+k+(el_not+k)*Na+(el_not-2*k)*NaNb]+(2*x[0]*lambda4/3.)*(alp[el_not+k+(el_not+k)*na+(el_not)*nanb+(el_not-2*k)*nanbnc]+2*alp[el_not+(el_not+k)*na+(el_not+k)*nanb+(el_not-2*k)*nanbnc]*double(1-2*(k%2))))*re[2*k];
                }

        for(int m=1;m<=el_not;m++) {
        ta=(2*lambda4/3.)*(re[m]*re[m]+im[m]*im[m]);
         hess[k+k*hs]+=(alp[el_not+k+(el_not+m)*na+(el_not-k)*nanb+(el_not-m)*nanbnc]+alp[el_not+k+(el_not-k)*na+(el_not+m)*nanb+(el_not-m)*nanbnc]*double(1-2*((k+m)%2))+alp[el_not+k+(el_not-m)*na+(el_not+m)*nanb+(el_not-k)*nanbnc])*ta;
 hess[k+el_not+(k+el_not)*hs]+=(alp[el_not+k+(el_not+m)*na+(el_not-k)*nanb+(el_not-m)*nanbnc]+alp[el_not+k+(el_not-k)*na+(el_not+m)*nanb+(el_not-m)*nanbnc]*double(1-2*((k+m)%2))+alp[el_not+k+(el_not-m)*na+(el_not+m)*nanb+(el_not-k)*nanbnc])*ta;

	if(2*k+m<=el_not) {
        ta=(lambda4/6.)*(re[m]*re[2*k+m]+im[m]*im[2*k+m]);
        hess[k+k*hs]+=(alp[el_not+k+(el_not+k)*na+(el_not+m)*nanb+(el_not-2*k-m)*nanbnc]*double(1-2*(m%2))+4*alp[el_not+k+(el_not+m)*na+(el_not+2*k+m)*nanb+(el_not-k)*nanbnc]*double(1-2*(k%2)))*ta;
hess[k+el_not+(k+el_not)*hs]-=(alp[el_not+k+(el_not+k)*na+(el_not+m)*nanb+(el_not-2*k-m)*nanbnc]*double(1-2*(m%2))+4*alp[el_not+k+(el_not+m)*na+(el_not+2*k+m)*nanb+(el_not-k)*nanbnc]*double(1-2*(k%2)))*ta;
        if(m>k){ hess[k+k*hs]+=4*double(1-2*(k%2))*alp[el_not+k+(el_not-m-2*k)*na+(el_not+m)*nanb+(el_not+k)*nanbnc]*ta;
         hess[k+el_not+(k+el_not)*hs]-=4*double(1-2*(k%2))*alp[el_not+k+(el_not-m-2*k)*na+(el_not+m)*nanb+(el_not+k)*nanbnc]*ta; }
        if(m<k) { hess[k+k*hs]+=4*double(1-2*(m%2))*alp[el_not+m+(el_not-m-2*k)*na+(el_not+k)*nanb+(el_not+k)*nanbnc]*ta;
          hess[k+el_not+(k+el_not)*hs]-=4*double(1-2*(m%2))*alp[el_not+m+(el_not-m-2*k)*na+(el_not+k)*nanb+(el_not+k)*nanbnc]*ta;}
       if(m==k) { hess[k+k*hs]+=4*double(1-2*(m%2))*alp[el_not+m+(el_not-m-2*m)*na+(el_not+m)*nanb+(el_not+m)*nanbnc]*ta;
          hess[k+el_not+(k+el_not)*hs]-=4*double(1-2*(m%2))*alp[el_not+m+(el_not-m-2*m)*na+(el_not+m)*nanb+(el_not+m)*nanbnc]*ta; }
                }

        if((2*k-m)<=el_not && (2*k-m)>=1) {
                        ta=(lambda4/6.)*(re[m]*re[2*k-m]-im[m]*im[2*k-m]);
                        if(m<k) { hess[k+k*hs]+=2*(alp[el_not+m+(el_not+2*k-m)*na+(el_not+k)*nanb+(el_not+k)*nanbnc]+2*alp[el_not+m+(el_not-k)*na+(el_not+2*k-m)*nanb+(el_not-k)*nanbnc]*double(1-2*((k+m)%2)))*ta;
                          hess[k+el_not+(k+el_not)*hs]-=2*(alp[el_not+m+(el_not+2*k-m)*na+(el_not+k)*nanb+(el_not+k)*nanbnc]+2*alp[el_not+m+(el_not-k)*na+(el_not+2*k-m)*nanb+(el_not-k)*nanbnc]*double(1-2*((k+m)%2)))*ta; }
                        if(m==k) { hess[k+k*hs]+=(alp[el_not+m+(el_not+2*k-m)*na+(el_not+k)*nanb+(el_not+k)*nanbnc]+2*alp[el_not+m+(el_not-k)*na+(el_not+2*k-m)*nanb+(el_not-k)*nanbnc]*double(1-2*((k+m)%2)))*ta;
                           hess[k+el_not+(k+el_not)*hs]-=(alp[el_not+m+(el_not+2*k-m)*na+(el_not+k)*nanb+(el_not+k)*nanbnc]+2*alp[el_not+m+(el_not-k)*na+(el_not+2*k-m)*nanb+(el_not-k)*nanbnc]*double(1-2*((k+m)%2)))*ta; }
                  if(m<=(2*k-1)) { hess[k+k*hs]+=(alp[el_not+k+(el_not+k)*na+(el_not+m)*nanb+(el_not+2*k-m)*nanbnc]+2*alp[el_not+k+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-2*k)*nanbnc]*double(1-2*((k+m)%2)))*ta;
                           hess[k+el_not+(k+el_not)*hs]-=(alp[el_not+k+(el_not+k)*na+(el_not+m)*nanb+(el_not+2*k-m)*nanbnc]+2*alp[el_not+k+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-2*k)*nanbnc]*double(1-2*((k+m)%2)))*ta; }
        }
        if((m-2*k)<=el_not && (m-2*k)>=1) {
                ta=(lambda4/6.)*(re[m]*re[m-2*k]+im[m]*im[m-2*k]);
                    if(m>=2*k+1) { hess[k+k*hs]+=(alp[el_not+k+(el_not+k)*na+(el_not+m)*nanb+(el_not+2*k-m)*nanbnc]*double(1-2*(m%2))+2*alp[el_not+k+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-2*k)*nanbnc]*double(1-2*(k%2)))*ta;
                           hess[k+el_not+(k+el_not)*hs]-=(alp[el_not+k+(el_not+k)*na+(el_not+m)*nanb+(el_not+2*k-m)*nanbnc]*double(1-2*(m%2))+2*alp[el_not+k+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-2*k)*nanbnc]*double(1-2*(k%2)))*ta; }

        }
                }
}
	// offdiagonal components dRkdc0 and dIkdc0

for (int k=1; k<=el_not; k++) {

        ta=2*lambda3*g3[el_not+k+(el_not-k)*Na+el_not*NaNb]*double(1-2*(k%2))+(2*lambda4/3)*(alp[el_not+el_not*na+(el_not+k)*nanb+(el_not-k)*nanbnc]*double(1-2*(k%2))+2*alp[el_not+(el_not+k)*na+el_not*nanb+(el_not-k)*nanbnc])*x[0];
        hess[0+k*hs]=ta*re[k];
        hess[k+0*hs]=ta*re[k];
      hess[0+(k+el_not)*hs]=ta*im[k];
      hess[(k+el_not)+0*hs]=ta*im[k];

        for (int m=1; m<=el_not; m++) {
        if (m+k<=el_not) {
        ta=(lambda4/3.)*(alp[el_not+(el_not+k)*na+(el_not+m)*nanb+(el_not-k-m)*nanbnc]*double(1-2*(m%2))+alp[el_not+k+(el_not+m)*na+(el_not)*nanb+(el_not-k-m)*nanbnc]+alp[el_not+(el_not+m)*na+(el_not+k)*nanb+(el_not-m-k)*nanbnc]*double(1-2*(k%2)));
        hess[0+k*hs]+=ta*(re[m]*re[m+k]+im[m]*im[m+k]);
        hess[k+0*hs]+=ta*(re[m]*re[m+k]+im[m]*im[m+k]);
      hess[0+(k+el_not)*hs]+=ta*(re[m]*im[m+k]-im[m]*re[m+k]);
      hess[(k+el_not)+0*hs]+=ta*(re[m]*im[m+k]-im[m]*re[m+k]);

        if(m>k) { ta=(lambda4/3.)*alp[el_not+k+(el_not-k-m)*na+(el_not+m)*nanb+el_not*nanbnc]*double(1-2*(k%2));
        hess[0+k*hs]+=ta*(re[m]*re[m+k]+im[m]*im[m+k]);
        hess[k+0*hs]+=ta*(re[m]*re[m+k]+im[m]*im[m+k]);
      hess[0+(k+el_not)*hs]+=ta*(re[m]*im[m+k]-im[m]*re[m+k]);
      hess[(k+el_not)+0*hs]+=ta*(re[m]*im[m+k]-im[m]*re[m+k]);
                }
        if(m<k) { ta=(lambda4/3.)*alp[el_not+m+(el_not-k-m)*na+(el_not+k)*nanb+el_not*nanbnc]*double(1-2*(m%2));
        hess[0+k*hs]+=ta*(re[m]*re[m+k]+im[m]*im[m+k]);
        hess[k+0*hs]+=ta*(re[m]*re[m+k]+im[m]*im[m+k]);
      hess[0+(k+el_not)*hs]+=ta*(re[m]*im[m+k]-im[m]*re[m+k]);
      hess[(k+el_not)+0*hs]+=ta*(re[m]*im[m+k]-im[m]*re[m+k]);
                }
       if(m==k) { ta=(lambda4/3.)*alp[el_not+m+(el_not-m-m)*na+(el_not+m)*nanb+el_not*nanbnc]*double(1-2*(m%2));
        hess[0+k*hs]+=ta*(re[m]*re[m+k]+im[m]*im[m+k]);
        hess[k+0*hs]+=ta*(re[m]*re[m+k]+im[m]*im[m+k]);
      hess[0+(k+el_not)*hs]+=ta*(re[m]*im[m+k]-im[m]*re[m+k]);
      hess[(k+el_not)+0*hs]+=ta*(re[m]*im[m+k]-im[m]*re[m+k]);
                }
        }
	if(m-k>=1) {
        ta=(lambda4/6.)*(alp[el_not+(el_not+k)*na+(el_not-m)*nanb+(el_not+m-k)*nanbnc]*double(1-2*((k+m)%2))+2*alp[el_not+(el_not+m)*na+(el_not-k)*nanb+(el_not+k-m)*nanbnc]+alp[el_not+k+(el_not-m)*na+(el_not)*nanb+(el_not+m-k)*nanbnc]*double(1-2*(k%2)));
        hess[0+k*hs]+=ta*(re[m]*re[m-k]+im[m]*im[m-k]);
        hess[k+0*hs]+=ta*(re[m]*re[m-k]+im[m]*im[m-k]);
     hess[0+(k+el_not)*hs]+=ta*(im[m]*re[m-k]-re[m]*im[m-k]);
     hess[(k+el_not)+0*hs]+=ta*(im[m]*re[m-k]-re[m]*im[m-k]);
        }

        if(1<=k-m) {
     ta=(lambda4/6.)*(2*alp[el_not+(el_not+m)*na+(el_not-k)*nanb+(el_not+k-m)*nanbnc]*double(1-2*((k+m)%2))+alp[el_not+(el_not+k)*na+(el_not-m)*nanb+(el_not+m-k)*nanbnc]+alp[el_not+m+(el_not-k)*na+(el_not)*nanb+(el_not+k-m)*nanbnc]*double(1-2*(m%2)));
        hess[0+k*hs]+=ta*(re[m]*re[k-m]-im[m]*im[k-m]);
        hess[k+0*hs]+=ta*(re[m]*re[k-m]-im[m]*im[k-m]);
     hess[0+(k+el_not)*hs]+=ta*(im[m]*re[k-m]+re[m]*im[k-m]);
     hess[(k+el_not)+0*hs]+=ta*(im[m]*re[k-m]+re[m]*im[k-m]);
        }

        if(2*m<=k) {
        ta=(lambda4/3.)*(alp[el_not+m+(el_not+k-m)*na+(el_not)*nanb+(el_not-k)*nanbnc]+alp[el_not+m+(el_not-k)*na+(el_not+k-m)*nanb+el_not*nanbnc]*double(1-2*(m%2)));
        if(2*m==k) ta=ta/2.;
        hess[0+k*hs]+=ta*(re[m]*re[k-m]-im[m]*im[k-m]);
        hess[k+0*hs]+=ta*(re[m]*re[k-m]-im[m]*im[k-m]);
     hess[0+(k+el_not)*hs]+=ta*(im[m]*re[k-m]+re[m]*im[k-m]);
     hess[(k+el_not)+0*hs]+=ta*(im[m]*re[k-m]+re[m]*im[k-m]);
                }

        }

        }

	//off diagonals /dRkdIg /dRgdIk

for (int k=1; k<=el_not; k++) {
        for (int g=1; g<=el_not; g++) {

        hess[g+(el_not+k)*hs]=0;
        hess[(el_not+g)+k*hs]=0;

        // g+k term
        if(g+k<=el_not) {
        ta=2*lambda3*g3[el_not+g+(el_not+k)*Na+(el_not-g-k)*NaNb]*double(1-2*((k+g)%2));
        ta+=(x[0]*lambda4/6.)*(3*alp[el_not+(el_not+k)*na+(el_not-k-g)*nanb+(el_not+g)*nanbnc]*double(1-2*(g%2))+4*alp[el_not+(el_not+k+g)*na+(el_not-k)*nanb+(el_not-g)*nanbnc]+3*alp[el_not+k+(el_not-k-g)*na+el_not*nanb+(el_not+g)*nanbnc]*double(1-2*(k%2)));

        if(g>k) ta+=(x[0]*lambda4/6.)*2*alp[el_not+k+(el_not-k-g)*na+(el_not+g)*nanb+el_not*nanbnc]*double(1-2*(k%2));
        if(g<k) ta+=(x[0]*lambda4/6.)*2*alp[el_not+g+(el_not-k-g)*na+(el_not+k)*nanb+el_not*nanbnc]*double(1-2*(g%2));
       if(g==k) ta+=(x[0]*lambda4/6.)*2*alp[el_not+g+(el_not-g-g)*na+(el_not+g)*nanb+el_not*nanbnc]*double(1-2*(g%2));

        hess[g+(el_not+k)*hs]+=ta*im[g+k];
        hess[(el_not+g)+k*hs]+=ta*im[g+k];
        }

        // k-g term
        if( k-g>=1 && k-g<=el_not ) {
        ta=2*lambda3*g3[el_not+k-g+(el_not+g)*Na+(el_not-k)*NaNb]*double(1-2*(k%2));
        hess[g+(el_not+k)*hs]+=ta*im[k-g];
        hess[(el_not+g)+k*hs]-=ta*im[k-g];

        ta=(x[0]*lambda4/6.)*(3*alp[el_not+(el_not+g)*na+(el_not-k)*nanb+(el_not+k-g)*nanbnc]*double(1-2*((g+k)%2))+4*alp[el_not+(el_not+k)*na+(el_not-g)*nanb+(el_not+g-k)*nanbnc]+3*alp[el_not+g+(el_not-k)*na+el_not*nanb+(el_not+k-g)*nanbnc]*double(1-2*(g%2)));

        if(2*g<k) ta+=(x[0]*lambda4/6.)*2*alp[el_not+g+(el_not-k)*na+(el_not+k-g)*nanb+el_not*nanbnc]*double(1-2*(g%2));
        if(2*g>k) ta+=(x[0]*lambda4/6.)*2*alp[el_not+k-g+(el_not-k)*na+(el_not+g)*nanb+el_not*nanbnc]*double(1-2*((k+g)%2));
       if(2*g==k) ta+=(x[0]*lambda4/6.)*2*alp[el_not+g+(el_not-k)*na+(el_not+g)*nanb+el_not*nanbnc]*double(1-2*(g%2));

        hess[g+(el_not+k)*hs]+=ta*im[k-g];
        hess[(el_not+g)+k*hs]-=ta*im[k-g];
        }
	// g-k term

        if ( k<g && g-k<=el_not) {

        ta=2*lambda3*g3[el_not+g-k+(el_not+k)*Na+(el_not-g)*NaNb]*double(1-2*(g%2));
        hess[g+(el_not+k)*hs]-=ta*im[g-k];
        hess[(el_not+g)+k*hs]+=ta*im[g-k];


        ta=(x[0]*lambda4/6.)*(3*alp[el_not+(el_not+k)*na+(el_not-g)*nanb+(el_not+g-k)*nanbnc]*double(1-2*((g+k)%2))+4*alp[el_not+(el_not+g)*na+(el_not-k)*nanb+(el_not+k-g)*nanbnc]+3*alp[el_not+k+(el_not-g)*na+el_not*nanb+(el_not+g-k)*nanbnc]*double(1-2*(k%2)));

        if(2*k<g) ta+=(x[0]*lambda4/6.)*2*alp[el_not+k+(el_not-g)*na+(el_not+g-k)*nanb+el_not*nanbnc]*double(1-2*(k%2));
        if(2*k>g) ta+=(x[0]*lambda4/6.)*2*alp[el_not+g-k+(el_not-g)*na+(el_not+k)*nanb+el_not*nanbnc]*double(1-2*((k+g)%2));
       if(2*k==g) ta+=(x[0]*lambda4/6.)*2*alp[el_not+k+(el_not-g)*na+(el_not+k)*nanb+el_not*nanbnc]*double(1-2*(k%2));

        hess[g+(el_not+k)*hs]-=ta*im[g-k];
        hess[(el_not+g)+k*hs]+=ta*im[g-k];

        }

        for(int m=1; m<=el_not; m++) {
        // Term 1 of 7
        if(k+g+m<=el_not) {
        ta=(lambda4/6.)*(alp[el_not+k+(el_not+g)*na+(el_not+m)*nanb+(el_not-m-g-k)*nanbnc]*double(1-2*(m%2))+2*alp[el_not+g+(el_not+m)*na+(el_not+k+g+m)*nanb+(el_not-k)*nanbnc]*double(1-2*(k%2))+2*alp[el_not+k+(el_not+m)*na+(el_not+g)*nanb+(el_not-k-m-g)*nanbnc]*double(1-2*(g%2)));

        if(m>g) ta+=(lambda4/6.)*2*alp[el_not+g+(el_not-g-m-k)*na+(el_not+m)*nanb+(el_not+k)*nanbnc]*double(1-2*(g%2));
        if(m<g) ta+=(lambda4/6.)*2*alp[el_not+m+(el_not-g-m-k)*na+(el_not+g)*nanb+(el_not+k)*nanbnc]*double(1-2*(m%2));
      if(m==g)  ta+=(lambda4/6.)*2*alp[el_not+g+(el_not-g-g-k)*na+(el_not+g)*nanb+(el_not+k)*nanbnc]*double(1-2*(g%2));
        if(m>k) ta+=(lambda4/6.)*2*alp[el_not+k+(el_not-g-m-k)*na+(el_not+m)*nanb+(el_not+g)*nanbnc]*double(1-2*(k%2));
        if(m<k) ta+=(lambda4/6.)*2*alp[el_not+m+(el_not-g-m-k)*na+(el_not+k)*nanb+(el_not+g)*nanbnc]*double(1-2*(m%2));
      if(m==k)  ta+=(lambda4/6.)*2*alp[el_not+k+(el_not-g-k-k)*na+(el_not+k)*nanb+(el_not+g)*nanbnc]*double(1-2*(k%2));

        hess[g+(el_not+k)*hs]+=ta*(re[m]*im[k+g+m]-im[m]*re[k+g+m]);
        hess[(el_not+g)+k*hs]+=ta*(re[m]*im[k+g+m]-im[m]*re[k+g+m]);
                }
	// Term 2 of 7
        if(g+m-k<=el_not && g+m-k>=1) {
        ta=(lambda4/3.)*(alp[el_not+g+(el_not+m)*na+(el_not+k)*nanb+(el_not+m+g-k)*nanbnc]+alp[el_not+g+(el_not-k)*na+(el_not+m)*nanb+(el_not+k-g-m)*nanbnc]*double(1-2*((k+m)%2))+alp[el_not+g+(el_not+k-g-m)*na+(el_not+m)*nanb+(el_not-k)*nanbnc]*double(1-2*((g+k)%2)));
        hess[g+(el_not+k)*hs]+=ta*(im[m]*re[g-k+m]-re[m]*im[g-k+m]);
        hess[(el_not+g)+k*hs]+=ta*(re[m]*im[g-k+m]-im[m]*re[g-k+m]);
        }
        // Term 3 of 7
        if(2*m<=k+g && g+k-m<=el_not) {
        ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+k+g-m)*na+(el_not+k)*nanb+(el_not+g)*nanbnc]+alp[el_not+m+(el_not-k)*na+(el_not+k+g-m)*nanb+(el_not-g)*nanbnc]*double(1-2*((g+m)%2))+alp[el_not+m+(el_not-g)*na+(el_not+k+g-m)*nanb+(el_not-k)*nanbnc]*double(1-2*((m+k)%2)));
        if(2*m==k+g) ta=ta/2.;
        hess[g+(el_not+k)*hs]+=ta*(re[m]*im[g+k-m]+im[m]*re[g+k-m]);
        hess[(el_not+g)+k*hs]+=ta*(re[m]*im[g+k-m]+im[m]*re[g+k-m]);
        }
        if(1<=k+g-m && g+k-m<=el_not) {
        if(k<g) ta=(lambda4/6.)*(alp[el_not+k+(el_not+g)*na+(el_not+m)*nanb+(el_not+k+g-m)*nanbnc]+2*alp[el_not+k+(el_not-m)*na+(el_not+g)*nanb+(el_not+m-k-g)*nanbnc]*double(1-2*((m+g)%2)));
        if(k>g) ta=(lambda4/6.)*(alp[el_not+k+(el_not+g)*na+(el_not+m)*nanb+(el_not+k+g-m)*nanbnc]+2*alp[el_not+g+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-k-g)*nanbnc]*double(1-2*((m+k)%2)));
       if(k==g) ta=(lambda4/6.)*(alp[el_not+k+(el_not+k)*na+(el_not+m)*nanb+(el_not+k+k-m)*nanbnc]+2*alp[el_not+k+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-k-k)*nanbnc]*double(1-2*((m+k)%2)));
        hess[g+(el_not+k)*hs]+=ta*(re[m]*im[g+k-m]+im[m]*re[g+k-m]);
        hess[(el_not+g)+k*hs]+=ta*(re[m]*im[g+k-m]+im[m]*re[g+k-m]);
        }
        // Term 4 of 7
        if(2*m<=k-g && k-g-m>=1) {
        ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+k-g-m)*na+(el_not+k)*nanb+(el_not-g)*nanbnc]*double(1-2*(g%2))+alp[el_not+m+(el_not-k)*na+(el_not+k-g-m)*nanb+(el_not+g)*nanbnc]*double(1-2*(m%2)));
        if(2*m==k-g) ta=ta/2.;
        hess[g+(el_not+k)*hs]+=ta*(re[m]*im[k-g-m]+im[m]*re[k-g-m]);
        hess[(el_not+g)+k*hs]-=ta*(re[m]*im[k-g-m]+im[m]*re[k-g-m]);
        }
        if(k-g-m<=el_not && k-g-m>=1) {
        if(m>g) ta=(lambda4/6.)*2*(alp[el_not+g+(el_not+m)*na+(el_not+k)*nanb+(el_not+g+m-k)*nanbnc]*double(1-2*((g+m+k)%2))+alp[el_not+g+(el_not-k)*na+(el_not+m)*nanb+(el_not+k-g-m)*nanbnc]*double(1-2*(g%2)));
        if(m<g) ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+g)*na+(el_not+k)*nanb+(el_not+m+g-k)*nanbnc]*double(1-2*((g+m+k)%2))+alp[el_not+m+(el_not-k)*na+(el_not+g)*nanb+(el_not+k-m-g)*nanbnc]*double(1-2*(m%2)));
       if(m==g) ta=(lambda4/6.)*2*(alp[el_not+g+(el_not+g)*na+(el_not+k)*nanb+(el_not+g+g-k)*nanbnc]*double(1-2*((g+g+k)%2))+alp[el_not+g+(el_not-k)*na+(el_not+g)*nanb+(el_not+k-g-g)*nanbnc]*double(1-2*(g%2)));
        hess[g+(el_not+k)*hs]+=ta*(re[m]*im[k-g-m]+im[m]*re[k-g-m]);
        hess[(el_not+g)+k*hs]-=ta*(re[m]*im[k-g-m]+im[m]*re[k-g-m]);
        }
	// Term 5 of 7
        if(k+m-g<=el_not && k+m-g>=1) {

        if(m>=max(k,g-k+1)) {
        ta=(lambda4/6.)*2*(alp[el_not+k+(el_not+m)*na+(el_not+g)*nanb+(el_not+k+m-g)*nanbnc]+alp[el_not+k+(el_not+g-k-m)*na+(el_not+m)*nanb+(el_not-g)*nanbnc]*double(1-2*((g+k)%2))+alp[el_not+k+(el_not-g)*na+(el_not+m)*nanb+(el_not+g-k-m)*nanbnc]*double(1-2*((m+g)%2)));
        if(m==k) ta=ta/2.;
        hess[g+(el_not+k)*hs]-=ta*(im[m]*re[k-g+m]-re[m]*im[k-g+m]);
        hess[(el_not+g)+k*hs]+=ta*(im[m]*re[k-g+m]-re[m]*im[k-g+m]);
                }
        if(m<=k && m>=max(1,g-k+1)) {
        ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+k)*na+(el_not+g)*nanb+(el_not+k+m-g)*nanbnc]+alp[el_not+m+(el_not+g-k-m)*na+(el_not+k)*nanb+(el_not-g)*nanbnc]*double(1-2*((g+m)%2))+alp[el_not+m+(el_not-g)*na+(el_not+k)*nanb+(el_not+g-k-m)*nanbnc]*double(1-2*((k+g)%2)));
        if(m==k) ta=ta/2.;
        hess[g+(el_not+k)*hs]-=ta*(im[m]*re[k-g+m]-re[m]*im[k-g+m]);
        hess[(el_not+g)+k*hs]+=ta*(im[m]*re[k-g+m]-re[m]*im[k-g+m]);
        }
        }
        // Term 6 of 7
        if(m-k-g>=1) {
        if(g>k) ta=(lambda4/6.)*(alp[el_not+k+(el_not+g)*na+(el_not+m)*nanb+(el_not+k+g-m)*nanbnc]*double(1-2*((k+g+m)%2))+2*alp[el_not+k+(el_not-m)*na+(el_not+g)*nanb+(el_not+m-k-g)*nanbnc]*double(1-2*(k%2)));
        if(g<k) ta=(lambda4/6.)*(alp[el_not+g+(el_not+k)*na+(el_not+m)*nanb+(el_not+k+g-m)*nanbnc]*double(1-2*((k+g+m)%2))+2*alp[el_not+g+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-k-g)*nanbnc]*double(1-2*(g%2)));
       if(g==k) ta=(lambda4/6.)*(alp[el_not+g+(el_not+g)*na+(el_not+m)*nanb+(el_not+g+g-m)*nanbnc]*double(1-2*((g+g+m)%2))+2*alp[el_not+g+(el_not-m)*na+(el_not+g)*nanb+(el_not+m-g-g)*nanbnc]*double(1-2*(g%2)));
        hess[g+(el_not+k)*hs]+=ta*(im[m]*re[m-k-g]-re[m]*im[m-k-g]);
        hess[(el_not+g)+k*hs]+=ta*(im[m]*re[m-k-g]-re[m]*im[m-g-k]);
        }
	// Term 7 of 7
        if(g-k-m>=1) {
        ta=(lambda4/6.)*alp[el_not+k+(el_not+m)*na+(el_not+g-k-m)*nanb+(el_not-g)*nanbnc]*double(1-2*((g+k+m)%2));
        hess[g+(el_not+k)*hs]-=ta*(im[m]*re[g-k-m]+re[m]*im[g-k-m]);
        hess[(el_not+g)+k*hs]+=ta*(im[m]*re[g-k-m]+re[m]*im[g-k-m]);

        if(2*m<=g-k) {
        ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+g-k-m)*na+(el_not+g)*nanb+(el_not-k)*nanbnc]*double(1-2*(k%2))+alp[el_not+m+(el_not-g)*na+(el_not+g-k-m)*nanb+(el_not+k)*nanbnc]*double(1-2*(m%2)));
        if(2*m==g-k) ta=ta/2.;
        hess[g+(el_not+k)*hs]-=ta*(im[m]*re[g-k-m]+re[m]*im[g-k-m]);
        hess[(el_not+g)+k*hs]+=ta*(im[m]*re[g-k-m]+re[m]*im[g-k-m]);
        }

        if(m>k) ta=(lambda4/6.)*(alp[el_not+k+(el_not+m)*na+(el_not+g)*nanb+(el_not+k+m-g)*nanbnc]*double(1-2*((m+g+k)%2))+2*alp[el_not+k+(el_not-g)*na+(el_not+m)*nanb+(el_not+g-k-m)*nanbnc]*double(1-2*(k%2)));
        if(m<k) ta=(lambda4/6.)*(alp[el_not+m+(el_not+k)*na+(el_not+g)*nanb+(el_not+m+k-g)*nanbnc]*double(1-2*((m+g+k)%2))+2*alp[el_not+m+(el_not-g)*na+(el_not+k)*nanb+(el_not+g-k-m)*nanbnc]*double(1-2*(m%2)));
       if(m==k) ta=(lambda4/6.)*(alp[el_not+k+(el_not+k)*na+(el_not+g)*nanb+(el_not+k+k-g)*nanbnc]*double(1-2*((k+g+k)%2))+2*alp[el_not+k+(el_not-g)*na+(el_not+k)*nanb+(el_not+g-k-k)*nanbnc]*double(1-2*(k%2)));
        hess[g+(el_not+k)*hs]-=ta*(im[m]*re[g-k-m]+re[m]*im[g-k-m]);
        hess[(el_not+g)+k*hs]+=ta*(im[m]*re[g-k-m]+re[m]*im[g-k-m]);

        }

        } //end of loop over m


        }       } //end of RkIg loops

// cross terms dRgdRk and dIgIk

        for (int k=1; k<=el_not; k++) {
                for (int g=k+1; g<=el_not; g++) {

        hess[g+k*hs]=0;
        hess[k+g*hs]=0;
        hess[k+el_not+(g+el_not)*hs]=0;
        hess[g+el_not+(k+el_not)*hs]=0;

	// g+k term

        if(g+k<=el_not) {
        ta=2*lambda3*g3[el_not+g+(el_not+k)*Na+(el_not-g-k)*NaNb]*double(1-2*((k+g)%2));
        ta+=(x[0]*lambda4/6.)*(3*alp[el_not+(el_not+k)*na+(el_not-k-g)*nanb+(el_not+g)*nanbnc]*double(1-2*(g%2))+4*alp[el_not+(el_not+k+g)*na+(el_not-k)*nanb+(el_not-g)*nanbnc]+3*alp[el_not+k+(el_not-k-g)*na+el_not*nanb+(el_not+g)*nanbnc]*double(1-2*(k%2)));
        ta+=(x[0]*lambda4/6.)*2*alp[el_not+k+(el_not-k-g)*na+(el_not+g)*nanb+el_not*nanbnc]*double(1-2*(k%2));

        hess[g+k*hs]+=ta*re[g+k];
        hess[k+g*hs]+=ta*re[g+k];
        hess[el_not+g+(k+el_not)*hs]-=ta*re[g+k];
        hess[el_not+k+(g+el_not)*hs]-=ta*re[g+k];
        }

        //g-k term

        ta=2*lambda3*g3[el_not+g-k+(el_not+k)*Na+(el_not-g)*NaNb]*double(1-2*(g%2));
        hess[g+k*hs]+=ta*re[g-k];
        hess[k+g*hs]+=ta*re[g-k];
        hess[g+el_not+(k+el_not)*hs]+=ta*re[g-k];
        hess[k+el_not+(g+el_not)*hs]+=ta*re[g-k];

        ta=(x[0]*lambda4/6.)*(3*alp[el_not+(el_not+k)*na+(el_not+g-k)*nanb+(el_not-g)*nanbnc]*double(1-2*((g+k)%2))+4*alp[el_not+k+(el_not+g-k)*na+(el_not)*nanb+(el_not-g)*nanbnc]+3*alp[el_not+(el_not+g-k)*na+(el_not+k)*nanb+(el_not-g)*nanbnc]*double(1-2*(k%2)));

        if(2*k<g) ta+=(x[0]*lambda4/6.)*2*alp[el_not+k+(el_not-g)*na+(el_not+g-k)*nanb+el_not*nanbnc]*double(1-2*(k%2));
        if(2*k>g) ta+=(x[0]*lambda4/6.)*2*alp[el_not+g-k+(el_not-g)*na+(el_not+k)*nanb+el_not*nanbnc]*double(1-2*((k+g)%2));
       if(2*k==g) ta+=(x[0]*lambda4/6.)*2*alp[el_not+k+(el_not-g)*na+(el_not+k)*nanb+el_not*nanbnc]*double(1-2*(k%2));
        hess[g+k*hs]+=ta*re[g-k];
        hess[k+g*hs]+=ta*re[g-k];
        hess[g+el_not+(k+el_not)*hs]+=ta*re[g-k];
        hess[k+el_not+(g+el_not)*hs]+=ta*re[g-k];



        for(int m=1; m<=el_not; m++) {
	// Term 1 of 7
        if(k+g+m<=el_not) {
        ta=(lambda4/6.)*(alp[el_not+k+(el_not+g)*na+(el_not+m)*nanb+(el_not-m-g-k)*nanbnc]*double(1-2*(m%2))+2*alp[el_not+g+(el_not+m)*na+(el_not+k+g+m)*nanb+(el_not-k)*nanbnc]*double(1-2*(k%2))+2*alp[el_not+k+(el_not+m)*na+(el_not+g)*nanb+(el_not-k-m-g)*nanbnc]*double(1-2*(g%2)));

        if(m>g) ta+=(lambda4/6.)*2*alp[el_not+g+(el_not-g-m-k)*na+(el_not+m)*nanb+(el_not+k)*nanbnc]*double(1-2*(g%2));
        if(m<g) ta+=(lambda4/6.)*2*alp[el_not+m+(el_not-g-m-k)*na+(el_not+g)*nanb+(el_not+k)*nanbnc]*double(1-2*(m%2));
      if(m==g)  ta+=(lambda4/6.)*2*alp[el_not+g+(el_not-g-g-k)*na+(el_not+g)*nanb+(el_not+k)*nanbnc]*double(1-2*(g%2));
        if(m>k) ta+=(lambda4/6.)*2*alp[el_not+k+(el_not-g-m-k)*na+(el_not+m)*nanb+(el_not+g)*nanbnc]*double(1-2*(k%2));
        if(m<k) ta+=(lambda4/6.)*2*alp[el_not+m+(el_not-g-m-k)*na+(el_not+k)*nanb+(el_not+g)*nanbnc]*double(1-2*(m%2));
      if(m==k)  ta+=(lambda4/6.)*2*alp[el_not+k+(el_not-g-k-k)*na+(el_not+k)*nanb+(el_not+g)*nanbnc]*double(1-2*(k%2));

        hess[g+k*hs]+=ta*(re[m]*re[k+g+m]+im[m]*im[k+g+m]);
        hess[k+g*hs]+=ta*(re[m]*re[k+g+m]+im[m]*im[k+g+m]);
        hess[g+el_not+(k+el_not)*hs]-=ta*(re[m]*re[k+g+m]+im[m]*im[k+g+m]);
        hess[k+el_not+(g+el_not)*hs]-=ta*(re[m]*re[k+g+m]+im[m]*im[k+g+m]);
                }
        // Term 2 of 7
        if(g+m-k<=el_not && g+m-k>=1) {
        ta=(lambda4/3.)*(alp[el_not+g+(el_not+m)*na+(el_not+k)*nanb+(el_not+m+g-k)*nanbnc]+alp[el_not+g+(el_not-k)*na+(el_not+m)*nanb+(el_not+k-g-m)*nanbnc]*double(1-2*((k+m)%2))+alp[el_not+g+(el_not+k-g-m)*na+(el_not+m)*nanb+(el_not-k)*nanbnc]*double(1-2*((g+k)%2)));
        hess[g+k*hs]+=ta*(re[m]*re[g-k+m]+im[m]*im[g-k+m]);
        hess[k+g*hs]+=ta*(re[m]*re[g-k+m]+im[m]*im[g-k+m]);
        hess[g+el_not+(k+el_not)*hs]+=ta*(re[m]*re[g-k+m]+im[m]*im[g-k+m]);
        hess[k+el_not+(g+el_not)*hs]+=ta*(re[m]*re[g-k+m]+im[m]*im[g-k+m]);
        }
        // Term 3 of 7
        if(2*m<=k+g && g+k-m<=el_not) {
        ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+k+g-m)*na+(el_not+k)*nanb+(el_not+g)*nanbnc]+alp[el_not+m+(el_not-k)*na+(el_not+k+g-m)*nanb+(el_not-g)*nanbnc]*double(1-2*((g+m)%2))+alp[el_not+m+(el_not-g)*na+(el_not+k+g-m)*nanb+(el_not-k)*nanbnc]*double(1-2*((m+k)%2)));
        if(2*m==k+g) ta=ta/2.;
        hess[g+k*hs]+=ta*(re[m]*re[g+k-m]-im[m]*im[g+k-m]);
        hess[k+g*hs]+=ta*(re[m]*re[g+k-m]-im[m]*im[g+k-m]);
        hess[(el_not+g)+(k+el_not)*hs]-=ta*(re[m]*re[g+k-m]-im[m]*im[g+k-m]);
        hess[(el_not+k)+(g+el_not)*hs]-=ta*(re[m]*re[g+k-m]-im[m]*im[g+k-m]);
        }
        if(m<=k+g-1 && k+g-m<=el_not) {
        if(k<g) ta=(lambda4/6.)*(alp[el_not+k+(el_not+g)*na+(el_not+m)*nanb+(el_not+k+g-m)*nanbnc]+2*alp[el_not+k+(el_not-m)*na+(el_not+g)*nanb+(el_not+m-k-g)*nanbnc]*double(1-2*((m+g)%2)));
        if(k>g) ta=(lambda4/6.)*(alp[el_not+k+(el_not+g)*na+(el_not+m)*nanb+(el_not+k+g-m)*nanbnc]+2*alp[el_not+g+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-k-g)*nanbnc]*double(1-2*((m+k)%2)));
       if(k==g) ta=(lambda4/6.)*(alp[el_not+k+(el_not+k)*na+(el_not+m)*nanb+(el_not+k+k-m)*nanbnc]+2*alp[el_not+k+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-k-k)*nanbnc]*double(1-2*((m+k)%2)));
        hess[g+k*hs]+=ta*(re[m]*re[g+k-m]-im[m]*im[g+k-m]);
        hess[k+g*hs]+=ta*(re[m]*re[g+k-m]-im[m]*im[g+k-m]);
        hess[(el_not+g)+(k+el_not)*hs]-=ta*(re[m]*re[g+k-m]-im[m]*im[g+k-m]);
        hess[(el_not+k)+(g+el_not)*hs]-=ta*(re[m]*re[g+k-m]-im[m]*im[g+k-m]);
        }
	// Term 4 of 7
        if(2*m<=k-g && k-g-m>=1) {
        ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+k-g-m)*na+(el_not+k)*nanb+(el_not-g)*nanbnc]*double(1-2*(g%2))+alp[el_not+m+(el_not-k)*na+(el_not+k-g-m)*nanb+(el_not+g)*nanbnc]*double(1-2*(m%2)));
        if(2*m==k-g) ta=ta/2.;
        hess[g+k*hs]+=ta*(re[m]*re[k-g-m]-im[m]*im[k-g-m]);
        hess[k+g*hs]+=ta*(re[m]*re[k-g-m]-im[m]*im[k-g-m]);
        hess[g+el_not+(k+el_not)*hs]+=ta*(re[m]*re[k-g-m]-im[m]*im[k-g-m]);
        hess[k+el_not+(g+el_not)*hs]+=ta*(re[m]*re[k-g-m]-im[m]*im[k-g-m]);
        }
        if(k-g-m<=el_not && k-g-m>=1) {
        if(m>g) ta=(lambda4/6.)*2*(alp[el_not+g+(el_not+m)*na+(el_not+k)*nanb+(el_not+g+m-k)*nanbnc]*double(1-2*((g+m+k)%2))+alp[el_not+g+(el_not-k)*na+(el_not+m)*nanb+(el_not+k-g-m)*nanbnc]*double(1-2*(g%2)));
        if(m<g) ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+g)*na+(el_not+k)*nanb+(el_not+m+g-k)*nanbnc]*double(1-2*((g+m+k)%2))+alp[el_not+m+(el_not-k)*na+(el_not+g)*nanb+(el_not+k-m-g)*nanbnc]*double(1-2*(m%2)));
       if(m==g) ta=(lambda4/6.)*2*(alp[el_not+g+(el_not+g)*na+(el_not+k)*nanb+(el_not+g+g-k)*nanbnc]*double(1-2*((g+g+k)%2))+alp[el_not+g+(el_not-k)*na+(el_not+g)*nanb+(el_not+k-g-g)*nanbnc]*double(1-2*(g%2)));
        hess[g+k*hs]+=ta*(re[m]*re[k-g-m]-im[m]*im[k-g-m]);
        hess[k+g*hs]+=ta*(re[m]*re[k-g-m]-im[m]*im[k-g-m]);
        hess[g+el_not+(k+el_not)*hs]+=ta*(re[m]*re[k-g-m]-im[m]*im[k-g-m]);
        hess[k+el_not+(g+el_not)*hs]+=ta*(re[m]*re[k-g-m]-im[m]*im[k-g-m]);
        }
        // Term 5 of 7
        if(k+m-g<=el_not && k+m-g>=1) {

        if(m>=k) {
        ta=(lambda4/6.)*2*(alp[el_not+k+(el_not+m)*na+(el_not+g)*nanb+(el_not+k+m-g)*nanbnc]+alp[el_not+k+(el_not+g-k-m)*na+(el_not+m)*nanb+(el_not-g)*nanbnc]*double(1-2*((g+k)%2))+alp[el_not+k+(el_not-g)*na+(el_not+m)*nanb+(el_not+g-k-m)*nanbnc]*double(1-2*((m+g)%2)));
        if(m==k) ta=ta/2.;
        hess[g+k*hs]+=ta*(re[m]*re[k-g+m]+im[m]*im[k-g+m]);
        hess[k+g*hs]+=ta*(re[m]*re[k-g+m]+im[m]*im[k-g+m]);
        hess[g+el_not+(k+el_not)*hs]+=ta*(re[m]*re[k-g+m]+im[m]*im[k-g+m]);
        hess[k+el_not+(g+el_not)*hs]+=ta*(re[m]*re[k-g+m]+im[m]*im[k-g+m]);
                }
        if(m<=k) {
        ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+k)*na+(el_not+g)*nanb+(el_not+k+m-g)*nanbnc]+alp[el_not+m+(el_not+g-k-m)*na+(el_not+k)*nanb+(el_not-g)*nanbnc]*double(1-2*((g+m)%2))+alp[el_not+m+(el_not-g)*na+(el_not+k)*nanb+(el_not+g-k-m)*nanbnc]*double(1-2*((k+g)%2)));
        if(m==k) ta=ta/2.;
        hess[g+k*hs]+=ta*(re[m]*re[k-g+m]+im[m]*im[k-g+m]);
        hess[k+g*hs]+=ta*(re[m]*re[k-g+m]+im[m]*im[k-g+m]);
        hess[g+el_not+(k+el_not)*hs]+=ta*(re[m]*re[k-g+m]+im[m]*im[k-g+m]);
        hess[k+el_not+(g+el_not)*hs]+=ta*(re[m]*re[k-g+m]+im[m]*im[k-g+m]);
        }
        }
	// Term 6 of 7
        if(m-k-g>=1) {
        if(g>k) ta=(lambda4/6.)*(alp[el_not+k+(el_not+g)*na+(el_not+m)*nanb+(el_not+k+g-m)*nanbnc]*double(1-2*((k+g+m)%2))+2*alp[el_not+k+(el_not-m)*na+(el_not+g)*nanb+(el_not+m-k-g)*nanbnc]*double(1-2*(k%2)));
        if(g<k) ta=(lambda4/6.)*(alp[el_not+g+(el_not+k)*na+(el_not+m)*nanb+(el_not+k+g-m)*nanbnc]*double(1-2*((k+g+m)%2))+2*alp[el_not+g+(el_not-m)*na+(el_not+k)*nanb+(el_not+m-k-g)*nanbnc]*double(1-2*(g%2)));
       if(g==k) ta=(lambda4/6.)*(alp[el_not+g+(el_not+g)*na+(el_not+m)*nanb+(el_not+g+g-m)*nanbnc]*double(1-2*((g+g+m)%2))+2*alp[el_not+g+(el_not-m)*na+(el_not+g)*nanb+(el_not+m-g-g)*nanbnc]*double(1-2*(g%2)));
        hess[g+k*hs]+=ta*(re[m]*re[m-k-g]+im[m]*im[m-k-g]);
        hess[k+g*hs]+=ta*(re[m]*re[m-k-g]+im[m]*im[m-k-g]);
        hess[g+el_not+(k+el_not)*hs]-=ta*(re[m]*re[m-k-g]+im[m]*im[m-k-g]);
        hess[k+el_not+(g+el_not)*hs]-=ta*(re[m]*re[m-k-g]+im[m]*im[m-k-g]);
        }
        // Term 7 of 7
        if(g-k-m>=1) {
        ta=(lambda4/6.)*alp[el_not+k+(el_not+m)*na+(el_not+g-k-m)*nanb+(el_not-g)*nanbnc]*double(1-2*((g+k+m)%2));
        hess[g+k*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        hess[k+g*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        hess[g+el_not+(k+el_not)*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        hess[k+el_not+(g+el_not)*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);

        if(2*m<=g-k) {
        ta=(lambda4/6.)*2*(alp[el_not+m+(el_not+g-k-m)*na+(el_not+g)*nanb+(el_not-k)*nanbnc]*double(1-2*(k%2))+alp[el_not+m+(el_not-g)*na+(el_not+g-k-m)*nanb+(el_not+k)*nanbnc]*double(1-2*(m%2)));
        if(2*m==g-k) ta=ta/2.;
        hess[g+k*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        hess[k+g*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        hess[g+el_not+(k+el_not)*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        hess[k+el_not+(g+el_not)*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        }

        if(m>k) ta=(lambda4/6.)*(alp[el_not+k+(el_not+m)*na+(el_not+g)*nanb+(el_not+k+m-g)*nanbnc]*double(1-2*((m+g+k)%2))+2*alp[el_not+k+(el_not-g)*na+(el_not+m)*nanb+(el_not+g-k-m)*nanbnc]*double(1-2*(k%2)));
        if(m<k) ta=(lambda4/6.)*(alp[el_not+m+(el_not+k)*na+(el_not+g)*nanb+(el_not+m+k-g)*nanbnc]*double(1-2*((m+g+k)%2))+2*alp[el_not+m+(el_not-g)*na+(el_not+k)*nanb+(el_not+g-k-m)*nanbnc]*double(1-2*(m%2)));
       if(m==k) ta=(lambda4/6.)*(alp[el_not+k+(el_not+k)*na+(el_not+g)*nanb+(el_not+k+k-g)*nanbnc]*double(1-2*((k+g+k)%2))+2*alp[el_not+k+(el_not-g)*na+(el_not+k)*nanb+(el_not+g-k-k)*nanbnc]*double(1-2*(k%2)));
        hess[g+k*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        hess[k+g*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        hess[g+el_not+(k+el_not)*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);
        hess[k+el_not+(g+el_not)*hs]+=ta*(re[m]*re[g-k-m]-im[m]*im[g-k-m]);

        }

        } //end of loop over m

        }
        } // end of RkRg IkIg loops



        return gg;
}


//FUNCTION EVALUATION
DP funk(Vec_I_DP &x, double* Gaunt_matrix, double* Gaunt_4th)
{
    double H_gradient = 0;
    double H_thirdorder=0;
    double H_fourthorder = 0;
    double H=0;

	double addedterm;
    
    complex<double> c[(int)(2*(el_not+1)+1)];
    c[0]=complex<double> (x[0],0);
    for(int i=1;i<el_not+1;i++)
        c[i]=complex<double> (x[i],x[el_not+i]);

   //Kinetic term in Hamiltonian 
   H_gradient += (1/2.0)*tau*(norm(c[0]));
   for(int m=1; m<=el_not; m++)
   {
       H_gradient += tau*(norm(c[m]));
   }
    
    
    ///3rd order term in Hamiltonian
    if((3*el_not)%2==0)
    {
        H_thirdorder +=(lambda3/6.0)*real(Gaunt_matrix[el_not + Na*el_not + Na*Nb*el_not])*pow(real(c[0]),3);
        for(int m2=1;m2<=el_not;m2++)
        {
     	    H_thirdorder+=4.0*(lambda3/6.0)*real(Gaunt_matrix[(m2+el_not) + Na*(-m2+el_not) + NaNb*el_not])*double(1-2*(m2%2))*real(c[0])*norm(c[m2]);
            for(int m1=m2;m1<=el_not;m1++)
            {
                if(abs(m1+m2)<=el_not) {
                   if(m1==m2) { H_thirdorder+=(lambda3/6.0)*real(Gaunt_matrix[(m1+el_not) + Na*(m2+el_not) + NaNb*(-m1-m2+el_not)])*double(1-2*((m1+m2)%2))*2*(real(c[m1]*c[m2]*conj(c[m1+m2]))); }
			else  { H_thirdorder+=2*(lambda3/6.0)*real(Gaunt_matrix[(m1+el_not) + Na*(m2+el_not) + NaNb*(-m1-m2+el_not)])*double(1-2*((m1+m2)%2))*2*(real(c[m1]*c[m2]*conj(c[m1+m2]))); }
     
			}
	     } 
            for(int m1=m2;m1<=m2;m1++)
            {
                    H_thirdorder+=(lambda3/6.0)*real(Gaunt_matrix[(m1+el_not) + Na*(-m2+el_not) + NaNb*(-m1+m2+el_not)])*double(1-2*(m2%2))*2*(real(c[m1]*conj(c[m2])*c[m2-m1]));
            }
            for(int m1=m2+1;m1<=el_not;m1++)
            {
                if(abs(m2-m1)<=el_not)
                    H_thirdorder+=2*(lambda3/6.0)*real(Gaunt_matrix[(m1+el_not) + Na*(-m2+el_not) + NaNb*(-m1+m2+el_not)])*double(1-2*(m1%2))*2*(real(c[m1]*conj(c[m2])*conj(c[m1-m2])));
            }
        }
     }
           

  //4th order term in Hamiltonian

		H_fourthorder=0;
		H_fourthorder+=(lambda4/24.0)*(Gaunt_4th[el_not + na*el_not + nanb*el_not + nanbnc*el_not]*pow(real(c[0]),4));

		for(int m3=1;m3<=el_not;m3++)
			H_fourthorder+=(lambda4/12.0)*(Gaunt_4th[el_not + na*el_not + nanb*(m3+el_not) + nanbnc*(el_not-m3)])*double(1-2*(m3%2))*norm(c[0])*norm(c[m3]);

		for(int m2=1;m2<=el_not;m2++)
		{

			H_fourthorder+=(lambda4/6.0)*(Gaunt_4th[el_not + na*(m2+el_not) + nanb*(el_not) + nanbnc*(el_not-m2)])*norm(c[0])*norm(c[m2]);

			for(int m3=1;m3<=m2;m3++)
			{
			   if(abs(m2+m3)<=el_not)
				H_fourthorder+=(lambda4/6.0)*(Gaunt_4th[el_not + na*(m2+el_not) + nanb*(m3+el_not) + nanbnc*(el_not-m2-m3)])*double(1-2*(m3%2))*(real(c[0]*c[m2]*c[m3]*conj(c[m2+m3])));
			   if(abs(m3-m2)<=el_not)
				H_fourthorder+=(lambda4/6.0)*(Gaunt_4th[el_not + na*(m2+el_not) + nanb*(el_not-m3) + nanbnc*(el_not+m3-m2)])*(real(c[0]*conj(c[m2])*c[m3]*c[m2-m3]));
			}

			for(int m3=m2+1;m3<=el_not;m3++)
			{
			   if(abs(m2+m3)<=el_not)
				H_fourthorder+=(lambda4/6.0)*(Gaunt_4th[el_not + na*(m2+el_not) + nanb*(m3+el_not) + nanbnc*(el_not-m2-m3)])*double(1-2*(m3%2))*(real(c[0]*c[m2]*c[m3]*conj(c[m2+m3])));
			   if(abs(m3-m2)<=el_not)
				H_fourthorder+=(lambda4/6.0)*(Gaunt_4th[el_not + na*(m2+el_not) + nanb*(el_not-m3) + nanbnc*(el_not+m3-m2)])*double(1-2*((m3+m2)%2))*(real(c[0]*conj(c[m2])*c[m3]*conj(c[m3-m2])));

			}

		}


		for(int m1=1;m1<=el_not;m1++)
		{

			for(int m2=m1;m2<=el_not;m2++)
			{

			addedterm=0;
		           if(m1+m2<=el_not)
				addedterm+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not+m2) + nanb*(el_not) + nanbnc*(el_not-m1-m2)])*(real(c[0]*c[m1]*c[m2]*conj(c[m1+m2])));
				for(int m3=1;m3<=m1+m2;m3++)
				{
				    if(abs(m1+m2+m3)<=el_not)
					addedterm+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not+m2) + nanb*(el_not+m3) + nanbnc*(el_not-m1-m2-m3)])*double(1-2*(m3%2))*(real(c[m1]*c[m2]*c[m3]*conj(c[m1+m2+m3])));
				    if(abs(-m1-m2+m3)<=el_not)
					addedterm+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not+m2) + nanb*(el_not-m3) + nanbnc*(el_not-m1-m2+m3)])*(real(conj(c[m1])*conj(c[m2])*c[m3]*c[m1+m2-m3]));
				}
				for(int m3=m1+m2+1;m3<=el_not;m3++)
				{
				    if(abs(m1+m2+m3)<=el_not)
					addedterm+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not+m2) + nanb*(el_not+m3) + nanbnc*(el_not-m1-m2-m3)])*double(1-2*(m3%2))*(real(c[m1]*c[m2]*c[m3]*conj(c[m1+m2+m3])));
				    if(abs(-m1-m2+m3)<=el_not)
					addedterm+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not+m2) + nanb*(el_not-m3) + nanbnc*(el_not-m1-m2+m3)])*double(1-2*((m1+m2+m3)%2))*(real(conj(c[m1])*conj(c[m2])*c[m3]*conj(c[-m1-m2+m3])));
				}
			
			if(m1==m2) {H_fourthorder+=addedterm; } else { H_fourthorder+=2*addedterm; }

			}
                  }


		for(int m1=1;m1<=el_not;m1++)
		{
                        for(int m2=m1;m2<=m1;m2++)
                        {
                                if(abs(-m1+m2)<=el_not)
                                        H_fourthorder+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not-m2)+ nanb*(el_not) + nanbnc*(el_not+m2-m1)])*double(1-2*(m2%2))*(real(c[0]*c[m1]*conj(c[m2])*conj(c[m1-m2])));
                        }
                        for(int m2=m1+1;m2<=el_not;m2++)
                        {
                                if(abs(-m1+m2)<=el_not)
                                        H_fourthorder+=2*(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not-m2) + nanb*(el_not) + nanbnc*(el_not+m2-m1)])*double(1-2*(m1%2))*(real(c[0]*c[m1]*conj(c[m2])*c[-m1+m2]));
                        }



			for(int m3=m1;m3<=m1;m3++)
			{
				for(int m2=1;m2<=min(m1+m3,el_not);m2++)
				{
					if(abs(-m1+m2-m3)<=el_not)
						H_fourthorder+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not-m2) + nanb*(el_not+m3) + nanbnc*(el_not+m2-m1-m3)])*double(1-2*((m2+m3)%2))*(real(c[m1]*conj(c[m2])*c[m3]*conj(c[m1+m3-m2])));
				}
				for(int m2=m1+m3+1;m2<=el_not;m2++)
				{
					if(abs(-m1+m2-m3)<=el_not)
						H_fourthorder+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not-m2) + nanb*(el_not+m3) + nanbnc*(el_not+m2-m1-m3)])*double(1-2*(m1%2))*(real(c[m1]*conj(c[m2])*c[m3]*c[-m1-m3+m2]));
				}
			}


			for(int m3=m1+1;m3<=el_not;m3++)
			{
				for(int m2=1;m2<=min(m1+m3,el_not);m2++)
				{
					if(abs(-m1+m2-m3)<=el_not)
						H_fourthorder+=2*(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not-m2) + nanb*(el_not+m3) + nanbnc*(el_not+m2-m1-m3)])*double(1-2*((m2+m3)%2))*(real(c[m1]*conj(c[m2])*c[m3]*conj(c[m1+m3-m2])));
				}
				for(int m2=m1+m3+1;m2<=el_not;m2++)
				{
					if(abs(-m1+m2-m3)<=el_not)
						H_fourthorder+=2*(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not-m2) + nanb*(el_not+m3) + nanbnc*(el_not+m2-m1-m3)])*double(1-2*(m1%2))*(real(c[m1]*conj(c[m2])*c[m3]*c[-m1-m3+m2]));
				}
			}

		}

		for(int m2=1;m2<=el_not;m2++)
		{
			for(int m3=m2;m3<=el_not;m3++)
			{

			addedterm=0;

				for(int m1=1;m1<=min(m2+m3,el_not);m1++)
				{
					if(abs(-m1+m2+m3)<=el_not)
						addedterm+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not-m2) + nanb*(el_not-m3) + nanbnc*(el_not+m3+m2-m1)])*double(1-2*((m1+m3)%2))*(real(c[m1]*conj(c[m2])*conj(c[m3])*c[-m1+m2+m3]));
				}
				for(int m1=m2+m3+1;m1<=el_not;m1++)
				{
					if(abs(-m1+m2+m3)<=el_not)
						addedterm+=(lambda4/12.0)*(Gaunt_4th[el_not+m1 + na*(el_not-m2) + nanb*(el_not-m3) + nanbnc*(el_not+m3+m2-m1)])*double(1-2*(m2%2))*(real(c[m1]*conj(c[m2])*conj(c[m3])*conj(c[m1-m2-m3])));
				}

			if(m2==m3) { H_fourthorder+=addedterm; } else { H_fourthorder+=2*addedterm; }

			}
		}

                
            

    H = H_gradient+H_thirdorder+H_fourthorder;
    return H;    
}



int counter=1;

DP tt;
DP NR::amotsa(Mat_IO_DP &p, Vec_O_DP &y, Vec_IO_DP &psum, Vec_O_DP &pb, DP &yb, DP funk(Vec_I_DP &, double* Gaunt_matrix, double* Gaunt_4th), const int ihi, DP &yhi, const DP fac)
{
    ofstream yfluc ("yfluc.txt", ios_base::app);
    ofstream y_try ("ytry.txt", ios_base::app);
    mt19937 gen(rd());

    int j;
    DP fac1,fac2,yflu,ytry;
    
    int ndim=p.ncols();
    Vec_DP ptry(ndim);
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for(j=0;j<ndim;j++)
        ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=funk(ptry,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order);

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

void NR::amebsa(Mat_IO_DP &p, Vec_IO_DP &y, Vec_O_DP &pb, DP &yb, const DP ftol, DP funk(Vec_I_DP &, double* Gaunt_matrix, double* Gaunt_4th), int &iter, const DP temptr)
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
                        y[i]=funk(psum,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order);
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
    ofstream T_dH,T_H,endcms,cms,parameters,time; 
    T_H.open("T_H.txt");
    T_dH.open("T_dH.txt");
    endcms.open("endcms.txt");
    cms.open("cms.txt");
    parameters.open("parameters.txt");
    time.open("time.txt");
    random_device rd;
    mt19937 gen(200);
    wig_table_init(2*100,3); //define using largest 2*3j value and largest wigner type (i.e. 3j, 6j, 9j) - however this particular 2*100, 9 is commonly used becuase it doesn't take a lot of memory and yet it works for most programs
    wig_temp_init(2*100);
    double t;

    double hess_array[(int)(2*el_not+1)*(2*el_not+1)];
    double hess_array2[(int)(2*el_not+1)*(2*el_not+1)];
   
    double g_array[(int)(2*el_not+1)];	
    for(int g=0;g<2*el_not+1;g++) g_array[g]=0;   
 
    parameters<<"el_not: "<<el_not<<endl;
    parameters<<"tau: "<<tau<<endl;
    parameters<<"lambda3: "<<lambda3<<endl;
    parameters<<"lambda4: "<<lambda4<<endl;
    parameters<<endl;

    //Gaunt matrix initialization for Hamiltonian 3rd order term
    if((3*el_not)%2==0)
    {
	Gaunt_matrix_3rd_order[el_not + Na*el_not + Na*Nb*el_not] = Gaunt(el_not,el_not,el_not,0,0,0);
        for(int m2=1;m2<=el_not;m2++)
        {
	    Gaunt_matrix_3rd_order[(m2+el_not) + Na*(-m2+el_not) + Na*Nb*el_not] = Gaunt(el_not,el_not,el_not,m2,-m2,0);
            for(int m1=1;m1<=el_not;m1++)
            {
                if(abs(m1+m2)<=el_not)
                    Gaunt_matrix_3rd_order[(m1+el_not) + Na*(m2+el_not) + Na*Nb*(-m1-m2+el_not)] = Gaunt(el_not,el_not,el_not,m1,m2,-m1-m2);
            }
            for(int m1=1;m1<=m2;m1++)
            {
                if(abs(m2-m1)<=el_not)
                    Gaunt_matrix_3rd_order[(m1+el_not) + Na*(-m2+el_not) + Na*Nb*(-m1+m2+el_not)] = Gaunt(el_not,el_not,el_not,m1,-m2,-m1+m2);
            }
            for(int m1=m2+1;m1<=el_not;m1++)
            {
                if(abs(m2-m1)<=el_not)
                    Gaunt_matrix_3rd_order[(m1+el_not) + Na*(-m2+el_not) + Na*Nb*(-m1+m2+el_not)] = Gaunt(el_not,el_not,el_not,m1,-m2,-m1+m2);
            }
        }
     }

   //Gaunt matrix initialization for Hamiltonian 4th order term
   double temp;
   int tm1, tm2, tm3, tm4;
   for(int i=0;i<na*nb*nc*nd; i++) {

	tm1=i%na;
	tm2=(i/na)%nb;
	tm3=(i/(na*nb))%nc;
	tm4=(i/(na*nb*nc))%nd;

	temp=0;
	for (int ll=0; ll<=2*el_not+2;ll++){ 
			temp+=Gaunt(el_not,el_not,ll,tm1-el_not,tm2-el_not,2*el_not-tm1-tm2)*Gaunt(el_not,el_not,ll,tm3-el_not,tm4-el_not,2*el_not-tm3-tm4);
					}

	Gaunt_matrix_4th_order[i]=temp;

	}


			
    //Define staring point values
    DP p_array[2*el_not+2][2*el_not+1];  

   
    parameters<<"Starting p matrix"<<endl;
//starting random cm's
    for(int a=0;a<2*el_not+2;a++)
    {
	for(int b=0;b<2*el_not+1;b++)
	{
	    p_array[a][b] = sqrt((-1)*tau/lambda4)*(double(gen())/double(gen.max()));
	    parameters<<p_array[a][b]<<" ";
	}
	parameters<<endl;
    }

    Mat_IO_DP p(*p_array,2*el_not+2,2*el_not+1);
    
    
//Set corresponding initial y-values
    double magg;
    DP psi_array[(int)(2*el_not+1)];
    DP y_array[(int)(2*el_not+2)];
    parameters<<"Starting y value: "<<endl;
    for(int xi=0; xi<2*el_not+2; xi++)
    {
        for(int xii=0; xii<2*el_not+1;xii++)
        {
            psi_array[xii]=p_array[xi][xii];
	}
        Vec_I_DP psi(psi_array,2*el_not+1);
        y_array[xi]=funk(psi,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order);
	parameters<<y_array[xi] << ", ";
    }
    parameters<<endl;


        for(int xii=0; xii<2*el_not+1;xii++)
        {
            psi_array[xii]=p_array[0][xii];
	}
        Vec_I_DP psi(psi_array,2*el_not+1);

    magg=grad(psi,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order);
    parameters<<"dH = " <<magg<< endl;

    Vec_IO_DP y(y_array,2*el_not+2);
   

 
    double pb_array[(int)(2*el_not+1)];
    Vec_O_DP pb(pb_array,2*el_not+1);
    const DP ftol = 1.0e-8;
    parameters<<"FTolerance: "<<ftol<<endl;
    
    DP yb = 10;
    parameters<<"yb starting: "<<yb<<endl; 
    int iter= 1000;
    parameters<<"Number of iterations: "<<iter<<endl;
    const DP temptr = 1;
    parameters<<"Starting temperature: "<<temptr<<endl;

   NR::amebsa(p,y,pb,yb,ftol,funk,iter,temptr);
    
   
    double tinterval = 0.1;
    parameters<<"Temperature interval: "<<tinterval<<endl;


//////ANNEALING SCHEDULE/////////////
	    for(t=temptr*0.5;t>temptr*0.1;t-=tinterval*0.05)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_dH<<t<<" "<<grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
		T_H<<t<<" "<<yb<< " " << grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
	       int zero = 0;
	       cms<<pb[0]<<" "<<zero<<endl;
	       for(int i=1;i<el_not+1;i++)
	       {
		   cms<<pb[i]<<" "<<pb[el_not+i]<<endl;
	       }
	       cms<<endl;
	    }

	    for(t=temptr*0.1;t>=temptr*0.01;t-=tinterval*0.01)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_dH<<t<<" "<<grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
		T_H<<t<<" "<<yb<< " " << grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;	
	       int zero = 0;
	       cms<<pb[0]<<" "<<zero<<endl;
	       for(int i=1;i<el_not+1;i++)
	       {
		   cms<<pb[i]<<" "<<pb[el_not+i]<<endl;
	       }
	       cms<<endl;
	    }

	    for(t=temptr*0.01;t>=0.001;t-=tinterval*0.001)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_dH<<t<<" "<<grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
		T_H<<t<<" "<<yb<< " " << grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
	       int zero = 0;
	       cms<<pb[0]<<" "<<zero<<endl;
	       for(int i=1;i<el_not+1;i++)
	       {
		   cms<<pb[i]<<" "<<pb[el_not+i]<<endl;
	       }
	       cms<<endl;
	    }

	    for(t=temptr*0.001;t>=0.0001;t-=tinterval*0.0001)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_dH<<t<<" "<<grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
		T_H<<t<<" "<<yb<< " " << grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
	       int zero = 0;
	       cms<<pb[0]<<" "<<zero<<endl;
	       for(int i=1;i<el_not+1;i++)
	       {
		   cms<<pb[i]<<" "<<pb[el_not+i]<<endl;
	       }
	       cms<<endl;
	    }

	    for(t=temptr*0.0001;t>=0.00001;t-=tinterval*0.00001)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_dH<<t<<" "<<grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
		T_H<<t<<" "<<yb<< " " << grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
	       int zero = 0;
	       cms<<pb[0]<<" "<<zero<<endl;
	       for(int i=1;i<el_not+1;i++)
	       {
		   cms<<pb[i]<<" "<<pb[el_not+i]<<endl;
	       }
	       cms<<endl;
	    }
	    
	    for(t=temptr*0.00001;t>=0.000001;t-=tinterval*0.000001)
	    {
		iter=10000;
		NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
		T_dH<<t<<" "<<grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
		T_H<<t<<" "<<yb<< " " << grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
	       int zero = 0;
	       cms<<pb[0]<<" "<<zero<<endl;
	       for(int i=1;i<el_not+1;i++)
	       {
		   cms<<pb[i]<<" "<<pb[el_not+i]<<endl;
	       }
	       cms<<endl;
	    }

	for(t=temptr*0.000001;t>=0.0;t-=tinterval*0.0000001)
    	{	
        	iter=10000;
        	NR::amebsa(p,y,pb,yb,ftol,funk,iter,t);
        	T_dH<<t<<" "<<grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
        	T_H<<t<<" "<<yb<< " " << grad(pb,g_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order) <<endl;
       		int zero = 0;
       		cms<<pb[0]<<" "<<zero<<endl;
       		for(int i=1;i<el_not+1;i++)
       		{
       		    cms<<pb[i]<<" "<<pb[el_not+i]<<endl;
       		}
       		cms<<endl;
	}

   int zero = 0;
   endcms<<pb[0]<<" "<<zero<<endl;
   for(int i=1;i<el_not+1;i++)
    {
	endcms<<pb[i]<<" "<<pb[el_not+i]<<endl;
    }

    magg=hessian(pb,hess_array,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order);
    magg=nhessian(pb,hess_array2,funk);

    Eigen::MatrixXd hhh(hs,hs);
    Eigen::MatrixXd hhh2(hs,hs);

    for (int m=0; m<hs*hs; m++) {
            int x = m%(2*el_not+1);
            int y = (m/hs)%hs;
    hhh(x,y)=hess_array[m];
    hhh2(x,y)=hess_array2[m];
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es(hhh);
    Eigen::EigenSolver<Eigen::MatrixXd> es2(hhh2);

    parameters << "Eigenvalues: " << endl << es.eigenvalues() << endl;
    parameters << "Eigenvalues 2: " << endl << es2.eigenvalues() << endl;
//      cout << "Eigenvectors: " << endl << es.eigenvectors() << endl;
   parameters<<"Lowest point derivative value: "<<grad(pb,g_array,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order)<<endl;
   parameters<<"Lowest_yb: "<<yb<<endl;

   // delete [] c;//is this where you delete it really?
   T_H.close();
   T_dH.close();
   endcms.close();
   cms.close();
   parameters.close();

   //free Gaunt matrices
   delete[] Gaunt_matrix_3rd_order;
   delete[] Gaunt_matrix_4th_order;
   delete[] Hamiltonian_values;
   delete[] Gradient_values;
   t = clock() - t;
   time<<"Computation Time: "<<t/CLOCKS_PER_SEC/60.0<<" minutes";
   time.close();

   return 0;    
}

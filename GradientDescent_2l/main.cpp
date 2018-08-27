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
#include "nrutil.h"
#include "nrtypes.h"
#include "nr.h"
#include "wigxjpf.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
const double pi = 3.1415926535897;
random_device rd;

const double kappa = 1.0;
const double tau = -1.0;
const double lambda4 = 1.0;

/////Set lambda3 to be any value between -1.0 and 1.0//////////
double lambda3 = 0.0;

////Set l_not to be any value between el_not and el_not +1//////////
double l_not = (2.0*el_not+1.0)/2.0;

int ncom;
DP (*nrfunc)(Vec_I_DP &, double*, double*, double*);
DP (*nrdfun)(DP(Vec_I_DP &, double*, double*, double*), const Vec_I_DP, Vec_O_DP &, const DP, DP &);
Vec_DP *pcom_p, *xicom_p;

namespace NR
{
	DP dfridr(DP func(const Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err);
	void mnbrak(DP &ax, DP &bx, DP &cx, DP &fa, DP &fb, DP &fc, DP f(const DP, double*, double*, double*));
	DP dbrent(const DP ax, const DP bx, const DP cx, DP f(const DP, double*, double*, double*), DP df(DP func(const DP, double*, double*, double*), const DP, const DP h, DP &err), const DP tol, DP &xmin, const DP h, DP &err);
	DP f1dim(const DP x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2);
        DP df1dim(DP func(const DP, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const DP x, const DP h, DP &err); 
	void dflinmin(Vec_IO_DP &p, Vec_IO_DP &xi, DP &fret, DP func(Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), DP dfridr(DP func(const Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err),const DP h, DP &err);
        void frprmn(Vec_IO_DP &p, const DP ftol, int &iter, DP &fret, DP func(Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), DP dfridr(DP func(const Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err),const DP h, DP &err);
}
namespace
{
	inline void shft3(DP &a, DP &b, DP &c, const DP d)
	{
		a = b;
		b = c;
		c = d;
	}
	inline void mov3(DP &a, DP &b, DP &c, const DP d, const DP e, const DP f)
	{
		a = d;
		b = e;
		c = f;
	}
}

double* Make6DArray(int arraySizeA, int arraySizeB, int arraySizeC, int arraySizeD, int arraySizeE, int arraySizeF)
{

    int size;
    size=arraySizeA*arraySizeB*arraySizeC*arraySizeD*arraySizeE*arraySizeF;
    double *Array = new double[size];
    return Array;
}

//Define size of Gaunt matrices
int Na = 2*el_not+3;
int Nb = 3*el_not+5;
int Nc = 3*el_not+5;
int Nd = 6*el_not+9;

int hs = (4*el_not+4);
int NaNb = (2*el_not+3)*(3*el_not+5);
int NaNbNc = (2*el_not+3)*(3*el_not+5)*(3*el_not+5);


double* Gaunt_matrix_3rd_order = Make6DArray(2,2,2,Na,Na,Na);
double* Gaunt_matrix_4th_order_1 = Make6DArray(2,2,Na,Nb,Nc,Nd);
double* Gaunt_matrix_4th_order_2 = Make6DArray(2,2,Na,Nb,Nc,Nd);

//Intialize tables for wigner3j coefficients
void wig_table_init(int max_two_j, int wigner_type);
void wig_temp_init(int max_two_j);        /* Single-threaded */
void wig_thread_temp_init(int max_two_j); /* Multi-threaded. */
double wig3jj(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3);
void wig_temp_free();  /* Per-thread when multi-threaded. */
void wig_table_free();

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

//FUNCTION EVALUTAION
DP func(Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2)
{
    double EPSILON = 1e-5;
    double H_gradient =0;
    double H_thirdorder = 0;
    double H_fourthorder = 0;
    double H = 0;

   
    complex<double>** c = Make2DArray(2, el_not+2);

    for(int q=0;q<=1;q++)
        c[q][0] = complex<double> (x[q],0);

    int count=2;

    for(int i=0;i<=1;i++)
    {
        for(int j=1;j<=el_not+i;j++)
        {
           c[i][j]=complex<double> (x[count],x[(2*el_not+1+count)]);
           count++;
        }
    }
    c[0][(el_not+1)]=complex<double> (0,0);

    //Gradient term
    for(int l=el_not; l<=el_not+1; l++)
    {
        H_gradient += (1/2.0)*(kappa*(pow(l-l_not,2))+tau)*real(c[(l-el_not)][0]*c[(l-el_not)][0]);
        for(int m=1; m<=l; m++)
        {
              H_gradient += ((kappa*(pow(l-l_not,2))+tau)*(norm(c[(l-el_not)][m])));
        }
    }


 ///3rd order term in Hamiltonian
    int FirstIndex,SecondIndex,ThirdIndex;

    for(int el_1=el_not;el_1<=el_not+1;el_1++)
    {
        FirstIndex = el_1-el_not;
        for(int el_2=el_not;el_2<=el_not+1;el_2++)
        {
            SecondIndex = el_2-el_not;
            for(int el_3=el_not;el_3<=min(el_not+1,el_1+el_2);el_3++)
            {
                ThirdIndex = el_3-el_not;
                if((el_3>=abs(el_1-el_2))&&((el_1+el_2+el_3)%2==0))
                {
                    H_thirdorder += (lambda3/6.0)*Gaunt_matrix[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*el_1 + 8*Na*el_2 + 8*Na*Na*el_3]*real(c[ThirdIndex][0]*c[SecondIndex][0]*c[FirstIndex][0]);
                    for(int m2=1;m2<=el_2;m2++)
                    {
                            H_thirdorder += 4.0*(lambda3/6.0)*Gaunt_matrix[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m2+el_1) + 8*Na*(el_2-m2) + 8*Na*Na*el_3]*double(1-2*(m2%2))*real(c[ThirdIndex][0]*c[FirstIndex][m2]*conj(c[SecondIndex][m2]));
                            for(int m1=1;m1<=el_1;m1++)
                            {
                                    if(abs(m1+m2)<=el_3)
                                            H_thirdorder += (lambda3/6.0)*Gaunt_matrix[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m1+el_1) + 8*Na*(m2+el_2) + 8*Na*Na*(-m1-m2+el_3)]*double(1-2*((m1+m2)%2))*2*(real(c[FirstIndex][m1]*c[SecondIndex][m2]*conj(c[ThirdIndex][m1+m2])));
                            }
                            for(int m1=1;m1<=m2;m1++)
                            {
                                    if((m2-m1)<=el_3)
                                            H_thirdorder += (lambda3/6.0)*Gaunt_matrix[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m1+el_1) + 8*Na*(-m2+el_2) + 8*Na*Na*(m2-m1+el_3)]*double(1-2*(m2%2))*2*(real(c[FirstIndex][m1]*conj(c[SecondIndex][m2])*c[ThirdIndex][m2-m1]));
                            }
                            for(int m1=m2+1;m1<=el_1;m1++)
                            {
                                    if((m2-m1)<=el_3)
                                            H_thirdorder += (lambda3/6.0)*Gaunt_matrix[FirstIndex + 2*SecondIndex + 4*ThirdIndex + 8*(m1+el_1) + 8*Na*(-m2+el_2) + 8*Na*Na*(m2-m1+el_3)]*double(1-2*(m1%2))*2*(real(c[FirstIndex][m1]*conj(c[SecondIndex][m2])*conj(c[ThirdIndex][m1-m2])));
                            }
                    }
                }
            }
        }
    }


//4th order term in Hamiltonian
  double FirstGauntTerm, SecondGauntTerm, ThirdGauntTerm, FourthGauntTerm, ZerothGauntTerm;
  int FirstIndex_1, SecondIndex_1, FourthIndex_1, FifthIndex_1, SixthIndex, FirstIndex_2, SecondIndex_2, FourthIndex_2, FifthIndex_2;

  for(int el_1=el_not;el_1<=el_not+1;el_1++)
  {
        FirstIndex_1 = el_1-el_not;
        FourthIndex_1 = el_1+el_not+2;
        for(int el_2=el_not;el_2<=el_not+1;el_2++)
        {
            SecondIndex_1 = el_2-el_not;
            FifthIndex_1 = el_2+el_not+2;
            for(int el_3=el_not;el_3<=el_not+1;el_3++)
            {
                FirstIndex_2= el_3-el_not;
                FourthIndex_2 = el_3+el_not+2;
                for(int el_4=el_not;el_4<=el_not+1;el_4++)
                {
                    SecondIndex_2 = el_4-el_not;
                    FifthIndex_2 = el_4+el_not+2;
                    for(int l=max(abs(el_1-el_2),abs(el_3-el_4));l<=min(el_1+el_2,el_3+el_4);l++)
                    {
                        if((el_1+el_2+l)%2==0 && (el_3+el_4+l)%2==0)
                        {
                                SixthIndex = l+2*el_not+2;
                                ZerothGauntTerm = Gaunt_matrix_4th_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*FourthIndex_1 + 4*NaNb*FifthIndex_1 + 4*NaNbNc*SixthIndex];

                                H_fourthorder += (lambda4/24.0)*ZerothGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*NaNb*FifthIndex_2 + 4*NaNbNc*SixthIndex]*real(c[FirstIndex_1][0]*c[SecondIndex_1][0]*c[FirstIndex_2][0]*c[SecondIndex_2][0]);

                                for(int m3=1;m3<=el_3;m3++)
                                        H_fourthorder += (lambda4/12.0)*ZerothGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*NaNb*(-m3+FourthIndex_2) + 4*NaNbNc*SixthIndex]*double(1-2*(m3%2))*real(c[FirstIndex_1][0]*c[SecondIndex_1][0]*c[FirstIndex_2][m3]*conj(c[SecondIndex_2][m3]));

                                for(int m2=1;m2<=el_2;m2++)
                                {
                                        FirstGauntTerm = Gaunt_matrix_4th_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*FourthIndex_1 + 4*NaNb*(m2+FifthIndex_1) + 4*NaNbNc*(-m2+SixthIndex)];

                                        H_fourthorder += (lambda4/6.0)*FirstGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*NaNb*(-m2+FifthIndex_2) + 4*NaNbNc*(m2+SixthIndex)]*real(c[FirstIndex_1][0]*c[FirstIndex_2][0]*c[SecondIndex_1][m2]*conj(c[SecondIndex_2][m2]));

                                        for(int m3=1;m3<=min(m2,el_3);m3++)
                                        {
                                           if(abs(m2+m3)<=el_4)
                                                H_fourthorder += (lambda4/12.0)*FirstGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*NaNb*(-m2-m3+FifthIndex_2)+4*NaNbNc*(m2+SixthIndex)]*double(1-2*(m3%2))*2*(real(c[FirstIndex_1][0]*c[SecondIndex_1][m2]*c[FirstIndex_2][m3]*conj(c[SecondIndex_2][m2+m3])));
 if(abs(m3-m2)<=el_4)
                                                H_fourthorder += (lambda4/12.0)*FirstGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*NaNb*(m3-m2+FifthIndex_2) + 4*NaNbNc*(m2+SixthIndex)]*2*(real(c[FirstIndex_1][0]*conj(c[SecondIndex_1][m2])*c[FirstIndex_2][m3]*c[SecondIndex_2][m2-m3]));
                                        }

                                        for(int m3=m2+1;m3<=el_3;m3++)
                                        {
                                           if(abs(m2+m3)<=el_4)
                                                H_fourthorder += (lambda4/12.0)*FirstGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*NaNb*(-m2-m3+FifthIndex_2) + 4*NaNbNc*(m2+SixthIndex)]*double(1-2*(m3%2))*2*(real(c[FirstIndex_1][0]*c[SecondIndex_1][m2]*c[FirstIndex_2][m3]*conj(c[SecondIndex_2][m2+m3])));

                                           if(abs(m3-m2)<=el_4)
                                                H_fourthorder += (lambda4/12.0)*FirstGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*NaNb*(m3-m2+FifthIndex_2) + 4*NaNbNc*(m2+SixthIndex)]*double(1-2*((m2+m3)%2))*2*(real(c[FirstIndex_1][0]*conj(c[SecondIndex_1][m2])*c[FirstIndex_2][m3]*conj(c[SecondIndex_2][m3-m2])));
                                        }
                                }

                                for(int m1=1;m1<=el_1;m1++)
                                {
                                        for(int m2=1;m2<=el_2;m2++)
                                        {
                                            if(m1+m2<=l)
                                            {
                                            ThirdGauntTerm = Gaunt_matrix_4th_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*NaNb*(m2+FifthIndex_1) + 4*NaNbNc*(-m1-m2+SixthIndex)];

                                                H_fourthorder += (lambda4/24.0)*ThirdGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*NaNb*(-m1-m2+FifthIndex_2) + 4*NaNbNc*(m1+m2+SixthIndex)]*2*(real(c[FirstIndex_2][0]*c[FirstIndex_1][m1]*c[SecondIndex_1][m2]*conj(c[SecondIndex_2][m1+m2])));
                                                for(int m3=1;m3<=min(m1+m2,el_3);m3++)
                                                {
                                                    if(abs(m1+m2+m3)<=el_4)
                                                        H_fourthorder += (lambda4/24.0)*ThirdGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*NaNb*(-m1-m2-m3+FifthIndex_2) + 4*NaNbNc*(m1+m2+SixthIndex)]*double(1-2*(m3%2))*2*(real(c[FirstIndex_1][m1]*c[SecondIndex_1][m2]*c[FirstIndex_2][m3]*conj(c[SecondIndex_2][m1+m2+m3])));

                                                    if(abs(-m1-m2+m3)<=el_4)
                                                        H_fourthorder += (lambda4/24.0)*ThirdGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*NaNb*(-m1-m2+m3+FifthIndex_2) + 4*NaNbNc*(m1+m2+SixthIndex)]*2*(real(conj(c[FirstIndex_1][m1])*conj(c[SecondIndex_1][m2])*c[FirstIndex_2][m3]*c[SecondIndex_2][m1+m2-m3]));
                                                }
for(int m3=m1+m2+1;m3<=el_3;m3++)
                                                {
                                                    if(abs(m1+m2+m3)<=el_4)
                                                        H_fourthorder += (lambda4/24.0)*ThirdGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*NaNb*(-m1-m2-m3+FifthIndex_2) + 4*NaNbNc*(m1+m2+SixthIndex)]*double(1-2*(m3%2))*2*(real(c[FirstIndex_1][m1]*c[SecondIndex_1][m2]*c[FirstIndex_2][m3]*conj(c[SecondIndex_2][m1+m2+m3])));

                                                    if(abs(-m1-m2+m3)<=el_4)
                                                        H_fourthorder += (lambda4/24.0)*ThirdGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*NaNb*(-m1-m2+m3+FifthIndex_2) + 4*NaNbNc*(m1+m2+SixthIndex)]*double(1-2*((m1+m2+m3)%2))*2*(real(conj(c[FirstIndex_1][m1])*conj(c[SecondIndex_1][m2])*c[FirstIndex_2][m3]*conj(c[SecondIndex_2][-m1-m2+m3])));
                                                }
                                            }
                                        }

                                        for(int m2=1;m2<=min(m1,el_2);m2++)
                                        {
                                                        if((abs(-m1+m2)<=el_4) && abs(m1-m2)<=l) {

                                                                FourthGauntTerm = Gaunt_matrix_4th_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*NaNb*(-m2+FifthIndex_1) + 4*NaNbNc*(m2-m1+SixthIndex)];
                                                                H_fourthorder += (lambda4/24.0)*FourthGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*Na*Nb*(-m1+m2+FifthIndex_2)+4*Na*Nb*Nc*(m1-m2+SixthIndex)]*double(1-2*(m2%2))*2*(real(c[FirstIndex_2][0]*c[FirstIndex_1][m1]*conj(c[SecondIndex_1][m2])*conj(c[SecondIndex_2][m1-m2])));
                                                                                }
                                        }

                                        for(int m2=m1+1;m2<=el_2;m2++)
                                        {
                                                if((abs(-m1+m2)<=el_4) && abs(m1-m2)<=l) {
                                                FourthGauntTerm = Gaunt_matrix_4th_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)];

                                                        H_fourthorder += (lambda4/24.0)*FourthGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*Na*Nb*(-m1+m2+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)]*double(1-2*(m1%2))*2*(real(c[FirstIndex_2][0]*c[FirstIndex_1][m1]*conj(c[SecondIndex_1][m2])*c[SecondIndex_2][-m1+m2]));

                                                                                        }
                                        }
for(int m3=1;m3<=el_3;m3++)
                                        {
                                                for(int m2=1;m2<=min(m1+m3,el_2);m2++)
                                                {
                                                                if((abs(-m1+m2-m3)<=el_4) && abs(m1-m2)<=l) {
                                                   FourthGauntTerm = Gaunt_matrix_4th_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)];
                                                                        H_fourthorder += (lambda4/24.0)*FourthGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*Na*Nb*(-m1+m2-m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)]*double(1-2*((m2+m3)%2))*2*(real(c[FirstIndex_1][m1]*conj(c[SecondIndex_1][m2])*c[FirstIndex_2][m3]*conj(c[SecondIndex_2][m1+m3-m2])));
                                                                                }
                                                }
for(int m2=m1+m3+1;m2<=el_2;m2++)
                                                {
                                                        if((abs(-m1+m2-m3)<=el_4) && abs(m1-m2)<=l) {
                                                        FourthGauntTerm = Gaunt_matrix_4th_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)];
                                                                H_fourthorder += (lambda4/24.0)*FourthGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2)+4*Na*Nb*(-m1+m2-m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)]*double(1-2*(m1%2))*2*(real(c[FirstIndex_1][m1]*conj(c[SecondIndex_1][m2])*c[FirstIndex_2][m3]*c[SecondIndex_2][-m1-m3+m2]));
                                                                        }
                                                }
                                        }
                                }
for(int m2=1;m2<=el_2;m2++)
                                {
                                        for(int m3=1;m3<=el_3;m3++)
                                        {
                                                for(int m1=1;m1<=min(m2+m3,el_1);m1++)
                                                {
                                                                if((abs(-m1+m2+m3)<=el_4) && abs(m1-m2)<=l) {
                                                                FourthGauntTerm = Gaunt_matrix_4th_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)];
                                                                        H_fourthorder += (lambda4/24.0)*FourthGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*Na*Nb*(-m1+m2+m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)]*double(1-2*((m1+m3)%2))*2*(real(c[FirstIndex_1][m1]*conj(c[SecondIndex_1][m2])*conj(c[FirstIndex_2][m3])*c[SecondIndex_2][-m1+m2+m3]));

                                                                                                        }
                                                }

                                                for(int m1=m2+m3+1;m1<=el_1;m1++)
                                                {
                                                        if((abs(-m1+m2+m3)<=el_4) && abs(m1-m2)<=l) {
                                                        FourthGauntTerm = Gaunt_matrix_4th_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)];
                                                                H_fourthorder += (lambda4/24.0)*FourthGauntTerm*Gaunt_matrix_4th_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*Na*Nb*(-m1+m2+m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)]*double(1-2*(m2%2))*2*(real(c[FirstIndex_1][m1]*conj(c[SecondIndex_1][m2])*conj(c[FirstIndex_2][m3])*conj(c[SecondIndex_2][m1-m2-m3])));

                                                        }
                                                }
                                        }
                                }
                           }
                        }
                }
            }
        }
    }
    H = H_gradient+H_thirdorder+H_fourthorder;

    for(int b = 0; b<2; b++)
        delete [] c[b];
    delete [] c;


   return H;
}

//DERIVATIVE CALCULATION - NUMERICAL APPROXIMATION
DP NR::dfridr(DP func(const Vec_I_DP &x, double* Gaunt_matrix, double* Gaunt_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err)
{
        const int NTAB = 10; //sets maximum size of tableau
        const DP CON = 1.4, CON2=CON*CON; //stepsize decreased by CON at each iteration
        const DP BIG = numeric_limits<DP>::max();
        const DP SAFE = 2.0; //return when error is SAFE worse than the best so far
        int i,j,m;
        const int dim=4*el_not+4;
        DP errt,fac,hh,ans;
        DP point1_array[dim];
        DP point2_array[dim];
        Mat_DP a(NTAB,NTAB);
        Vec_O_DP pder(dim);
        DP dH = 0.0;

        for(int index=0;index<dim;index++)
        {

                if(h==0.0)
                        nrerror("h must be nonzero in dfridr.");
                hh=h;
                for(int m=0;m<dim;m++)
                {
                        point1_array[m]=x[m];
                        point2_array[m]=x[m];

                }
                point1_array[index]+=hh;
                point2_array[index]-=hh;
                Vec_IO_DP point1(point1_array,dim);
                Vec_IO_DP point2(point2_array,dim);
  
                a[0][0] = (func(point1, Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2)-func(point2,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2))/(2.0*hh);

                err=BIG;
                for(int i=1;i<NTAB;i++)
                {
                        point1[index]-=hh;
                        point2[index]+=hh;
                        hh/=CON;
			point1[index]+=hh;
                        point2[index]-=hh;
                        a[0][i]=(func(point1,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2)-func(point2,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2))/(2.0*hh);
                        fac=CON2;
                        for(int j=1;j<=i;j++) //computes extrapolations of various orders
                        {
                                a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
                                fac=CON2*fac;
                                errt=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
                                if(errt<=err)
                                {
                                        err=errt;
                                        ans=a[j][i];
                                }

                        }
                        if(fabs(a[i][i]-a[i-1][i-1])>=SAFE*err)
                        {
                                break;
                        }
                }
                pd[index] = ans;
                dH += pd[index]*pd[index];
        }
        return sqrt(dH);
}

//CALCULATES HESSIAN - NUMERICAL APPROXIMATION
double nhessian(Vec_I_DP &x, double *hess,  DP funk(Vec_I_DP &, double* Gaunt_matrix, double* Gaunt_4th_1, double *Gaunt_matrix_4th_2))
{

        DP ptempl[(4*el_not+4)];
        DP ptempr[(4*el_not+4)];
        DP ptempc[(4*el_not+4)];
        DP ptempc2[(4*el_not+4)];
        for (int i=0; i<=4*el_not+3; i++)
        {
                ptempl[i]=x[i];
                ptempr[i]=x[i];
                ptempc[i]=x[i];
                ptempc2[i]=x[i];
        }
        double del=1e-5;

        for(int k=0; k<=4*el_not+3; k++)
        {
                ptempr[k]=ptempr[k]+del;
                ptempl[k]=ptempl[k]-del;

                Vec_I_DP tempr(ptempr,4*el_not+4);
                Vec_I_DP templ(ptempl,4*el_not+4);
                Vec_I_DP tempc(ptempc,4*el_not+4);
                hess[k+k*hs]=(funk(tempr,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2)+funk(templ,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2)-2*funk(tempc,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2))/(del*del);
                ptempr[k]=x[k];
                ptempl[k]=x[k];

        }

        for(int k=0; k<=4*el_not+3; k++) {
         for(int g=k+1; g<=4*el_not+3; g++) {
                ptempr[k]=ptempr[k]+del;
                ptempr[g]=ptempr[g]+del;
                ptempl[k]=ptempl[k]-del;
                ptempl[g]=ptempl[g]-del;
                ptempc[k]=ptempc[k]+del;
                ptempc[g]=ptempc[g]-del;
                ptempc2[k]=ptempc2[k]-del;
                ptempc2[g]=ptempc2[g]+del;

		Vec_I_DP tempr(ptempr,4*el_not+4);
                Vec_I_DP templ(ptempl,4*el_not+4);
                Vec_I_DP tempc(ptempc,4*el_not+4);
                Vec_I_DP tempc2(ptempc2,4*el_not+4);
        hess[k+g*hs]=(funk(tempr,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2)+funk(templ,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2)-funk(tempc,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2)-funk(tempc2,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2))/(4*del*del);
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

void NR::mnbrak(DP &ax, DP &bx, DP &cx, DP &fa, DP &fb, DP &fc, DP f(const DP, double*, double*, double*))
{
	const DP GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
	DP ulim,u,r,q,fu;
	
	fa = f(ax,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
	fb = f(bx,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
	if(fb > fa)
	{
		SWAP(ax,bx);
		SWAP(fb,fa);
	}
	cx = bx+GOLD*(bx-ax);
	fc = f(cx,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
	while(fb > fc)
	{
		r = (bx-ax)*(fb-fc);
		q = (bx-cx)*(fb-fa);
		u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim = bx+GLIMIT*(cx-bx);
		if((bx-u)*(u-cx) > 0.0)
		{
			fu = f(u,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
			if(fu < fc)
			{
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
				return;
			}
			else if (fu > fb)
			{
				cx = u;
				fc = fu;
				return;
			}
			u = cx+GOLD*(cx-bx);
			fu = f(u,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
		}
		else if((cx-u)*(u-ulim) > 0.0)
		{
			fu = f(u,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
			if(fu < fc)
			{
				shft3(bx,cx,u,cx+GOLD*(u-bx));
				shft3(fb,fc,fu,f(u,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2));
			}
		}
		else if((u-ulim)*(ulim-cx)>=0.0)
		{
			u = ulim;
			fu = f(u,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
		}
		else
		{
			u = cx+GOLD*(cx-bx);
			fu = f(u,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
		}
		shft3(ax,bx,cx,u);
		shft3(fa,fb,fc,fu);
	}
}
DP NR::dbrent(const DP ax, const DP bx, const DP cx, DP f(const DP, double*, double*, double*), DP df(DP func(const DP, double*, double*, double*), const DP, const DP, DP &), const DP tol, DP &xmin, const DP h, DP &err)
{
	const int ITMAX = 200; //200 works normally
	const DP ZEPS = numeric_limits<DP>::epsilon()*1.0e-3;
	bool ok1,ok2;
	int iter;
	DP a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
	DP fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	a = (ax < cx ? ax:cx);
	b = (ax > cx ? ax:cx);
	x=w=v=bx;
	fw=fv=fx=f(x,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
	dw=dv=dx=df(f,x,h,err);
	for(iter=0;iter<ITMAX;iter++)
	{
		xm = 0.5*(a+b);
		tol1 = tol*fabs(x)+ZEPS;
		tol2 = 2.0*tol1;
		if(fabs(x-xm) <= (tol2-0.5*(b-a)))
		{
			xmin = x;
			return fx;
		}
		if(fabs(e) > tol1)
		{
			d1 = 2.0*(b-a);
			d2 = d1;
			if(dw!=dx)
				d1 = (w-x)*dx/(dx-dw);
			if(dv!=dx)
				d2 = (v-x)*dx/(dx-dv);
			u1 = x+d1;
			u2 = x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde = e;
			e = d;
			if(ok1||ok2)
			{
				if(ok1&&ok2)
					d=(fabs(d1) < fabs(d2) ? d1:d2);
				else if (ok1)
					d = d1;
				else
					d = d2;
				if(fabs(d) <= fabs(0.5*olde))
				{
					u = x+d;
					if(u-a < tol2 || b-u < tol2)
						d = SIGN(tol1,xm-x);
				}
				else
				{
					d = 0.5*(e = (dx >= 0.0 ? a-x : b-x));
				}
			}
			else
			{
				d = 0.5*(e = (dx >= 0.0 ? a-x : b-x));
			}
		}
		else
		{
			d = 0.5*(e = (dx >= 0.0 ? a-x : b-x));
		}
		if(fabs(d) >= tol1)
		{
			u = x+d;
			fu = f(u,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
		}
		else
		{
			u = x+SIGN(tol1,d);
			fu = f(u,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
			if(fu > fx)
			{
				xmin=x;
				return fx;
			}
		}
		du = df(f,u,h,err);
		if (fu <=fx)
		{
			if (u>=x)
				a=x;
			else
				b=x;
			mov3(v,fv,dv,w,fw,dw);
			mov3(w,fw,dw,u,fu,du);
			mov3(x,fx,dx,u,fu,du);
		}
		else
		{
			if(u<x)
				a=u;
			else
				b=u;
			if(fu<=fw || w==x)
			{
				mov3(v,fv,dv,w,fw,dw);
				mov3(w,fw,dw,u,fu,du);
			}
			else if (fu < fv || v==x || v==w)
			{
				mov3(v,fv,dv,u,fu,du);
			}
		}
	}
	nrerror("Too many interations in routine dbrent");
	xmin = x;
	return fx;
}

DP NR::f1dim(const DP x, double* Gaunt_matrix, double* Gaunt_4th_1, double *Gaunt_matrix_4th_2)
{
	int j;
	
	Vec_DP xt(ncom);
	Vec_DP &pcom = *pcom_p;
	Vec_DP &xicom = *xicom_p;
	for(j = 0; j<ncom; j++)
		xt[j] = pcom[j] + x*xicom[j];
	return nrfunc(xt,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);
}
DP NR::df1dim(DP func(const DP x, double* Gaunt_matrix, double* Gaunt_4th_1, double *Gaunt_matrix_4th_2), const DP x, const DP h, DP &err)
{
        int j;
        DP df1 = 0.0;
        Vec_DP xt(ncom);
	Vec_DP df(ncom);

        Vec_DP &pcom = *pcom_p;
	Vec_DP &xicom = *xicom_p;
        for(j=0;j<ncom;j++)
                xt[j] = pcom[j]+x*xicom[j];

	nrdfun(nrfunc,xt,df,h,err);
	for(j=0;j<ncom;j++)
                df1+=df[j]*xicom[j];
	return df1;
}
void NR::dlinmin(Vec_IO_DP &p, Vec_IO_DP &xi, DP &fret, DP func(Vec_I_DP &x, double* Gaunt_matrix, double* Gaunt_4th_1, double *Gaunt_matrix_4th_2), DP dfridr(DP func(const Vec_I_DP &x, double* Gaunt_matrix, double* Gaunt_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err), const DP h, DP &err)
{
        const DP TOL=2.0e-8;
        int j;
        DP xx,xmin,fx,fb,fa,bx,ax;

        int n=p.size();
        ncom=n;
        pcom_p=new Vec_DP(n);
        xicom_p = new Vec_DP(n);
        nrfunc = func;
        nrdfun = dfridr;
        Vec_DP &pcom = *pcom_p,&xicom=*xicom_p;
        for(j=0;j<n;j++)
        {
                pcom[j] = p[j];
                xicom[j] = xi[j];
        }
        ax = 0.0;
        xx = 1.0;
        mnbrak(ax,xx,bx,fa,fx,fb,f1dim);
        fret = dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin,h,err);
	for(j=0;j<n;j++)
        {
                xi[j] *= xmin;
                p[j] += xi[j];
        }
	int endcount=2;
        delete xicom_p;
        delete pcom_p;
}

void NR::frprmn(Vec_IO_DP &p, const DP ftol, int &iter, DP &fret, DP func(Vec_I_DP &x, double* Gaunt_matrix, double* Gaunt_4th_1, double *Gaunt_matrix_4th_2), DP dfridr(DP func(const Vec_I_DP &x, double* Gaunt_matrix, double* Gaunt_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err), const DP h, DP &err)
{
        const int ITMAX=200; //max allowed number of iterations
        const DP EPS=1.0e-18; //EPS is small number to rectify the special case of converging to exactly zero function value
        int j,its;
        DP gg,gam,fp,dgg;
        int n=p.size();
        Vec_DP g(n),hp(n),xi(n);
        fp=func(p,Gaunt_matrix_3rd_order,Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2);     //initializatios
	dfridr(func,p,xi,h,err);
	for(j=0;j<n;j++)
        {
                g[j] = -xi[j];
                xi[j] = hp[j] = g[j];
        }
        for(its=0;its<ITMAX;its++) //loop over iterations
        {
      		iter=its;
                dlinmin(p,xi,fret,func,dfridr,h,err);    //next statement is the normal return
   		if(2.0*fabs(fret-fp) <= ftol*(fabs(fret)+fabs(fp)+EPS))
		{
			return;
		}
                fp=fret;
		dfridr(func,p,xi,h,err);

                dgg=gg=0.0;
                for(j=0;j<n;j++)
                {
                        gg += g[j]*g[j];
                        dgg += (xi[j]+g[j])*xi[j];  
                }
                if(gg==0.0)            
                        return;
                gam=dgg/gg;
                for(j=0;j<n;j++)
                {
                        g[j] = -xi[j];
                        xi[j]=hp[j]=g[j]+gam*hp[j];
                }
        }
        nrerror("Too many iterations in frprmn");

}

int main()
{
    clock_t t1;
    t1 = clock();
    ofstream endcms,startingcms,cms,parameters,time; 
    endcms.open("endcms.txt");
    startingcms.open("startingcms.txt");
    parameters.open("parameters.txt");
    time.open("time.txt");
    random_device rd;
    mt19937 gen(300);
    wig_table_init(2*100,3); //define using largest 2*3j value and largest wigner type (i.e. 3j, 6j, 9j) - however this particular 2*100, 9 is commonly used becuase it doesn't take a lot of memory and yet it works for most programs
    wig_temp_init(2*100);
    double hess_array[(int)(4*el_not+4)*(4*el_not+4)];
    double t;
    DP fret,der;
    double g_array[(int)(4*el_not+4)];
    for(int g=0;g<4*el_not+4;g++) g_array[g]=0; 
    Vec_IO_DP g(g_array,4*el_not+4);
    parameters<<"l: "<<el_not<<" , "<<el_not+1<<endl;
    parameters<<"el_not: "<<l_not<<endl;
    parameters<<"tau: "<<tau<<endl;
    parameters<<"lambda3: "<<lambda3<<endl;
    parameters<<"lambda4: "<<lambda4<<endl;
    parameters<<endl;

//3rd order Gaunt Matrix
    int FirstIndex,SecondIndex,ThirdIndex;
    for(int el_1=el_not;el_1<=el_not+1;el_1++)
    {
        FirstIndex = el_1-el_not;
        for(int el_2=el_not;el_2<=el_not+1;el_2++)
        {
            SecondIndex = el_2-el_not;
            for(int el_3=el_not;el_3<=el_not+1;el_3++)
            {
                ThirdIndex = el_3-el_not;
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

//Fourth Order Gaunt Matrix
  int FirstIndex_1, SecondIndex_1, FourthIndex_1, FifthIndex_1, SixthIndex, FirstIndex_2, SecondIndex_2, FourthIndex_2, FifthIndex_2;
  for(int el_1=el_not;el_1<=el_not+1;el_1++)
  {
        FirstIndex_1 = el_1-el_not;

        FourthIndex_1 = el_1+el_not+2;
        for(int el_2=el_not;el_2<=el_not+1;el_2++)
        {
            SecondIndex_1 = el_2-el_not;
            FifthIndex_1 = el_2+el_not+2;
            for(int el_3=el_not;el_3<=el_not+1;el_3++)
            {
                FirstIndex_2 = el_3-el_not;
                FourthIndex_2 = el_3+el_not+2;
                for(int el_4=el_not;el_4<=el_not+1;el_4++)
                {
                    SecondIndex_2 = el_4-el_not;
                    FifthIndex_2 = el_4+el_not+2;
                    for(int l=0;l<=(2*el_not)+2;l++)
                    {
                        if((el_1+el_2+l)%2==0 && (el_3+el_4+l)%2==0)
                        {
                                SixthIndex = l+2*el_not+2;

                                Gaunt_matrix_4th_order_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*FourthIndex_1 + 4*Na*Nb*FifthIndex_1 + 4*Na*Nb*Nc*SixthIndex] = Gaunt(el_1,el_2,l,0,0,0);
                                Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*Na*Nb*FifthIndex_2 + 4*Na*Nb*Nc*SixthIndex] = Gaunt(el_3,el_4,l,0,0,0);

                                for(int m3=1;m3<=el_3;m3++)
                                {
                                     Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*Na*Nb*(FourthIndex_2-m3) + 4*Na*Nb*Nc*SixthIndex] = Gaunt(el_3,el_4,l,m3,-m3,0);
                                }
                                for(int m2=1;m2<=el_2;m2++)
                                {
                                        Gaunt_matrix_4th_order_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*FourthIndex_1 + 4*Na*Nb*(m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m2+SixthIndex)] =  Gaunt(el_1,el_2,l,0,m2,-m2);

                                        Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*Na*Nb*(-m2+FifthIndex_2) + 4*Na*Nb*Nc*(m2+SixthIndex)] = Gaunt(el_3,el_4,l,0,-m2,m2); 

for(int m3=1;m3<=m2;m3++)
                                        {
                                           if(abs(m2+m3)<=el_4 && m3<=el_3)
                                           {
Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*Na*Nb*(-m2-m3+FifthIndex_2)+4*Na*Nb*Nc*(m2+SixthIndex)] = Gaunt(el_3,el_4,l,m3,-m2-m3,m2);
                                           }
                                           if(abs(m3-m2)<=el_4 && m3<=el_3)
                                           {
                                                Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*Na*Nb*(m3-m2+FifthIndex_2) + 4*Na*Nb*Nc*(m2+SixthIndex)] = Gaunt(el_3,el_4,l,-m3,m3-m2,m2);
                                           }
                                        }

                                        for(int m3=m2+1;m3<=el_3;m3++)
                                        {
                                           if(abs(m2+m3)<=el_4)
                                           {
                                                Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*Na*Nb*(-m2-m3+FifthIndex_2) + 4*Na*Nb*Nc*(m2+SixthIndex)] = Gaunt(el_3,el_4,l,m3,-m2-m3,m2);
                                           }
                                           if(abs(m3-m2)<=el_4)
                                           {
                                                Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*Na*Nb*(m3-m2+FifthIndex_2) + 4*Na*Nb*Nc*(m2+SixthIndex)] = Gaunt(el_3,el_4,l,-m3,m3-m2,m2);
                                           }
                                        }
                                }
for(int m1=1;m1<=el_1;m1++)
                                {
                                        for(int m2=1;m2<=el_2;m2++)
                                        {
                                            Gaunt_matrix_4th_order_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1-m2+SixthIndex)] = Gaunt(el_1,el_2,l,m1,m2,-m1-m2);
                                            if(m1+m2<=l)
                                            {
                                                Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*Na*Nb*(-m1-m2+FifthIndex_2) + 4*Na*Nb*Nc*(m1+m2+SixthIndex)] = Gaunt(el_3,el_4,l,0,-m1-m2,m1+m2);
                                                for(int m3=1;m3<=m1+m2;m3++)
                                                {
                                                    if(abs(m1+m2+m3)<=el_4 && m3<=el_3)
                                                    {
                                                        Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*Na*Nb*(-m1-m2-m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1+m2+SixthIndex)] = Gaunt(el_3,el_4,l,m3,-m1-m2-m3,m1+m2);
                                                    }
if(abs(-m1-m2+m3)<=el_4 && m3<=el_3)
                                                    {
                                                        Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*Na*Nb*(-m1-m2+m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1+m2+SixthIndex)] = Gaunt(el_3,el_4,l,-m3,-m1-m2+m3,m1+m2);
                                                    }
                                                }
                                                for(int m3=m1+m2+1;m3<=el_3;m3++)
                                                {
                                                    if(abs(m1+m2+m3)<=el_4)
                                                    {
                                                        Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*Na*Nb*(-m1-m2-m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1+m2+SixthIndex)] = Gaunt(el_3,el_4,l,m3,-m1-m2-m3,m1+m2);
                                                    }
                                                    if(abs(-m1-m2+m3)<=el_4)
                                                    {
                                                        Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*Na*Nb*(-m1-m2+m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1+m2+SixthIndex)] = Gaunt(el_3,el_4,l,-m3,-m1-m2+m3,m1+m2);
                                                     }
                                                }
                                            }
                                        }
for(int m2=1;m2<=m1;m2++)
                                        {
                                                Gaunt_matrix_4th_order_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(m2-m1+SixthIndex)] = Gaunt(el_1,el_2,l,m1,-m2,-m1+m2);
                                                if((abs(-m1+m2)<=el_4) && abs(m1-m2)<=l)
                                                {
                                                        Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*Na*Nb*(-m1+m2+FifthIndex_2)+4*Na*Nb*Nc*(m1-m2+SixthIndex)] = Gaunt(el_3,el_4,l,0,-m1+m2,m1-m2);
                                                }
                                        }for(int m2=m1+1;m2<=el_2;m2++)
                                        {
                                                Gaunt_matrix_4th_order_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)] = Gaunt(el_1,el_2,l,m1,-m2,-m1+m2);
                                                if((abs(-m1+m2)<=el_4) && abs(m1-m2)<=l)
                                                {
                                                        Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*FourthIndex_2 + 4*Na*Nb*(-m1+m2+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)] = Gaunt(el_3,el_4,l,0,-m1+m2,m1-m2);
                                                }
                                        }

                                        for(int m3=1;m3<=el_3;m3++)
                                        {
                                                for(int m2=1;m2<=m1+m3;m2++)
                                                {
                                                        if(m2<=el_2) 
                                                        {
                                                                Gaunt_matrix_4th_order_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)] = Gaunt(el_1,el_2,l,m1,-m2,-m1+m2);
                                                                if((abs(-m1+m2-m3)<=el_4) && abs(m1-m2)<=l)
                                                                {
                                                                        Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2) + 4*Na*Nb*(-m1+m2-m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)] = Gaunt(el_3,el_4,l,m3,-m1+m2-m3,m1-m2);
                                                                }
                                                        }
                                                }
for(int m2=m1+m3+1;m2<=el_2;m2++)
                                                {
                                                        Gaunt_matrix_4th_order_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)] = Gaunt(el_1,el_2,l,m1,-m2,-m1+m2);
                                                if((abs(-m1+m2-m3)<=el_4) && abs(m1-m2)<=l)
                                                        {
                                                                Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(m3+FourthIndex_2)+4*Na*Nb*(-m1+m2-m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)] = Gaunt(el_3,el_4,l,m3,-m1+m2-m3,m1-m2);
                                                        }
                                                }
                                        }
                                }
                                for(int m2=1;m2<=el_2;m2++)
                                {
                                        for(int m3=1;m3<=el_3;m3++)
                                        {
                                                for(int m1=1;m1<=m2+m3;m1++)
                                                {
                                                        if(m1<=el_1) 
                                                        {
                                                                Gaunt_matrix_4th_order_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)] = Gaunt(el_1,el_2,l,m1,-m2,-m1+m2);
                                                                if((abs(-m1+m2+m3)<=el_4) && abs(m1-m2)<=l)
                                                                {
                                                                        Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*Na*Nb*(-m1+m2+m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)] = Gaunt(el_3,el_4,l,-m3,-m1+m2+m3,m1-m2);
                                                                }
                                                        }
                                                }
                                                for(int m1=m2+m3+1;m1<=el_1;m1++)
                                                {
                                                        Gaunt_matrix_4th_order_1[FirstIndex_1 + 2*SecondIndex_1 + 4*l + 4*Na*(m1+FourthIndex_1) + 4*Na*Nb*(-m2+FifthIndex_1) + 4*Na*Nb*Nc*(-m1+m2+SixthIndex)] = Gaunt(el_1,el_2,l,m1,-m2,-m1+m2);
                                                        if((abs(-m1+m2+m3)<=el_4) && abs(m1-m2)<=l)
                                                        {
                                                                Gaunt_matrix_4th_order_2[FirstIndex_2 + 2*SecondIndex_2 + 4*l + 4*Na*(-m3+FourthIndex_2) + 4*Na*Nb*(-m1+m2+m3+FifthIndex_2) + 4*Na*Nb*Nc*(m1-m2+SixthIndex)] = Gaunt(el_3,el_4,l,-m3,-m1+m2+m3,m1-m2);
                                                        }
                                                }
                                        }
                                }
                           }
                        }
                  }
            }
        }
   }



/*Required definitons for gradient descent:
 * starting p point, ftol, iter, fret*/
//Defines starting point, p:    
    DP p_array[4*el_not+4];
    parameters<<"Starting cms:"<<endl;

    for(int a=0;a<4*el_not+4;a++)
    {
        p_array[a] = 1.0*sqrt((-1)*tau/lambda4)*(double(gen())/double(gen.max()));
        parameters<<p_array[a]<<" ";
    }

   Vec_IO_DP p(p_array,4*el_not+4);

    startingcms<<p[0]<<" 0"<<endl;
    for(int i=2;i<=el_not+1;i++)
        startingcms<<p[i]<<" "<<p[i+2*el_not+1]<<endl;

    startingcms<<p[1]<<" 0"<<endl;
    for(int j=el_not+2;j<=2*el_not+2;j++)
        startingcms<<p[j]<<" "<<p[j+2*el_not+1]<<endl;
       

    parameters<<endl;    
    double magg;
    const DP h = 1.0; 
    DP pd_array[4*el_not+4];
    for(int i=0;i<4*el_not+4;i++)
        pd_array[i]=0.0;
    Vec_O_DP pd(pd_array,4*el_not+4);
    DP err = 1.0;

    parameters<<"Function evaluation at starting cms: "<<func(p,Gaunt_matrix_3rd_order, Gaunt_matrix_4th_order_1, Gaunt_matrix_4th_order_2)<<endl;
    parameters<<"Starting derivative: "<<NR::dfridr(func,p,pd,h,err)<<endl;
    const DP ftol = 1.0e-6; //1.0e-8;
    parameters<<"FTolerance: "<<ftol<<endl;
     
    int iter= 10000;
    parameters<<"Number of iterations: "<<iter<<endl;
    NR::frprmn(p,ftol,iter,fret,func,NR::dfridr,h,err); 
   
    endcms<<p[0]<<" 0"<<endl;
    int endcount_2=2;
    for(int i=0;i<el_not;i++)
    {
        endcms<<p[endcount_2]<<" "<<p[2*el_not+1+endcount_2]<<endl;
        endcount_2++;
    }

    endcms<<p[1]<<" 0"<<endl;
    for(int j=0;j<=el_not;j++)
    {
        endcms<<p[endcount_2]<<" "<<p[2*el_not+1+endcount_2]<<endl;
        endcount_2++;
    }


   magg=nhessian(p,hess_array,func);
   Eigen::MatrixXd hhh(hs,hs);
   for (int m=0; m<hs*hs; m++) {
                int x = m%(4*el_not+4); 
                int y = (m/hs)%hs;
        hhh(x,y)=hess_array[m];
        }

    Eigen::EigenSolver<Eigen::MatrixXd> es(hhh);

    parameters << "Eigenvalues: " << endl << es.eigenvalues() << endl;
   parameters<<"Lowest point derivative value: "<<NR::dfridr(func,p,pd,h,err)<<endl;

   parameters<<"Lowest_Hamiltonian_value: "<<fret<<endl;
   endcms.close();
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

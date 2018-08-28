#ifndef _NR_H_
#define _NR_H_
#include <fstream>
#include <complex>
#include <cmath>
#include "nrutil.h"
#include "nrtypes.h"
#define l0 3
using namespace std;


namespace NR
{
    DP dfridr(DP funk(Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err);
    void amebsa(Mat_IO_DP &p, Vec_IO_DP &y, Vec_O_DP &pb, DP &yb, const DP ftol, DP funk(Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), int &iter, const DP temptr);
    DP amotsa(Mat_IO_DP &p, Vec_O_DP &y, Vec_IO_DP &psum, Vec_IO_DP &pb, DP &yb, DP funk(Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const int ihi, DP &yhi, const DP fac);

	//DP ran1(int &idum);
}
#endif /*_NR_H_*/

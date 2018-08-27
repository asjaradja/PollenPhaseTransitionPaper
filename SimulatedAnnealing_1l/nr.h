#ifndef _NR_H_
#define _NR_H_
#include <fstream>
#include <complex>
#include <cmath>
#include "nrutil.h"
#include "nrtypes.h"
#define el_not 4

using namespace std;


namespace NR

{
   void amebsa(Mat_IO_DP &p, Vec_IO_DP &y, Vec_O_DP &pb, DP &yb, const DP ftol, DP func(Vec_I_DP &x, double *Gaunt_matrix, double* Gaunt_4th), int &iter, const DP temptr);
   DP amotsa(Mat_IO_DP &p, Vec_O_DP &y, Vec_IO_DP &psum, Vec_IO_DP &pb, DP &yb, DP func(Vec_I_DP &x, double *Gaunt_matrix, double* Gaunt_4th), const int ihi, DP &yhi, const DP fac);
   DP dfridr(DP func(const Vec_I_DP &x, double *Gaunt_matrix, double* Gaunt_4th), const Vec_I_DP x, const DP h, DP &err, int index);
}
#endif /*_NR_H_*/

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
	DP dfridr(DP func(const Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err);
        void mnbrak(DP &ax, DP &bx, DP &cx, DP &fa, DP &fb, DP &fc, DP f(const DP, double*, double*, double*));
        DP dbrent(const DP ax, const DP bx, const DP cx, DP f(const DP, double*, double*, double*), DP df(DP func(const DP, double*, double*, double*), const DP, const DP h, DP &err), const DP tol, DP &xmin, const DP h, DP &err);
        DP f1dim(const DP x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2);
        DP df1dim(DP func(const DP, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const DP x, const DP h, DP &err);
        void dlinmin(Vec_IO_DP &p, Vec_IO_DP &xi, DP &fret, DP func(Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), DP dfridr(DP func(const Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err),const DP h, DP &err);
        void frprmn(Vec_IO_DP &p, const DP ftol, int &iter, DP &fret, DP func(Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), DP dfridr(DP func(const Vec_I_DP &x, double *Gaunt_matrix, double *Gaunt_matrix_4th_1, double *Gaunt_matrix_4th_2), const Vec_I_DP x, Vec_O_DP &pd, const DP h, DP &err),const DP h, DP &err);
}
#endif /*_NR_H_*/

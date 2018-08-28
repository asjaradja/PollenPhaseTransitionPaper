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
	void mnbrak(DP &ax, DP &bx, DP &cx, DP &fa, DP &fb, DP &fc, DP f(const DP, double*, double*));
        DP dbrent(const DP ax, const DP bx, const DP cx, DP f(const DP, double*, double*), DP df(const DP, double*, double*), const DP tol, DP &xmin);
	DP f1dim(const DP x, double* Gaunt_matrix, double* Gaunt_4th);
        DP df1dim(const DP x, double* Gaunt_matrix, double* alp);        
        void dlinmin(Vec_IO_DP &p, Vec_IO_DP &xi, DP &fret, DP func(Vec_I_DP &x, double* Gaunt_matrix, double* Gaunt_4th), DP dfunc(Vec_I_DP &x, Vec_O_DP &g, double *Gaunt_matrix, double* alp));
        void frprmn(Vec_IO_DP &p, const DP ftol, int &iter, DP &fret, DP func(Vec_I_DP &x, double* Gaunt_matrix, double* Gaunt_4th), DP dfunc(Vec_I_DP &x, Vec_O_DP &g, double *Gaunt_matrix, double* alp));	
}
#endif /*_NR_H_*/

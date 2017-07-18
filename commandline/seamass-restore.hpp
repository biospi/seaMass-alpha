//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> bristol.ac.uk>
//
// Copyright (C) 2016  biospi Laboratory, University of Bristol, UK
//
// This file is part of seaMass.
//
// seaMass is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// seaMass is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with seaMass.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef SEAMASS_SEAMASS_PEAK_HPP
#define SEAMASS_SEAMASS_PEAK_HPP

#include "../kernel/intel/types.hpp"
#include <iostream>
#include <BasisBsplineMz.hpp>


template <typename T>
T* alcMat(T *M, int m, int n)
{
    M = (T *)mkl_malloc(m * n * sizeof(T), 64);
    if (M == NULL)
    {
        cout << "ERROR: Can't allocate memory for matrices. Aborting... \n\n"<<endl;
        mkl_free(M);
        return NULL;
    }
    return M;
};

template  <typename T>
void delMat(T *M) { mkl_free(M); };

void matDmul(fp *A, fp *B, fp *C, int  m, int k, int n)
{
    /*
     * Computes real matrix C=alpha*A*B+beta*C using
     * Intel(R) MKL function dgemm, where A, B, and  C are matrices
     * alpha and beta are single precision scalars
     * Input Mat A and B, output C.
     */

    fp alpha = 1.0;
    fp beta = 0.0;

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
             m, n, k, alpha, A, k, B, n, beta, C, n);
};

void matDmul(double *A, double *B, double *C, int  m, int k, int n)
{
    /*
     * Computes real matrix C=alpha*A*B+beta*C using
     * Intel(R) MKL function dgemm, where A, B, and  C are matrices
     * alpha and beta are single precision scalars
     * Input Mat A and B, output C.
     */

    double alpha = 1.0;
    double beta = 0.0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
             m, n, k, alpha, A, k, B, n, beta, C, n);
};

void genMZAxis(vector<double> &mz,Seamass::ControlPoints &cpts, int n, int res)
{
    mz.resize(n);
    double offset = double(cpts.offset[0]+1);
    /*
    //double dmz=1.0/(res -1);
    double dmz=1.0/(res);
    double ppbmz = 1.0033548378 / (pow(2.0, cpts.scale[0]) * 60.0);
    for (lli i = 0; i < mz.size(); ++i) {
        mz[i] = (offset + i*dmz) * ppbmz;
    }
    */
    double dmz=1.0/(res);
    for (lli i = 0; i < mz.size(); ++i)
    {
        mz[i] = pow(2.0, (offset + i*dmz) / double(1L << cpts.scale[0])) + BasisBsplineMz::PROTON_MASS;
    }
}
#endif //SEAMASS_SEAMASS_PEAK_HPP

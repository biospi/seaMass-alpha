//
// Original author: Andrew Dowsey <ranjeet.bhameber <a.t> bristol.ac.uk>
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


template <typename T>
T* alcMat(T *M, int m, int n) { return M = (T *)mkl_malloc(m * n * sizeof(T), 64); };


template  <typename T>
void delMat(T *M) { mkl_free(M); };


template <typename T>
void matDmul(T *A, T *B, T *C, int  m, int k, int n)
{
    /*
     * Computes real matrix C=alpha*A*B+beta*C using
     * Intel(R) MKL function dgemm, where A, B, and  C are matrices
     * alpha and beta are single precision scalars
     * Input Mat A and B, output C.
     */

    int i;
    double alpha = 1.0;
    double beta = 0.0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
             m, n, k, alpha, A, k, B, n, beta, C, n);

    /*
    int i;
    double alpha, beta;

    printf ("\n This example computes real matrix C=alpha*A*B+beta*C using \n"
                    " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
                    " alpha and beta are double precision scalars\n\n");

    printf (" Initializing data for matrix multiplication C=A*B for matrix \n"
                    " A(%ix%i) and matrix B(%ix%i)\n\n", m, p, p, n);
    alpha = 1.0; beta = 0.0;

    printf (" Allocating memory for matrices aligned on 64-byte boundary for better \n"
                    " performance \n\n");
//    A = (double *)mkl_malloc( m*p*sizeof( double ), 64 );
//    B = (double *)mkl_malloc( p*n*sizeof( double ), 64 );
//    C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );
//    if (A == NULL || B == NULL || C == NULL) {
//        printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
//        mkl_free(C);
//    }

    for (i = 0; i < (m*n); i++) {
        C[i] = 0.0;
    }

    printf (" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, p, alpha, A, p, B, n, beta, C, n);
    printf ("\n Computations completed.\n\n");

    printf (" Example completed. \n\n");*/
};

void genMZAxis(vector<double> &mz,Seamass::ControlPoints &cpts, int n, int res)
{
    mz.resize(n);
    double offset = double(cpts.offset[0]+1);
    double dmz=1.0/(res -1);
    double ppbmz = 1.0033548378 / (pow(2.0, cpts.scale[0]) * 60.0);
    for (lli i = 0; i < mz.size(); ++i) {
        mz[i] = (offset + i*dmz) * ppbmz;
    }
}
#endif //SEAMASS_SEAMASS_PEAK_HPP

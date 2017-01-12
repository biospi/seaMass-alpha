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


#ifndef _SEAMASS_MATH_MATRIXSPARSE_HPP_
#define _SEAMASS_MATH_MATRIXSPARSE_HPP_


#include <iostream>
#include <vector>


#define MKL_ILP64 // use 64 bit addressing (comment out for 32 bit)
#include <mkl.h>
#include <mkl_spblas.h>


typedef float fp; // fp is the selected floating point precision (float or double)
typedef MKL_INT ii; // ii is the selected addressing (32 or 64 bit)
typedef MKL_INT64 li; // li is always 64 bit


void printNumThreads();
double getWallTime();


class MatrixSparse
{
public:
	MatrixSparse();
	~MatrixSparse();
    
    void free();

	void init(ii m, ii n); // create empty matrix
	void init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind); // create from COO matrix
    void set(fp v); // set all non-zeros to constant 'v'
    bool isTransposed() const;
    
    ii m() const;
    ii n() const;
    li size() const;
    ii nnz() const;
    li mem() const;
  
    enum class Operation { NONE, TRANSPOSE, PACK_ROWS, UNPACK_ROWS };
    void copy(const MatrixSparse& a, Operation operation = Operation::NONE);
    void copy(const MatrixSparse& a, fp pruneThreshold);
    void output(fp* vs) const;

    void deleteRows(const MatrixSparse& a);

    // generalised sparse matrix multiplication
	void mul(bool transposeA, const MatrixSparse& a, const MatrixSparse& b, bool accumulate, bool transpose);

    // elementwise operations
	void elementwiseAdd(fp beta);
	void elementwiseMul(fp beta);
    void elementwiseMul(const MatrixSparse& a);
    void elementwiseDiv(const MatrixSparse& a);
	void elementwiseSqr();
	void elementwiseSqrt();
	void elementwiseLn();
	void elementwisePow(fp power);

    // aggregate operations
	fp sum() const;
	fp sumSqrs() const;
    fp sumSqrDiffs(const MatrixSparse& a) const;
    
    // note: these elementwise operations ONLY considers the non-zero elements of THIS matrix
    void subsetElementwiseCopy(const MatrixSparse& a);
    //void subsetElementwiseMul(const MatrixSparse& a);
    void subsetElementwiseDiv(const MatrixSparse& a); // needed for error calculation atm
    
    // note: this aggregate operation ONLY considers the non-zero elements of THIS matrix
	//double subsetSumSqrDiffs(const MatrixSparse& a) const;

private:
	ii m_;   // number of rows
	ii n_;   // number of columns
    bool isEmpty_; // MKL doesn't like empty sparse matrices
    bool isTransposed_;
    
    ii* is0_; ii* is1_; ii* js_; fp* vs_; // pointers to CSR array
    sparse_matrix_t mat_;
    bool isMklData_; // true if CSR arrays owned by mat_, false is owned by this object

    sparse_status_t status_; // last MKL function status
    
    friend std::ostream& operator<<(std::ostream& os, const MatrixSparse& mat);
};

std::ostream& operator<<(std::ostream& os, const MatrixSparse& mat);


#endif


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


#ifndef _SEAMASS_MATH_MATRIXSPARSEMKL_HPP_
#define _SEAMASS_MATH_MATRIXSPARSEMKL_HPP_


#include <iostream>
#include <vector>


#define MKL_ILP64 // use 64 bit addressing (comment out for 32 bit)
#include <mkl.h>
#include <mkl_spblas.h>


typedef float fp; // fp is the selected floating point precision (float or double)
typedef MKL_INT ii; // ii is the selected addressing (32 or 64 bit)
typedef MKL_INT64 li; // li is always 64 bit


std::string getThreadInfo();

void resetElapsedTime();
double getElapsedTime();

li getUsedMemory();

std::string getTimeStamp();

void setDebugLevel(int debugLevel);
int getDebugLevel();


class MatrixSparseMKL
{
public:
	MatrixSparseMKL();
	~MatrixSparseMKL();
    
    void init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind); // create from COO matrix
    void free();
    
    ii m() const;
    ii n() const;
    li size() const;
    ii nnz() const;
    ii nnzActual() const;
  
    enum class Operation { NONE, TRANSPOSE, PACK_ROWS, UNPACK_ROWS };
    void copy(const MatrixSparseMKL& a, Operation operation = Operation::NONE);
    void prune(const MatrixSparseMKL& a, fp pruneThreshold);
    void output(fp* vs) const;
    void write(const std::string& filename) const;

    void zeroRowsOfZeroColumns(const MatrixSparseMKL& a, const MatrixSparseMKL& x);

    void add(fp alpha, bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b);
	void matmul(bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b, bool accumulate);
    void mul(fp beta);
    void mul(const MatrixSparseMKL& a);
    void sqr();
    void sqrt();

    // operate only on non-zero elements
    void setNonzeros(fp v);
	void addNonzeros(fp beta);
	void lnNonzeros();
    void expNonzeros();
	void powNonzeros(fp power);
    void divCorrespondingNonzeros(const MatrixSparseMKL& a);

    // aggregate operations
	fp sum() const;
	fp sumSqrs() const;
    fp sumSqrDiffsCorrespondingNonzeros(const MatrixSparseMKL& a) const;
    
    // note: these elementwise operations ONLY considers the non-zero elements of THIS matrix
    void subsetElementwiseCopy(const MatrixSparseMKL& a);

public:
	ii m_;   // number of rows
	ii n_;   // number of columns
    
    ii* is0_; ii* is1_; ii* js_; fp* vs_; // pointers to CSR array
    sparse_matrix_t mat_;
    bool isMklData_; // true if CSR arrays owned by mat_, false is owned by this object

    sparse_status_t status_; // last MKL function status
    
    friend std::ostream& operator<<(std::ostream& os, const MatrixSparseMKL& mat);
};

std::ostream& operator<<(std::ostream& os, const MatrixSparseMKL& mat);


#endif


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


#ifndef SEAMASS_KERNEL_INTEL_MATRIXSPARSE_HPP
#define SEAMASS_KERNEL_INTEL_MATRIXSPARSE_HPP


#include "Matrix.hpp"
#include <vector>


class MatrixSparseView;


class MatrixSparse
{
public:
    MatrixSparse(ii m = 0, ii n = 0);
    ~MatrixSparse();

    void init(ii m = 0, ii n = 0);
    void free();
    void swap(MatrixSparse& a);

    // accessors
    ii m() const;
    ii n() const;
    li size() const;
    ii nnz() const;
    ii nnzActual() const;
    fp* vs() const;

    // these functions allocate memory
    void copy(const MatrixSparse& a, bool transpose = false);
    void copy(ii m, ii n, ii nnz, const ii* rowind, const ii* colind, const fp* acoo); // create from COO matrix
    void copy(const Matrix& a); // create from dense matrix a
    void copy(ii m, ii n, fp v); // create from dense matrix of constant value
    void copyConcatenate(const std::vector<MatrixSparse>& xs); // the xs must be row vectors
    void copySubset(const MatrixSparse& a); // only non-zero elements of this matrix are overwritten by corresponding elements in a
    void copySubset(const MatrixSparse& a, const MatrixSparse& b); // only non-zero elements of b are copied from a to this matrix
    ii copyPrune(const MatrixSparse &a, fp threshold); // prune values under threshold
    ii copyPruneRows(const MatrixSparse& a, const MatrixSparse& b, bool bRows, fp threshold); // prune rows of this matrix when rows or columns of a are empty

    // exports
    void exportTo(ii* rowind, ii* colind, fp* acoo) const; // export as COO matrix
    void exportTo(fp *vs) const; // export as dense matrix

    // elementwise operations
    void add(fp alpha, bool transposeA, const MatrixSparse& a, const MatrixSparse& b);
    void matmul(bool transposeA, const MatrixSparse& a, const MatrixSparse& b, bool accumulate, bool denseOutput = false);
    void mul(fp beta);
    void mul(const MatrixSparse& a);
    void sqr();
    void sqr(const MatrixSparse& a);
    void sqrt();
    void pow(fp power);

    // elementwise operations only operating on non-zero elements
    void addNonzeros(fp beta);
    void addNonzeros(const MatrixSparse& a);
    void lnNonzeros();
    void lnNonzeros(const MatrixSparse& a);
    void expNonzeros();
    void divNonzeros(const MatrixSparse& b); // a is denominator
    void divNonzeros(const MatrixSparse& a, const MatrixSparse& b); // a/b
    void div2Nonzeros(const MatrixSparse& a); // a is numerator
    void div2Nonzeros(const Matrix& a); // a is numerator & must be dense

    // aggregate operations
    fp sum() const;
    fp sumSqrs() const;
    fp sumSqrDiffsNonzeros(const MatrixSparse& a) const;

    static double sortElapsed_;

protected:
    void sort() const;

    ii m_; // number of rows
    ii n_; // number of columns
    
    ii* is0_; ii* is1_; ii* js_; fp* vs_; // pointers to CSR array
    sparse_matrix_t mat_; // opaque MKL sparse matrix object
    bool isSorted_; // true if we definately know the sparse matrix is sorted
    bool isOwned_; // true if data arrays owned by this object (false is owned by MKL or by a parent matrix)

    sparse_status_t status_; // last MKL function status

    friend MatrixSparseView;
    friend std::ostream& operator<<(std::ostream& os, const MatrixSparse& a);
};

std::ostream& operator<<(std::ostream& os, const MatrixSparse& a);


class MatrixSparseView : public MatrixSparse
{
public:
    MatrixSparseView(const MatrixSparse &a, ii row);
    ~MatrixSparseView();

private:
    bool isOwned_;
};


#endif


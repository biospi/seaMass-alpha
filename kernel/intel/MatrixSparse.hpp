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
#include "../SubjectMatrixSparse.hpp"
#include <vector>


class MatrixSparseView;


class MatrixSparse : public SubjectMatrixSparse
{
public:
    MatrixSparse();
    ~MatrixSparse();

    void clear();
    void init(ii m, ii n);
    void swap(MatrixSparse& a);

    // accessors
    ii m() const;
    ii n() const;
    li size() const;
    ii nnz() const;
    ii nnzActual() const;

public:
    // these functions allocate memory
    void importFromCoo(ii a_m, ii a_n, ii a_nnz, const ii *a_is, const ii *a_js, const fp *a_vs); // create from COO matrix
    void importFromMatrix(const Matrix &a); // create from dense matrix a
    void importFromMatrix(ii m, ii n, fp v); // create from dense matrix of constant value
    void concatenateSparseVectors(const std::vector<MatrixSparse> &xs); // the xs must be row vectors
    void copySubset(const MatrixSparse& a, const MatrixSparse& b); // only non-zero elements of b are copied from a to this matrix
    ii copyPrune(const MatrixSparse &a, fp threshold = 0.0); // prune values under threshold
    ii pruneRows(const MatrixSparse &a, const MatrixSparse &b, bool bRows, fp threshold); // prune rows of this matrix when rows or columns of a are empty

    void copy(const MatrixSparse& a);
    void transpose(const MatrixSparse& a);

    // exports
    void exportTo(ii* is, ii* js, fp* vs) const; // export as COO matrix
    void exportTo(fp *vs) const; // export as dense matrix

    // matrix multiplication
    void matmul(bool transposeA, const MatrixSparse& a, const MatrixSparse& b, bool accumulate);
    void matmulDense(bool transposeA, const MatrixSparse &a, const MatrixSparse &b, bool accumulate);

    // elementwise inplace operations
    void mul(fp beta);

    // elementwise copy operations
    void add(fp alpha, bool transposeA, const MatrixSparse& a, const MatrixSparse& b);

    // elementwise inplace and copy operations
    void mul(const MatrixSparse& a, const MatrixSparse& b);
    void sqr(const MatrixSparse& a);
    void sqrt(const MatrixSparse& a);
    void pow(const MatrixSparse& a, fp power);
    void censorLeft(const MatrixSparse& a, fp threshold);

    // elementwise operations only operating on non-zero elements
    void addNonzeros(fp beta);
    void addNonzeros(const MatrixSparse& a, const MatrixSparse& b);
    void lnNonzeros(const MatrixSparse& a);
    void divNonzeros(const MatrixSparse& b); // a is denominator
    void divNonzeros(const MatrixSparse& a, const MatrixSparse& b); // a/b
    void div2Nonzeros(const MatrixSparse& a); // a is numerator
    void div2(const Matrix &a); // a is numerator & must be dense

    // aggregate operations
    fp sum() const;
    fp sumSqrs() const;
    fp sumSqrDiffsNonzeros(const MatrixSparse& a) const;

    static double sortElapsed_;

protected:
    bool initCsr(ii m, ii n, ii a_nnz);
    void commitCsr(bool isSorted);
    void initMkl(ii m, ii n);
    void commitMkl(bool isSorted);
    void sort() const;

    ii m_; // number of rows
    ii n_; // number of columns

    //! Pointer to opaque MKL sparse matrix output if an MKL Inspector-executor Sparse BLAS routien has been run
    //!
    sparse_matrix_t mat_;

    //! Pointers to CSR array if we had to create of use matrix outside of MKL Inspector-executor Sparse BLAS routines
    //! If isOwned_, this object is responsible for freeing the CSR array
    //! If both mat_ and is1_ are 0, the matrix has no non-zeros
    //!
    ii* ijs_;
    ii* ijs1_;
    ii* js_;
    fp* vs_;
    bool isCsrOwned_; // true if data arrays owned by this object (false is owned by MKL or by a parent matrix)

    //! isSorted_ should be true if we definately know column indices are sorted within each row
    bool isSorted_;

    sparse_status_t status_; // last MKL function status

    friend ObserverMatrixSparse;
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


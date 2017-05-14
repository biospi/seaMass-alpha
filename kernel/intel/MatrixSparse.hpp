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
    void swap(MatrixSparse& a);

    // accessors
    ii m() const;
    ii n() const;
    li size() const;
    ii nnz() const;
    ii nnzActual() const;

    // these free the existing contents and create new contents
    void free();
    void init(ii m, ii n);
    void initWithSharedIndexes(MatrixSparse& x, const MatrixSparse& a);
    ii pruneCells(fp threshold = 0.0); // prune values under threshold

    // elementwise operators: inputs must have shared indexes
    void divNonzeros(MatrixSparse& a, MatrixSparse& b);

    // todo:

    void copy(const MatrixSparse& a);

    void transpose(const MatrixSparse& a);
    void importFromCoo(ii a_m, ii a_n, ii a_nnz, const ii *a_is, const ii *a_js, const fp *a_vs); // create from COO matrix
    void importFromMatrix(const Matrix &a); // create from dense matrix a
    void initDense(ii m, ii n, fp v); // create from dense matrix of constant value
    void concatenateRows(const std::vector<MatrixSparse> &xs); // the xs must be row vectors
    ii pruneRows(const MatrixSparse &a, const MatrixSparse &b, bool bRows, fp threshold); // prune rows of this matrix when rows or columns of a are empty

    // constructors with channel ops

    // exports
    void exportToCoo(ii *is, ii *js, fp *vs) const; // export as COO matrix
    void exportToDense(fp *vs) const; // export as dense matrix

    // matrix multiplication
    void matmul(bool transposeA, const MatrixSparse& a, const MatrixSparse& b, bool accumulate);
    void matmulDense(bool transposeA, const MatrixSparse &a, const MatrixSparse &b);

    // elementwise channel operations
    void mul(fp beta);
    //void divNonzeros(ii a, ii b, ii c = 0); // a/b -> c

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
    void div2Dense(const Matrix &a); // a is numerator & must be dense

    // aggregate operations
    fp sum() const;
    fp sumSqrs() const;
    fp sumSqrDiffsNonzeros(const MatrixSparse& a) const;

    static double sortElapsed_;

public:
    bool createCsr(ii m, ii n, ii a_nnz);
    bool createCsrWithSharedIndexes(MatrixSparse& x);
    bool initCsr(bool notEmpty) const;
    void commitCsr(bool isSorted) const;

    bool createMkl(ii m, ii n, bool notEmpty);
    bool initMkl(bool notEmpty) const;
    void commitMkl(bool isSorted) const;

    ii createChannel(ii nnz);

    void sort() const;

    //! number of rows
    //!
    ii m_;

    //! number of columns
    //!
    ii n_;

    //! number of nnz (this is a convenience for debugging only, it is also available in ijs_[m_])
    //!
    ii nnz_;

    //! CSR row indexes of size m_+1. They point to where in js_ the row starts.
    //!
    ii* ijs_;

    //! CSR column indexes of size nnz_
    //!
    ii* js_;

    //! CSR cell values of size nnz_.
    fp* vs_;

    //! True if ijs_ and js_ were malloc'ed by us, False if allocated by _mat or shared with another MatrixSparse.
    //
    bool isCsrIndexesOwned_;

    //! True if vs_ was malloc'ed by us, False if allocated by _mat.
    //
    bool isCsrValuesOwned_;

    //! Pointer to opaque MKL sparse matrix output if an MKL Inspector-executor Sparse BLAS routine created this matrix.
    //!
    sparse_matrix_t mat_;

    //! isSorted_ should be true only if we definitely know column indices are sorted within each row
    //!
    bool isSorted_;

    //! pointer to MatrixSparse that holds the ijs_ and js_
    //!
    MatrixSparse* parent_;

    //! pointers to MatrixSparse's that share this ijs_ and js_
    //!
    std::vector<MatrixSparse*> children_;

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


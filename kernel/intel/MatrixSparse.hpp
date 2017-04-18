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


class MatrixSparse
{
public:
	MatrixSparse();
	~MatrixSparse();

    // accessors
    ii m() const;
    ii n() const;
    li size() const;
    ii nnz() const;
    ii nnzActual() const;
    ii* js() const;
    fp* vs() const;

	// shallow inits
    void init(ii m = 0, ii n = 0);
	void init(const MatrixSparse &a, ii row);

    // deep copies
    void copy(const MatrixSparse& a, bool transpose = false);
    void copy(const std::vector<MatrixSparse> &xs); // the xs must be row vectors
    void copy(ii m, ii n, const std::vector<ii> &is, const std::vector<ii> &js, const std::vector<fp> &vs); // create from COO matrix
    void copy(ii m, ii n, fp v); // create from dense matrix of constant value
    void copy(const Matrix &a); // create from dense matrix a
    void copySubset(const MatrixSparse &a); // only non-zero elements of this matrix are overwritten by corresponding elements in a

    // exports
    void exportTo(std::vector<ii> &is, std::vector<ii> &js, std::vector<fp> &vs) const; // export as COO matrix
    void exportTo(fp *vs) const; // export as dense matrix

    // management
    void clear();
    void sort();
    ii prune(const MatrixSparse& a, fp pruneThreshold);
    ii pruneRows(const MatrixSparse& a, ii aNnzRows, const MatrixSparse& b, bool bRows, fp threshold);

    // elementwise operations
    void add(fp alpha, bool transposeA, const MatrixSparse& a, const MatrixSparse& b);
	void matmul(bool transposeA, const MatrixSparse& a, const MatrixSparse& b, bool accumulate, bool denseOutput = false);
    void mul(fp beta);
    void mul(const fp* a_vs);
    void sqr();
    void sqrt();
    void pow(fp power);

    // elementwise operations only operating on non-zero elements
	void addNonzeros(fp beta);
	void addNonzeros(const fp* a_vs);
	void lnNonzeros();
    void expNonzeros();
    void divNonzeros(const fp* a_vs); // a_vs is denominator
    void div2Nonzeros(const fp* a_vs); // a_vs is numerator

    // aggregate operations
	fp sum() const;
	fp sumSqrs() const;
    fp sumSqrDiffsNonzeros(const fp *a_vs) const;

private:
	ii m_; // number of rows
	ii n_; // number of columns
    
    ii* is0_; ii* is1_; ii* js_; fp* vs_; // pointers to CSR array
    sparse_matrix_t mat_; // opaque MKL sparse matrix object
    bool isOwned_; // true if data arrays owned by this object (false is owned by MKL or by a parent matrix)

    sparse_status_t status_; // last MKL function status
    
    friend std::ostream& operator<<(std::ostream& os, const MatrixSparse& a);
};

std::ostream& operator<<(std::ostream& os, const MatrixSparse& a);


#endif


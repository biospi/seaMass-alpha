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


#include "MatrixMKL.hpp"
#include <vector>


class MatrixSparseMKL
{
public:
	MatrixSparseMKL();
	~MatrixSparseMKL();
    
    // deep copies
    void init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind); // create from COO matrix
    void init(ii m, ii n, fp v);
    void init(const std::vector<MatrixSparseMKL>& as); // stack these rows
    void init(const MatrixMKL& a);
    
    // shallow copy
    void wrap(const MatrixSparseMKL& a, ii row);
    
    void clear();
    void free();
    
    ii m() const;
    ii n() const;
    li size() const;
    ii nnz() const;
    ii nnzActual() const;
    ii* js() const;
    fp* vs() const;
  
    void sort();
    ii prune(const MatrixSparseMKL& a, fp pruneThreshold);
    ii pruneRows(const MatrixSparseMKL& a, ii aNnzRows, const MatrixSparseMKL& b, bool bRows, fp threshold);
    void output(fp* vs) const;
    void write(const std::string& filename) const;

    void copy(const MatrixSparseMKL& a, bool transpose = false);
    void add(fp alpha, bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b);
	void matmul(bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b, bool accumulate, bool denseOutput = false);
    void mul(fp beta);
    void mul(const fp* a_vs);
    void sqr();
    void sqrt();
    void pow(fp power);

    // operate only on non-zero elements
    void setNonzeros(fp v);
	void addNonzeros(fp beta);
	void addNonzeros(const fp* a_vs);
	void lnNonzeros();
    void expNonzeros();
    void divNonzeros(const fp* a_vs);
    void div2Nonzeros(const fp* a_vs);

    // aggregate operations
	fp sum() const;
	fp sumSqrs() const;
    fp sumSqrDiffsNonzeros(const fp *a_vs) const;
    
    // note: these elementwise operations ONLY considers the non-zero elements of THIS matrix
    void subsetCopy(const MatrixSparseMKL& a);

private:
	ii m_;   // number of rows
	ii n_;   // number of columns
    
    ii* is0_; ii* is1_; ii* js_; fp* vs_; // pointers to CSR array
    sparse_matrix_t mat_;
    bool isOwned_; // true if data arrays owned by this object (false is owned by MKL or by a parent matrix)

    sparse_status_t status_; // last MKL function status
    
    class MyComparator;
    
    friend std::ostream& operator<<(std::ostream& os, const MatrixSparseMKL& mat);
};

std::ostream& operator<<(std::ostream& os, const MatrixSparseMKL& mat);


#endif


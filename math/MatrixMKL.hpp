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


#ifndef _SEAMASS_MATH_MATRIXMKL_HPP_
#define _SEAMASS_MATH_MATRIXMKL_HPP_


#include "MatrixSparseMKL.hpp"


/*class MatrixMKL
{
public:
	MatrixMKL();
	~MatrixMKL();
    
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
    void copy(const MatrixSparseMKL& a, Operation operation = Operation::NONE);
    void copy(const MatrixSparseMKL& a, fp pruneThreshold);
    void output(fp* vs) const;
    
    void deleteRows(const MatrixSparseMKL& a);
    
    // generalised sparse matrix multiplication
    void mul(bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b, bool accumulate, bool transpose);
    
    // elementwise operations
    void elementwiseAdd(fp beta);
    void elementwiseMul(fp beta);
    void elementwiseMul(const MatrixSparseMKL& a);
    void elementwiseDiv(const MatrixSparseMKL& a);
    void elementwiseSqr();
    void elementwiseSqrt();
    void elementwiseLn();
    void elementwisePow(fp power);
    
    // aggregate operations
    fp sum() const;
    fp sumSqrs() const;
    fp sumSqrDiffs(const MatrixSparseMKL& a) const;
    
    // note: these elementwise operations ONLY considers the non-zero elements of THIS matrix
    void subsetElementwiseCopy(const MatrixSparseMKL& a);
    void subsetElementwiseDiv(const MatrixSparseMKL& a); // needed for error calculation atm
    
 
private:
	li m_; // rows
	ii n_; // columns
	fp* vs_; // data
	bool isOwned_; // does this object own its data?
};

std::ostream& operator<<(std::ostream& os, const MatrixMKL& mat);
*/

#endif


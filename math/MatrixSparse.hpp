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

#include <mkl.h>
#include <mkl_spblas.h>


typedef float fp; // fp is the selected floating point precision (float or double)
typedef MKL_INT ii; // ii is the selected indexing integer size (32 or 64 bit)
typedef long long li;

class MatrixSparse
{
public:
	MatrixSparse();
	~MatrixSparse();

	void init(const MatrixSparse& a);
	void init(const MatrixSparse& a, fp pruneThreshold);
	void init(ii m, ii n); // create empty csr
	void init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind);
	void init(ii m, ii n, fp v);
	void convertFromDense(ii m, ii n, const fp* vs);
	void convertToDense(fp* vs);
	void free();

	ii m() const;
	ii n() const;
	li size() const;

	bool operator!() const;
	ii nnz(bool actual = false) const;
	li mem() const;

	enum class Transpose { NO, YES };
	enum class Accumulate { NO, YES };
	void mul(const MatrixSparse& a, const Transpose transpose, const MatrixSparse& x, const Accumulate accumulate);

	void elementwiseAdd(fp beta);
	void elementwiseMul(fp beta);
	void elementwiseMul(const MatrixSparse& a);
	void elementwiseDiv(const MatrixSparse& a);
	void elementwiseSqr();
	void elementwiseSqrt();
	void elementwiseLn();
	void elementwisePow(fp power);

	double sum() const;
	double sumSqrs() const;
	double sumSqrDiffs(const MatrixSparse& a) const;

	const fp* getVs() const;

private:
	ii m_;
	ii n_;
	ii* is0_;
	ii* is1_;
	ii* js_;
	fp* vs_;
	bool isOwned_;

	sparse_matrix_t mat_;

	friend class Matrix;
};

std::ostream& operator<<(std::ostream& os, const MatrixSparse& mat);


#endif


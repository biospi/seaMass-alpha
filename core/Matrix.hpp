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


#ifndef _SEAMASS_CORE_MATRIX_HPP_
#define _SEAMASS_CORE_MATRIX_HPP_


#include "MatrixSparse.hpp"


class Matrix
{
public:
	Matrix();
	~Matrix();

	void init(li m, ii n); // create matrix with storage but do not initialise
	void init(li m, ii n, fp v); // create matrix with data all set to 'v'  
	void init(li m, ii n, ii s, fp* vs); // create matrix using external data 'vs'
	void init(const Matrix& a, li i, ii j, li m, ii n, ii s = 0); // define submatrix of 'a' 
	void free();

	li m() const;
	ii n() const;
	li nnz() const;
	li size() const;

	void copy(const Matrix& a);
	void mul(const MatrixSparse& a, const Matrix& x, bool accumulate, bool transposeA);
	void elementwiseMul(const Matrix& a, const Matrix& b);
	void elementwiseMul(fp scale, const Matrix& a);
	void elementwiseDiv(const Matrix& n, const Matrix& d);
	void elementwiseSqrt(const Matrix& a);
	void shrinkage(const Matrix& cE, const Matrix& c0, const Matrix& l1, fp lambda);

	void prune(fp threshold);

	double sum() const;
	double sumSqrs() const;
	double sumSqrDiffs(const Matrix& a) const;

	bool operator!() const;
	li mem() const;

	const fp* getVs() const;

private:
	li m_; // rows
	ii n_; // columns
	ii s_; // stride
	fp* vs_; // data
	bool isOwned_; // does this object own its data?
};

std::ostream& operator<<(std::ostream& os, const Matrix& mat);


#endif


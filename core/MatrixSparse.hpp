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


#ifndef _SEAMASS_CORE_MATRIXSPARSE_HPP_
#define _SEAMASS_CORE_MATRIXSPARSE_HPP_


#include <iostream>
#include <vector>

#include <mkl.h>


typedef float fp; // fp is the selected floating point precision (float or double)
typedef MKL_INT ii; // ii is the selected indexing integer size (32 or 64 bit)
typedef long long li;


class MatrixSparse
{
public:
	MatrixSparse();
	//MatrixSparse(const MatrixSparse& a);
	~MatrixSparse();

	void init(ii m, ii n, ii nnz);
	void init(const MatrixSparse& a);
	void init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind);
	void free();

	ii m() const;
	ii n() const;
	ii nnz() const;
	li size() const;

	bool operator!() const;
	li mem() const;

	void elementwiseSqr(const MatrixSparse& a);

public:
	ii m_;
	ii n_;
	ii nnz_;
	ii* is_;
	ii* js_;
	fp* vs_;
	bool isIsJsOwned_;
	bool isVsOwned_;
};

std::ostream& operator<<(std::ostream& os, const MatrixSparse& mat);


#endif


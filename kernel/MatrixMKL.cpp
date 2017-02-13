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


#include "MatrixMKL.hpp"

#include "MatrixSparseMKL.hpp"

#include <iomanip>
#include <iostream>
#include <cmath>

using namespace std;


MatrixMKL::MatrixMKL()
	: m_(0), n_(0), vs_(0)
{
}


MatrixMKL::~MatrixMKL()
{
	free();
}


void MatrixMKL::init(ii m, ii n, const fp* vs)
{
	free();

	m_ = m;
	n_ = n;
    vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * m_ * n_, 64));
    memcpy(vs_, vs, sizeof(fp) * m_ * n_);
}


void MatrixMKL::free()
{
	m_ = 0;
	n_ = 0;
    mkl_free(vs_);
}


li MatrixMKL:: size() const
{
	return (li)m_ * n_;
}


ii MatrixMKL::m() const
{
	return m_;
}


ii MatrixMKL::n() const
{
	return n_;
}


fp* MatrixMKL::vs() const
{
	return vs_;
}


ostream& operator<<(ostream& os, const MatrixMKL& a)
{
	if (a.m() == 0)
	{
		os << "[]";
	}
	else
	{
		os << "[" << a.m() << "," << a.n() << "]";
	}

	return  os;
}

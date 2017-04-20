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


#include "Matrix.hpp"
#include <cstring>
using namespace std;


Matrix::Matrix()
    : m_(0), n_(0), vs_(0)
{
}


Matrix::~Matrix()
{
    free();
}


void Matrix::alloc(ii m, ii n)
{
    free();

    m_ = m;
    n_ = n;
    vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * m_ * n_, 64));
}


void Matrix::copy(ii m, ii n, const fp *vs)
{
    alloc(m, n);
    memcpy(vs_, vs, sizeof(fp) * m_ * n_);
}


void Matrix::free()
{
    m_ = 0;
    n_ = 0;

    if (vs_)
    {
        mkl_free(vs_);
        vs_ = 0;
    }
}


li Matrix:: size() const
{
    return (li)m_ * n_;
}


ii Matrix::m() const
{
    return m_;
}


ii Matrix::n() const
{
    return n_;
}


fp* Matrix::vs() const
{
    return vs_;
}


ostream& operator<<(ostream& os, const Matrix& a)
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

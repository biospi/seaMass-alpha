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
#include "kernel.hpp"
#include <iomanip>
#include <sstream>
#include <ippcore.h>
#include <ipps.h>
using namespace std;
using namespace kernel;


Matrix::Matrix()
    : m_(0), n_(0), vs_(0)
{
}


Matrix::~Matrix()
{
    free();
}


void Matrix::init(ii m, ii n)
{
    free();

    m_ = m;
    n_ = n;
    vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * m_ * n_, 64));
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


void Matrix::importFromArray(ii m, ii n, const fp *vs)
{
    init(m, n);

    for (ii i = 0; i < m_; i++)
        ippsCopy_32f(&vs[i * n_], &vs_[i * n_], n_);
}


void Matrix::exportToArray(fp *vs) const
{
    for (ii i = 0; i < m_; i++)
        ippsCopy_32f(&vs_[i * n_], &vs[i * n_], n_);
}


fp Matrix::sum() const
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       sum(X" << *this << ") := ...";
        info(oss.str());
    }


    fp sum = 0.0;
    for (ii i = 0; i < m_; i++)
    {
        fp sumRow;
        ippsSum_32f(&vs_[i * n_], n_, &sumRow, ippAlgHintFast);
        sum += sumRow;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << sum;
        info(oss.str(), this);
    }

    return sum;
}


li Matrix:: size() const
{
    return li(m_) * n_;
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

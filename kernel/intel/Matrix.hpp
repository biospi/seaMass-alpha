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


#ifndef SEAMASS_KERNEL_INTEL_MATRIX_HPP
#define SEAMASS_KERNEL_INTEL_MATRIX_HPP


#include "types.hpp"
#include "../SubjectMatrix.hpp"
#include <iostream>


class Matrix : public SubjectMatrix
{
public:
    Matrix();
    ~Matrix();

    void init(ii m, ii n);
    void free();

    ii m() const;
    ii n() const;
    li size() const;
    fp* vs() const;

    void copy(ii m, ii n, const fp *vs);
    fp sum() const;

private:
    ii m_; // rows
    ii n_; // columns
    fp* vs_; // data
};

std::ostream& operator<<(std::ostream& os, const Matrix& mat);


#endif


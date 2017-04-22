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


#ifndef SEAMASS_ASRL_BASISMATRIX_HPP
#define SEAMASS_ASRL_BASISMATRIX_HPP


#include "Basis.hpp"


class BasisMatrix : public Basis
{
public:
    BasisMatrix(std::vector<Basis*>& bases, std::vector<MatrixSparse>& aT, std::vector<MatrixSparse>* gT, bool transient);
    virtual ~BasisMatrix();

    virtual void synthesise(std::vector<MatrixSparse> &f, const std::vector<MatrixSparse> &x, bool accumulate);
    virtual void analyse(std::vector<MatrixSparse> &xE, const std::vector<MatrixSparse> &fE, bool sqrA = false) const;

    virtual const std::vector<MatrixSparse>* getGroups(bool transpose) const;

private:
    std::vector<MatrixSparse>& aTs_;
    std::vector<ii> aTnnzRows_;
    std::vector<MatrixSparse> as_;

    const std::vector<MatrixSparse>* gT_;
    std::vector<MatrixSparse>* g_;
};


#endif


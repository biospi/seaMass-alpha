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


#ifndef SEAMASS_ASRL_BASIS_HPP
#define SEAMASS_ASRL_BASIS_HPP


#include "../kernel/Subject.hpp"
#include <MatrixSparse.hpp>
#include <vector>


class Basis : public Subject
{
public:
    Basis(std::vector<Basis*>& bases, bool transient, int parentIndex = -1);
    virtual ~Basis();

    virtual void synthesize(std::vector<MatrixSparse> &f, const std::vector<MatrixSparse> &x, bool accumulate) = 0;
    virtual void analyze(std::vector<MatrixSparse> &xE, const std::vector<MatrixSparse> &fE, bool sqrA) = 0;

    virtual void synthesizeGroups(std::vector<MatrixSparse> &g, const std::vector<MatrixSparse> &x, bool accumulate);
    virtual const std::vector<MatrixSparse> * getColGroups(bool transpose) const;

    int getIndex() const;
    int getParentIndex() const;
    bool isTransient() const;

private:
    int index_;       // index of this basis in the serialised tree
    int parentIndex_; // parent node
    bool transient_; // if transient, coefficients not part of fitting
};


#endif


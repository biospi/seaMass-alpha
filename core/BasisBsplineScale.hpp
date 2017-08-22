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


#ifndef SEAMASS_CORE_BASISBSPLINESCALE_HPP
#define SEAMASS_CORE_BASISBSPLINESCALE_HPP


#include "BasisBspline.hpp"


class BasisBsplineScale : public BasisBspline
{ 
public:
    BasisBsplineScale(std::vector<Basis*>& bases, int parentIndex, short dimension0, short dimension1,
                      bool transient);
    virtual ~BasisBsplineScale();

    virtual void synthesize(std::vector<MatrixSparse> &f, const std::vector<MatrixSparse> &x, bool accumulate);
    virtual void analyze(std::vector<MatrixSparse> &xE, const std::vector<MatrixSparse> &fE, bool sqrA = false);

    virtual const std::vector<MatrixSparse> * getColGroups(bool transpose) const;

private:
    MatrixSparse aT_;
    MatrixSparse a_;

    short dimension0_;
    short dimension1_;

    std::vector<MatrixSparse> gTs_;
    std::vector<MatrixSparse> gs_;
};


#endif


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


#ifndef SEAMASS_ASRL_OPTIMIZER_HPP
#define SEAMASS_ASRL_OPTIMIZER_HPP


#include "Basis.hpp"


class Optimizer
{
public:    
    Optimizer();
    virtual ~Optimizer();
    
    virtual void setLambda(fp lambda, fp lambdaGroup = fp(0.0)) = 0;
    virtual const std::vector<Basis*>& getBases() const = 0;
    virtual ii getIteration() const = 0;

    virtual fp step() = 0;

    virtual void synthesize(std::vector<MatrixSparse>& f, std::vector< std::vector<MatrixSparse> >& xEs, ii basis = -1) = 0;
    virtual void analyze(std::vector< std::vector<MatrixSparse> > &xEs, std::vector<MatrixSparse> &fE, bool l2, bool l2Normalize = true) const = 0;

    virtual std::vector< std::vector<MatrixSparse> >& xs() = 0;
    virtual std::vector< std::vector<MatrixSparse> >& l2s() = 0;
    virtual std::vector< std::vector<MatrixSparse> >& l1l2s() = 0;
};


#endif


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


#ifndef SEAMASS_ASRL_ASRL_HPP
#define SEAMASS_ASRL_ASRL_HPP


#include "Optimizer.hpp"


/**
* SeamassAsrl performs Accelerated Sparse Richardson Lucy optimisation on the input.
*/
class Asrl
{
public:
    static void notice();

    struct Input
    {
        std::vector<MatrixSparse> aT; // transpose of A from Ax = b
        std::vector<MatrixSparse> xT; // transpose of x from Ax = b (leave x empty to autogenerate seed)
        std::vector<Matrix> bT;       // transpose of b from Ax = b
        std::vector<MatrixSparse> gT; // transpose of Group indicator matrix
    };

    struct Output
    {
        std::vector<MatrixSparse> xT;  // from Ax = b
        std::vector<Matrix> xTaT;       // transpose of Ax
        std::vector<MatrixSparse> xTgT; // transpose of Gx
    };

    Asrl(Input &input, fp lambda, fp lambdaGroup, bool taperShrinkage, fp tolerance);
    virtual ~Asrl();

    bool step();
    ii getIteration() const;

    void getOutput(Output& output) const;

private:
    std::vector<Basis*> bases_;
    const std::vector<Matrix>& bT_;

    Optimizer* innerOptimizer_;
    Optimizer* optimizer_;

    fp lambda_;
    fp lambdaGroup_;
    fp lambdaGroupStart_;
    bool taperShrinkage_;
    fp tolerance_;
    int iteration_;
};


#endif


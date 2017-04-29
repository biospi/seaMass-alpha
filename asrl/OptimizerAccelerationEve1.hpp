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


#ifndef SEAMASS_ASRL_OPTIMIZERACCELERATIONEVE1_HPP
#define SEAMASS_ASRL_OPTIMIZERACCELERATIONEVE1_HPP


#include "Optimizer.hpp"


/**
* OptimizerAccelerationEve1 accelerates a wrapped Optimizer with 1st Order Exponential Vector Extrapolation (EVE1)
*
* Notes:
*  i) This class implements Section 3.6.8 of the PhD Thesis of David DC Biggs (1998) available from
*   https://researchspace.auckland.ac.nz/handle/2292/1760 . Note that this is an exponentiated version of the algorithm
*   presented in [Biggs & Andrews, Applied Optics 1997] which does not need a Poisson non-negativety constraint when used
*   with Richardson-Lucy optimization. It is interesting to note that Biggs does not present results of the
*   non-exponentiated version in Richardson-Lucy in his thesis, despite it being the topic of his Applied Optics paper. 
*  ii) [Wang & Miller, IEEE Trans Imag Proc 2014] provide a Scaled Heavy-Ball method with a convergence rate proof and
    show Vector Extrapolation is a special case. However, again this is the non-exponentiated version.
*/
class OptimizerAccelerationEve1 : public Optimizer
{
public:    
    OptimizerAccelerationEve1(Optimizer* optimizer);
    virtual ~OptimizerAccelerationEve1();
    
    virtual void setLambda(fp lambda, fp lambdaGroup = fp(0.0));
    virtual const std::vector<Basis*>& getBases() const;
    virtual ii getIteration() const;

    virtual fp step();

    virtual void synthesize(std::vector<MatrixSparse>& f, std::vector< std::vector<MatrixSparse> >& xEs, ii basis = -1) const;
    virtual void analyze(std::vector< std::vector<MatrixSparse> > &xEs, std::vector<MatrixSparse> &fE, bool l2, bool l2Normalize = true) const;

    std::vector< std::vector<MatrixSparse> >& xs();
    std::vector< std::vector<MatrixSparse> >& l2s();
    std::vector< std::vector<MatrixSparse> >& l1l2s();

private:
    Optimizer* optimizer_;

    std::vector< std::vector<MatrixSparse> > x0s_;
    std::vector< std::vector<MatrixSparse> > y0s_;
    std::vector< std::vector<MatrixSparse> > u0s_;
    
    double accelerationDuration_;
};


#endif


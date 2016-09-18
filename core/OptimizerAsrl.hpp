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


#ifndef _SEAMASS_CORE_OPTIMIZERASRL_HPP_
#define _SEAMASS_CORE_OPTIMIZERASRL_HPP_


#include "Basis.hpp"


class OptimizerAsrl
{
public:    
	OptimizerAsrl(const std::vector<Basis*>& bases, const Matrix& g, ii accelleration = 2);
	virtual ~OptimizerAsrl();
    
	double step(fp lambda);

    //void threshold(fp threshold);
	//std::vector<Matrix>& getCoeffs() { return cs_; }
	//const std::vector<Matrix>& getL1Norms() const { return l1s_; }
	//const std::vector<Matrix>& getL2Norms() const { return l2s_; }

	void synthesis(Matrix& f, ii basis = -1) const;

protected:
	void error(Matrix& f) const;
	void analysis(std::vector<Matrix>& cE, const Matrix& fE) const;
	void shrinkage(std::vector<Matrix>& cE, fp lambda);
	double acceleration(std::vector<Matrix>& cE, fp lambda);

private:
	const std::vector<Basis*>& bases_;
	const Matrix& g_;

	std::vector<Matrix> cs_;
	std::vector<Matrix> l1s_;
	std::vector<Matrix> l2s_;

	ii accelleration_;
	std::vector<Matrix> c0s_;
	std::vector<Matrix> u0s_;
	std::vector<Matrix> q0s_;

	ii iteration_;
};


#endif


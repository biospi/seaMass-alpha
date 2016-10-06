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


#ifndef _SEAMASS_CORE_OPTIMIZERSRL_HPP_
#define _SEAMASS_CORE_OPTIMIZERSRL_HPP_


#include "Basis.hpp"


class OptimizerSrl
{
public:    
	OptimizerSrl(const std::vector<Basis*>& bases, const Matrix& g, fp pruneThreshold);
	virtual ~OptimizerSrl();
    
	void init(fp lamba);
	double step();

	void synthesis(Matrix& f, ii basis = -1) const;
	const std::vector<Matrix>& getCs() const;

protected:
	virtual void update(std::vector<Matrix>& c, std::vector<Matrix>& c1, ii iteration);

private:
	const std::vector<Basis*>& bases_;
	const Matrix& g_;
	fp pruneThreshold_;

	std::vector<Matrix> cs_;
	std::vector<Matrix> l2s_;
	std::vector<Matrix> l1l2s_;

	fp lambda_;
	ii iteration_;
};


#endif


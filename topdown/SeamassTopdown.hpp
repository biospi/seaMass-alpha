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


#ifndef _SEAMASS_CORE_SEAMASSTOPDOWN_HPP_
#define _SEAMASS_CORE_SEAMASSTOPDOWN_HPP_


#include "../core/Basis.hpp"
#include "../core/Optimizer.hpp"


/**
* SeamassTopdown deconvolution of the input spectrum.
*/
class SeamassTopdown
{
public:
	static void notice();

	struct Input {
		std::vector<fp> binCounts;
		ii scale;
		ii offset;
	};

	struct Output {
		std::vector<fp> weights;
		std::vector< std::vector<ii> > scales;
		std::vector< std::vector<li> > offsets;
		std::vector<ii> baselineScale;
		std::vector<ii> baselineOffset;
		std::vector<ii> baselineExtent;
	};

	SeamassTopdown(Input& input, ii maxMass, ii binsPerDalton, double shrinkage, double tolerance, ii debugLevel = 0);
	virtual ~SeamassTopdown();

	bool step();
	ii getIteration() const;

private:
	void init(Input& input, ii maxMass, ii binsPerDalton);

	MatrixSparse b_;
	ii dimensions_;
	std::vector<Basis*> bases_;
	Optimizer* inner_optimizer_;
	Optimizer* optimizer_;
	double shrinkage_;
	double tolerance_;
	ii iteration_;
	ii debugLevel_;
};


#endif // _SEAMASS_HPP_


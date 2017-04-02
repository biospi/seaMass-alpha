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


#ifndef SEAMASS_ASRL_SEAMASSASRL_HPP
#define SEAMASS_ASRL_SEAMASSASRL_HPP


#include "Optimizer.hpp"


/**
* SeamassAsrl performs Accelerated Sparse Richardson Lucy optimisation on the input.
*/
class SeamassAsrl
{
public:
	static void notice();

	struct Input {
		// Ax = b
		ii aM, aN;
		std::vector<fp> aVs;
		std::vector<ii> aIs;
		std::vector<ii> aJs;

		std::vector<fp> xs; // if not zero size, used as starting estimate
		std::vector<fp> bs;

		// GROUPS
		ii gM, gN; // if gM = 0 then we do not use group sparse inference
		std::vector<fp> gVs;
		std::vector<ii> gIs;
		std::vector<ii> gJs;
	};

	struct Output {
		std::vector<fp> xs;
		std::vector<fp> ys; // Gx = y
	};

	SeamassAsrl(Input& input, double shrinkage, double tolerance);
	virtual ~SeamassAsrl();

	bool step();
	ii getIteration() const;

	void getOutput(Output& output) const;

private:
	std::vector<Basis*> bases_;
	Optimizer* innerOptimizer_;
	Optimizer* optimizer_;
	double shrinkage_;
	double tolerance_;
	int iteration_;
};


#endif


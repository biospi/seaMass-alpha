//
// $Id$
//
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


#ifndef _SEAMASS_HPP_
#define _SEAMASS_HPP_


#include <vector>
#include "core.hpp"
#include "BasisFunctions.hpp"
#include "OptimiserASRL.hpp"


/**
* seaMass fitting of a 1-d curve or multi-dimensional surface to the input spectr(um|a).
*/
class seaMass
{
public:
	static void notice();

	struct Input {
		std::vector<fp> bin_counts;
		std::vector<li> spectrum_index;
		std::vector<double> bin_edges;
		std::vector<double> start_times;
		std::vector<double> finish_times;
		std::vector<fp> exposures;
	};

	struct Output {
		std::vector<fp> weights;
		std::vector< std::vector<ii> > scales;
		std::vector< std::vector<li> > offsets;
	};

	struct ControlPoints {
		std::vector<fp> coeffs;
		std::vector<ii> size;
		std::vector<ii> scale;
		std::vector<ii> offset;
	};

	struct ControlPoints1D {
		std::vector< std::vector<fp> > coeffs;
		ii scale;
		std::vector<ii> offsets;
	};

	seaMass(Input& input, const std::vector<ii>& scales);
	seaMass(Input& input, const Output& seed);

	virtual bool step(double shrinkage, double tolerance);
	void get_output(Output& output) const;

	void get_output_residuals(std::vector<fp>& residuals) const;
	void get_output_control_points(ControlPoints& control_points) const;
	void get_output_control_points_1d(ControlPoints1D& control_points) const;

private:
	void init(Input& input, const std::vector<ii>& scales);

	short dimensions;
	std::vector<Basis*> bases;
	OptimiserASRL* optimiser;
};


#endif // _SEAMASS_HPP_


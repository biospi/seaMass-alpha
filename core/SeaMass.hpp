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


#ifndef _SEAMASS_CORE_SEAMASS_HPP_
#define _SEAMASS_CORE_SEAMASS_HPP_


#include "Basis.hpp"
#include "OptimizerSrl.hpp"


//void remove_zeros(std::vector< std::vector<fp> >& mzs, std::vector< std::vector<fp> >& intensities);
//void merge_bins(std::vector< std::vector<fp> >& mzs, std::vector< std::vector<fp> >& intensities, double width);


/**
* SeaMass fitting of a 1-d curve or multi-dimensional surface to the input spectr(um|a).
*/
class SeaMass
{
public:
	static void notice();

	struct Input {
		std::vector<fp> binCounts;
		std::vector<li> spectrumIndex;
		std::vector<double> binEdges;
		std::vector<double> startTimes;
		std::vector<double> finishTimes;
		std::vector<fp> exposures;
	};

	struct Output {
		std::vector<fp> weights;
		std::vector< std::vector<ii> > scales;
		std::vector< std::vector<li> > offsets;
		std::vector<ii> baselineScale;
		std::vector<ii> baselineOffset;
		std::vector<ii> baselineExtent;
	};

	struct ControlPoints {
		std::vector<fp> coeffs;
		std::vector<ii> scale;
		std::vector<ii> offset;
		std::vector<ii> extent;
	};

	struct ControlPoints1D {
		std::vector< std::vector<fp> > coeffs;
		ii scale;
		std::vector<ii> offsets;
	};

	SeaMass(Input& input, const std::vector<ii>& scales, double shrinkage, double tolerance, ii debugLevel = 0);
	SeaMass(Input& input, const Output& seed, ii debugLevel = 0);
	virtual ~SeaMass();

	bool step();
	ii getIteration() const;

	void getOutput(Output& output) const;
	void getOutputBinCounts(std::vector<fp>& binCounts) const;
	void getOutputControlPoints(ControlPoints& controlPoints) const;
	void getOutputControlPoints1d(ControlPoints1D& controlPoints) const;

private:
	void init(Input& input, const std::vector<ii>& scales);

	Matrix b_;
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


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


#ifndef _SEAMASS_CORE_SEAMASSCORE_HPP_
#define _SEAMASS_CORE_SEAMASSCORE_HPP_


#include "../asrl/Basis.hpp"
#include "../asrl/OptimizerSrl.hpp"


//void remove_zeros(std::vector< std::vector<fp> >& mzs, std::vector< std::vector<fp> >& intensities);
//void merge_bins(std::vector< std::vector<fp> >& mzs, std::vector< std::vector<fp> >& intensities, double width);


/**
* SeamassCore fitting of a 1-d curve or multi-dimensional surface to the input spectr(um|a).
*/
class SeamassCore
{
public:
	static void notice();

	struct Input {
		std::vector<fp> binCounts;
		std::vector<li> binCountsIndex;
		std::vector<double> binEdges;
		std::vector<double> startTimes;
		std::vector<double> finishTimes;
		std::vector<fp> exposures;
	};

	struct Output {
		std::vector<fp> weights;
		std::vector< std::vector<short> > scales;
		std::vector< std::vector<li> > offsets;
		std::vector<short> baselineScale;
		std::vector<ii> baselineOffset;
		std::vector<ii> baselineExtent;
	};

	struct ControlPoints {
		std::vector<fp> coeffs;
		std::vector<short> scale;
		std::vector<ii> offset;
		std::vector<ii> extent;
	};

	SeamassCore(Input& input, const std::vector<short>& scales, double shrinkage, double tolerance);
	SeamassCore(Input& input, const Output& seed);
	virtual ~SeamassCore();

	bool step();
	ii getIteration() const;

	// get seaMass output (for smv file)
	void getOutput(Output& output) const;

	// get restored bin counts derived from seaMass output
	void getOutputBinCounts(std::vector<fp>& binCounts) const;

	// get restored 1D control points (i.e. per spectra) derived from seaMass output
	void getOutputControlPoints1d(ControlPoints& controlPoints) const;

	// get restored control points with dimension depending on input (i.e. 1D or 2D)
	void getOutputControlPoints(ControlPoints& controlPoints) const;

private:
	void init(Input& input, const std::vector<short>& scales);

	short dimensions_;
	std::vector<Basis*> bases_;
	Optimizer* inner_optimizer_;
	Optimizer* optimizer_;
	double shrinkage_;
	double tolerance_;
	int iteration_;
};


#endif // _SEAMASS_HPP_


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


#include "SeamassCore.hpp"

#include "../asrl/OptimizerAccelerationEve1.hpp"

#include "BasisBsplineMz.hpp"
#include "BasisBsplineScale.hpp"
#include "BasisBsplineScantime.hpp"

#include <iostream>
#include <iomanip>


using namespace std;


SeamassCore::SeamassCore(Input& input, const std::vector<short>& scales, double shrinkage, double tolerance) : shrinkage_(shrinkage), tolerance_(tolerance), iteration_(0)
{
	init(input, scales);
}


SeamassCore::SeamassCore(Input& input, const Output& seed) : iteration_(0)
{
	vector<short> scales(seed.scales.size());
	init(input, scales);

    throw runtime_error("not implemented yet");

    // todo
	/*cout << seed.weights.size() << endl;
	for (ii i = 0; i < seed.weights.size(); i++)
	{
		//cout << seed.scales[0][i] << "," << seed.scales[1][i] << ":" << seed.offsets[0][i] << "," << seed.offsets[1][i] << ":" << seed.weights[i] << endl;
	}*/
}


SeamassCore::~SeamassCore()
{
	delete optimizer_;
	delete inner_optimizer_;

	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		delete bases_[i];
	}
}


void SeamassCore::notice()
{
	cout << "seaMass - Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
	cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
	cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
}


void SeamassCore::init(Input& input, const std::vector<short>& scales)
{
	// for speed only, merge bins if rc_mz is set more than 8 times higher than the bin width
	// this is conservative, 4 times might be ok, but 2 times isn't enough
	//merge_bins(mzs, intensities, 0.125 / rc_mz);
    
	// INIT BASIS FUNCTIONS
	// Create our tree of bases
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "  Initialising overcomplete tree of basis functions ..." << endl;
    }
	if (input.binCountsIndex.size() <= 2)
	{
		dimensions_ = 1;

        new BasisBsplineMz(bases_, input.binCounts, input.binCountsIndex, input.binEdges, scales[0], Basis::Transient::NO);
        
		while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().scale[0] > -6)
		{
            new BasisBsplineScale(bases_, bases_.back()->getIndex(), 0, Basis::Transient::NO);
		}
	}
	else
	{
		dimensions_ = 2;

        new BasisBsplineMz(bases_, input.binCounts, input.binCountsIndex, input.binEdges, scales[0], Basis::Transient::YES);
        Basis* previousBasis = new BasisBsplineScantime(bases_, bases_.back()->getIndex(), input.startTimes, input.finishTimes, input.exposures, scales[1], Basis::Transient::NO);
 
        for (ii i = 0; static_cast<BasisBspline*>(bases_.back())->getGridInfo().scale[0] > -6; i++)
        {
            if (i > 0)
            {
                previousBasis = new BasisBsplineScale(bases_, previousBasis->getIndex(), 0, Basis::Transient::NO);
            }
            
            while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().extent[1] > 4)
            {
                new BasisBsplineScale(bases_, bases_.back()->getIndex(), 1, Basis::Transient::NO);
            }
        }
	}
    
    // INIT OPTIMISER
	inner_optimizer_ = new OptimizerSrl(bases_, input.binCounts, input.binCountsIndex);
	optimizer_ = new OptimizerAccelerationEve1(inner_optimizer_);
	optimizer_->init((fp)shrinkage_);
}


bool SeamassCore::step()
{
	if (iteration_ == 0 && getDebugLevel() % 10 >= 1)
	{
		li nnz = 0;
		li nx = 0;
		for (ii j = 0; j < (ii)bases_.size(); j++)
		{
            if (static_cast<BasisBspline*>(bases_[j])->getTransient() == Basis::Transient::NO)
			{
                for (size_t k = 0; k < optimizer_->xs()[j].size(); k++)
                {
                    nnz += optimizer_->xs()[j][k].nnz();
                    nx += optimizer_->xs()[j][k].size();
                }
			}
		}

        cout << getTimeStamp();
        cout << "   it:     0 nx: " << setw(10) << nx << " nnz: " << setw(10) << nnz;
        cout << " tol:  " << fixed << setprecision(8) << setw(10) << tolerance_ << endl;
	}

	iteration_++;
	double grad = optimizer_->step();

	if (getDebugLevel() % 10 >= 1)
	{
		li nnz = 0;
		for (ii j = 0; j < (ii)bases_.size(); j++)
		{
			if (static_cast<BasisBspline*>(bases_[j])->getTransient() == Basis::Transient::NO)
			{
                for (size_t k = 0; k < optimizer_->xs()[j].size(); k++)
                {
                    nnz += optimizer_->xs()[j][k].nnz();
                }
			}
		}
        cout << getTimeStamp();
		cout << "   it: " << setw(5) << iteration_;
		cout << " shrink: ";
		cout.unsetf(ios::floatfield); 
		cout << setprecision(4) << setw(6) << shrinkage_;
		cout << " nnz: " << setw(10) << nnz;
		cout << " grad: " << fixed << setprecision(8) << setw(10) << grad << endl;
	}

	if (grad <= tolerance_)
	{
		if (shrinkage_ == 0)
		{
			if (getDebugLevel() % 10 == 0) cout << "o" << endl;
			return false;
		}
		else
		{
			if (getDebugLevel() % 10 == 0) cout << "o" << flush;
			shrinkage_ *= (shrinkage_ > 0.0625 ? 0.5 : 0.0);
			optimizer_->init((fp)shrinkage_);
		}
	}
	else
	{
		if (getDebugLevel() % 10 == 0) cout << "." << flush;
	}
    
    if (grad != grad)
    {
        cout << "ARGGH!" << endl;
        return false;
    }
    
	return true;
}


ii SeamassCore::getIteration() const
{
	return iteration_;
}


void SeamassCore::getOutput(Output& output) const
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "  Getting output ..." << endl;
    }

    throw runtime_error("not implemented yet");

	/*output.scales.resize(dimensions_);
	output.offsets.resize(dimensions_);
	output.baselineScale.resize(dimensions_);
	output.baselineOffset.resize(dimensions_);
	output.baselineExtent.resize(dimensions_);

	for (ii i = 0; i < dimensions_; i++)
	{
		const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[dimensions_ - 1])->getGridInfo();

		output.baselineScale[i] = meshInfo.scale[i];
		output.baselineOffset[i] = meshInfo.offset[i];
		output.baselineExtent[i] = meshInfo.extent[i];
	}

	for (ii j = 0; j < (ii)bases_.size(); j++)
	{
        if (bases_[j]->getTransient() == Basis::Transient::NO)
		{
			const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[j])->getGridInfo();

			/*for (ii nz = 0; nz < optimizer_->xs()[j].nnz(); nz++)
			{
                fp v = 0.0;//optimizer_->xs()[j].getVs()[nz];
				if (v != 0.0)
				{
					output.weights.push_back(v);

					ii index = 0;
					for (ii d = 0; d < dimensions_; d++)
					{
						output.scales[d].push_back(meshInfo.scale[d]);
						output.offsets[d].push_back(meshInfo.offset[d] + index % meshInfo.extent[d]);
						index /= meshInfo.extent[d];
					}
				}
			}
		}
	}*/
}


void SeamassCore::getOutputBinCounts(std::vector<fp>& binCounts) const
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "  Deriving restored bin counts ..." << endl;

	vector<MatrixSparse> f;
	optimizer_->synthesis(f);
	f[0].output(binCounts.data());
}


void SeamassCore::getOutputControlPoints(ControlPoints& controlPoints) const
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "  Deriving control points ..." << endl;

	const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[dimensions_ - 1])->getGridInfo();

	vector<MatrixSparse> c(1);
	optimizer_->synthesis(c, dimensions_ - 1);
    vector<fp>(meshInfo.size()).swap(controlPoints.coeffs);
	c[0].output(controlPoints.coeffs.data());

	controlPoints.scale = meshInfo.scale;
	controlPoints.offset = meshInfo.offset;
	controlPoints.extent = meshInfo.extent;
}


void SeamassCore::getOutputControlPoints1d(ControlPoints& controlPoints) const
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "  Deriving 1D control points ..." << endl;

    const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[0])->getGridInfo();

    vector<MatrixSparse> c(1);
    optimizer_->synthesis(c, 0);
    vector<fp>(meshInfo.size()).swap(controlPoints.coeffs);
    c[0].output(controlPoints.coeffs.data());

    controlPoints.scale = meshInfo.scale;
    controlPoints.offset = meshInfo.offset;
    controlPoints.extent = meshInfo.extent;
    controlPoints.extent.push_back(meshInfo.count);
}

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

#include "OptimizerSrl.hpp"
#include "OptimizerAccelerationEve1.hpp"

#include "BasisBsplineMz.hpp"
#include "BasisBsplineScale.hpp"
#include "BasisBsplineScantime.hpp"

#include <iostream>
#include <iomanip>


using namespace std;


SeamassCore::SeamassCore(Input& input, const std::vector<short>& scales, double shrinkage, double tolerance, int debugLevel) : shrinkage_(shrinkage), tolerance_(tolerance), iteration_(0), debugLevel_(debugLevel)
{
	init(input, scales);
}


SeamassCore::SeamassCore(Input& input, const Output& seed, int debugLevel) : iteration_(0), debugLevel_(debugLevel)
{
	vector<short> scales(seed.scales.size());
	init(input, scales);

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
	//delete inner_optimizer_;

	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		delete bases_[i];
	}
}


void SeamassCore::notice()
{
	cout << endl;
	cout << "seaMass - Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
	cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
	cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
	cout << endl;
}


void SeamassCore::init(Input& input, const std::vector<short>& scales)
{
	// for speed only, merge bins if rc_mz is set more than 8 times higher than the bin width
	// this is conservative, 4 times might be ok, but 2 times isn't enough
	//merge_bins(mzs, intensities, 0.125 / rc_mz);

	////////////////////////////////////////////////////////////////////////////////////
	// INIT BASIS FUNCTIONS

	// Create our tree of bases
	if (input.spectrumIndex.size() <= 2)
	{
		dimensions_ = 1;

		new BasisBsplineMz(bases_, input.binCounts, input.spectrumIndex, input.binEdges, scales[0], Basis::Transient::NO);
		while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().extent[0] > 4)
		{
            new BasisBsplineScale(bases_, bases_.back()->getIndex(), 0, Basis::Transient::NO);
		}
	}
	/*else
	{
		dimensions_ = 2;

		BasisBsplineMz* basisMz = new BasisBsplineMz(bases_, input.binCounts, input.spectrumIndex, input.binEdges, scales[0], order, true);
		Basis* previousBasis = new BasisBsplineScantime(bases_, basisMz->getIndex(), input.startTimes, input.finishTimes, input.exposures, scales[1], order);
		while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().extent[0] > order + 1)
		{
			new BasisBsplineScale(bases_, bases_.back()->getIndex(), 0, order);
		}

		while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().extent[1] > order + 1)
		{
			previousBasis = new BasisBsplineScale(bases_, previousBasis->getIndex(), 1, order);
			while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().extent[0] > order + 1)
			{
				new BasisBsplineScale(bases_, bases_.back()->getIndex(), 0, order);
			}
		}
	}*/

	// Initialize the optimizer
	b_.init(1, (ii)input.binCounts.size(), input.binCounts.data());
	
	//inner_optimizer_ = new OptimizerSrl(bases_, b_, debugLevel_);
	optimizer_ = new OptimizerSrl(bases_, b_, debugLevel_);
	//optimizer_ = new OptimizerAccelerationEve1(inner_optimizer_, debugLevel_);
	optimizer_->init((fp)shrinkage_);
}


bool SeamassCore::step()
{
	if (iteration_ == 0 && debugLevel_ > 0)
	{
		li nnz = 0;
		li nx = 0;
		for (ii j = 0; j < (ii)bases_.size(); j++)
		{
            if (static_cast<BasisBspline*>(bases_[j])->getTransient() == Basis::Transient::NO)
			{
				nnz += optimizer_->xs()[j].nnz();
				nx += optimizer_->xs()[j].size();
			}
		}

		cout << " it:     0 nx: " << setw(10) << nx << " nnz: " << setw(10) << nnz << " tol:  " << fixed << setprecision(8) << setw(10) << tolerance_ << endl;
	}

	iteration_++;
	double grad = optimizer_->step();

	if (debugLevel_ > 0)
	{
		li nnz = 0;
		for (ii j = 0; j < (ii)bases_.size(); j++)
		{
			if (static_cast<BasisBspline*>(bases_[j])->getTransient() == Basis::Transient::NO)
			{
				nnz += optimizer_->xs()[j].nnz();

				//vector<fp> ts(optimizer_->xs()[j].m() * optimizer_->xs()[j].n());
				//optimizer_->xs()[j].convertToDense(ts.data());
				//for (ii i = 0; i < ts.size(); i++) cout << j << ":" << i << ":" << ts[i] << endl;
			}
		}
		cout << " it: " << setw(5) << iteration_;
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
			if (debugLevel_ == 0) cout << "o" << endl;
			return false;
		}
		else
		{
			if (debugLevel_ == 0) cout << "o" << flush;
			shrinkage_ *= (shrinkage_ > 0.0625 ? 0.5 : 0.0);
			optimizer_->init((fp)shrinkage_);
		}
	}
	else
	{
		if (debugLevel_ == 0) cout << "." << flush;
	}
    
	return true;
}


ii SeamassCore::getIteration() const
{
	return iteration_;
}


void SeamassCore::getOutput(Output& output) const
{
#ifndef NDEBUG
	cout << iteration_ << " getOutput" << endl;
#endif

	output.scales.resize(dimensions_);
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

			for (ii nz = 0; nz < optimizer_->xs()[j].nnz(); nz++)
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
	}
}


void SeamassCore::getOutputBinCounts(std::vector<fp>& binCounts) const
{
#ifndef NDEBUG
	cout << iteration_ << " getOutputBinCounts" << endl;
#endif

	MatrixSparse f;
	optimizer_->synthesis(f);
	binCounts.resize(b_.size());
	f.output(binCounts.data());
}


void SeamassCore::getOutputControlPoints(ControlPoints& controlPoints) const
{
#ifndef NDEBUG
	cout << iteration_ << " getOutputControlPoints" << endl;
#endif

	const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[dimensions_ - 1])->getGridInfo();

	MatrixSparse c;
	optimizer_->synthesis(c, dimensions_ - 1);
	controlPoints.coeffs.resize(meshInfo.size());
	c.output(controlPoints.coeffs.data());

	controlPoints.scale = meshInfo.scale;
	controlPoints.offset = meshInfo.offset;
	controlPoints.extent = meshInfo.extent;
}


/*double CatmullRomInterpolate(
	double y0, double y1,
	double y2, double y3,
	double mu)
{
	double a0, a1, a2, a3, mu2;

	mu2 = mu*mu;
	a0 = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3;
	a1 = y0 - 2.5*y1 + 2 * y2 - 0.5*y3;
	a2 = -0.5*y0 + 0.5*y2;
	a3 = y1;

	return a0*mu*mu2 + a1*mu2 + a2*mu + a3;
}


void remove_zeros(vector< vector<double> >& mzs, vector< vector<double> >& intensities)
{
	for (ii j = 0; j < intensities.size(); j++)
	{
		ii shift = 0;
		for (ii i = 0; i < intensities[j].size(); i++)
		{
			if (intensities[j][i] <= 0.0)
			{
				shift++;
			}
			else if (shift > 0)
			{
				intensities[j][i - shift] = intensities[j][i];
				mzs[j][i - shift] = mzs[j][i];
			}
		}
		intensities[j].resize(intensities[j].size() - shift);
		mzs[j].resize(mzs[j].size() - shift);
	}
}


void merge_bins(vector< vector<fp> >& mzs,
	vector< vector<fp> >& intensities,
	double width)
{
#pragma omp parallel for
	for (ii j = 0; j < mzs.size(); j++)
	{
		double w = 0;
		double v = 0;
		ii n = 0;
		ii k = 1;
		for (ii i = 1; i < mzs[j].size();)
		{
			w += mzs[j][i] - mzs[j][i - 1];
			n++;

			if (w > width || i == mzs[j].size() - 1)
			{
				if (n == 1)
				{
					mzs[j][k] = mzs[j][i];
					intensities[j][k - 1] = intensities[j][i - 1];
					i++;
				}
				else
				{
					mzs[j][k] = mzs[j][i - 1];
					intensities[j][k - 1] = (fp)v;
				}

				w = 0;
				v = 0;
				n = 0;
				k++;
			}
			else
			{
				i++;
			}

			if (i < intensities[j].size()) v += intensities[j][i - 1];
		}

		if (k != 1) mzs[j].resize(k);
		intensities[j].resize(k - 1);
	}
}
*/

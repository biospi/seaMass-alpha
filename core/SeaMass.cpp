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


#include "SeaMass.hpp"

#include "OptimizerSrl.hpp"
#include "OptimizerAccelerationEve1.hpp"

#include "BasisBsplineMz.hpp"
#include "BasisBsplineScale.hpp"
#include "BasisBsplineScantime.hpp"

#include <iostream>
#include <iomanip>


using namespace std;


SeaMass::SeaMass(Input& input, const std::vector<ii>& scales, double shrinkage, double tolerance, ii debugLevel) : shrinkage_(shrinkage), tolerance_(tolerance), iteration_(0), debugLevel_(debugLevel)
{
	init(input, scales);
}


SeaMass::SeaMass(Input& input, const Output& seed, ii debugLevel) : iteration_(0), debugLevel_(debugLevel)
{
	vector<ii> scales(seed.scales.size());
	init(input, scales);

	// todo
	/*cout << seed.weights.size() << endl;
	for (ii i = 0; i < seed.weights.size(); i++)
	{
		//cout << seed.scales[0][i] << "," << seed.scales[1][i] << ":" << seed.offsets[0][i] << "," << seed.offsets[1][i] << ":" << seed.weights[i] << endl;
	}*/
}


SeaMass::~SeaMass()
{
	delete optimizer_;
	delete inner_optimizer_;

	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		delete bases_[i];
	}
}


void SeaMass::notice()
{
	cout << endl;
	cout << "seaMass - Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
	cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
	cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
	cout << endl;
}


void SeaMass::init(Input& input, const std::vector<ii>& scales)
{
	// for speed only, merge bins if rc_mz is set more than 8 times higher than the bin width
	// this is conservative, 4 times might be ok, but 2 times isn't enough
	//merge_bins(mzs, intensities, 0.125 / rc_mz);

	////////////////////////////////////////////////////////////////////////////////////
	// INIT BASIS FUNCTIONS

	// Create our tree of bases
	ii order = 3; // B-spline order
	if (input.spectrumIndex.size() <= 2)
	{
		dimensions_ = 1;

		BasisBsplineMz* basisMz = new BasisBsplineMz(bases_, input.binCounts, input.spectrumIndex, input.binEdges, scales[0], order);
		while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().extent[0] > order + 1)
		{
			new BasisBsplineScale(bases_, bases_.back()->getIndex(), 0, order);
		}
	}
	else
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
	}

	// Initialize the optimizer
	/*vector<fp>* binCounts = new vector<fp>(input.binCounts.size(), 0.0);
	for (li i = 1; i < (li)input.binCounts.size() - 1; i++)
	{
		(*binCounts)[i - 1] += 0.25 * input.binCounts[i];
		(*binCounts)[i] += 0.5 * input.binCounts[i];
		(*binCounts)[i + 1] += 0.25 * input.binCounts[i];
	}
	b_.init((li)binCounts->size(), 1, binCounts->data());*/
	b_.init((li)input.binCounts.size(), 1, input.binCounts.data());
	inner_optimizer_ = new OptimizerSrl(bases_, b_, debugLevel_);
	//optimizer_ = new OptimizerSrl(bases_, b_);
	optimizer_ = new OptimizerAccelerationEve1(inner_optimizer_, debugLevel_);
	optimizer_->init((fp)shrinkage_);
}


bool SeaMass::step()
{
	if (iteration_ == 0 && debugLevel_ > 0)
	{
		li nnz = 0;
		li nx = 0;
		for (ii j = 0; j < (ii)bases_.size(); j++)
		{
			if (!static_cast<BasisBspline*>(bases_[j])->isTransient())
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
			if (!static_cast<BasisBspline*>(bases_[j])->isTransient())
			{
				nnz += optimizer_->xs()[j].nnz();
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


ii SeaMass::getIteration() const
{
	return iteration_;
}


void SeaMass::getOutput(Output& output) const
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
		if (!bases_[j]->isTransient())
		{
			const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[j])->getGridInfo();

			for (ii i = 0; i < optimizer_->xs()[j].size(); i++)
			{
				fp x = optimizer_->xs()[j].getVs()[i];
				if (x > 0.0)
				{
					output.weights.push_back(x);

					ii index = i;
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


void SeaMass::getOutputBinCounts(std::vector<fp>& binCounts) const
{
#ifndef NDEBUG
	cout << iteration_ << " getOutputBinCounts" << endl;
#endif

	binCounts.assign(b_.size(), 0.0);
	Matrix f; f.init(b_.m(), b_.n(), binCounts.data());
	optimizer_->synthesis(f);
}


void SeaMass::getOutputControlPoints(ControlPoints& controlPoints) const
{
#ifndef NDEBUG
	cout << iteration_ << " getOutputControlPoints" << endl;
#endif

	const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[dimensions_ - 1])->getGridInfo();

	controlPoints.coeffs.assign(meshInfo.size(), 0.0);
	Matrix c; c.init(meshInfo.m(), meshInfo.n,  controlPoints.coeffs.data());
	optimizer_->synthesis(c, dimensions_ - 1);

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
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

#include "BasisBsplineMz.hpp"
#include "BasisBsplineScale.hpp"

#include <iostream>
#include <iomanip>


using namespace std;


SeaMass::SeaMass(Input& input, const std::vector<ii>& scales) : iteration_(0)
{
	init(input, scales);
}


SeaMass::SeaMass(Input& input, const Output& seed) : iteration_(0)
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

	cout << "Init Basis Functions..." << endl;

	// Create our tree of bases
	ii order = 3; // B-spline order

	if (input.spectrumIndex.size() <= 2)
	{
		dimensions_ = 1;

		BasisBsplineMz* basisMz = new BasisBsplineMz(bases_, input.binCounts, input.spectrumIndex, input.binEdges, scales[0], order);
		while (static_cast<BasisBspline*>(bases_.back())->getMeshInfo().extent[0] > order + 1)
		{
			new BasisBsplineScale(bases_, static_cast<BasisBspline*>(bases_.back())->getIndex(), 0, order);
		}
	}
	else
	{
		/*dimensions_ = 2;

		BasisResampleMZ* bResampleMZ = new BasisResampleMZ(bases_, input.bin_counts, input.spectrum_index, input.bin_edges, scales[0], order, true);
		BasisResampleRT* bResampleRT = new BasisResampleRT(bases_, bResampleMZ, input.start_times, input.finish_times, input.exposures, scales[1], order);
		Basis* last = bResampleRT;
		while (last->get_cm().n[1] > order + 1)
		{
			last = new BasisDyadicScale(bases_, last, 1, order);
			while (bases_.back()->get_cm().n[0] > order + 1)
			{
				new BasisDyadicScale(bases_, bases_.back(), 0, order);
			}
		}*/
	}

	g_.init((li)input.binCounts.size(), 1, 1, input.binCounts.data());
	optimiser_ = new OptimizerAsrl(bases_, g_, 2);
}


bool SeaMass::step(double shrinkage, double tolerance)
{
	if (iteration_ == 0)
	{
		li nc = 0;
		for (ii j = 0; j < (ii)bases_.size(); j++)
		{
			if (!static_cast<BasisBspline*>(bases_[j])->isTransient())
			{
				nc += static_cast<BasisBspline*>(bases_[j])->getMeshInfo().size();
			}
		}
		cout << "L1 nc=" << nc << " shrinkage=" << setprecision(3) << defaultfloat << shrinkage << " tolerance=" << tolerance << endl;
	}
	iteration_++;

	double grad = optimiser_->step((fp)shrinkage);

	cout << "  f: " << fixed << setw(5) << iteration_;
	cout << "  grad: " << fixed << setprecision(6) << setw(8) << grad << endl;

	return grad > tolerance;
}


ii SeaMass::getIteration() const
{
	return iteration_;
}


void SeaMass::getOutput(Output& output) const
{
	/*output.scales.resize(dimensions_);
	output.offsets.resize(dimensions_);
	output.baselineScale.resize(dimensions_);
	output.baselineOffset.resize(dimensions_);
	output.baselineExtent.resize(dimensions_);

	for (ii i = 0; i < dimensions_; i++)
	{
		const BasisBspline::MeshInfo& meshInfo = static_cast<BasisBspline*>(bases_[dimensions_ - 1])->getMeshInfo();

		output.baselineScale[i] = meshInfo.scale[i];
		output.baselineOffset[i] = meshInfo.offset[i];
		output.baselineExtent[i] = meshInfo.extent[i];
	}

	for (ii j = 0; j < (ii)bases_.size(); j++)
	{
		if (!bases_[j]->isTransient())
		{
			const BasisBspline::MeshInfo& meshInfo = static_cast<BasisBspline*>(bases_[j])->getMeshInfo();

			for (ii i = 0; i < (ii)bases_[j]->getSize(); i++)
			{
				fp c = optimiser_->getCoeffs()[j].vs_[i];
				if (c > 0.0)
				{
					output.weights.push_back(c);

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
	}*/
}


void SeaMass::getOutputBinCounts(std::vector<fp>& binCounts) const
{
	binCounts.assign(g_.size(), 0.0);
	Matrix f; f.init(g_.m(), g_.n(), g_.n(), binCounts.data());
	optimiser_->synthesis(f);
}


void SeaMass::getOutputControlPoints(ControlPoints& controlPoints) const
{
	const BasisBspline::MeshInfo& meshInfo = static_cast<BasisBspline*>(bases_[dimensions_ - 1])->getMeshInfo();

	controlPoints.coeffs.assign(meshInfo.size(), 0.0);
	Matrix c; c.init(meshInfo.m(), meshInfo.n, meshInfo.n, controlPoints.coeffs.data());
	optimiser_->synthesis(c, dimensions_ - 1);

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
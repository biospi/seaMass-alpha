//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  biospi Laboratory, EEE, University of Liverpool, UK
//
// This file is part of seaMass-TD.
//
// seaMass-TD is free software: you can redistribute it and/or modify
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
// along with seaMass-TD.  If not, see <http://www.gnu.org/licenses/>.
//


#include "BasisChargeDistribution.hpp"
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;


/*double
BasisChargeDistribution::
protonMass_ = 1.007276466879; // in Daltons


BasisChargeDistribution::
BasisChargeDistribution(std::vector<Basis*>& bases, const std::vector<fp>& binCounts, ii scale, ii offset, ii maxMass, ii binsPerDalton, bool transient)
	: Basis(bases, transient)
{
	// min and max mz of input stream
	double binWidth = pow(2.0, scale) / 60.0 / 1.0033548378;
	double mz0 = offset * binWidth;
	double mz1 = (offset + binCounts.size()) * binWidth;

	// min mass of output
	double minMass = mz0 - protonMass_;

	// offset and extent of output bins
	li yOffset = (li)floor(minMass * binsPerDalton);
	li yExtent = maxMass * binsPerDalton - yOffset;

	// output bin edges (mass in Daltons)
	vector<double> binEdges(yExtent + 1);
	binEdges[0] = minMass;
	for (li y = 1; y < (li)binEdges.size(); y++)
	{
		binEdges[y] = (yOffset + y) / (double) binsPerDalton;
	}

	// minimum and maximum charge state z that could appear at each bin edge
	groupZOffsets_.resize(binEdges.size());
	vector<short> groupZOffsets1(binEdges.size());
	for (li y = 0; y < (li)binEdges.size(); y++)
	{
		groupZOffsets_[y] = (short)ceil(binEdges[y] / (mz1 - protonMass_));
		groupZOffsets1[y] = (short)floor(binEdges[y] / (mz0 - protonMass_));
	}

	/*for (li y = 0; y < 5; y++)
	{
		cout << y << ":" << setprecision(20) << binEdges[y] << " -> " << groupZOffsets1[y] << ":" << groupZOffsets_[y] << endl;
	}
	cout << "..." << endl;
	for (li y = binEdges.size() - 5; y < (li)binEdges.size(); y++)
	{
		cout << y << ":" << setprecision(20) << binEdges[y] << " -> " << groupZOffsets1[y] << ":" << groupZOffsets_[y] << endl;
	}*/

	// groupXOffsets_ contains start indices to each group of coefficients, followed by total number of coefficients
	// groupZOffsets_ contains the first charge state z in each group (the last value can be ignored)
	/*groupXOffsets_.resize(groupZOffsets_.size());
	groupXOffsets_[0] = 0;
	for (li y = 1; y < (li)groupXOffsets_.size(); y++)
	{
		groupXOffsets_[y] = (groupZOffsets1[y] - groupZOffsets_[y - 1] + 1) + groupXOffsets_[y - 1];
	}*/

	/*cout << endl;
	for (li y = 0; y < 5; y++)
	{
		cout << y << ":" << setprecision(20) << groupXOffsets_[y] << "," << groupZOffsets_[y] << endl;
	}
	cout << "..." << endl;
	for (li y = groupXOffsets_.size() - 5; y < (li)groupXOffsets_.size(); y++)
	{
		cout << y << ":" << setprecision(20) << groupXOffsets_[y] << "," << groupZOffsets_[y] << endl;
	}*/

	// populate coo matrix
	/*vector<fp> acoo;
	vector<ii> rowind;
	vector<ii> colind;
	//Bspline bspline(order, 65536); // bspline basis function lookup table
	ii minRow = 1000000000, maxRow = 0, minCol = 1000000000, maxCol = 0;
	for (li y = 0; y < (li)groupXOffsets_.size() - 1; y++)
	{
		for (li x = groupXOffsets_[y], z = groupZOffsets_[y]; x < groupXOffsets_[y + 1]; x++, z++)
		{
			ii i = (ii)(((binEdges[y] / z + protonMass_) - mz0) / binWidth);

			acoo.push_back(1.0);
			rowind.push_back(i);
			colind.push_back(x);

			minRow = minRow < i ? minRow : i;
			maxRow = maxRow > i ? maxRow : i;
			minCol = minCol < x ? minCol : x;
			maxCol = maxCol > x ? maxCol : x;
		}
	}
	cout << minRow << "," << maxRow << endl;
	cout << minCol << "," << maxCol << endl;
	cout << binCounts.size() << "," << groupXOffsets_[groupXOffsets_.size() - 1] << endl;
	a_.init(binCounts.size(), groupXOffsets_[groupXOffsets_.size() - 1], (ii)acoo.size(), acoo.data(), rowind.data(), colind.data());

//#ifndef NDEBUG
	cout << " " << getIndex() << " BasisChargeDistribution";
	if (isTransient()) cout << " (t)";
	cout << "  A" << a_ << " (";
	cout.unsetf(ios::floatfield);
	cout << setprecision(2) << a_.mem() / 1024.0 / 1024.0 << "Mb)" << endl;
//#endif
}


BasisChargeDistribution::~BasisChargeDistribution()
{
}


ii BasisChargeDistribution::getM() const
{
	return 0;
}


ii BasisChargeDistribution::getN() const
{
	return groupXOffsets_[groupXOffsets_.size() - 1];
}


void
BasisChargeDistribution::
synthesis(Matrix& f, const Matrix& x, bool accumulate) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisChargeDistribution::synthesis" << endl;
#endif

	f.mul(a_, x, accumulate, false, false);
}


void
BasisChargeDistribution::
analysis(Matrix& xE, const Matrix& fE, bool sqrA) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisChargeDistribution::analysis" << endl;
#endif

	if (sqrA)
	{
		MatrixSparse aSqrd;
		aSqrd.elementwiseSqr(a_);
		xE.mul(aSqrd, fE, false, true, false);
	}
	else
	{
		xE.mul(a_, fE, false, true, false);
	}
}*/


/*void
BasisChargeDistribution::
shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage) const
{
	// GROUP-WISE SHRINKAGE! (note - intuitive implemention, not mathematically verified yet)
	// also horrible gather operation
	vector<fp> gcs(ci1s.back() - ci0s.front() + 1, 0.0);

	// sum up coefficients per group
	for (ii z = 0; z < cos.size() - 1; z++)
	for (ii i = cos[z]; i < cos[z + 1]; i++)
	for (ii j = 0; j < as.size(); j++)
	if (es[j*cos.back() + i] > 0.0)
	{
		gcs[i - cos[z] + ci0s[z] - ci0s.front()] += cs[j*cos.back() + i];
	}

	// scale the shrinkage to be proportional to the contribution of this coefficient to the group total
	#pragma omp parallel for
	for (ii z = 0; z < cos.size(); z++)
	for (ii i = cos[z]; i < cos[z + 1]; i++)
	for (ii j = 0; j < as.size(); j++)
	if (es[j*cos.back() + i] > 0.0)
	{
		if (cs[j*cos.back() + i] > 0.0)
		{
			double scale = cs[j*cos.back() + i] / gcs[i - cos[z] + ci0s[z] - ci0s.front()];
			es[j*cos.back() + i] *= cs[j*cos.back() + i] / (scale * shrinkage * l2[j*cos.back() + i] + wcs[j*cos.back() + i]);
		}
		else
		{
			es[j*cos.back() + i] = 0.0;
		}
	}
}*/


/*void
BasisChargeDistribution::
write_cs(const std::vector<fp>& cs) const
{
	ostringstream oss; oss << "profile" << get_index() << ".csv";
	ofstream ofs(oss.str().c_str());

	// sum up coefficients per group
	vector<fp> gcs(ci1s.back() - ci0s.front() + 1, 0.0);
	for (ii z = 0; z < cos.size() - 1; z++)
	for (ii i = cos[z]; i < cos[z + 1]; i++)
	for (li j = 0; j < as.size(); j++)
	{
		gcs[i - cos[z] + ci0s[z] - ci0s.front()] += cs[j*cos.back() + i];
	}

	ofs << "mass,intensity" << setprecision(10) << endl;
	for (ii i = 0; i < gcs.size(); ++i)
	{
		if (gcs[i] > 0.0 || i > 0 && gcs[i - 1] > 0.0 || i < gcs.size() - 1 && gcs[i + 1] > 0.0)
		{
			ofs << (ci0s.front() + i) * mass_interval << "," << gcs[i] << endl;
		}
	}
}*/


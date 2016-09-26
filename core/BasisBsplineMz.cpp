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


#include "BasisBsplineMZ.hpp"

#include "Bspline.hpp"

#include <limits>
#include <iomanip>


using namespace std;


BasisBsplineMz::BasisBsplineMz(std::vector<Basis*>& bases, const std::vector<fp>& binCounts, const std::vector<li>& spectrumIndex,
	                           const std::vector<double>& binEdges, short resolution, ii order, bool transient)
	: BasisBspline(bases, 1, transient)
{
	if (spectrumIndex.size() > 0)
	{
		is = spectrumIndex;
	}
	else
	{
		is.push_back(0);
		is.push_back(binCounts.size());
	}

	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	// init As
	as.resize(is.size() - 1);

	// find min and max m/z across spectra, and m for each A
	mzMin = numeric_limits<double>::max();
	mzMax = 0.0;
	double mzDiff = numeric_limits<double>::max();
	vector<ii> ms(as.size());
	for (ii k = 0; k < (ii)as.size(); k++)
	{
		ms[k] = (ii)(is[k + 1] - is[k]);
		mzMin = binEdges[is[k] + k] < mzMin ? binEdges[is[k] + k] : mzMin;
		mzMax = binEdges[is[k + 1] + k] > mzMax ? binEdges[is[k + 1] + k] : mzMax;

		for (ii i = 0; i < ms[k]; i++)
		{
			double diff = binEdges[is[k] + k + i + 1] - binEdges[is[k] + k + i];
			mzDiff = diff < mzDiff ? diff : mzDiff;
		}
	}

	ii resolutionAuto = (ii) floor(log2(1.0 / mzDiff / 60.0 / 1.0033548378));
	if (resolution == numeric_limits<short>::max())
	{
		resolution = resolutionAuto;
	}

	// Bases per 1.0033548378Th (difference between carbon12 and carbon13)
	double bpi = pow(2.0, (double)resolution) * 60 / 1.0033548378;

	// fill in b-spline mesh info
	meshInfo().n = (ii)as.size();
	meshInfo().scale[0] = resolution;
	meshInfo().offset[0] = (ii)floor(mzMin * bpi);
	meshInfo().extent[0] = ((ii)ceil(mzMax * bpi)) + order - meshInfo().offset[0];

	// populate coo matrix
	ii done = 0;
	Bspline bspline(order, 65536); // bspline basis function lookup table
	#pragma omp parallel for
	for (ii k = 0; k < (ii)as.size(); k++)
	{
		vector<fp> acoo;
		vector<ii> rowind;
		vector<ii> colind;

		for (ii i = 0; i < ms[k]; i++)
		{
			if (binCounts[is[k] + i] >= 0.0)
			{
				double cfMin = binEdges[is[k] + k + i] * bpi;
				double cfMax = binEdges[is[k] + k + i + 1] * bpi;

				ii cMin = (ii)floor(cfMin);
				ii cMax = ((ii)ceil(cfMax)) + order;

				// work out basis coefficients
				for (ii c = cMin; c < cMax; c++)
				{
					double bfMin = (double)(c - order);
					double bfMax = (double)(c + 1);

					// intersection of bin and basis, between 0 and order+1
					double bMin = cfMin > bfMin ? cfMin - bfMin : 0.0;
					double bMax = cfMax < bfMax ? cfMax - bfMin : bfMax - bfMin;

					// basis coefficient b is _integral_ of area under b-spline basis
					fp b = (fp)(bspline.ibasis(bMax) - bspline.ibasis(bMin));

					acoo.push_back(b);
					rowind.push_back(i);
					colind.push_back(c - meshInfo().offset[0]);
				}
			}
		}

		// create A
		as[k].init(ms[k], meshInfo().extent[0], (ii)acoo.size(), acoo.data(), rowind.data(), colind.data());

		// display progress update
		#pragma omp critical
		{
			done++;
			if (done % 100 == 0)
			{
				for (int i = 0; i < 256; ++i) cout << '\b';
				cout << getIndex() << " BasisBsplineMz " << setw(1 + (int)(log10((float)as.size()))) << done << "/" << as.size() << " " << flush;
			}
		}
	}
	for (int i = 0; i < 256; ++i) cout << '\b';

#ifndef NDEBUG
	li m = 0; for (ii k = 0; k < (ii)as.size(); k++) m += as[k].m();
	li n = 0; for (ii k = 0; k < (ii)as.size(); k++) n += as[k].n();
	li nnz = 0; for (ii k = 0; k < (ii)as.size(); k++) nnz += as[k].nnz();
	li mem = 0; for (ii k = 0; k < (ii)as.size(); k++) mem += as[k].mem();
	cout << " " << getIndex() << " BasisBsplineMz";
	if (isTransient()) cout << " (t)";
	cout << " range=" << setprecision(3) << mzMin << ":" << defaultfloat << mzDiff << ":" << fixed << mzMax << "Th";
	cout << " resolution=" << fixed << setprecision(1) << resolution << " (" << bpi << " bases per 1.0033548378Th)";
	cout << " " << meshInfo() << endl;
	cout << "  A{" << m << "," << n << "}:" << nnz << "/" << m * n << "=" << defaultfloat << setprecision(2) << nnz / (double) (m * n) << "% (" << defaultfloat << setprecision(2) << (2 * mem) / 1024.0 / 1024.0 << "Mb)" << endl;
#endif

	if (resolutionAuto != resolution)
	{
		cerr << endl << "WARNING: resolution is not the suggested value of " << resolutionAuto << ". Continue at your own risk!" << endl << endl;
	}
}


BasisBsplineMz::~BasisBsplineMz()
{
}


void BasisBsplineMz::synthesis(Matrix& f, const Matrix& c, bool accumulate) const
{
	if (!f) f.init(is.back(), 1);

	# pragma omp parallel for
	for (ii k = 0; k < (ii)as.size(); k++)
	{
		Matrix f_sub; f_sub.init(f, is[k], 0, as[k].m(), 1);
		Matrix c_sub; c_sub.init(c, k * as[k].n(), 0, as[k].n(), 1, 1);

#ifndef NDEBUG
		cout << " " << getIndex() << " BasisBsplineMz::synthesis[" << k << "]" << endl;
#endif

		f_sub.mul(as[k], c_sub, accumulate, false);
	}
}


void BasisBsplineMz::analysis(Matrix& cE, const Matrix& fE, bool sqrA) const
{
	if (!cE) cE.init(getMeshInfo().m(), getMeshInfo().n);

	# pragma omp parallel for
	for (ii k = 0; k < (ii)as.size(); k++)
	{
		Matrix cESub; cESub.init(cE, k * as[k].n(), 0, as[k].n(), 1, 1);
		Matrix fESub; fESub.init(fE, is[k], 0, as[k].m(), 1);

#ifndef NDEBUG
		cout << " " << getIndex() << " BasisBsplineMz::analysis[" << k << "]" << endl;
#endif

		if (sqrA)
		{
			MatrixSparse a_sqrd; a_sqrd.elementwiseSqr(as[k]);
			cESub.mul(a_sqrd, fESub, false, true);
		}
		else
		{
			cESub.mul(as[k], fESub, false, true);
		}
	}
}



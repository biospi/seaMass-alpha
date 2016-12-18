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


#include "BasisBsplineMz.hpp"

#include "Bspline.hpp"

#include <limits>
#include <iomanip>
#include <cmath>
#include <iostream>


using namespace std;


BasisBsplineMz::BasisBsplineMz(std::vector<Basis*>& bases, const std::vector<fp>& binCounts, const std::vector<li>& spectrumIndex,
	                           const std::vector<
                               double>& binEdges, short scale, Transient transient, int order)
	: BasisBspline(bases, 1, transient), nnzRows_(0)
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineMz";
    if (getTransient() == Transient::YES) cout << " (t)";
	cout << endl;
#endif

	if (spectrumIndex.size() > 0)
	{
		is_ = spectrumIndex;
	}
	else
	{
		is_.push_back(0);
		is_.push_back(binCounts.size());
	}

	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	// init As
	as_.resize(is_.size() - 1);
    aTs_.resize(is_.size() - 1);

	// find min and max m/z across spectra, and m for each A
	double mzMin = numeric_limits<double>::max();
	double mzMax = 0.0;
	double mzDiff = numeric_limits<double>::max();
	vector<ii> ms(as_.size());
	for (ii k = 0; k < (ii)as_.size(); k++)
	{
		ms[k] = (ii)(is_[k + 1] - is_[k]);
		mzMin = binEdges[is_[k] + k] < mzMin ? binEdges[is_[k] + k] : mzMin;
		mzMax = binEdges[is_[k + 1] + k] > mzMax ? binEdges[is_[k + 1] + k] : mzMax;

		for (ii i = 0; i < ms[k]; i++)
		{
			double diff = binEdges[is_[k] + k + i + 1] - binEdges[is_[k] + k + i];
			mzDiff = diff < mzDiff ? diff : mzDiff;
		}
	}

	ii scaleAuto = (ii) floor(log2(1.0 / mzDiff / 60.0 / 1.0033548378));
	if (scale == numeric_limits<short>::max())
	{
		scale = scaleAuto;
		cout << "Autodetected mz_scale=" << scale << endl;
	}

	// Bases per 1.0033548378Th (difference between carbon12 and carbon13)
	double bpi = pow(2.0, (double)scale) * 60 / 1.0033548378;

	// fill in b-spline grid info
	gridInfo().n = (ii)as_.size();
	gridInfo().scale[0] = scale;
	gridInfo().offset[0] = (ii)floor(mzMin * bpi);
	gridInfo().extent[0] = ((ii)ceil(mzMax * bpi)) + order - gridInfo().offset[0];

	// populate coo matrix
	ii done = 0;
	Bspline bspline(order, 65536); // bspline basis function lookup table
	#pragma omp parallel for
	for (ii k = 0; k < (ii)as_.size(); k++)
	{
		vector<fp> acoo;
		vector<ii> rowind;
		vector<ii> colind;

		for (ii i = 0; i < ms[k]; i++)
		{
			if (binCounts[is_[k] + i] >= 0.0)
			{
				double xfMin = binEdges[is_[k] + k + i] * bpi;
				double xfMax = binEdges[is_[k] + k + i + 1] * bpi;

				ii xMin = (ii)floor(xfMin);
				ii xMax = ((ii)ceil(xfMax)) + order;

				// work out basis coefficients
				for (ii x = xMin; x < xMax; x++)
				{
					double bfMin = (double)(x - order);
					double bfMax = (double)(x + 1);

					// intersection of bin and basis, between 0 and order+1
					double bMin = xfMin > bfMin ? xfMin - bfMin : 0.0;
					double bMax = xfMax < bfMax ? xfMax - bfMin : bfMax - bfMin;

					// basis coefficient b is _integral_ of area under b-spline basis
					fp b = (fp)(bspline.ibasis(bMax) - bspline.ibasis(bMin));

					acoo.push_back(b);
					rowind.push_back(i);
					colind.push_back(x - gridInfo().offset[0]);
				}
			}
		}

		// create A
		as_[k].init(ms[k], gridInfo().extent[0], (ii)acoo.size(), acoo.data(), rowind.data(), colind.data());
        aTs_[k].init(gridInfo().extent[0], ms[k], (ii)acoo.size(), acoo.data(), colind.data(), rowind.data());

		// display progress update
		#pragma omp critical
		{
			done++;
			if (done % 100 == 0)
			{
				for (int i = 0; i < 256; ++i) cout << '\b';
				cout << getIndex() << " BasisBsplineMz " << setw(1 + (int)(log10((float)as_.size()))) << done << "/" << as_.size() << " " << flush;
			}
		}
	}
	for (int i = 0; i < 256; ++i) cout << '\b';
    
    for (ii k = 0; k < (ii)as_.size(); k++) nnzRows_ += as_[k].n();

#ifndef NDEBUG
	li m = 0; for (ii k = 0; k < (ii)as_.size(); k++) m += as_[k].m();
	li n = 0; for (ii k = 0; k < (ii)as_.size(); k++) n += as_[k].n();
	li nnz = 0; for (ii k = 0; k < (ii)as_.size(); k++) nnz += as_[k].nnz();
	li mem = 0; for (ii k = 0; k < (ii)as_.size(); k++) mem += as_[k].mem();
	cout << "  range=" << setprecision(3) << mzMin << ":";
    cout.unsetf(std::ios::floatfield);
    cout << mzDiff << ":" << fixed << mzMax << "Th";
	cout << " scale=" << fixed << setprecision(1) << scale << " (" << bpi << " bases per 1.0033548378Th)";
    cout.unsetf(ios::floatfield);
	cout << " " << gridInfo();
	cout << " mem=" << fixed << setprecision(2) << (2 * mem) / 1024.0 / 1024.0 << "Mb" << endl;
	cout << "  A{" << m << "," << n << "}:" << nnz << "/" << m * n << "=";
	cout << setprecision(2) << nnz / (double)(m * n) << "%";
	cout.unsetf(ios::floatfield); 
	cout << endl;
#endif

	if (scaleAuto != scale)
	{
		cerr << endl << "WARNING: scale is not the suggested value of " << scaleAuto << ". Continue at your own risk!" << endl << endl;
	}
}


BasisBsplineMz::~BasisBsplineMz()
{
}


void BasisBsplineMz::synthesis(MatrixSparse& f, const MatrixSparse& x, bool accumulate) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineMz::synthesis[" << 0 << "]" << endl;
#endif
    
	f.mul(accumulate ? MatrixSparse::Accumulate::YES : MatrixSparse::Accumulate::NO, x, MatrixSparse::Transpose::NO, aTs_[0]);
    
	/*if (!f) f.init(is_.back(), 1);

#ifdef NDEBUG
	# pragma omp parallel for
#endif
	for (ii k = 0; k < (ii)as_.size(); k++)
	{
		Matrix fSub; fSub.init(f, is_[k], 0, as_[k].m(), 1);
		Matrix xSub; xSub.init(x, k * as_[k].n(), 0, as_[k].n(), 1);

#ifndef NDEBUG
		cout << " " << getIndex() << " BasisBsplineMz::synthesis[" << k << "]" << endl;
#endif

		fSub.mul(as_[k], xSub, accumulate, false, false);
	}*/
}


void BasisBsplineMz::analysis(MatrixSparse& xE, const MatrixSparse& fE, bool sqrA) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineMz::analysis[" << 0 << "]" << endl;
#endif

	if (sqrA)
	{
		MatrixSparse t;
		t.copy(as_[0]);
		t.elementwiseSqr();
		xE.mul(MatrixSparse::Accumulate::NO, fE, MatrixSparse::Transpose::NO, t);
	}
	else
	{
		xE.mul(MatrixSparse::Accumulate::NO, fE, MatrixSparse::Transpose::NO, as_[0]);
	}

	/*if (!xE)
	{
		xE.init(getGridInfo().n, getGridInfo().m());
	}

#ifdef NDEBUG
	# pragma omp parallel for
#endif
	for (ii k = 0; k < (ii)as_.size(); k++)
	{
		Matrix xESub; xESub.init(xE, k * as_[k].n(), 0, as_[k].n(), 1);
		Matrix fESub; fESub.init(fE, is_[k], 0, as_[k].m(), 1);

#ifndef NDEBUG
		cout << " " << getIndex() << " BasisBsplineMz::analysis[" << k << "]" << endl;
#endif

		if (sqrA)
		{
			MatrixSparse aSqrd;
			aSqrd.elementwiseSqr(as_[k]);
			xESub.mul(aSqrd, fESub, false, true, false);
		}
		else
		{
			xESub.mul(as_[k], fESub, false, true, false);
		}
	}*/
}


void BasisBsplineMz::deleteRows(const MatrixSparse& x, ii threshold)
{
    if(nnzRows_ - x.nnz() >= threshold)
    {
        // delete rows in aTs we don't need anymore
        aTs_[0].deleteRows(x);
        as_[0].copy(aTs_[0], MatrixSparse::Transpose::YES);
        
        nnzRows_ = x.nnz();
    }
}



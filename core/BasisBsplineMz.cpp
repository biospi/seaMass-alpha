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
                               const std::vector<double>& binEdges, short scale, Transient transient, int order)
	: BasisBspline(bases, 1, transient), nnzBasisFunctions_(0)
{
    if (getDebugLevel() % 10 >= 2)
    {
        cout << getTimeStamp() << " " << getIndex() << " BasisBsplineMz";
        if (getTransient() == Transient::YES) cout << " (transient)";
        cout << endl;
    }
    
	if (spectrumIndex.size() > 0)
	{
		js_ = spectrumIndex;
	}
	else
	{
		js_.push_back(0);
		js_.push_back(binCounts.size());
	}

	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	// find min and max m/z across spectra, and m for each A
	double mzMin = numeric_limits<double>::max();
	double mzMax = 0.0;
	double mzDiff = numeric_limits<double>::max();
	for (ii j = 0; j < (ii)js_.size() - 1; j++)
	{
		mzMin = binEdges[js_[j] + j] < mzMin ? binEdges[js_[j] + j] : mzMin;
		mzMax = binEdges[js_[j + 1] + j] > mzMax ? binEdges[js_[j + 1] + j] : mzMax;

		for (ii i = js_[j]; i < js_[j + 1]; i++)
		{
			double diff = binEdges[i + j + 1] - binEdges[i + j];
			mzDiff = diff < mzDiff ? diff : mzDiff;
		}
	}
    
	ii scaleAuto = (ii) floor(log2(1.0 / mzDiff / 60.0 / 1.0033548378));
	if (scale == numeric_limits<short>::max())
	{
		scale = scaleAuto;
        
        if (getDebugLevel() % 10 >= 1)
        {
            cout << getTimeStamp() << "   autodetected_mz_scale=" << fixed << setprecision(1) << scale << endl;
        }
	}

	// Bases per 1.0033548378Th (difference between carbon12 and carbon13)
	double bpi = pow(2.0, (double)scale) * 60 / 1.0033548378;
    
    // fill in b-spline grid info
    gridInfo().count = (ii)js_.size() - 1;
    gridInfo().scale[0] = scale;
    gridInfo().offset[0] = (ii)floor(mzMin * bpi);
    gridInfo().extent[0] = ((ii)ceil(mzMax * bpi)) + order - gridInfo().offset[0];
    
    if (getDebugLevel() % 10 >= 2)
    {
        cout << getTimeStamp() << "   range=" << fixed << setprecision(3) << mzMin << ":"; cout.unsetf(std::ios::floatfield); cout << mzDiff << ":" << fixed << mzMax << "Th" << endl;
        cout << getTimeStamp() << "   scale=" << fixed << setprecision(1) << scale << " (" << bpi << " bases per 1.0033548378Th)" << endl;
        cout << getTimeStamp() << "   " << gridInfo() << endl;
    }
   
	// populate coo matrix
    vector<fp> acoo;
    vector<ii> rowind;
    vector<ii> colind;
    
	ii done = 0;
	Bspline bspline(order, 65536); // bspline basis function lookup table
	for (ii j = 0; j < (ii)js_.size() - 1; j++)
	{
        for (ii i = js_[j]; i < js_[j + 1]; i++)
		{
			if (binCounts[i] >= 0.0)
			{
				double xfMin = binEdges[i + j] * bpi;
				double xfMax = binEdges[i + j + 1] * bpi;

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
					colind.push_back(j * gridInfo().extent[0] + (x - gridInfo().offset[0]));
				}
			}
		}

		// display progress update
        done++;
        if (done % 100 == 0)
        {
            for (int i = 0; i < 256; ++i) cout << '\b';
            cout << getIndex() << " BasisBsplineMz " << setw(1 + (int)(log10((float)js_.size() - 1))) << done << "/" << (js_.size() - 1) << " " << flush;
        }
	}
	for (int i = 0; i < 256; ++i) cout << '\b';

    
    // create A
    a_ = new MatrixSparse();
    a_->init(js_[js_.size() - 1], getGridInfo().n() * getGridInfo().count, (ii)acoo.size(), acoo.data(), rowind.data(), colind.data());
    aT_ = new MatrixSparse();
    aT_->copy(*a_, MatrixSparse::Operation::TRANSPOSE);
    nnzBasisFunctions_ = getGridInfo().n() * getGridInfo().count;
    
	if (scaleAuto != scale)
	{
		cerr << "WARNING: scale is not the suggested value of " << scaleAuto << ". Continue at your own risk!" << endl;
	}
}


BasisBsplineMz::~BasisBsplineMz()
{
    delete a_;
    delete aT_;
}


void BasisBsplineMz::synthesis(MatrixSparse& f, const MatrixSparse& x, bool accumulate) const
{
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "   " << getIndex() << " BasisBsplineMz::synthesis" << endl;
    }
    
    MatrixSparse t;
    t.copy(x, MatrixSparse::Operation::UNPACK_ROWS);
    
	f.matmul(false, t, *aT_, accumulate);
    
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "   " << getIndex() << "   " << f << endl;
    }
}


void BasisBsplineMz::analysis(MatrixSparse& xE, const MatrixSparse& fE, bool sqrA) const
{
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "   " << getIndex() << " BasisBsplineMz::analysis" << endl;
    }

    MatrixSparse t;
    
	if (sqrA)
	{
		MatrixSparse aSqr;
		aSqr.copy(*a_);
		aSqr.sqr();
		t.matmul(false, fE, aSqr, false);
	}
	else
	{
		t.matmul(false, fE, *a_, false);
	}
    
    xE.copy(t, MatrixSparse::Operation::PACK_ROWS);
    
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "   " << getIndex() << "   " << xE << endl;
    }
}


// delete basis functions we don't need anymore
void BasisBsplineMz::deleteBasisFunctions(const MatrixSparse& x, ii threshold)
{
    
    
    
    /*if(x.nnz() / (double) nnzBasisFunctions_ <= 0.5)
    {
        cout << "deleting " << nnzBasisFunctions_ - x.nnz() << " basis functions" << endl;
        
        delete a_;
        
        MatrixSparse* aT = new MatrixSparse();
        aT->zeroRowsOfZeroColumns(*aT_, x);
        delete aT_;
        aT_ = aT;
        
        a_ = new MatrixSparse();
        a_->copy(*aT_, MatrixSparse::Operation::TRANSPOSE);
        
        nnzBasisFunctions_ = x.nnz();
    }*/
}



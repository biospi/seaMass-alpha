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
#include <assert.h>


using namespace std;
using namespace kernel;


BasisBsplineMz::BasisBsplineMz(std::vector<Basis*>& bases, const std::vector<fp>& binCounts, const std::vector<li>& binCountsIndex,
                               const std::vector<double>& binEdges, char scale, Transient transient, int order)
	: BasisBspline(bases, 1, transient)
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp();
        if (getDebugLevel() % 10 >= 2)
            cout << "   " << getIndex() << " BasisBsplineMz";
        else
            cout << "   BasisBsplineMz";
        if (getTransient() == Transient::YES) cout << " (transient)";
        cout << " ..." << endl;
    }
    
    std::vector<li> bci;
    std::vector<li> bei;
	if (binCountsIndex.size() > 0)
	{
		bci = binCountsIndex;
        bei = binCountsIndex;
        for (ii i = 0; i < bei.size(); i++) bei[i] += i;
	}
	else
	{
		bci.push_back(0);
		bci.push_back(binCounts.size());
        bei.push_back(0);
        bei.push_back(binEdges.size());
	}

	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	// find min and max m/z across spectra, and m for each A
	double mzMin = numeric_limits<double>::max();
	double mzMax = 0.0;
	double mzDiff = numeric_limits<double>::max();
	for (ii k = 0; k < (ii)bei.size() - 1; k++)
	{
		mzMin = binEdges[bei[k]] < mzMin ? binEdges[bei[k]] : mzMin;
		mzMax = binEdges[bei[k + 1] - 1] > mzMax ? binEdges[bei[k + 1] - 1] : mzMax;

 		for (ii i = bei[k]; i < bei[k + 1] - 1; i++)
		{
			double diff = binEdges[i + 1] - binEdges[i];
 			mzDiff = diff < mzDiff ? diff : mzDiff;
		}
	}
    
	ii scaleAuto = (ii) floor(log2(1.0 / mzDiff / 60.0 / 1.0033548378));
	if (scale == numeric_limits<char>::max())
	{
		scale = scaleAuto;
        
        if (getDebugLevel() % 10 >= 1)
        {
            cout << getTimeStamp() << "     autodetected_mz_scale=" << fixed << setprecision(1) << (int) scale << endl;
        }
	}

	// Bases per 1.0033548378Th (difference between carbon12 and carbon13)
	double bpi = pow(2.0, (double)scale) * 60 / 1.0033548378;
    
    // fill in b-spline grid info
    gridInfo().count = (ii)bei.size() - 1;
    gridInfo().scale[0] = scale;
    gridInfo().offset[0] = (ii)floor(mzMin * bpi);
    gridInfo().extent[0] = ((ii)ceil(mzMax * bpi)) + order - gridInfo().offset[0];
    
    if (getDebugLevel() % 10 >= 2)
    {
        cout << getTimeStamp() << "     range=" << fixed << setprecision(3) << mzMin << ":"; cout.unsetf(std::ios::floatfield);
        cout << mzDiff << ":" << fixed << mzMax << "Th" << endl;
        cout << getTimeStamp() << "     scale=" << fixed << setprecision(1) << (int) scale << " (" << bpi << " bases per 1.0033548378Th)" << endl;
        cout << getTimeStamp() << "     " << gridInfo() << endl;
    }
   
    aTs_.resize(bei.size() - 1);
    aTnnzRows_.resize(bei.size() - 1);
    
    as_.resize(bei.size() - 1);
    
	Bspline bspline(order, 65536); // bspline basis function lookup table
	for (ii k = 0; k < (ii)bei.size() - 1; k++)
	{
        vector<ii> is;
        vector<ii> js;
        vector<fp> vs;

        for (ii i = bei[k]; i < bei[k + 1] - 1; i++)
		{
			if (binCounts[i - bei[k] + bci[k]] >= 0.0)
			{
				double xfMin = binEdges[i] * bpi;
				double xfMax = binEdges[i + 1] * bpi;

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

					vs.push_back(b);
					is.push_back(i - bei[k]);
					js.push_back(x - gridInfo().offset[0]);
				}
			}
		}
        
        // create A
        aTs_[k].copy(getGridInfo().n(), bci[k + 1] - bci[k], js, is, vs);
        aTnnzRows_[k] = getGridInfo().n();
        
        as_[k].copy(aTs_[k], true);

		// display progress update
        if (getDebugLevel() % 10 >= 1)
        {
            if ((k + 1) % 100 == 0 || k == (ii)bei.size() - 2)
            {
                cout << getTimeStamp() << "     " << setw(1 + (int)(log10((float)bci.size() - 1))) << (k + 1) << "/" << (bci.size() - 1) << endl;
            }
        }
	}

	if (scaleAuto != scale)
	{
		cerr << "WARNING: mz_scale is not the suggested value of " << scaleAuto << ". Continue at your own risk!" << endl;
	}
}


BasisBsplineMz::~BasisBsplineMz()
{
}


void BasisBsplineMz::synthesis(vector<MatrixSparse>& f, const vector<MatrixSparse>& x, bool accumulate) const
{
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "     " << getIndex() << " BasisBsplineMz::synthesis" << endl;
    }
    
    if (!f.size()) f.resize(aTs_.size());
    for (size_t k = 0; k < f.size(); k++)
    {
        MatrixSparse t;
        t.init(x[0], k);
        
        f[k].matmul(false, t, aTs_[k], accumulate, true);
        
        if (getDebugLevel() % 10 >= 3)
        {
            cout << getTimeStamp() << "       " << f[k] << endl;
        }
    }
}


void BasisBsplineMz::analysis(vector<MatrixSparse>& xE, const vector<MatrixSparse>& fE, bool sqrA) const
{
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "     " << getIndex() << " BasisBsplineMz::analysis" << endl;
    }
    
    vector<MatrixSparse> xEs(aTs_.size());
    for (size_t k = 0; k < aTs_.size(); k++)
    {
        if (sqrA)
        {
            MatrixSparse aSqr;
            aSqr.copy(as_[k]);
            aSqr.sqr();
            xEs[k].matmul(false, fE[k], aSqr, false);
        }
        else
        {
            xEs[k].matmul(false, fE[k], as_[k], false);
        }
    }
    
    xE.resize(1);
    xE[0].copy(xEs);
    
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "       " << xE[0] << endl;
    }
}


// delete basis functions we don't need anymore
void BasisBsplineMz::deleteBasisFunctions(const vector<MatrixSparse>& x, fp threshold)
{
    for (size_t k = 0; k < x.size(); k++)
    {
        ii aTnnzRows = aTs_[k].pruneRows(aTs_[k], aTnnzRows_[k], x[k], false, threshold);
        
        if (aTnnzRows < aTnnzRows_[k])
        {
            as_[k].copy(aTs_[k], true);
            
            if (getDebugLevel() % 10 >= 3)
            {
                cout << getTimeStamp() << "     " << getIndex() << " BasisBsplineMz::deleteBasisFunctions " << aTnnzRows_[k] - aTnnzRows << endl;
            }
            
            aTnnzRows_[k] = aTnnzRows;
        }
    }
}



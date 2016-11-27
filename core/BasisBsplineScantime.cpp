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


#include "BasisBsplineScantime.hpp"

#include "Bspline.hpp"

#include <limits>
#include <iomanip>
#include <cmath>


using namespace std;


BasisBsplineScantime::BasisBsplineScantime(std::vector<Basis*>& bases, ii parentIndex, const std::vector<double>& startTimes, const std::vector<double>& finishTimes,
	const std::vector<fp>& exposures, short resolution, ii order, bool transient) : BasisBspline(bases, 2, transient, parentIndex)
{
	double scantimeMin = startTimes.front();
	double scantimeMax = finishTimes.back();

	double scantimeDiff = numeric_limits<double>::max();
	for (ii j = 0; j < (ii)startTimes.size() - 1; j++)
	{
		double diff = 0.5 * (startTimes[j + 1] + finishTimes[j + 1]) - 0.5 * (startTimes[j] + finishTimes[j]);
		scantimeDiff = diff < scantimeDiff ? diff : scantimeDiff;
	}

	ii resolutionAuto = (ii)floor(log2(1.0 / scantimeDiff));
	if (resolution == numeric_limits<short>::max())
	{
		resolution = resolutionAuto;
	}

	// Bases per second
	double bpi = pow(2.0, (double)resolution);

	// fill in b-spline grid info
	const GridInfo parentGridInfo = static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo();
	gridInfo().n = 1;
	gridInfo().scale[0] = parentGridInfo.scale[0];
	gridInfo().offset[0] = parentGridInfo.offset[0];
	gridInfo().extent[0] = parentGridInfo.extent[0];
	gridInfo().scale[1] = resolution;
	gridInfo().offset[1] = (ii)floor(scantimeMin * bpi);
	gridInfo().extent[1] = ((ii)ceil(scantimeMax * bpi)) + order - gridInfo().offset[1];
	ii m = parentGridInfo.n;
	ii n = gridInfo().extent[1];
    
    // populate coo matrix
    vector<fp> acoo;
    vector<ii> rowind;
    vector<ii> colind;
	Bspline bspline(order, 65536); // bspline basis function lookup table
	for (ii i = 0; i < startTimes.size(); i++)
    {
		double xfMin = startTimes[i] * bpi;
		double xfMax = finishTimes[i] * bpi;

        ii xMin = (ii) floor(xfMin);
        ii xMax = ((ii) ceil(xfMax)) + order;
        
		// work out basis coefficients
		for (int x = xMin; x < xMax; x++)
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
			colind.push_back(x - gridInfo().offset[1]);
        }
    }

	// create 'a' and 'aT'
	a_.init(m, n, (ii)acoo.size(), acoo.data(), rowind.data(), colind.data());
	aT_.init(n, m, (ii)acoo.size(), acoo.data(), colind.data(), rowind.data());

#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineScantime";
	if (isTransient()) cout << " (t)";
	cout << " parent=" << getParentIndex();
	cout << " range=" << setprecision(3) << fixed << scantimeMin << ":";
	cout.unsetf(ios::floatfield); 
	cout << scantimeDiff << ":" << fixed << scantimeMax << "s";
	cout << " resolution=" << fixed << setprecision(1) << resolution << " (" << bpi << " bases per second)";
	cout << " " << gridInfo() << endl;
	cout << "  A" << a_ << " (";
	cout.unsetf(ios::floatfield); 
	cout << setprecision(2) << (a_.mem() + aT_.mem()) / 1024.0 / 1024.0 << "Mb)" << endl;
#endif
	
	if (resolutionAuto != resolution)
	{
		cerr << endl << "WARNING: resolution is not the suggested value of " << resolutionAuto << ". Continue at your own risk!" << endl << endl;
	}
}


BasisBsplineScantime::~BasisBsplineScantime()
{
}


void
BasisBsplineScantime::
synthesis(Matrix& f, const Matrix& x, bool accumulate) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineScantime::synthesis" << endl;
#endif

	f.mul(a_, x, accumulate, false, false);
}


void
BasisBsplineScantime::
analysis(Matrix& xE, const Matrix& fE, bool sqrA) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineScantime::analysis" << endl;
#endif

	if (sqrA)
	{
		MatrixSparse aTSqrd;
		aTSqrd.elementwiseSqr(aT_);
		xE.mul(aTSqrd, fE, false, false, false);
	}
	else
	{
		xE.mul(aT_, fE, false, false, false);
	}
}
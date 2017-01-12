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


// todo: support ion mobility
BasisBsplineScantime::BasisBsplineScantime(std::vector<Basis*>& bases, ii parentIndex, const std::vector<double>& startTimes, const std::vector<double>& finishTimes,
                                           const std::vector<fp>& exposures, short scale, Transient transient, ii order) : BasisBspline(bases, 2, transient, parentIndex)
{
#ifndef NDEBUG
    cout << " " << getIndex() << " BasisBsplineScantime";
    if (getTransient() == Basis::Transient::YES) cout << " (transient)";
    cout << endl;
#endif
    
	double scantimeMin = startTimes.front();
	double scantimeMax = finishTimes.back();

	double scantimeDiff = numeric_limits<double>::max();
	for (ii j = 0; j < (ii)startTimes.size() - 1; j++)
	{
		double diff = 0.5 * (startTimes[j + 1] + finishTimes[j + 1]) - 0.5 * (startTimes[j] + finishTimes[j]);
		scantimeDiff = diff < scantimeDiff ? diff : scantimeDiff;
	}

	ii scaleAuto = (ii)floor(log2(1.0 / scantimeDiff));
	if (scale == numeric_limits<short>::max())
	{
		scale = scaleAuto;
		cout << "Autodetected st_scale=" << scale << endl;
	}

	// Bases per second
	double bpi = pow(2.0, (double)scale);

	// fill in b-spline grid info
	const GridInfo parentGridInfo = static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo();
	gridInfo().count = 1;
    gridInfo().scale[0] = parentGridInfo.scale[0];
    gridInfo().offset[0] = parentGridInfo.offset[0];
    gridInfo().extent[0] = parentGridInfo.extent[0];
	gridInfo().scale[1] = scale;
	gridInfo().offset[1] = (ii)floor(scantimeMin * bpi);
	gridInfo().extent[1] = ((ii)ceil(scantimeMax * bpi)) + order - gridInfo().offset[1];
    
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

	// create transformation matrix 'a' and its transpose 'aT'
	a_.init(parentGridInfo.m(), getGridInfo().m(), (ii)acoo.size(), acoo.data(), rowind.data(), colind.data());
	aT_.init(getGridInfo().m(), parentGridInfo.m(), (ii)acoo.size(), acoo.data(), colind.data(), rowind.data());
 
#ifndef NDEBUG
    cout << "  parent=" << getParentIndex() << " range=" << fixed << setprecision(3) << scantimeMin << ":";
    cout.unsetf(std::ios::floatfield);
    cout << scantimeDiff << ":" << fixed << scantimeMax << "Th";
    cout << " scale=" << fixed << setprecision(1) << scale << " (" << bpi << " bases per second) ";
    cout.unsetf(ios::floatfield);
    cout << gridInfo() << " mem=";
    cout << fixed << setprecision(2) << (a_.mem() + aT_.mem()) / 1024.0 / 1024.0 << "Mb" << endl;
#endif

	if (scaleAuto != scale)
	{
		cerr << endl << "WARNING: st_scale is not the suggested value of " << scaleAuto << ". Continue at your own risk!" << endl << endl;
	}
}


BasisBsplineScantime::~BasisBsplineScantime()
{
}


void BasisBsplineScantime::synthesis(MatrixSparse& f, const MatrixSparse& x, bool accumulate) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineScantime::synthesis" << endl;
#endif

    f.mul(true, x, aT_, accumulate, true);
}


void BasisBsplineScantime::analysis(MatrixSparse& xE, const MatrixSparse& fE, bool sqrA) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineScantime::analysis" << endl;
#endif
    
    if (sqrA)
    {
        MatrixSparse t;
        t.copy(a_);
        t.elementwiseSqr();
        xE.mul(true, fE, t, false, true);
    }
    else
    {
        xE.mul(true, fE, a_, false, true);
    }
}


void BasisBsplineScantime::deleteRows(const MatrixSparse& x, ii threshold)
{
    /*if(nnzRows_ - x.nnz() >= threshold)
     {
     // delete rows in aTs we don't need anymore
     aT_.deleteRows(x);
     a_.copy(aT_, MatrixSparse::Operation::TRANSPOSE);
     
     nnzRows_ = x.nnz();
     }*/
}

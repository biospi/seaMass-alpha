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
using namespace kernel;


// TODO: support ion mobility
BasisBsplineScantime::BasisBsplineScantime(std::vector<Basis*>& bases, ii parentIndex, const std::vector<double>& startTimes, const std::vector<double>& finishTimes,
                                           const std::vector<fp>& exposures, char scale, bool transient, ii order) : BasisBspline(bases, 2, transient, parentIndex)
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp();
        if (getDebugLevel() % 10 >= 2)
            cout << "   " << getIndex() << " BasisBsplineScantime";
        else
            cout << "   BasisBsplineScantime";
        if (isTransient()) cout << " (transient)";
        cout << " ..." << endl;
    }
    
    double scantimeMin = startTimes.front();
    double scantimeMax = finishTimes.back();

    double scantimeDiff = numeric_limits<double>::max();
    for (ii j = 0; j < (ii)startTimes.size() - 1; j++)
    {
        double diff = 0.5 * (startTimes[j + 1] + finishTimes[j + 1]) - 0.5 * (startTimes[j] + finishTimes[j]);
        scantimeDiff = diff < scantimeDiff ? diff : scantimeDiff;
    }

    ii scaleAuto = (ii)floor(log2(1.0 / scantimeDiff));
    if (scale == numeric_limits<char>::max())
    {
        scale = scaleAuto;
        
        if (getDebugLevel() % 10 >= 1)
            cout << getTimeStamp() << "     autodetected_st_scale=" << fixed << setprecision(1) << (int) scale << endl;
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
    
    if (getDebugLevel() % 10 >= 2)
    {
        cout << getTimeStamp() << "     parent=" << getParentIndex() << endl;
        cout << getTimeStamp() << "     range=" << fixed << setprecision(3) << scantimeMin << ":";
        cout.unsetf(std::ios::floatfield);
        cout << scantimeDiff << ":" << fixed << scantimeMax << "seconds" << endl;
        cout << getTimeStamp() << "     scale=" << fixed << setprecision(1) << (int) scale << " (" << bpi << " bases per second)" << endl;
        cout << getTimeStamp() << "     " << gridInfo() << endl;
    }

    // populate coo matrix
    vector<ii> rowind;
    vector<ii> colind;
    vector<fp> acoo;
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

            rowind.push_back(i);
            colind.push_back(x - gridInfo().offset[1]);
            acoo.push_back(b);
        }
    }

    // create transformation matrix 'a'
    aT_.copy(getGridInfo().m(), parentGridInfo.m(), acoo.size(), colind.data(), rowind.data(), acoo.data());

    if (scaleAuto != scale)
        cerr << "WARNING: st_scale is not the suggested value of " << scaleAuto << ". Continue at your own risk!" << endl;
}


BasisBsplineScantime::~BasisBsplineScantime()
{
}


void BasisBsplineScantime::synthesize(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "     " << getIndex() << " BasisBsplineScantime::synthesise" << endl;

    if (!f.size())
        f.resize(1);

    // zero basis functions that are no longer needed
    MatrixSparse t;
    ii rowsPruned = t.copyPruneRows(aT_, x[0], true, 0.75);
    if (rowsPruned > 0)
    {
        aT_.swap(t);

        if (getDebugLevel() % 10 >= 2)
            cout << getTimeStamp() << "      " << getIndex() << " pruned " << rowsPruned << " basis functions" << endl;
    }

    // synthesise
    f[0].matmul(true, aT_, x[0], accumulate);
        
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << f[0] << endl;
}


void BasisBsplineScantime::analyze(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA)
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "     " << getIndex() << " BasisBsplineScantime::analyse" << endl;

    if (!xE.size())
        xE.resize(1);

    if (sqrA)
    {
        MatrixSparse t;
        t.copy(aT_);
        t.sqr();
        xE[0].matmul(false, t, fE[0], false);
    }
    else
    {
        xE[0].matmul(false, aT_, fE[0], false);
    }
    
    if (getDebugLevel() % 10 >= 3)
         cout << getTimeStamp() << "       " << xE[0] << endl;
}


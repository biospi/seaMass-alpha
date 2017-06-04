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
#include <sstream>
#include <cmath>
using namespace std;
using namespace kernel;


// TODO: support ion mobility
BasisBsplineScantime::BasisBsplineScantime(std::vector<Basis*>& bases, ii parentIndex,
                                           const std::vector<double>& startTimes,
                                           const std::vector<double>& finishTimes,
                                           const std::vector<fp>& exposures,
                                           char scale, bool transient)
        : BasisBspline(bases, 2, transient, parentIndex)
{
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp();
        if (getDebugLevel() % 10 >= 2)
            oss << "   " << getIndex() << " BasisBsplineScantime";
        else
            oss << "   BasisBsplineScantime";
        if (isTransient()) oss << " (transient)";
        oss << " ...";
        info(oss.str());
    }

    double scantimeMin = startTimes.front();
    double scantimeMax = finishTimes.back();

    double scantimeDiff = 0.0;
    for (ii j = 0; j < ii(startTimes.size()) - 1; j++)
        scantimeDiff += 0.5 * (startTimes[j + 1] + finishTimes[j + 1]) - 0.5 * (startTimes[j] + finishTimes[j]);
    scantimeDiff /= ii(startTimes.size()) - 1;

    ii scaleAuto = ii(ceil(log2(1.0 / scantimeDiff))) + 1;
    if (scale == numeric_limits<char>::max())
    {
        scale = scaleAuto;
        
        if (getDebugLevel() % 10 >= 1)
        {
            ostringstream oss;
            oss << getTimeStamp() << "     autodetected_st_scale=" << fixed << setprecision(1) << int(scale);
            info(oss.str());
        }
    }
    
    // fill in b-spline grid info
    const GridInfo parentGridInfo = static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo();
    gridInfo().count = 1;
    gridInfo().scale[0] = parentGridInfo.scale[0];
    gridInfo().offset[0] = parentGridInfo.offset[0];
    gridInfo().extent[0] = parentGridInfo.extent[0];
    gridInfo().scale[1] = scale;
    gridInfo().offset[1] = ii(floor(scantimeMin * (1L << scale))) - 1;
    gridInfo().extent[1] = (ii(ceil(scantimeMax * (1L << scale))) + 1)
                           - gridInfo().offset[1] + 1;
    
    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     parent=" << getParentIndex();
        info(oss.str());
        ostringstream oss2;
        oss2 << getTimeStamp() << "     range=" << fixed << setprecision(3) << scantimeMin << ":";
        oss2.unsetf(std::ios::floatfield);
        oss2 << scantimeDiff << ":" << fixed << scantimeMax << "seconds";
        info(oss2.str());
        ostringstream oss3;
        oss3 << getTimeStamp() << "     " << gridInfo();
        info(oss3.str());
    }

    // populate coo matrix
    vector<ii> rowind;
    vector<ii> colind;
    vector<fp> acoo;
    Bspline bspline(3, 65536); // bspline basis function lookup table
    for (ii i = 0; i < startTimes.size(); i++)
    {
        double xfMin = startTimes[i] * (1L << scale);
        double xfMax = finishTimes[i] * (1L << scale);

        ii xMin = ii(floor(xfMin)) - 1;
        ii xMax = ii(ceil(xfMax)) + 1;
        
        // work out basis coefficients
        for (int x = xMin; x <= xMax; x++)
        {
            double bfMin = (double)(x - 2);
            double bfMax = (double)(x + 2);
            
            // intersection of bin and basis, between 0.0 and 4.0
            double bMin = xfMin > bfMin ? xfMin - bfMin : 0.0;
            double bMax = xfMax < bfMax ? xfMax - bfMin : bfMax - bfMin;

            // basis coefficient b is _integral_ of area under b-spline basis
            fp b = exposures[i] * fp(bspline.ibasis(bMax) - bspline.ibasis(bMin));
            if (b > 0.0)
            {
                acoo.push_back(b);
                rowind.push_back(i);
                colind.push_back(x - gridInfo().offset[1]);
            }
       }
    }

    // create transformation matrix 'a'
    aT_.importFromCoo(getGridInfo().m(), parentGridInfo.m(), acoo.size(), colind.data(), rowind.data(), acoo.data());

    if (scaleAuto != scale && getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << "WARNING: st_scale is not the suggested value of " << scaleAuto << ". Continue at your own risk!";
        warning(oss.str());
    }
}


BasisBsplineScantime::~BasisBsplineScantime()
{
}


void BasisBsplineScantime::synthesize(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineScantime::synthesise";
        info(oss.str());
    }

    if (!f.size())
        f.resize(1);

    // zero basis functions that are no longer needed
    MatrixSparse t;
    ii rowsPruned = t.pruneRows(aT_, x[0], true, 0.75);
    if (rowsPruned > 0)
    {
        aT_.swap(t);

        if (getDebugLevel() % 10 >= 2)
        {
            ostringstream oss;
            oss << getTimeStamp() << "      " << getIndex() << " pruned " << rowsPruned << " basis functions";
            info(oss.str());
        }
    }

    // synthesise
    f[0].matmul(true, aT_, x[0], accumulate);
        
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << f[0];
        info(oss.str());
    }
}


void BasisBsplineScantime::analyze(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineScantime::analyse";
        info(oss.str());
    }

    if (!xE.size())
        xE.resize(1);

    if (sqrA)
    {
        MatrixSparse t;
        t.sqr(aT_);
        xE[0].matmul(false, t, fE[0], false);
    }
    else
    {
        xE[0].matmul(false, aT_, fE[0], false);
    }
    
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << xE[0];
        info(oss.str());
    }
}


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


#include "BasisBsplineCharge.hpp"
#include "Bspline.hpp"
#include <limits>
#include <iomanip>
#include <sstream>
#include <cmath>
using namespace std;
using namespace kernel;


BasisBsplineCharge::BasisBsplineCharge(std::vector<Basis*>& bases, ii parentIndex, bool transient)
        : BasisBspline(bases, 2, transient, parentIndex)
{
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp();
        if (getDebugLevel() % 10 >= 2)
            oss << "   " << getIndex() << " BasisBsplineCharge";
        else
            oss << "   BasisBsplineCharge";
        if (isTransient()) oss << " (transient)";
        oss << " ...";
        info(oss.str());
    }

    // fill in b-spline grid info
    vector<fp> hs(4);
    hs[0] = Bspline::im(1.0, 4) - Bspline::im(0.0, 4);
    hs[1] = Bspline::im(2.0, 4) - Bspline::im(1.0, 4);
    hs[2] = Bspline::im(3.0, 4) - Bspline::im(2.0, 4);
    hs[3] = Bspline::im(4.0, 4) - Bspline::im(3.0, 4);

    // fill in b-spline grid info
    const GridInfo parentGridInfo = static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo();
    gridInfo() = parentGridInfo;
    gridInfo().colExtent[1] += ii(hs.size() - 1);
    
    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     parent=" << getParentIndex();
        info(oss.str());
        ostringstream oss3;
        oss3 << getTimeStamp() << "     " << gridInfo();
        info(oss3.str());
    }

    // create A as a temporary COO matrix
    ii m = parentGridInfo.colExtent[1];
    ii n = gridInfo().colExtent[1];
    vector<ii> is;
    vector<ii> js;
    vector<fp> vs;

    for (ii j = 0; j < n; j++)
    {
        for (ii k = 0; k < ii(hs.size()); k++)
        {
            ii i = j + k - ii(hs.size() - 1);
            if (i >= 0 && i < m)
            {
                is.push_back(i);
                js.push_back(j);
                vs.push_back(hs[k]);
            }
        }
    }

    // create A
    aT_.importFromCoo(n, m, vs.size(), js.data(), is.data(), vs.data());
}


BasisBsplineCharge::~BasisBsplineCharge()
{
}


void BasisBsplineCharge::synthesize(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineCharge::synthesise";
        info(oss.str());
    }

    if (!f.size())
        f.resize(1);

    // zero basis functions that are no longer needed
    /*MatrixSparse t;
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
    }*/

    // synthesise
    f[0].matmul(true, aT_, x[0], accumulate);
        
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << f[0];
        info(oss.str());
    }
}


void BasisBsplineCharge::analyze(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineCharge::analyse";
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


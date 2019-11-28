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


#include "BasisBsplinePeak.hpp"
#include "Bspline.hpp"
#include <limits>
#include <iomanip>
#include <cmath>
#include <sstream>
using namespace std;
using namespace kernel;


BasisBsplinePeak::BasisBsplinePeak(std::vector<Basis*>& bases, int parentIndex, double fwhm, bool transient) :
        BasisBspline(bases,
                     static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo().rowDimensions(),
                     static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo().colDimensions(),
                     transient, parentIndex)
{
    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp();
        if (getDebugLevel() % 10 >= 2)
            oss << "   " << getIndex() << " BasisBsplinePeak";
        else
            oss << "   BasisBsplinePeak";
        if (isTransient()) oss << " (transient)";
        oss << " ...";
        info(oss.str());
    }

    // create our kernel
    Bspline bspline(3, 65536); // bspline basis function lookup table
    ii nh = 2 * ii(ceil(2.0 * fwhm)) + 1;
    vector<fp> hs(nh, 0.0);
    for (ii k = 0; k < nh; k++)
    {
        ii l = k - (nh/2);

        double low = (l - 0.5) / fwhm;
        low = low > -2.0 ? low : -2.0;
        low = low < 2.0 ? low : 2.0;

        double high = (l + 0.5) / fwhm;
        high = high > -2.0 ? high : -2.0;
        high = high < 2.0 ? high : 2.0;

        hs[k] = Bspline::im(high + 2.0, 4) - Bspline::im(low + 2.0, 4);
    }

    const GridInfo parentGridInfo = static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo();
    gridInfo() = parentGridInfo;
    gridInfo().colOffset[0] = parentGridInfo.colOffset[0] - (nh/2);
    gridInfo().colExtent[0] = parentGridInfo.colExtent[0] + 2 * (nh/2);

    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     parent=" << getParentIndex();
        info(oss.str());

        for (ii k = 0; k < nh; k++)
        {
            ostringstream oss2;
            oss2 << getTimeStamp() << "     kernel=" << hs[k];
            info(oss2.str());
        }

        ostringstream oss3;
        oss3 << getTimeStamp() << "     " << gridInfo();
        info(oss3.str());
    }

    // create A as a temporary COO matrix
    ii m = parentGridInfo.colExtent[0];
    ii n = gridInfo().colExtent[0];
    vector<ii> is;
    vector<ii> js;
    vector<fp> vs;

    for (ii j = 0; j < n; j++)
    {
        for (ii k = 0; k < nh; k++)
        {
            ii i = j + k - 2 * (nh/2);
            if (i < 0 || i >= m) continue;

            if (hs[k] > 0.0)
            {
                is.push_back(i);
                js.push_back(j);
                vs.push_back(hs[k]);
            }
        }
    }

    // create A
    aT_.importFromCoo(n, m, vs.size(), js.data(), is.data(), vs.data());
    a_.transpose(aT_);
}


BasisBsplinePeak::~BasisBsplinePeak()
{
}


void
BasisBsplinePeak::
synthesize(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplinePeak::synthesise";
        info(oss.str());
    }

    if (!f.size())
        f.resize(1);

    // zero basis functions that are no longer needed
    /*MatrixSparse t;
    ii rowsPruned = t.pruneRows(aT_, x[0], false, 0.75);
    if (rowsPruned > 0)
    {
        aT_.swap(t);
        a_.transpose(aT_);

        if (getDebugLevel() % 10 >= 3)
        {
            ostringstream oss;
            oss << getTimeStamp() << "      " << getIndex() << " pruned " << rowsPruned << " basis functions";
            info(oss.str());
        }
    }*/

    // synthesise
    f[0].matmul(false, x[0], aT_, accumulate);

    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << f[0];
        info(oss.str());
    }
}


void BasisBsplinePeak::analyze(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplinePeak::analyse";
        info(oss.str());
    }

    if (!xE.size())
        xE.resize(1);

    if (sqrA)
    {
        MatrixSparse t;
        t.sqr(a_);

        xE[0].matmul(false, fE[0], t, false);
    }
    else
    {
        xE[0].matmul(false, fE[0], a_, false);
    }

    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << xE[0];
        info(oss.str());
    }
}




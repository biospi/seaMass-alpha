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
#include <sstream>
using namespace std;
using namespace kernel;


double BasisBsplineMz::PROTON_MASS = 1.007276466879;


BasisBsplineMz::BasisBsplineMz(std::vector<Basis*>& bases, const std::vector<fp>& binCounts,
                               const std::vector<li>& binCountsIndex_, const std::vector<double>& binEdges,
                               char scale, bool transient) : BasisBspline(bases, 1, transient)
{
    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp();
        if (getDebugLevel() % 10 >= 2 % 10 >= 2)
            oss << "   " << getIndex() << " BasisBsplineMz";
        else
            oss << "   BasisBsplineMz";
        if (isTransient()) oss << " (transient)";
        oss << " ...";
        info(oss.str());
    }
    
    std::vector<li> binCountsIndex;
    std::vector<li> binEdgesIndex;
    if (binCountsIndex_.size() > 0)
    {
        binCountsIndex = binCountsIndex_;
        binEdgesIndex = binCountsIndex_;
        for (ii i = 0; i < ii(binEdgesIndex.size()); i++) binEdgesIndex[i] += i;
    }
    else
    {
        binCountsIndex.push_back(0);
        binCountsIndex.push_back(binCounts.size());
        binEdgesIndex.push_back(0);
        binEdgesIndex.push_back(binEdges.size());
    }

    // find min and max m/z across spectra, m for each A
    double mzMin = numeric_limits<double>::max();
    double mzMax = 0.0;
    double xDiff = 0.0;
    li n = 0;
    for (ii k = 0; k < ii(binEdgesIndex.size()) - 1; k++)
    {
        mzMin = binEdges[binEdgesIndex[k]] < mzMin ? binEdges[binEdgesIndex[k]] : mzMin;
        mzMax = binEdges[binEdgesIndex[k + 1] - 1] > mzMax ? binEdges[binEdgesIndex[k + 1] - 1] : mzMax;

        // find mean difference in index between edges, ignoring first, last and zeros
        for (ii i = 1; i < binCountsIndex[k + 1] - binCountsIndex[k] - 1; i++)
        {
            if (binCounts[binCountsIndex[k] + i] != 0)
            {
                ii ei = binEdgesIndex[k] + i;
                xDiff += log2(binEdges[ei + 1] - PROTON_MASS) - log2(binEdges[ei] - PROTON_MASS);
                n++;
            }
        }
    }
    xDiff /= double(n);

    ii scaleAuto = ii(ceil(log2(1.0 / xDiff))) + 1;
    if (scale == numeric_limits<char>::max())
    {
        scale = scaleAuto;

        if (getDebugLevel() % 10 >= 1)
        {
            ostringstream oss;
            oss << getTimeStamp() << "     autodetected_mz_scale=" << fixed << setprecision(1) << int(scale);
            info(oss.str());
        }
    }

    // fill in b-spline grid info
    gridInfo().count = ii(binCountsIndex.size()) - 1;
    gridInfo().scale[0] = scale;
    gridInfo().offset[0] = ii(floor(log2(mzMin - PROTON_MASS) * (1L << scale))) - 1;
    gridInfo().extent[0] = (ii(ceil(log2(mzMax - PROTON_MASS) * (1L << scale))) + 1)
                           - gridInfo().offset[0] + 1;

    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     range=" << fixed << setprecision(3) << mzMin << ":" << mzMax << "Th";
        info(oss.str());
        ostringstream oss2;
        oss2 << getTimeStamp() << "     " << gridInfo();
        info(oss2.str());
    }

    // create A matrices
    Bspline bspline(3, 65536); // bspline basis function lookup table
    aTs_.resize(binCountsIndex.size() - 1);
    as_.resize(aTs_.size());
    for (ii k = 0; k < ii(aTs_.size()); k++)
    {
        vector<ii> rowind;
        vector<ii> colind;
        vector<fp> acoo;

        for (ii i = 0; i < binEdgesIndex[k + 1] - binEdgesIndex[k] - 1; i++)
        {
            if (binCounts[binCountsIndex[k] + i] >= 0.0)
            {
                ii ei = binEdgesIndex[k] + i;
                double xfMin = log2(binEdges[ei] - PROTON_MASS) * (1L << scale);
                double xfMax = log2(binEdges[ei + 1] - PROTON_MASS) * (1L << scale);

                ii xMin = ii(floor(xfMin)) - 1;
                ii xMax = ii(ceil(xfMax)) + 1;

                // work out basis coefficients
                for (ii x = xMin; x <= xMax; x++)
                {
                    double bfMin = double(x - 2);
                    double bfMax = double(x + 2);

                    // intersection of bin and basis, between 0.0 and 4.0
                    double bMin = xfMin > bfMin ? xfMin - bfMin : 0.0;
                    double bMax = xfMax < bfMax ? xfMax - bfMin : bfMax - bfMin;

                    // basis coefficient b is _integral_ of area under b-spline basis
                    fp b = fp(bspline.ibasis(bMax) - bspline.ibasis(bMin));
                    if (b > 0.0)
                    {
                        acoo.push_back(b);
                        rowind.push_back(i);
                        colind.push_back(x - gridInfo().offset[0]);
                    }
                }
            }
        }
        
        // create A
        aTs_[k].importFromCoo(getGridInfo().n(), binCountsIndex[k + 1] - binCountsIndex[k], acoo.size(),
                              colind.data(), rowind.data(), acoo.data());
        as_[k].transpose(aTs_[k]);

        // display progress update
        if (getDebugLevel() % 10 >= 2)
        {
            if ((k + 1) % 100 == 0 || k == (ii)binEdgesIndex.size() - 2)
            {
                ostringstream oss;
                oss << getTimeStamp() << "     " << setw(1 + int(log10((float)binCountsIndex.size() - 1)))
                    << (k + 1) << "/" << (binCountsIndex.size() - 1);
                info(oss.str());
            }
        }
    }

    if (scaleAuto != scale && getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << "WARNING: mz_scale is not the suggested value of " << scaleAuto << ". Continue at your own risk!";
        warning(oss.str());
    }
}


BasisBsplineMz::~BasisBsplineMz()
{
}


void BasisBsplineMz::synthesize(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineMz::synthesise";
        info(oss.str());
    }

    if (!f.size())
        f.resize(aTs_.size());

    for (ii k = 0; k < ii(aTs_.size()); k++)
    {
        MatrixSparseView row(x[0], k);

        // prune basis functions that are no longer needed
        MatrixSparse t;
        ii rowsPruned = t.pruneRows(aTs_[k], row, false, 0.75);
        if (rowsPruned > 0)
        {
            aTs_[k].swap(t);
            as_[k].transpose(aTs_[k]);

            if (getDebugLevel() % 10 >= 3)
            {
                ostringstream oss;
                oss << getTimeStamp() << "      " << getIndex() << " pruned " << rowsPruned << " basis functions";
                info(oss.str());
            }
        }

        // synthesise with dense result
        f[k].matmulDense(false, row, aTs_[k]);

        if (getDebugLevel() % 10 >= 3)
        {
            ostringstream oss;
            oss << getTimeStamp() << "       " << f[k];
            info(oss.str());
        }
    }
}


void BasisBsplineMz::analyze(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineMz::analyse";
        info(oss.str());
    }

    vector<MatrixSparse> xEs(aTs_.size());
    for (ii k = 0; k < ii(aTs_.size()); k++)
    {
        if (sqrA)
        {
            MatrixSparse aSqr;
            aSqr.sqr(as_[k]);
            xEs[k].matmul(false, fE[k], aSqr, false);
        }
        else
        {
            xEs[k].matmul(false, fE[k], as_[k], false);
        }
    }
    
    xE.resize(1);
    xE[0].concatenateRows(xEs);

    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << xE[0];
        info(oss.str());
    }
}




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


#include "BasisBsplineScale.hpp"
#include "Bspline.hpp"
#include <sstream>
#include <cmath>
using namespace std;
using namespace kernel;


BasisBsplineScale::
BasisBsplineScale(vector<Basis*>& bases, int parentIndex, short dimension0, short dimension1, bool group,
                  bool transient) :
        BasisBspline(bases,
                     static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo().rowDimensions(),
                     static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo().colDimensions(),
                     transient, parentIndex), dimension0_(dimension0), dimension1_(dimension1)
{
    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "   " << getIndex() << " BasisBsplineScale";
        if (isTransient()) oss << " (transient)";
        //oss;
        info(oss.str());
    }

    ii order = 3; // b-spline order
    ii count, m, n, offset;

    // todo: support stride for non-major dimension!!
    const GridInfo parentGridInfo = static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo();
    gridInfo() = parentGridInfo;
    if (dimension0_ == 0)
    {
        gridInfo().rowScale[dimension1_] = parentGridInfo.rowScale[dimension1_] - 1;
        gridInfo().rowOffset[dimension1_] = parentGridInfo.rowOffset[dimension1_] / 2;
        gridInfo().rowExtent[dimension1_] = (parentGridInfo.rowOffset[dimension1_] +
            parentGridInfo.rowExtent[dimension1_]) / 2 + 1 - gridInfo().rowOffset[dimension1_];

        count = 1;
        for (ii i = 0; i < dimension1_; i++)
            count *= parentGridInfo.rowExtent[i];

        m = parentGridInfo.rowExtent[dimension1_];
        n = gridInfo().rowExtent[dimension1_];

        offset = order + ((parentGridInfo.rowOffset[dimension1_] + 1) % 2);
    }
    else
    {
        gridInfo().colScale[dimension1_] = parentGridInfo.colScale[dimension1_] - 1;
        gridInfo().colOffset[dimension1_] = parentGridInfo.colOffset[dimension1_] / 2;
        gridInfo().colExtent[dimension1_] = (parentGridInfo.colOffset[dimension1_] +
            parentGridInfo.colExtent[dimension1_]) / 2 +1 - gridInfo().colOffset[dimension1_];

        count = 1;
        for (ii i = 0; i < dimension1_; i++)
            count *= parentGridInfo.colExtent[i];

        m = parentGridInfo.colExtent[dimension1_];
        n = gridInfo().colExtent[dimension1_];

        offset = order + ((parentGridInfo.colOffset[dimension1_] + 1) % 2);
    }

    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     parent=" << getParentIndex();
        info(oss.str());
        ostringstream oss2;
        oss2 << getTimeStamp() << "     dimension=" << dimension0_ << ":" << dimension1_;
        info(oss2.str());
        ostringstream oss3;
        oss3 << getTimeStamp() << "     " << gridInfo();
        info(oss3.str());
    }

    // create our kernel
    ii nh = order + 2;
    vector<fp> hs(nh);
    double sum = 0.0;
    for (ii k = 0; k < nh; k++)
    {
        hs[k] = fp(1.0 / pow(2.0, double(order)) * Bspline::factorial(order + 1)
                   / double(Bspline::factorial(k)* Bspline::factorial(order + 1 - k)));
        sum += hs[k];
    }
    for (ii i = 0; i < nh; i++)
        hs[i] /= (fp) sum;

    // create A as a temporary COO matrix
    vector<ii> is;
    vector<ii> js;
    vector<fp> vs;

    for (ii l = 0; l < count; l++)
    {
        for (ii j1 = 0; j1 < n; j1++)
        {
            for (ii k = 0; k < nh; k++)
            {
                ii i1 = 2 * j1 + k - offset;
                if (i1 < 0 || i1 >= m) continue;

                is.push_back(i1 + l*m);
                js.push_back(j1 + l*n);
                vs.push_back(hs[k]);
            }
        }
    }

    // create A
    aT_.importFromCoo(count * n, count * m, vs.size(), js.data(), is.data(), vs.data());

    if (dimension0)
        a_.transpose(aT_);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Gt = m x n matrix where m are the coefficients and n are the groups (monoisotope centroid mass).
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (dimension0 == 1 && group)
    {
        a_.transpose(aT_);

        vector<ii> is;
        vector<ii> js;
        vector<fp> vs;

        ii m = getGridInfo().colExtent[0] * getGridInfo().colExtent[1];
        ii n = getGridInfo().colExtent[1] +  ii(round(log2(double(getGridInfo().colExtent[0])) *
                                                              (1L << getGridInfo().colScale[1])));

        vector<ii> gSizes(n, 0);
        for (ii z = 0; z < gridInfo().colExtent[0]; z++)
        {
            auto g0 = ii(round(log2(double(z + 1)) * (1L << getGridInfo().colScale[1])));

            for (ii x = 0; x < gridInfo().colExtent[1]; x++)
            {
                ii g = g0 + x;
                gSizes[g]++;

                //double mass = pow(2.0, (gridInfo().colOffset[1] + g) / double(1L << gridInfo().colScale[1]));
                //cout << mass << endl;

                is.push_back(x + z * gridInfo().colExtent[1]);
                js.push_back(g);
                vs.push_back(1.0);
                //vs.push_back(1.0 / sqrt(mass); // this does not work
                //vs.push_back(1.0 / pow(300.0*mass, 1.0/4.0));
                //vs.push_back(1.0 / pow(6.0*mass, 1.0/3.0)); //vs.push_back(1.0 / sqrt(pow(6.0*mass, 2.0/3.0)));
            }
        }

        /*for (ii nz = 0; nz < ii(vs.size()); nz++)
        {
            vs[nz] /= sqrt(fp(gSizes[js[nz]]));
        }*/

        gTs_.resize(1);
        gs_.resize(1);

        gTs_[0].importFromCoo(m, n, vs.size(), is.data(), js.data(), vs.data());
        gs_[0].transpose(gTs_[0]);
    }
}


BasisBsplineScale::~BasisBsplineScale()
{
}


void
BasisBsplineScale::
synthesize(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineScale::synthesise";
        info(oss.str());
    }

    if (!f.size())
        f.resize(1);

    // zero basis functions that are no longer needed
    MatrixSparse t;
    ii rowsPruned = t.pruneRows(aT_, x[0], dimension0_ == 0, 0.75);
    if (rowsPruned > 0)
    {
        aT_.swap(t);

        if (dimension0_ == 1)
            a_.transpose(aT_);

        if (getDebugLevel() % 10 >= 3)
        {
            ostringstream oss;
            oss << getTimeStamp() << "      " << getIndex() << " pruned " << rowsPruned << " basis functions";
            info(oss.str());
        }
    }

    // synthesise
    if (dimension0_ == 0)
        f[0].matmul(true, aT_, x[0], accumulate);
    else
        f[0].matmul(false, x[0], aT_, accumulate);

    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << f[0];
        info(oss.str());
    }
}


void BasisBsplineScale::analyze(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineScale::analyse";
        info(oss.str());
    }

    if (!xE.size())
        xE.resize(1);

    if (sqrA)
    {
        if (dimension0_ == 0)
        {
            MatrixSparse t;
            t.sqr(aT_);

            xE[0].matmul(false, t, fE[0], false);
        }
        else
        {
            MatrixSparse t;
            t.sqr(a_);

            xE[0].matmul(false, fE[0], t, false);
        }

    }
    else
    {
        if (dimension0_ == 0)
            xE[0].matmul(false, aT_, fE[0], false);
        else
            xE[0].matmul(false, fE[0], a_, false);
    }

    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << xE[0];
        info(oss.str());
    }
}

const vector<MatrixSparse> * BasisBsplineScale::getColGroups(bool transpose) const
{
    if (dimension0_ == 1)
    {
        if (transpose)
            return &gTs_;
        else
            return &gs_;
    }
    else
    {
        return 0;
    }
}

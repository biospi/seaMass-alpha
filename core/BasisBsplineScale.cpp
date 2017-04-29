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
#include <iomanip>
#include <cmath>
using namespace std;
using namespace kernel;


BasisBsplineScale::
BasisBsplineScale(vector<Basis*>& bases, int parentIndex, char dimension, bool transient, int order)
    : BasisBspline(bases, static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo().dimensions, transient, parentIndex), dimension_(dimension)
{
    if (getDebugLevel() % 10 >= 2)
    {
        cout << getTimeStamp() << "   " << getIndex() << " BasisBsplineScale";
        if (isTransient()) cout << " (transient)";
        cout << endl;
    }

    const GridInfo parentGridInfo = static_cast<BasisBspline*>(bases[parentIndex])->getGridInfo();
    gridInfo() = parentGridInfo;
    gridInfo().scale[dimension_] = parentGridInfo.scale[dimension_] - 1;
    gridInfo().offset[dimension_] = parentGridInfo.offset[dimension_] / 2;
    gridInfo().extent[dimension_] = (parentGridInfo.offset[dimension_] + parentGridInfo.extent[dimension_]) / 2 + 1 - gridInfo().offset[dimension_];
    
    if (getDebugLevel() % 10 >= 2)
    {
        cout << getTimeStamp() << "     parent=" << getParentIndex() << endl;
        cout << getTimeStamp() << "     dimension=" << (int) dimension_ << endl;
        cout << getTimeStamp() << "     " << gridInfo() << endl;
    }
    
    ii stride = 1;
    for (ii j = 0; j < dimension_; j++) stride *= gridInfo().extent[j];

    // create our kernel
    ii nh = order + 2;
    vector<fp> hs(nh);
    double sum = 0.0;
    for (ii i = 0; i < nh; i++)
    {
        hs[i] = (fp) (1.0 / pow(2.0, (double)order) * Bspline::factorial(order + 1) / (double)(Bspline::factorial(i)*Bspline::factorial(order + 1 - i)));
        sum += hs[i];
    }
    for (ii i = 0; i < nh; i++)
        hs[i] /= (fp) sum;

    // create A as a temporary COO matrix
    ii m = parentGridInfo.extent[dimension_];
    ii n = gridInfo().extent[dimension_];
    vector<ii> rowind(nh * n);
    vector<ii> colind(nh * n);
    vector<fp> acoo(nh * n);

    ii nnz = 0;
    ii offset = order + ((parentGridInfo.offset[dimension_] + 1) % 2);
    for (ii j = 0; j < n; j++)
    {
        for (ii i = 0; i < nh; i++)
        {
            rowind[nnz] = 2 * j + i - offset;
            if (rowind[nnz] < 0 || rowind[nnz] >= m) continue;
            acoo[nnz] = hs[i];
            colind[nnz] = j;

            nnz++;
        }
    }

    // create A
    aT_.copy(n, m, acoo.size(), colind.data(), rowind.data(), acoo.data());

    if (dimension == 0)
        a_.copy(aT_, true);
}


BasisBsplineScale::~BasisBsplineScale()
{
}


void
BasisBsplineScale::
synthesize(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "     " << getIndex() << " BasisBsplineScale::synthesise" << endl;

    if (!f.size())
        f.resize(1);

    // zero basis functions that are no longer needed
    MatrixSparse t;
    ii rowsPruned = t.copyPruneRows(aT_, x[0], dimension_ > 0, 0.75);
    if (rowsPruned > 0)
    {
        aT_.swap(t);

        if (dimension_ == 0)
            a_.copy(aT_, true);

        if (getDebugLevel() % 10 >= 3)
            cout << getTimeStamp() << "      " << getIndex() << " pruned " << rowsPruned << " basis functions" << endl;
    }

    // synthesise
    if (dimension_ == 0)
        f[0].matmul(false, x[0], aT_, accumulate);
    else
        f[0].matmul(true, aT_, x[0], accumulate);

    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << f[0] << endl;
}


void BasisBsplineScale::analyze(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA)
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "     " << getIndex() << " BasisBsplineScale::analyse" << endl;

    if (!xE.size())
        xE.resize(1);
    
    if (sqrA)
    {
        if (dimension_ == 0)
        {
            MatrixSparse t;
            t.copy(a_);
            t.sqr();
            
            xE[0].matmul(false, fE[0], t, false);
        }
        else
        {
            MatrixSparse t;
            t.copy(aT_);
            t.sqr();
            
            xE[0].matmul(false, t, fE[0], false);
        }
    }
    else
    {
        if (dimension_ == 0)
            xE[0].matmul(false, fE[0], a_, false);
        else
            xE[0].matmul(false, aT_, fE[0], false);
    }
    
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << xE[0] << endl;
}



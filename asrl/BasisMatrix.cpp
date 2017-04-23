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


#include "BasisMatrix.hpp"
using namespace std;
using namespace kernel;


BasisMatrix::BasisMatrix(std::vector<Basis*>& bases, std::vector<MatrixSparse>& aT, std::vector<MatrixSparse>* gT, bool transient) : Basis(bases, transient), aTs_(aT), gT_(gT)
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "   BasisMatrix";
        if (isTransient()) cout << " (transient)";
        cout << " ..." << endl;
    }

    aTs_ = aT;
    aTnnzRows_.resize(aTs_.size());
    for (ii i = 0; i < ii(aTs_.size()); i++)
        aTnnzRows_[i] = aTs_[i].m();

    as_.resize(aTs_.size());
    for (ii i = 0; i < ii(as_.size()); i++)
        as_[i].copy(aTs_[i], true);

    if (gT_)
    {
        g_ = new vector<MatrixSparse>(gT_->size());

        for (ii i = 0; i < ii(g_->size()); i++)
            (*g_)[i].copy((*gT_)[i], true);
    }
}


BasisMatrix::~BasisMatrix()
{
    if (gT_)
         delete g_;
}


void BasisMatrix::synthesise(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "      BasisMatrix::synthesise" << endl;

    if (!f.size())
        f.resize(x.size());

    for (ii k = 0; k < ii(as_.size()); k++)
    {
        /*//cout << aTnnzRows_[k] << ",";
        // zero basis functions that are no longer needed
        ii aTnnzRows = aTs_[k].pruneRows(aTs_[k], aTnnzRows_[k], x[0], false, 0.75);
        //cout << aTnnzRows << endl;
        if (aTnnzRows < aTnnzRows_[k])
        {
            as_[k].copy(aTs_[k], true);

            if (getDebugLevel() % 10 >= 3)
                cout << getTimeStamp() << "      " << getIndex() << " zeroed " << aTnnzRows_[k] - aTnnzRows << " basis functions" << endl;

            aTnnzRows_[k] = aTnnzRows;
        }*/

        // synthesise
        f[k].matmul(false, x[k], aTs_[k], accumulate);
    }

    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << f[0] << endl;
}


void BasisMatrix::analyse(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA) const
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "      BasisMatrix::analyse" << endl;

    if (!xE.size())
        xE.resize(fE.size());

    for (ii i = 0; i < ii(as_.size()); i++)
    {
        if (sqrA)
        {
            MatrixSparse t;
            t.copy(as_[i]);
            t.sqr();

            xE[i].matmul(false, fE[i], t, false);
        }
        else
            xE[i].matmul(false, fE[i], as_[i], false);
    }



    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << xE[0] << endl;
}


const std::vector<MatrixSparse>* BasisMatrix::getGroups(bool transpose) const
{
    if (transpose)
        return gT_;
    else
        return g_;
}

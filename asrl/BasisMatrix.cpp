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


BasisMatrix::BasisMatrix(std::vector<Basis*>& bases, std::vector<MatrixSparse>& a, std::vector<MatrixSparse>* g, bool transient) : Basis(bases, transient), a_(a), g_(g)
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "   BasisMatrix";
        if (isTransient()) cout << " (transient)";
        cout << " ..." << endl;
    }

    a_ = a;
    aT_.resize(a_.size());
    aTnnzRows_.resize(a_.size());
    for (ii i = 0; i < ii(a_.size()); i++)
    {
        aT_[i].copy(a_[i], true);
        aTnnzRows_[i] = aT_[i].m();
    }

    if (g_)
    {
        gT_ = new vector<MatrixSparse>(g_->size());
        for (ii i = 0; i < ii(g_->size()); i++)
            (*gT_)[i].copy((*g_)[i], true);
    }
}


BasisMatrix::~BasisMatrix()
{
    if (g_)
    {
         delete gT_;
    }
}


void BasisMatrix::synthesise(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate) const
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "      BasisMatrix::synthesise" << endl;

    if (!f.size())
        f.resize(x.size());

    for (ii i = 0; i < ii(a_.size()); i++)
        f[i].matmul(false, x[i], aT_[i], accumulate);

    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << f[0] << endl;
}


void BasisMatrix::analyse(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA) const
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "      BasisMatrix::analyse" << endl;

    if (!xE.size())
        xE.resize(fE.size());

    for (ii i = 0; i < ii(a_.size()); i++)
    {
        if (sqrA)
        {
            MatrixSparse t;
            t.copy(a_[i]);
            t.sqr();

            xE[i].matmul(false, fE[i], t, false);
        }
        else
            xE[i].matmul(false, fE[i], a_[i], false);
    }



    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << xE[0] << endl;
}


void BasisMatrix::deleteBasisFunctions(const vector<MatrixSparse>& x, fp threshold)
{
    for (ii i = 0; i < ii(a_.size()); i++)
    {
        ii aTnnzRows = aT_[i].pruneRows(aT_[i], aTnnzRows_[i], x[0], false, threshold);

        if (aTnnzRows < aTnnzRows_[i])
        {
            a_[i].copy(aT_[i], true);

            if (getDebugLevel() % 10 >= 3)
                cout << getTimeStamp() << "      BasisMatrix::deleteBasisFunctions " << aTnnzRows_[i] - aTnnzRows << endl;

            aTnnzRows_[i] = aTnnzRows;
        }
    }
}


const std::vector<MatrixSparse>* BasisMatrix::getGroups(bool transpose) const
{
    if (transpose)
        return gT_;
    else
        return g_;
}

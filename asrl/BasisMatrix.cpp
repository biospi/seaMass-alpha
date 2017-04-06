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
#include <iostream>
using namespace std;


BasisMatrix::BasisMatrix(std::vector<Basis*>& bases, ii m, ii n, std::vector<fp>& aV, std::vector<ii>& aI, std::vector<ii>& aJ, Transient transient) : Basis(bases, transient)
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "   BasisMatrix";
        if (getTransient() == Basis::Transient::YES) cout << " (transient)";
        cout << " ..." << endl;
    }

    aT_.init(n, m, (ii)aV.size(), aV.data(), aJ.data(), aI.data());
    aTnnzRows_ = m;
    a_.copy(aT_, true);
}


BasisMatrix::~BasisMatrix()
{
}


void BasisMatrix::synthesis(vector<MatrixSparse>& f, const vector<MatrixSparse>& x, bool accumulate) const
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "      BasisMatrix::synthesis" << endl;

    if (!f.size()) f.resize(1);

    f[0].matmul(false, x[0], aT_, accumulate);

    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << f[0] << endl;
}


void BasisMatrix::analysis(vector<MatrixSparse>& xE, const vector<MatrixSparse>& fE, bool sqrA) const
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "      BasisMatrix::analysis" << endl;

    if (!xE.size())
        xE.resize(1);

    if (sqrA)
    {
        MatrixSparse t;
        t.copy(a_);
        t.sqr();

        xE[0].matmul(false, fE[0], t, false);
    }
    else
        xE[0].matmul(false, fE[0], a_, false);

    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << xE[0] << endl;
}


void BasisMatrix::deleteBasisFunctions(const vector<MatrixSparse>& x, fp threshold)
{
    ii aTnnzRows = aT_.pruneRows(aT_, aTnnzRows_, x[0], false, threshold);

    if (aTnnzRows < aTnnzRows_)
    {
        a_.copy(aT_, true);

        if (getDebugLevel() % 10 >= 3)
            cout << getTimeStamp() << "      BasisMatrix::deleteBasisFunctions " << aTnnzRows_ - aTnnzRows << endl;

        aTnnzRows_ = aTnnzRows;
    }
}


ii BasisMatrix::getM() const
{
    return aT_.n();
}


ii BasisMatrix::getN() const
{
    return aT_.m();
}
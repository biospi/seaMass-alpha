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


#include "BasisMatrixGroup.hpp"
#include <iostream>
using namespace std;


BasisMatrixGroup::BasisMatrixGroup(std::vector<Basis*>& bases, ii aM, ii aN, std::vector<fp>& aV, std::vector<ii>& aI, std::vector<ii>& aJ,
                                   ii gM, std::vector<fp>& gV, std::vector<ii>& gI, std::vector<ii>& gJ, Transient transient) : BasisMatrix(bases, aM, aN, aV, aI, aJ, transient)
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "   BasisMatrixGroup";
        if (getTransient() == Basis::Transient::YES) cout << " (transient)";
        cout << " ..." << endl;
    }

    g_.init(gM, aN, (ii)aV.size(), gV.data(), gI.data(), gJ.data());
    gT_.copy(g_, true);
}


BasisMatrixGroup::~BasisMatrixGroup()
{
}


void BasisMatrixGroup::groupSynthesis(std::vector<MatrixSparse>& f, const std::vector<MatrixSparse>& x, bool accumulate) const
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "      BasisMatrixGroup::groupSynthesis" << endl;

    if (!f.size()) f.resize(1);

    f[0].matmul(false, x[0], gT_, accumulate);

    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "       " << f[0] << endl;
}


void BasisMatrixGroup::shrinkage(std::vector<MatrixSparse>& y, const std::vector<MatrixSparse>& x, const std::vector<MatrixSparse>& xE, const std::vector<MatrixSparse>& l1l2, fp lambda) const
{
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "     " << getIndex() << " BasisMatrixGroup::shrinkage" << endl;
    }

    if (!y.size())
        y.resize(x.size());

    for (size_t k = 0; k < y.size(); k++)
    {
        // y = groupNorm(x)
        MatrixSparse t;
        t.copy(x[k]);
        t.sqr();
        y[k].matmul(false, t, gT_, false);
        t.matmul(false, y[k], g_, false);
        t.sort();
        y[k].copy(x[k]);
        y[k].subsetCopy(t);
        y[k].sqrt();

        // y = x * groupNorm(x)^-1)
        y[k].div2Nonzeros(x[k].vs());

        // y = lambda * x * groupNorm(x)^-1
        y[k].mul(lambda);

        // y = l1l2 + lambda * x * groupNorm(x)^-1
        y[k].addNonzeros(l1l2[k].vs());

        // y = x / (l1l2 + lambda * x * groupNorm(x)^-1)
        y[k].div2Nonzeros(x[k].vs());

        // y = xE * x / (l1l2 + lambda * x * groupNorm(x)^-1)
        y[k].mul(xE[k].vs());
    }
}


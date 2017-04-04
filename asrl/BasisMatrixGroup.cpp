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

    /*vector<fp> t = gV;

    for (ii i = 0; i < gM; i++)
    {
        ii count = 0;
        for (ii nz = 0; nz < gV.size(); nz++)
        {
            if (gI[nz] == i) count++;
        }
        for (ii nz = 0; nz < gV.size(); nz++)
        {
            if (gI[nz] == i) gV[nz] /= count;
        }
    }

    for (ii i = 0; i < gV.size(); i++)
    {
        cout << t[i] << "," << gV[i] << "," << gI[i] << "," << gJ[i] << endl;
    }*/

    //exit(0);

    /*vector<fp> tV(aN);
    vector<ii> tI(aN);
    vector<ii> tJ(aN);
    for (ii i = 0; i < aN; i++)
    {
        tV[i] = 1.0;
        tI[i] = i;
        tJ[i] = i;
    }
    gT_.init(aN, aN, (ii)tV.size(), tV.data(), tJ.data(), tI.data());
    g_.copy(gT_, true);*/

    gT_.init(aN, gM, (ii)aV.size(), gV.data(), gJ.data(), gI.data());
    g_.copy(gT_, true);
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


void BasisMatrixGroup::shrinkage(std::vector<MatrixSparse>& y, std::vector<MatrixSparse>& x, const std::vector<MatrixSparse>& xE, const std::vector<MatrixSparse>& l1l2, fp lambda) const
{
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "     " << getIndex() << " BasisMatrixGroup::shrinkage" << endl;
    }

    if (!y.size())
        y.resize(x.size());

    for (size_t k = 0; k < y.size(); k++)
    {
        // Gx
        MatrixSparse gX;
        gX.matmul(false, x[k], gT_, false);

        //for (ii nz = 0; nz < gT_.nnz(); nz++)
        //{
        //    cout << gT_.vs()[nz] << endl;
        //}
        //for (ii nz = 0; nz < gXgT.nnz(); nz++)
        //{
        //    cout << gXgT.vs()[nz] << endl;
        //}
        //exit(0);

        // GxGt
        MatrixSparse t;
        t.matmul(false, gX, g_, false);
        t.sort();
        MatrixSparse& gXgT = gX;
        gXgT.copy(x[k]);
        gXgT.subsetCopy(t);

        //for (ii nz = 0; nz < x[k].nnz(); nz++)
        //{
        //    cout << x[k].vs()[nz] << "," << gXgT.vs()[nz] << endl;
        //}
        //exit(0);

        // x./GxGt
        MatrixSparse& l1l2PlusLambdas = t;
        l1l2PlusLambdas.copy(x[k]);
        l1l2PlusLambdas.divNonzeros(gXgT.vs());

        // lambda * x./GxGt + l1l2
        l1l2PlusLambdas.mul(lambda);
        l1l2PlusLambdas.addNonzeros(l1l2[k].vs());

        // y = x * xE / (lambda * x./GxGt + l1l2)
        y[k].copy(x[k]);
        y[k].divNonzeros(l1l2PlusLambdas.vs());
        y[k].mul(xE[k].vs());
    }
}


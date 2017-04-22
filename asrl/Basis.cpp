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


#include "Basis.hpp"
using namespace std;
using namespace kernel;


Basis::Basis(vector<Basis*>& bases, bool transient, int parentIndex) : parentIndex_(parentIndex), transient_(transient)
{
    index_ = (ii) bases.size();
    bases.push_back(this);
}


Basis::~Basis()
{
}


const std::vector<MatrixSparse>* Basis::getGroups(bool transpose) const
{
    return 0;
}


void Basis::synthesiseGroups(std::vector<MatrixSparse> &f, const std::vector<MatrixSparse> &x, bool accumulate)
{
    const std::vector<MatrixSparse>* gT = getGroups(true);
    if (gT)
    {
        if (!f.size())
            f.resize(x.size());

        for (ii k = 0; k < ii(f.size()); k++)
            f[k].matmul(false, x[k], (*gT)[k], accumulate);
    }
    else
    {
        synthesise(f, x, accumulate);
    }
}


int Basis::getIndex() const
{
    return index_;
}


int Basis::getParentIndex() const
{
    return parentIndex_;
}


bool Basis::isTransient() const
{
    return transient_;
}

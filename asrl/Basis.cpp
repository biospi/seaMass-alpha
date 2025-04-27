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
#include <kernel.hpp>
using namespace std;
using namespace kernel;


Basis::Basis(vector<Basis*>& bases, bool transient, int parentIndex, std::vector<MatrixSparse>* gT) : parentIndex_(parentIndex), transient_(transient), gT_(gT), g_(0)
{
    index_ = (ii) bases.size();
    bases.push_back(this);

    if (gT_)
    {
        g_ = new vector<MatrixSparse>(gT_->size());

        for (ii i = 0; i < ii(g_->size()); i++)
            (*g_)[i].transpose((*gT_)[i]);
    }
}


Basis::~Basis()
{
    if (gT_)
        delete g_;
}


void Basis::synthesizeGroups(std::vector<MatrixSparse> &g, const vector<MatrixSparse> &x, bool accumulate, bool prune)
{
    const std::vector<MatrixSparse>* gT = getColGroups(true);
    if (gT)
    {
        if (!g.size())
            g.resize(x.size());

        for (ii k = 0; k < ii(g.size()); k++)
            g[k].matmul(false, x[k], (*gT)[k], accumulate);
    }
    else
    {
        synthesize(g, x, accumulate, prune);
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


const vector<MatrixSparse>* Basis::getColGroups(bool transpose) const
{
    if (transpose)
        return gT_;
    else
        return g_;
}
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


#include <iostream>


using namespace std;


Basis::Basis(vector<Basis*>& bases, Transient transient, int parentIndex)
	: parentIndex_(parentIndex), transient_(transient)
{
	index_ = (ii) bases.size();
	bases.push_back(this);
}


Basis::~Basis()
{
}


void Basis::shrinkage(vector<MatrixSparse>& y, vector<MatrixSparse>& x, const vector<MatrixSparse>& xE, const vector<MatrixSparse>& l1l2PlusLambda) const
{
    if (getDebugLevel() % 10 >= 3)
    {
        cout << getTimeStamp() << "   " << getIndex() << " Basis::shrinkage" << endl;
    }

    if (!y.size()) y.resize(x.size());
    for (size_t k = 0; k < y.size(); k++)
    {
        y[k].copy(x[k]);
        y[k].divNonzeros(l1l2PlusLambda[k].vs());
        y[k].mul(xE[k]);
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


Basis::Transient Basis::getTransient() const
{
	return transient_;
}

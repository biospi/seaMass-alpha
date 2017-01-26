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


Basis::Basis(vector<Basis*>& bases, Transient transient, int parentIndex)
	: parentIndex_(parentIndex), transient_(transient)
{
	index_ = (ii) bases.size();
	bases.push_back(this);
}


Basis::~Basis()
{
}


void Basis::shrinkage(MatrixSparse& y, MatrixSparse& x, const MatrixSparse& xE, const MatrixSparse& l1l2PlusLambda) const
{
#ifndef NDEBUG
	cout << getTimeStamp() << " " << getIndex() << " Basis::shrinkage" << endl;
#endif
    
	y.copy(x);
	y.elementwiseDiv(l1l2PlusLambda);
	y.elementwiseMul(xE);
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

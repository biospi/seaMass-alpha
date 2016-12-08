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


#ifndef _SEAMASS_CORE_BASIS_HPP_
#define _SEAMASS_CORE_BASIS_HPP_


#include "../math/Matrix.hpp"


class Basis
{
public:
	Basis(std::vector<Basis*>& bases, bool isTransient = false, ii parentIndex = -1);
	virtual ~Basis();

	virtual void synthesis(Matrix& f, const Matrix& x, bool accumulate) const = 0;
	virtual void analysis(Matrix& xE, const Matrix& fE, bool sqrA) const = 0;
	virtual void shrinkage(Matrix& x, const Matrix& xE, const Matrix& x0, const Matrix& l1, fp lambda) const;

	virtual ii getM() const = 0;
	virtual ii getN() const = 0;

	ii getIndex() const;
	ii getParentIndex() const;
	bool isTransient() const;

private:
	ii index_;       // index of this basis in the serialised tree
	ii parentIndex_; // parent node
	bool isTransient_; // if transient, coefficients not part of fitting
};


#endif


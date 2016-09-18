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


#ifndef _SEAMASS_CORE_BASISBSPLINESCALE_HPP_
#define _SEAMASS_CORE_BASISBSPLINESCALE_HPP_


#include "BasisBspline.hpp"


class BasisBsplineScale : public BasisBspline
{ 
public:
	BasisBsplineScale(std::vector<Basis*>& bases, ii parentIndex, ii dim, ii order = 3, bool isTransient = false);
	virtual ~BasisBsplineScale();

	void synthesis(Matrix& f, const Matrix& c, bool accumulate) const;
	void analysis(Matrix& cE, const Matrix& fE, bool sqrA = false) const;

private:
	MatrixSparse a;
	MatrixSparse aT;

	ii dimension;
};


#endif

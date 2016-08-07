//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  biospi Laboratory, EEE, University of Liverpool, UK
//
// This file is part of seaMass-TD.
//
// seaMass-TD is free software: you can redistribute it and/or modify
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
// along with seaMass-TD.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef _SEAMASS_CORE_BSPLINE_HPP_
#define _SEAMASS_CORE_BSPLINE_HPP_


#include "seaMass.hpp"


class BSpline
{
protected:
	ii order, n;
	std::vector<double> lookup;

public:
	BSpline(ii order, ii n);
	double ibasis(double x);

	static double m(double x, ii k, ii i, std::vector<fp>& ks);
	static double m(double x, ii k, ii i);
	static double im(double x, ii k);
	static ii factorial(ii n);
};


#endif


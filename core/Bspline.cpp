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


#include "Bspline.hpp"

#include <iomanip>


using namespace std;


Bspline::Bspline(ii order, ii n)
	: order_(order), n_(n), lookup_(n)
{
	for (ii i = 0; i < n; i++)
	{
		lookup_[i] = im(i / (double)(n - 1) * (order + 1), order + 1);
 	}
}


double Bspline::ibasis(double x)
{
	if (x >= order_ + 1)
	{
		return 1.0;
	}
	else if (x <= 0.0)
	{
		return 0.0;
	}
	else
	{
		double f = x / (order_ + 1) * (n_ - 1);
		ii i = (ii)f;
		f = f - i;
		return (1 - f) * lookup_[i] + f * lookup_[i + 1];
	}
}


double Bspline::m(double x, ii k, ii i, vector<fp>& ks)
{
	if (k == 1)
	{
		if (ks[i] <= x && x < ks[i + 1])
		{
			return 1.0 / (ks[i + 1] - ks[i]);
		}
		else
		{
			return 0.0;
		}
	}
	else
	{
		if (ks[i + k] - ks[i] == 0.0)
		{
			return 0.0;
		}
		else
		{
			return (k*((x - ks[i]) * m(x, k - 1, i, ks) +
				(ks[i + k] - x) * m(x, k - 1, i + 1, ks))) / ((k - 1)*(ks[i + k] - ks[i]));
		}
	}
}


double Bspline::m(double x, ii k, ii i)
{
	if (k == 1)
	{
		if (i <= x && x < i + 1)
		{
			return 1.0;
		}
		else
		{
			return 0.0;
		}
	}
	else
	{
		return (k*((x - i) * m(x, k - 1, i) +
			((i + k) - x) * m(x, k - 1, i + 1))) / ((k - 1)*k);
	}
}


double Bspline::im(double x, ii k)
{
	double v = 0.0;
	for (ii i = 0; i < k + 1; ++i)
	{
		v += m(x, k + 1, i);
	}
	return v;
}


ii Bspline::factorial(ii n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

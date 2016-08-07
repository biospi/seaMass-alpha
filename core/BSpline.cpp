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


#include "BSpline.hpp"
using namespace std;


BSpline::
BSpline(ii _order, ii _n) :
	order(_order),
	n(_n),
	lookup(n+1)
{
	for (ii i = 0; i < n; i++)
	{
		lookup[i] = im(i / (double)(n - 1) * (order + 1), order + 1);
 	}
	lookup[n] = 1.0;
}


double
BSpline::
ibasis(double x)
{
	double f = x / (order + 1) * n;
	ii i = (ii) f;
	f = f - i;
	return (1-f)*lookup[i] + f*lookup[i+1];
}


double
BSpline::
m(double x, ii k, ii i, vector<fp>& ks)
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


double
BSpline::
m(double x, ii k, ii i)
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


double
BSpline::
im(double x, ii k)
{
	double v = 0.0;
	for (ii i = 0; i < k + 1; ++i)
	{
		v += m(x, k + 1, i);
	}
	return v;
}


ii
BSpline::
factorial(ii n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

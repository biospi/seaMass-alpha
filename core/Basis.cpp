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


Basis::Basis(vector<Basis*>& bases, bool isTransient, ii parentIndex)
	: parentIndex_(parentIndex), isTransient_(isTransient)
{
	index_ = (ii) bases.size();
	bases.push_back(this);
	if (parentIndex_ >= 0) bases[parentIndex_]->childCount_++;
}


Basis::~Basis()
{
}


Basis::ErrorInfo Basis::error(Matrix& fE, const Matrix& f, const Matrix& g) const
{
	ErrorInfo info;

	ii size_d = 0;
	double discrepancy = 0.0;
	double error = 0.0;
	double sum_g = 0.0;
	double sum_f = 0.0;
	info.maxError = 0;
	for (li i = 0; i < g.size(); i++)
	{
		double v = fabs(g.vs_[i] - f.vs_[i]);
		info.maxError = info.maxError > v ? info.maxError : v;
	}

	// bug in this openmp section at present for err and discrep
	//#pragma omp parallel for simd reduction(+:dis,err,size_d,sum_g,sum_f)
	for (li i = 0; i < g.size(); i++)
	{
		sum_g += g.vs_[i];
		sum_f += f.vs_[i];

		if (f.vs_[i] > 0.0 && g.vs_[i] >= 0.0)
		{
			if (g.vs_[i] > 0.0)
			{
				discrepancy += ((g.vs_[i] - ceil(f.vs_[i]))*(g.vs_[i] - ceil(f.vs_[i]))) / ceil(f.vs_[i]);
				error += fabs(g.vs_[i] - f.vs_[i]);
				size_d++;
			}

			fE.vs_[i] = g.vs_[i] / f.vs_[i];
		}
		else
		{
			fE.vs_[i] = 1.0;
		}
	}
	info.volume = sum_f / sum_g;
	info.discrepancy = discrepancy / size_d;
	info.error = error / size_d;

	return info;
}


void Basis::shrink(Matrix& c, const Matrix& cE, const Matrix& c0, const Matrix& l1, const Matrix& l2, fp shrinkage) const
{
	// Standard shrinkage
	#pragma omp parallel for
	for (li i = 0; i < cE.size(); i++)
	{
		if (cE.vs_[i] >= FLT_MIN && c.vs_[i] >= FLT_MIN)
		{
			c.vs_[i] = cE.vs_[i] * c0.vs_[i] / (shrinkage * l2.vs_[i] + l1.vs_[i]);
		}
		else
		{
			c.vs_[i] = 0.0;
		}
	}
}


ii Basis::getIndex() const
{
	return index_;
}


ii Basis::getParentIndex() const
{
	return parentIndex_;
}


bool Basis::isTransient() const
{
	return isTransient_;
}

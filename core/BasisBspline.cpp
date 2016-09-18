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


#include "BasisBspline.hpp"
using namespace std;


BasisBspline::BasisBspline(std::vector<Basis*>& bases, ii dims, bool isTransient, ii parentIndex)
	: Basis(bases, isTransient, parentIndex), meshInfo_(dims)
{
}


BasisBspline::~BasisBspline()
{
}


BasisBspline::MeshInfo::MeshInfo(ii dimensions_)
	: dimensions(dimensions_), scale(dimensions_), offset(dimensions_), extent(dimensions_), n(0)
{
}


void BasisBspline::MeshInfo::operator=(const BasisBspline::MeshInfo& mi)
{
	dimensions = mi.dimensions;
	scale = mi.scale;
	offset = mi.offset;
	extent = mi.extent;
	n = mi.n;
}


BasisBspline::MeshInfo::~MeshInfo()
{
}


ii BasisBspline::MeshInfo::m() const
{
	ii m = 1;
	for (ii i = 0; i < dimensions; i++)
	{
		m *= extent[i];
	}
	return m;
}


li BasisBspline::MeshInfo::size() const
{
	li size = n;
	for (ii i = 0; i < dimensions; i++)
	{
		size *= extent[i];
	}
	return size;
}


ii BasisBspline::getM() const
{
	return meshInfo_.m();
}


ii BasisBspline::getN() const
{
	return meshInfo_.n;
}


const BasisBspline::MeshInfo& BasisBspline::getMeshInfo() const
{
	return meshInfo_;
}


BasisBspline::MeshInfo& BasisBspline::meshInfo()
{
	return meshInfo_;
}


ostream&
operator<<(ostream& os, const BasisBspline::MeshInfo& meshInfo)
{
	os << "meshInfo=(n=" << meshInfo.n;

	os << ",scale=[";
	for (ii i = 0; i < meshInfo.dimensions; i++)
	{
		os << meshInfo.scale[i];
		if (i < meshInfo.dimensions - 1) os << ",";
	}
	os << "]";

	os << ",offset=[";
	for (ii i = 0; i < meshInfo.dimensions; i++)
	{
		os << meshInfo.offset[i];
		if (i < meshInfo.dimensions - 1) os << ",";
	}
	os << "]";

	os << ",extent=[";
	for (ii i = 0; i < meshInfo.dimensions; i++)
	{
		os << meshInfo.extent[i];
		if (i < meshInfo.dimensions - 1) os << ",";
	}
	os << "])";

	return os;
}





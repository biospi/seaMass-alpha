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

#include <iostream>


using namespace std;


BasisBspline::BasisBspline(std::vector<Basis*>& bases, char dimensions, bool transient, int parentIndex)
    : Basis(bases, transient, parentIndex), gridInfo_(dimensions)
{
}


BasisBspline::~BasisBspline()
{
}


BasisBspline::GridInfo::GridInfo(char dimensions_)
    : dimensions(dimensions_), scale(dimensions_), offset(dimensions_), extent(dimensions_), count(0)
{
}


void BasisBspline::GridInfo::operator=(const BasisBspline::GridInfo& mi)
{
    dimensions = mi.dimensions;
    scale = mi.scale;
    offset = mi.offset;
    extent = mi.extent;
    count = mi.count;
}


BasisBspline::GridInfo::~GridInfo()
{
}


ii BasisBspline::GridInfo::m() const
{
    ii m = count;
    for (ii i = 1; i < dimensions; i++)
    {
        m *= extent[i];
    }
    return m;
}


ii BasisBspline::GridInfo::n() const
{
    return extent[0];
}


li BasisBspline::GridInfo::size() const
{
    li size = count;
    for (ii i = 0; i < dimensions; i++)
    {
        size *= extent[i];
    }
    return size;
}


const BasisBspline::GridInfo& BasisBspline::getGridInfo() const
{
    return gridInfo_;
}


BasisBspline::GridInfo& BasisBspline::gridInfo()
{
    return gridInfo_;
}


ostream&
operator<<(ostream& os, const BasisBspline::GridInfo& gridInfo)
{
    os << "gridInfo=[" << gridInfo.m() << "," << gridInfo.n() << "]:(count=" << gridInfo.count;

    os << ",scale=[";
    for (ii i = 0; i < gridInfo.dimensions; i++)
    {
        os << (int) gridInfo.scale[i];
        if (i < gridInfo.dimensions - 1) os << ",";
    }
    os << "]";

    os << ",offset=[";
    for (ii i = 0; i < gridInfo.dimensions; i++)
    {
        os << gridInfo.offset[i];
        if (i < gridInfo.dimensions - 1) os << ",";
    }
    os << "]";

    os << ",extent=[";
    for (ii i = 0; i < gridInfo.dimensions; i++)
    {
        os << gridInfo.extent[i];
        if (i < gridInfo.dimensions - 1) os << ",";
    }
    os << "])";

    return os;
}





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

#include <limits>


using namespace std;


BasisBspline::BasisBspline(std::vector<Basis*>& bases, short rowDimensions, short colDimensions,
                           bool transient, int parentIndex)
    : Basis(bases, transient, parentIndex), gridInfo_(rowDimensions, colDimensions)
{
}


BasisBspline::~BasisBspline()
{
}


BasisBspline::GridInfo::GridInfo(short rowDimensions, short colDimensions)
    : rowScale(rowDimensions), rowOffset(rowDimensions), rowExtent(rowDimensions),
      colScale(colDimensions), colOffset(colDimensions), colExtent(colDimensions)
{
}

short BasisBspline::GridInfo::rowDimensions() const
{
    return rowScale.size();
}


short BasisBspline::GridInfo::colDimensions() const
{
    return colScale.size();
}


/*void BasisBspline::GridInfo::operator=(const BasisBspline::GridInfo& mi)
{
    rowScale = mi.colScale;
    colOffset = mi.colOffset;
    colExtent = mi.colExtent;

    colScale = mi.colScale;
    colOffset = mi.colOffset;
    colExtent = mi.colExtent;
}*/


BasisBspline::GridInfo::~GridInfo()
{
}


ii BasisBspline::GridInfo::m() const
{
    ii m = 1;
    for (short i = 0; i < rowDimensions(); i++)
    {
        m *= rowExtent[i];
    }
    return m;
}


ii BasisBspline::GridInfo::n() const
{
    ii n = 1;
    for (short i = 0; i < colDimensions(); i++)
    {
        n *= colExtent[i];
    }
    return n;
}


li BasisBspline::GridInfo::size() const
{
    return li(m()) * li(n());
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
    os << "gridInfo=[" << gridInfo.m() << "," << gridInfo.n() << "]";

    os << ",extent=[[";
    for (short i = 0; i < gridInfo.rowDimensions(); i++)
    {
        os << gridInfo.rowExtent[i];
        if (i < gridInfo.rowDimensions() - 1)
            os << ",";
    }
    os << "],[";
    for (short i = 0; i < gridInfo.colDimensions(); i++)
    {
        os << gridInfo.colExtent[i];
        if (i < gridInfo.colDimensions() - 1)
            os << ",";
    }
    os << "]]";

    os << ",offset=[[";
    for (short i = 0; i < gridInfo.rowDimensions(); i++)
    {
        os << gridInfo.rowOffset[i];
        if (i < gridInfo.rowDimensions() - 1)
            os << ",";
    }
    os << "],[";
    for (short i = 0; i < gridInfo.colDimensions(); i++)
    {
        os << gridInfo.colOffset[i];
        if (i < gridInfo.colDimensions() - 1)
            os << ",";
    }
    os << "]]";

    os << ",scale=[[";
    for (short i = 0; i < gridInfo.rowDimensions(); i++)
    {
        if (gridInfo.rowScale[i] == numeric_limits<short>::min())
            os << "NA";
        else
            os << gridInfo.rowScale[i];

        if (i < gridInfo.rowDimensions() - 1)
            os << ",";
    }
    os << "],[";
    for (short i = 0; i < gridInfo.colDimensions(); i++)
    {
        if (gridInfo.colScale[i] == numeric_limits<short>::min())
            os << "NA";
        else
            os << gridInfo.colScale[i];

        if (i < gridInfo.colDimensions() - 1)
            os << ",";
    }
    os << "]]";



    return os;
}





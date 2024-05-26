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


#include "BasisGrid.hpp"

#include <limits>


using namespace std;


BasisGrid::BasisGrid(std::vector<Basis*>& bases, const GridInfo& parentGridInfo, bool transient)
    : Basis(bases, transient, parentGridInfo.index), gridInfo_(parentGridInfo)
{
    gridInfo_.index = getIndex();
}


BasisGrid::~BasisGrid()
{
}


BasisGrid::GridInfo::GridInfo(short rowDimensions, short colDimensions)
    : index(-1),
      rowScale(rowDimensions), rowOffset(rowDimensions), rowExtent(rowDimensions),
      colScale(colDimensions), colOffset(colDimensions), colExtent(colDimensions)
{
}

short BasisGrid::GridInfo::rowDimensions() const
{
    return rowScale.size();
}


short BasisGrid::GridInfo::colDimensions() const
{
    return colScale.size();
}


BasisGrid::GridInfo::~GridInfo()
{
}


ii BasisGrid::GridInfo::m() const
{
    ii m = 1;
    for (short i = 0; i < rowDimensions(); i++)
    {
        m *= rowExtent[i];
    }
    return m;
}


ii BasisGrid::GridInfo::n() const
{
    ii n = 1;
    for (short i = 0; i < colDimensions(); i++)
    {
        n *= colExtent[i];
    }
    return n;
}


li BasisGrid::GridInfo::size() const
{
    return li(m()) * li(n());
}


const BasisGrid::GridInfo& BasisGrid::getGridInfo() const
{
    return gridInfo_;
}

const MatrixSparse& BasisGrid::getAt() const
{
    return aT_;
}

const std::vector<ii>& BasisGrid::getIDs() const
{
    return ids_;
}

BasisGrid::GridInfo& BasisGrid::gridInfo()
{
    return gridInfo_;
}


const std::string& BasisGrid::getType() const
{
    return type_;
}


std::string& BasisGrid::type()
{
    return type_;
}

ostream&
operator<<(ostream& os, const BasisGrid::GridInfo& gridInfo)
{
    os << "A[" << gridInfo.m() << ", " << gridInfo.n() << "] extent=[[";
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
    os << "]] offset=[[";
    for (short i = 0; i < gridInfo.rowDimensions(); i++)
    {
        if (gridInfo.rowOffset[i] == -9999)
            os << "NA";
        else
            os << gridInfo.rowOffset[i];
        if (i < gridInfo.rowDimensions() - 1)
            os << ",";
    }
    os << "],[";
    for (short i = 0; i < gridInfo.colDimensions(); i++)
    {
        if (gridInfo.colOffset[i] == -9999)
            os << "NA";
        else
            os << gridInfo.colOffset[i];
        if (i < gridInfo.colDimensions() - 1)
            os << ",";
    }
    os << "]] scale=[[";
    for (short i = 0; i < gridInfo.rowDimensions(); i++)
    {
        if (gridInfo.rowScale[i] == -9999)
            os << "NA";
        else
            os << gridInfo.rowScale[i];

        if (i < gridInfo.rowDimensions() - 1)
            os << ",";
    }
    os << "],[";
    for (short i = 0; i < gridInfo.colDimensions(); i++)
    {
        if (gridInfo.colScale[i] == -9999)
            os << "NA";
        else
            os << gridInfo.colScale[i];

        if (i < gridInfo.colDimensions() - 1)
            os << ",";
    }
    os << "]]";



    return os;
}





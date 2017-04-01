//
// $Id$
//
//
// Original author: Ranjeet Bhamber <ranjeet.bhamber <a.t> bristol.ac.uk>
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

#include "NetcdfWriter.hpp"


using namespace std;


NetcdfWriter::NetcdfWriter(const string& fileName) :
        fileName_(fileName)
{
    file_.open(fileName_,NC_NETCDF4);
}


NetcdfWriter::~NetcdfWriter()
{
    file_.close();
}


void NetcdfWriter::writeSmi(SeamassCore::Input& input)
{
    if (getDebugLevel() % 10 >= 2)
        cout << getTimeStamp() << "  Writing " << fileName_ << " ..." << endl;

    if (input.startTimes.size() > 0) {
        if (getDebugLevel() % 10 >= 3) cout << getTimeStamp() << "    startTimes ..." << endl;
        write("startTimes", input.startTimes);
     }

    if (input.finishTimes.size() > 0) {
        if (getDebugLevel() % 10 >= 3) cout << getTimeStamp() << "    finishTimes ..." << endl;
        write("finishTimes", input.finishTimes);
    }

    if (input.exposures.size() > 0) {
        if (getDebugLevel() % 10 >= 3) cout << getTimeStamp() << "    exposures ..." << endl;
        write("exposures", input.exposures);
    }

    if (input.binCountsIndex.size() > 0) {
        if (getDebugLevel() % 10 >= 3) cout << getTimeStamp() << "    binCountsIndex ..." << endl;
        write("binCountsIndex", input.binCountsIndex);
    }

    if (getDebugLevel() % 10 >= 3) cout << getTimeStamp() << "    binCounts" << endl;
    write("binCounts", input.binCounts);

    if (getDebugLevel() % 10 >= 3) cout << getTimeStamp() << "    binEdges" << endl;
    write("binEdges", input.binEdges);
}

void NetcdfWriter::writeSmv(SeamassCore::Output& output, ii shrinkage, ii tolerance, ii page_size)
{
    vector<ii> shrink(1,shrinkage);
    vector<ii> tol(1,tolerance);

    file_.write_AttNC("seamassIndex","extent",output.baselineExtent,NC_INT);
    file_.write_AttNC("seamassIndex","scale",output.baselineScale,NC_SHORT);
    file_.write_AttNC("seamassIndex","offset",output.baselineOffset,NC_INT);
    file_.write_AttNC("seamassIndex","shrinkage",shrink,NC_INT);
    file_.write_AttNC("seamassIndex","tolerance",tol,NC_INT);
}


void NetcdfWriter::writeSmo(SeamassCore::ControlPoints& controlPoints)
{
    if (getDebugLevel() % 10 >= 2)
        cout << getTimeStamp() << "  Writing " << fileName_ << " ..." << endl;

    VecMat<fp> cpMat(controlPoints.extent[1], controlPoints.extent[0], controlPoints.coeffs);
    file_.write_MatNC("controlPoints", cpMat, sizeof(controlPoints.coeffs[0]) == 4 ? NC_FLOAT : NC_DOUBLE);

    file_.write_AttNC("controlPoints","scale",controlPoints.scale,NC_SHORT);
    file_.write_AttNC("controlPoints","offset",controlPoints.offset,sizeof(controlPoints.offset[0]) == 4 ? NC_INT : NC_INT64);
}


void NetcdfWriter::write(const string& objectname, vector<short>& cdata)
{
    write(objectname, cdata, NC_SHORT);
}


void NetcdfWriter::write(const string& objectname, vector<float>& cdata)
{
    write(objectname, cdata, NC_FLOAT);
}


void NetcdfWriter::write(const string& objectname, vector<double>& cdata)
{
    write(objectname, cdata, NC_DOUBLE);
}


void  NetcdfWriter::write(const string& objectname,vector<long>& cdata)
{
    write(objectname, cdata, NC_INT);
}


void NetcdfWriter::write(const string& objectname, vector<long long>& cdata)
{
    write(objectname, cdata, NC_INT64);
}


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

NetcdfWriter::NetcdfWriter(const string& _filename) :
        filename(_filename)
{
    fileout.open(filename,NC_NETCDF4);
}


NetcdfWriter::~NetcdfWriter()
{
    fileout.close();
}


void NetcdfWriter::write_input(SeamassCore::Input& input)
{
    write("binEdges", input.binEdges);
    write("binCounts", input.binCounts);
    if (input.spectrumIndex.size() > 0) write("spectrumIndex", input.spectrumIndex);
    if (input.startTimes.size() > 0) write("startTimes", input.startTimes);
    if (input.finishTimes.size() > 0) write("finishTimes", input.finishTimes);
    if (input.exposures.size() > 0) write("exposures", input.exposures);
}

void NetcdfWriter::write_output(SeamassCore::Output& output, ii shrinkage, ii tolerance, ii page_size)
{
    vector<ii> shrik(1,shrinkage);
    vector<ii> tol(1,tolerance);

    fileout.write_AttNC("seamassIndex","extent",output.baselineExtent,NC_INT);
    fileout.write_AttNC("seamassIndex","scale",output.baselineScale,NC_SHORT);
    fileout.write_AttNC("seamassIndex","offset",output.baselineOffset,NC_INT);
    fileout.write_AttNC("seamassIndex","shrinkage",shrik,NC_INT);
    fileout.write_AttNC("seamassIndex","tolerance",tol,NC_INT);
}


void NetcdfWriter::write_output_control_points(SeamassCore::ControlPoints& controlPoints)
{

    VecMat<float> cpMat(uli(controlPoints.extent[0]),uli(controlPoints.extent.size()),controlPoints.coeffs);

    fileout.write_MatNC("controlPoints",cpMat,NC_FLOAT);
    fileout.write_AttNC("controlPoints","offset",controlPoints.offset,NC_INT);
    fileout.write_AttNC("controlPoints","scale",controlPoints.scale,NC_INT);
}


void NetcdfWriter::write(const string& objectname, vector<unsigned char>& cdata)
{
    write(objectname, cdata, NC_UBYTE);
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


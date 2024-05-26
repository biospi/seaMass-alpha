//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
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


#include "DatasetSeamass.hpp"
#include <kernel.hpp>
#include <iomanip>
using namespace kernel;


DatasetSeamass::DatasetSeamass(const std::string& filePathIn, const std::string& filePathStemOut, Dataset::WriteType writeType) : fileIn_(0), fileOut_(0), finished_(false)
{
    if (!filePathIn.empty())
        fileIn_ = new FileNetcdf(filePathIn);

    if (!filePathStemOut.empty())
        fileOut_ = new FileNetcdf(filePathStemOut + (writeType == Dataset::WriteType::InputOutput ? ".smv" : ".smb"), NC_NETCDF4);
}


DatasetSeamass::~DatasetSeamass()
{
    if (fileIn_)
        delete fileIn_;

    if (fileOut_)
        delete fileOut_;
}


bool DatasetSeamass::read(Seamass::Input &input, std::string &id)
{
    input = Seamass::Input();

    if(finished_ == true)
        return false;

    vector<short> pols;
    fileIn_->readAttribute(pols, "polarity", "");
    if (pols.size() == 1)
        input.polarity = pols[0];
    else
        input.polarity = 0;

    if (fileIn_->exists("counts"))
        fileIn_->readVector(input.counts, "counts");
    else
        throw runtime_error("ERROR: 'counts' dataset not found in smb file");

    if (fileIn_->exists("binLocations"))
    {
        fileIn_->readVector(input.locations, "binLocations");
        input.type = Seamass::Input::Type::Binned;
    }
    else if (fileIn_->exists("sampleLocations"))
    {
        fileIn_->readVector(input.locations, "sampleLocations");
        input.type = Seamass::Input::Type::Sampled;
    }
    else if (fileIn_->exists("centroidLocations"))
    {
        fileIn_->readVector(input.locations, "centroidLocations");
        input.type = Seamass::Input::Type::Centroided;
    }
    else
        throw runtime_error("ERROR: one dataset called 'binLocations', 'sampleLocations' or 'centroidLocations' is needed in smb file");

    if (fileIn_->exists("countsIndex"))
        fileIn_->readVector(input.countsIndex, "countsIndex");

    if (fileIn_->exists("startTimes"))
        fileIn_->readVector(input.startTimes, "startTimes");

    if (fileIn_->exists("finishTimes"))
        fileIn_->readVector(input.finishTimes, "finishTimes");

    if (fileIn_->exists("exposures"))
        fileIn_->readVector(input.exposures, "exposures");

    id = "";

    return finished_ = true;
}


void DatasetSeamass::write(const Seamass::Input &input, const std::string &id)
{
    fileOut_->writeAttribute(input.polarity, "polarity", "");

    if (input.startTimes.size() > 0)
        fileOut_->writeVector(input.startTimes, "startTimes");

    if (input.finishTimes.size() > 0)
        fileOut_->writeVector(input.finishTimes, "finishTimes");

    if (input.exposures.size() > 0)
        fileOut_->writeVector(input.exposures, "exposures");

    if (input.countsIndex.size() > 0)
        fileOut_->writeVector(input.countsIndex, "countsIndex");

    if (input.counts.size() > 0)
        fileOut_->writeVector(input.counts, "counts");

    if (input.locations.size() > 0)
    {
        switch (input.type)
        {
            case Seamass::Input::Type::Binned:
                fileOut_->writeVector(input.locations, "binLocations");
                break;
            case Seamass::Input::Type::Sampled:
                fileOut_->writeVector(input.locations, "sampleLocations");
                break;
            case Seamass::Input::Type::Centroided:
                fileOut_->writeVector(input.locations, "centroidLocations");
                break;
            default:
                throw runtime_error("BUG: input has no type");
        }
    }

    fileOut_->flush();
}


bool DatasetSeamass::read(Seamass::Input &input, Seamass::Output &output, std::string &id)
{
    output = Seamass::Output();

    if (!read(input, id))
        return false;

    int groupId = fileIn_->openGroup("seamass");

    fileIn_->readAttribute(output.scale, "scale", "", groupId);
    output.lambda = fileIn_->readAttribute<double>("lambda", "", groupId);
    output.lambdaGroup = fileIn_->readAttribute<double>("lambdaGroup", "", groupId);
    output.tolerance = fileIn_->readAttribute<double>("tolerance", "", groupId);
    output.peakFwhm = fileIn_->readAttribute<double>("peakFwhm", "", groupId);
    output.chargeStates = fileIn_->readAttribute<short>("chargeStates", "", groupId);
    fileIn_->readAttribute(output.dbFilename, "dbFilename", "", groupId);

    {
        ii n = 0;
        for (;; n++)
        {
            ostringstream oss1; oss1 << setw(4) << setfill('0') << n << " X";
            if (fileIn_->openGroup(oss1.str(), groupId) == -1)
                break;
        }

        output.xs.resize(n);
        for (ii k = 0; k < n; k++)
        {
            ostringstream oss1; oss1 << setw(4) << setfill('0') << k << " X";
            fileIn_->readMatrixSparseCsr(output.xs[k], oss1.str(), groupId);
        }
    }

    {
        ii n = 0;
        for (;; n++)
        {
            ostringstream oss1; oss1 << setw(4) << setfill('0') << n << " L2";
            if (fileIn_->openGroup(oss1.str(), groupId) == -1)
                break;
        }

        output.l2s.resize(n);
        for (ii k = 0; k < n; k++)
        {
            ostringstream oss1; oss1 << setw(4) << setfill('0') << k << " L2";
            fileIn_->readMatrixSparseCsr(output.l2s[k], oss1.str(), groupId);
        }
    }

    {
        ii n = 0;
        for (;; n++)
        {
            ostringstream oss1; oss1 << setw(4) << setfill('0') << n << " L1L2";
            if (fileIn_->openGroup(oss1.str(), groupId) == -1)
                break;
        }

        output.l1l2s.resize(n);
        for (ii k = 0; k < n; k++)
        {
            ostringstream oss1; oss1 << setw(4) << setfill('0') << k << " L1L2";
            fileIn_->readMatrixSparseCsr(output.l1l2s[k], oss1.str(), groupId);
        }
    }

    id = "";
    return finished_ = true;
}


void DatasetSeamass::write(const Seamass::Input &input, const Seamass::Output &output, const std::string &id)
{
    write(input, id);

    int groupId = fileOut_->createGroup("seamass");

    fileOut_->writeAttribute(output.scale, "scale", "", groupId);
    fileOut_->writeAttribute(output.lambda, "lambda", "", groupId);
    fileOut_->writeAttribute(output.lambdaGroup, "lambdaGroup", "", groupId);
    fileOut_->writeAttribute(output.tolerance, "tolerance", "", groupId);
    fileOut_->writeAttribute(output.peakFwhm, "peakFwhm", "", groupId);
    fileOut_->writeAttribute(output.chargeStates, "chargeStates", "", groupId);
    fileOut_->writeAttribute(output.dbFilename, "dbFilename", "", groupId);

    if (output.b.size() > 0)
    {
        int matrixId = fileOut_->writeMatrixSparseCsr(output.b, "B", groupId);

        fileOut_->writeAttribute(output.bGridInfo.rowScale, "gridInfo.rowScale", "", matrixId);
        fileOut_->writeAttribute(output.bGridInfo.rowOffset, "gridInfo.rowOffset", "", matrixId);
        fileOut_->writeAttribute(output.bGridInfo.rowExtent, "gridInfo.rowExtent", "", matrixId);
        fileOut_->writeAttribute(output.bGridInfo.colScale, "gridInfo.colScale", "", matrixId);
        fileOut_->writeAttribute(output.bGridInfo.colOffset, "gridInfo.colOffset", "", matrixId);
        fileOut_->writeAttribute(output.bGridInfo.colExtent, "gridInfo.colExtent", "", matrixId);

        fileOut_->flush();
    }

    li nnz = 0;
    ii nk = ii(output.xs.size());
    for (ii k = 0; k < nk; k++)
    {
        ostringstream oss;
        oss << setw(5) << setfill('0') << k;
        int groupId2 = fileOut_->createGroup(oss.str(), groupId);
 
        fileOut_->writeAttribute(output.types[k], "basis.type", "", groupId2);
        fileOut_->writeAttribute(output.parents[k], "basis.parentIndex", "", groupId2);
        fileOut_->writeAttribute(output.gridInfos[k].rowScale, "gridInfo.rowScale", "", groupId2);
        fileOut_->writeAttribute(output.gridInfos[k].rowOffset, "gridInfo.rowOffset", "", groupId2);
        fileOut_->writeAttribute(output.gridInfos[k].rowExtent, "gridInfo.rowExtent", "", groupId2);
        fileOut_->writeAttribute(output.gridInfos[k].colScale, "gridInfo.colScale", "", groupId2);
        fileOut_->writeAttribute(output.gridInfos[k].colOffset, "gridInfo.colOffset", "", groupId2);
        fileOut_->writeAttribute(output.gridInfos[k].colExtent, "gridInfo.colExtent", "", groupId2);

        if (output.xs[k].size() > 0) {
            fileOut_->writeMatrixSparseCsr(output.xs[k], "X", groupId2);
            nnz += output.xs[k].nnz();
        }
        if (output.l2s[k].size() > 0) {
            fileOut_->writeMatrixSparseCsr(output.l2s[k], "L2", groupId2);
            nnz += output.l2s[k].nnz();
        }
        if (output.l1l2s[k].size() > 0) {
            fileOut_->writeMatrixSparseCsr(output.l1l2s[k], "L1L2", groupId2);
            nnz += output.l1l2s[k].nnz();
        }
        if (output.aTs[k]->size() > 0) {
            fileOut_->writeMatrixSparseCsr(*output.aTs[k], "At", groupId2);
            nnz += output.aTs[k]->nnz();
        }
        if (output.ids[k]->size() > 0) {
            fileOut_->writeVector<ii>(*output.ids[k], "ids", groupId2);
            nnz += output.ids[k]->size();
        }

        if (nnz > 1048576) {           
            fileOut_->flush();
            nnz = 0;

            if (getDebugLevel() % 10 >= 1) {
                ostringstream oss;
                oss << getTimeStamp() << "   Saved " << setw(1 + (int)(log10((float)nk))) << (k + 1) << " / " << nk;
                info(oss.str());
            }
            else {
                oss << "." << flush;
            }
        }
    }
}

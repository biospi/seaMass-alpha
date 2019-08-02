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
    
    // THIS IS ALL HARDCODED TO IMAGE2SML OUTPUT ATM
    fileIn_->readMatrixSparseCsr(input.b, "xScale=0");

    id = "";

    return finished_ = true;
}


void DatasetSeamass::write(const Seamass::Input &input, const std::string &id)
{
    /*if (input.startTimes.size() > 0)
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
    }*/
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
    fileIn_->readAttribute(output.isotopesFilename, "isotopesFilename", "", groupId);

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

    /*for (ii k = 0; k < (ii)output.xs.size(); k++)
    {
        if (output.xs[k])
        {
            ostringstream oss; oss << "xs[" << k << "]";
            fileOut_->write(*output.xs[k], oss.str(), grpid);
        }
    }

    for (ii k = 0; k < (ii)output.l2s.size(); k++)
    {
        if (output.l2s[k])
        {
            ostringstream oss; oss << "l2s[" << k << "]";
            fileOut_->write(*output.l2s[k], oss.str(), grpid);
        }
    }

    for (ii k = 0; k < (ii)output.l1l2s.size(); k++)
    {
        if (output.l1l2s[k])
        {
            ostringstream oss; oss << "l1l2s[" << k << "]";
            fileOut_->write(*output.l1l2s[k], oss.str(), grpid);
        }
    }*/

    /*fileIn_->read_AttNC("baselineScale", NC_GLOBAL, output.baselineScale, grpid);
    fileIn_->read_AttNC("baselineOffset", NC_GLOBAL, output.baselineOffset, grpid);
    fileIn_->read_AttNC("baselineExtent", NC_GLOBAL, output.baselineExtent, grpid);

    vector<double> shrinkage;
    fileIn_->read_AttNC("shrinkage", NC_GLOBAL, shrinkage, grpid);
    output.shrinkage = shrinkage[0];

    vector<double> tolerance;
    fileIn_->read_AttNC("tolerance", NC_GLOBAL, tolerance, grpid);
    output.tolerance = tolerance[0];

    if (fileIn_->read_VarIDNC("weights", grpid) != -1)
        fileIn_->read_VecNC("weights", output.weights, grpid);

    output.scales.resize(output.baselineExtent.size());
    output.offsets.resize(output.baselineExtent.size());
    for (ii d = 0; d < output.baselineExtent.size(); d++)
    {
        ostringstream oss1; oss1 << "scales[" << d << "]";
        fileIn_->read_VecNC(oss1.str(), output.scales[d], grpid);

        ostringstream oss2; oss2 << "offsets[" << d << "]";
        fileIn_->readVector(oss2.str(), output.offsets[d], grpid);
    }*/

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
    fileOut_->writeAttribute(output.isotopesFilename, "isotopesFilename", "", groupId);

    for (ii k = 0; k < ii(output.xs.size()); k++)
    {
        ostringstream oss;
        oss << setw(4) << setfill('0') << k << " X";
        int matrixId = fileOut_->writeMatrixSparseCsr(output.xs[k], oss.str(), groupId);

        fileOut_->writeAttribute(output.gridInfos[k].rowScale, "gridInfo.rowScale", "", matrixId);
        fileOut_->writeAttribute(output.gridInfos[k].rowOffset, "gridInfo.rowOffset", "", matrixId);
        fileOut_->writeAttribute(output.gridInfos[k].rowExtent, "gridInfo.rowExtent", "", matrixId);
        fileOut_->writeAttribute(output.gridInfos[k].colScale, "gridInfo.colScale", "", matrixId);
        fileOut_->writeAttribute(output.gridInfos[k].colOffset, "gridInfo.colOffset", "", matrixId);
        fileOut_->writeAttribute(output.gridInfos[k].colExtent, "gridInfo.colExtent", "", matrixId);
    }

    for (ii k = 0; k < ii(output.l2s.size()); k++)
    {
        ostringstream oss;
        oss << setw(4) << setfill('0') << k << " L2";
        fileOut_->writeMatrixSparseCsr(output.l2s[k], oss.str(), groupId);
    }

    for (ii k = 0; k < ii(output.l1l2s.size()); k++)
    {
        ostringstream oss;
        oss << setw(4) << setfill('0') << k << " L1L2";
        fileOut_->writeMatrixSparseCsr(output.l1l2s[k], oss.str(), groupId);
    }

    /*fileOut_->write_AttNC("", "baselineScale", output.baselineScale, NC_BYTE, grpid);

    vector<double> shrinkage(1); shrinkage[0] = output.shrinkage;
    fileOut_->write_AttNC("", "shrinkage", shrinkage, NC_DOUBLE, grpid);

    vector<double> tolerance(1); tolerance[0] = output.tolerance;
    fileOut_->write_AttNC("", "tolerance", tolerance, NC_DOUBLE, grpid);

    for (ii k = 0; k < output.baselineExtent.size(); d++)
    {
        ostringstream oss1; oss1 << "scales[" << d << "]";
        fileOut_->write_VecNC(oss1.str(), output.scales[d], NC_BYTE, grpid);

        ostringstream oss2; oss2 << "offsets[" << d << "]";
        fileOut_->write_VecNC(oss2.str(), output.offsets[d], sizeof(output.offsets[0]) == 4 ? NC_INT : NC_INT64, grpid);
    }*/

    /*fileOut_->write_AttNC("", "baselineOffset", output.baselineOffset, sizeof(output.baselineOffset[0]) == 4 ? NC_INT : NC_INT64, grpid);
    fileOut_->write_AttNC("", "baselineExtent", output.baselineExtent, sizeof(output.baselineExtent[0]) == 4 ? NC_INT : NC_INT64, grpid);

    if (output.weights.size() > 0)
    {
        fileOut_->write_VecNC("weights", output.weights, sizeof(output.weights[0]) == 4 ? NC_FLOAT : NC_DOUBLE, grpid);

        // ought to be a compound type rather than multiple datasets!
        for (ii d = 0; d < output.baselineExtent.size(); d++)
        {
            ostringstream oss1; oss1 << "scales[" << d << "]";
            fileOut_->write_VecNC(oss1.str(), output.scales[d], NC_BYTE, grpid);

            ostringstream oss2; oss2 << "offsets[" << d << "]";
            fileOut_->write_VecNC(oss2.str(), output.offsets[d], sizeof(output.offsets[0]) == 4 ? NC_INT : NC_INT64, grpid);
        }
    }*/
}





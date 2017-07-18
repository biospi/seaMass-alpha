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
using namespace kernel;


DatasetSeamass::DatasetSeamass(const std::string filePathIn, const std::string filePathStemOut, Dataset::WriteType writeType) : fileIn_(0), fileOut_(0), finished_(false)
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

    if (fileIn_->read_VarIDNC("countsIndex") != -1)
        fileIn_->read_VecNC("counts", input.counts);
    else
        throw runtime_error("ERROR: 'counts' dataset not found in smb file");

    if (fileIn_->read_VarIDNC("binLocations") != -1)
    {
        fileIn_->read_VecNC("binLocations", input.locations);
        input.type = Seamass::Input::Type::Binned;
    }
    else if (fileIn_->read_VarIDNC("sampleLocations") != -1)
    {
        fileIn_->read_VecNC("sampleLocations", input.locations);
        input.type = Seamass::Input::Type::Sampled;
    }
    else if (fileIn_->read_VarIDNC("centroidLocations") != -1)
    {
        fileIn_->read_VecNC("centroidLocations", input.locations);
        input.type = Seamass::Input::Type::Centroided;
    }
    else
        throw runtime_error("ERROR: one dataset called 'binLocations', 'sampleLocations' or 'centroidLocations' is needed in smb file");

    if (fileIn_->read_VarIDNC("countsIndex") != -1)
        fileIn_->read_VecNC("countsIndex", input.countsIndex);

    if (fileIn_->read_VarIDNC("startTimes") != -1)
        fileIn_->read_VecNC("startTimes", input.startTimes);

    if (fileIn_->read_VarIDNC("finishTimes") != -1)
        fileIn_->read_VecNC("finishTimes", input.finishTimes);

    if (fileIn_->read_VarIDNC("exposures") != -1)
        fileIn_->read_VecNC("exposures", input.exposures);

    id = "";

    return finished_ = true;
}


void DatasetSeamass::write(const Seamass::Input &input, const std::string &id)
{
    if (input.startTimes.size() > 0)
        fileOut_->write_VecNC("startTimes", input.startTimes, NC_DOUBLE);

    if (input.finishTimes.size() > 0)
        fileOut_->write_VecNC("finishTimes", input.finishTimes, NC_DOUBLE);

    if (input.exposures.size() > 0)
        fileOut_->write_VecNC("exposures", input.exposures, sizeof(input.counts[0]) == 4 ? NC_FLOAT : NC_DOUBLE);

    if (input.countsIndex.size() > 0)
        fileOut_->write_VecNC("countsIndex", input.countsIndex, sizeof(input.countsIndex[0]) == 4 ? NC_INT : NC_INT64);

    if (input.counts.size() > 0)
        fileOut_->write_VecNC("counts", input.counts, sizeof(input.counts[0]) == 4 ? NC_FLOAT : NC_DOUBLE);

    if (input.locations.size() > 0)
    {
        switch (input.type)
        {
            case Seamass::Input::Type::Binned:
                fileOut_->write_VecNC("binLocations", input.locations, NC_DOUBLE);
                break;
            case Seamass::Input::Type::Sampled:
                fileOut_->write_VecNC("sampleLocations", input.locations, NC_DOUBLE);
                break;
            case Seamass::Input::Type::Centroided:
                fileOut_->write_VecNC("centroidLocations", input.locations, NC_DOUBLE);
                break;
            default:
                throw runtime_error("BUG: input has no type");
        }
    }
}


bool DatasetSeamass::read(Seamass::Input &input, Seamass::Output &output, std::string &id)
{
    output = Seamass::Output();

    if (!read(input, id))
        return false;

    int grpid = fileIn_->open_Group("seamass");

    fileIn_->read_AttNC("scale", NC_GLOBAL, output.scale, grpid);

    vector<double> shrinkage;
    fileIn_->read_AttNC("shrinkage", NC_GLOBAL, shrinkage, grpid);
    output.shrinkage = shrinkage[0];

    vector<double> tolerance;
    fileIn_->read_AttNC("tolerance", NC_GLOBAL, tolerance, grpid);
    output.tolerance = tolerance[0];

    vector<double> peakFwhm;
    fileIn_->read_AttNC("peakFwhm", NC_GLOBAL, peakFwhm, grpid);
    output.peakFwhm = peakFwhm[0];

    {
        ii n = 0;
        for (;; n++)
        {
            ostringstream oss1; oss1 << "xs[" << n << "]";
            if (fileIn_->open_Group(oss1.str(), grpid) == -1)
                break;
        }

        output.xs.resize(n);
        for (ii k = 0; k < n; k++)
        {
            ostringstream oss1; oss1 << "xs[" << k << "]";
            fileIn_->read(output.xs[k], oss1.str(), grpid);
        }
    }

    {
        ii n = 0;
        for (;; n++)
        {
            ostringstream oss1; oss1 << "l2s[" << n << "]";
            if (fileIn_->open_Group(oss1.str(), grpid) == -1)
                break;
        }

        output.l2s.resize(n);
        for (ii k = 0; k < n; k++)
        {
            ostringstream oss1; oss1 << "l2s[" << k << "]";
            fileIn_->read(output.l2s[k], oss1.str(), grpid);
        }
    }

    {
        ii n = 0;
        for (;; n++)
        {
            ostringstream oss1; oss1 << "l1l2s[" << n << "]";
            if (fileIn_->open_Group(oss1.str(), grpid) == -1)
                break;
        }

        output.l1l2s.resize(n);
        for (ii k = 0; k < n; k++)
        {
            ostringstream oss1; oss1 << "l1l2s[" << k << "]";
            fileIn_->read(output.l1l2s[k], oss1.str(), grpid);
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
        fileIn_->read_VecNC(oss2.str(), output.offsets[d], grpid);
    }*/

    id = "";

    return finished_ = true;
}


void DatasetSeamass::write(const Seamass::Input &input, const Seamass::Output &output, const std::string &id)
{
    write(input, id);

    int grpid = fileOut_->create_Group("seamass");

    fileOut_->write_AttNC("", "scale", output.scale, NC_BYTE, grpid);

    vector<double> shrinkage(1); shrinkage[0] = output.shrinkage;
    fileOut_->write_AttNC("", "shrinkage", shrinkage, NC_DOUBLE, grpid);

    vector<double> tolerance(1); tolerance[0] = output.tolerance;
    fileOut_->write_AttNC("", "tolerance", tolerance, NC_DOUBLE, grpid);

    vector<double> peakFwhm(1); peakFwhm[0] = output.peakFwhm;
    fileOut_->write_AttNC("", "peakFwhm", peakFwhm, NC_DOUBLE, grpid);

    for (ii k = 0; k < (ii)output.xs.size(); k++)
    {
        ostringstream oss; oss << "xs[" << k << "]";
        fileOut_->write(output.xs[k], oss.str(), grpid);
    }

    for (ii k = 0; k < (ii)output.l2s.size(); k++)
    {
        ostringstream oss; oss << "l2s[" << k << "]";
        fileOut_->write(output.l2s[k], oss.str(), grpid);
    }

    for (ii k = 0; k < (ii)output.l1l2s.size(); k++)
    {
        ostringstream oss; oss << "l1l2s[" << k << "]";
        fileOut_->write(output.l1l2s[k], oss.str(), grpid);
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





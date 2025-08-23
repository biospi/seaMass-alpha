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


#include "DatasetLibrary.hpp"
#include <kernel.hpp>
#include <iomanip>
using namespace kernel;


DatasetLibrary::DatasetLibrary(const std::string& filePathIn, const std::string& filePathStemOut, Dataset::WriteType writeType) : fileIn_(0), fileOut_(0), finished_(false)
{
    if (!filePathIn.empty())
        fileIn_ = new FileNetcdf(filePathIn);

    if (!filePathStemOut.empty())
        fileOut_ = new FileNetcdf(filePathStemOut + (writeType == Dataset::WriteType::InputOutput ? ".sml" : ".sml"), NC_NETCDF4);
}


DatasetLibrary::~DatasetLibrary()
{
    if (fileIn_)
        delete fileIn_;

    if (fileOut_)
        delete fileOut_;
}


bool DatasetLibrary::read(Seamass::Input &input, std::string &id)
{
    if (finished_ == true)
        return false;

    input = Seamass::Input();

    ostringstream oss2;
    oss2 << "m" << setfill('0') << setw(2) << 10;

    // read in b
    fileIn_->readMatrixSparseCsr(input.b, oss2.str());

    //  read in offset
    input.offset = fileIn_->readAttribute<ii>("offset", "", fileIn_->openGroup(oss2.str()));

    id = "m10";

    return finished_ = true;
}


void DatasetLibrary::write(const Seamass::Input &input, const std::string &id)
{
}


bool DatasetLibrary::read(Seamass::Input &input, Seamass::Output &output, std::string &id)
{
    return false;
}


void DatasetLibrary::write(const Seamass::Input &input, const Seamass::Output &output, const std::string &id)
{
}

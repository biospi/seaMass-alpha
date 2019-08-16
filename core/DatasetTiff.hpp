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


#ifndef SEAMASS_DATASETTIFF_HPP
#define SEAMASS_DATASETTIFF_HPP


#include "Dataset.hpp"
#include "../io/FileNetcdf.hpp"
#include <tiffio.h>


class DatasetTiff: public Dataset
{
public:
    DatasetTiff(const std::string& filePathIn, const std::string& filePathStemOut, Dataset::WriteType writeType = Dataset::WriteType::InputOutput);
    virtual ~DatasetTiff();

    virtual bool read(std::string& filePathSml, std::string &id);
    virtual void write(const std::string& filePathSml, const std::string &id);

    virtual bool read(std::string& filePathSml, Seamass::Output &output, std::string &id);
    virtual void write(const std::string& filePathSml, const Seamass::Output &output, const std::string &id);

private:
    string filePathSml_;
    TIFF* fileIn_;
    TIFF* fileOut_;
    bool finished_;
    //uint32 width;
    //uint32 height;
};


#endif //SEAMASS_DATASETTIFF_HPP

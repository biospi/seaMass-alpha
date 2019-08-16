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

#ifndef SEAMASS_DATASET_HPP
#define SEAMASS_DATASET_HPP


#include "Seamass.hpp"


class Dataset
{
public:
    virtual ~Dataset();

    enum class WriteType { InputOutput, Input };

    virtual bool read(std::string& smlFilename, std::string &id) = 0;
    virtual void write(const std::string& smlFilename, const std::string &id) = 0;

    virtual bool read(std::string& smlFilename, Seamass::Output &output, std::string &id) = 0;
    virtual void write(const std::string& smlFilename, const Seamass::Output &output, const std::string &id) = 0;
};


class FileFactory
{
public:
    static Dataset* createFileObj(const std::string filePathIn, const std::string filePathStemOut, Dataset::WriteType writeType = Dataset::WriteType::InputOutput);
};


#endif

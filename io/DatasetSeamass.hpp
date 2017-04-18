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

#ifndef SEAMASS_DATASETSEAMASSCORE_HPP
#define SEAMASS_DATASETSEAMASSCORE_HPP


#include "Dataset.hpp"
#include "../kernel/FileNetcdf.hpp"


class DatasetSeamass: public Dataset
{
public:
	DatasetSeamass(std::string &fileName, bool onlyWrite = false);
	virtual ~DatasetSeamass();

	virtual bool read(Seamass::Input &input, std::string &id);
	virtual bool read(Seamass::Input &input, Seamass::Output &output, std::string &id);

    virtual void write(const Seamass::Input& input, const std::string& id);
    virtual void write(const Seamass::Input& input, const Seamass::Output& output, const std::string& id);

private:
    std::string fileName_;
    FileNetcdf* fileIn_;
	FileNetcdf* fileOut_;
 	bool processed_;
};


#endif

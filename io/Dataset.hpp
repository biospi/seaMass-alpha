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

#ifndef _SEAMASS_INPUTFILE_HPP_
#define _SEAMASS_INPUTFILE_HPP_


#include "../core/SeamassCore.hpp"


class Dataset;

class FileFactory
{
public:
    static Dataset* createFileObj(std::string fName);
};


class Dataset
{
public:
    virtual ~Dataset() {};

    virtual bool next(SeamassCore::Input& output, std::string& id) = 0;
	virtual void writeData(SeamassCore& sm_, SeamassCore::Input& input_, bool centriod_, double threshold_) = 0;
};


#endif // _SEAMASS_INPUTFILE_HPP_

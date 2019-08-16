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

#include "Dataset.hpp"
#include "DatasetMzmlb.hpp"
#include "DatasetSeamass.hpp"
#include "DatasetTiff.hpp"
using namespace std;


Dataset* FileFactory::createFileObj(const std::string filePathIn, const std::string filePathStemOut, Dataset::WriteType writeType)
{
    size_t pos = filePathIn.find_last_of(".");
    if (pos != string::npos)
    {
        string ext = filePathIn.substr(pos);

        if(ext == ".smb" || ext == ".smv")
        {
            return new DatasetSeamass(filePathIn, filePathStemOut, writeType);
        }
        else if(ext == ".mzMLb" || ext == ".mzMLv")
        {
            return new DatasetMzmlb(filePathIn, filePathStemOut, writeType);
        }
        else if(ext == ".tiff" || ext == ".tiff")
        {
            return new DatasetTiff(filePathIn, filePathStemOut, writeType);
        }
    }

    return 0;
}


Dataset::~Dataset()
{
}




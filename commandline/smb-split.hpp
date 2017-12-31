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

#ifndef SEAMASS_SMB_SPLIT_HPP
#define SEAMASS_SMB_SPLIT_HPP

#include <boost/filesystem.hpp>

void createSmbTiles(int sections, const boost::filesystem::path &fileIn);

struct SmbTile
{
    SmbTile(boost::filesystem::path _file, vector<double> _vec, vector<li> cs):fileName(_file),
        startTime(_vec),countsIndex(cs),id(0),begIdx(0),endIdx(0),offset(0) {};

    boost::filesystem::path fileName;
    vector<double> startTime;
    vector<li> countsIndex;
    vector<fp> exposures;
    vector<double> finishTimes;
    int id;
    li begIdx;
    li endIdx;
    li len;
    li csBeg;
    li csEnd;
    li csLen;
    li csIdxLen;
    li mzBeg;
    li mzEnd;
    li mzLen;
    li offset;
};

template<typename T>
T getSmbIndex(string fileIdx);

template<typename T>
T getSmbIndex(string fileIdx)
{
    T index;
    size_t idx=fileIdx.find_last_of(".");

    istringstream(fileIdx.substr(idx+1))>>index;

    return index;
};

#endif //SEAMASS_SMB_SPLIT_HPP

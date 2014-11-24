//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> manchester.ac.uk>
//
// Copyright (C) 2013  CADET Bioinformatics Laboratory, University of Manchester, UK
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


#ifndef _SEAMASSRESTORATION_SMOFILE_HPP_
#define _SEAMASSRESTORATION_SMOFILE_HPP_


#include "core.hpp"


class SMOFile
{
protected:
    string filename;
    int file;
    
public:
	SMOFile(const string& filename);
	~SMOFile();
    
    void write_cs(const string& objectname,
                  const CoeffsMetadata& cm,
                  const vector<fp>& cs) const;
    
    void write_fs(const string& objectname,
                  const vector<fp>& fs,
                  const vector<li>& is,
                  const vector<ii>& js) const;
};


#endif // _SEAMASSRESTORATION_SMOFILE_HPP_


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


#ifndef _SEAMASSRESTORATION_HPP_
#define _SEAMASSRESTORATION_HPP_


#include <vector>
#include "HDF5File.hpp"

class seaMassRestoration
{
private:
    const std::string& id;
    HDF5File h5out;
    
    int instrument_type;
    int debug;

public:

    seaMassRestoration(const std::string& id,
                       int instrument_type,
                       int debug);
    
    void process(const std::string& config_id,
                 std::vector<double>& rts,
                 std::vector< std::vector<double> >& mzs,
                 std::vector< std::vector<double> >& intensities,
                 int rc0_mz, int rc1_mz,
                 int rc0_rt, int rc1_rt,
                 int shrinkage0, int shrinkage1,
                 int tolerance0, int tolerance1);
    
};


#endif // _SEAMASSRESTORATION_HPP_


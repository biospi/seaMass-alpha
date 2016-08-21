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


#ifndef _SEAMASSRESTORATION_AUX_HPP_
#define _SEAMASSRESTORATION_AUX_HPP_


#include <vector>
#include <mkl.h>

typedef float fp; // fp is the selected floating point precision
typedef MKL_INT ii; // ii is the selected indexing integer size
typedef long long li;


void remove_zeros(std::vector< std::vector<fp> >& mzs, std::vector< std::vector<fp> >& intensities);
void merge_bins(std::vector< std::vector<fp> >& mzs, std::vector< std::vector<fp> >& intensities, double width);


#endif // _SEAMASSRESTORATION_AUX_HPP_


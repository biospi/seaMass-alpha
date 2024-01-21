//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> bristol.ac.uk>
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


#ifndef SEAMASS_KERNEL_INTEL_TYPES_HPP
#define SEAMASS_KERNEL_INTEL_TYPES_HPP


// ilp64 64 bit addressing NOT IMPLEMENTED in seamass_kernel/intel yet (and might never be as it increases memory usage)
// Note the lp64 or ilp64 define is handled in CMakeLists as it has to link in the correct library too
#include <mkl.h>
#include <mkl_spblas.h>
typedef float fp; // fp is the selected floating point precision (float or double)
typedef MKL_INT ii; // ii is the selected addressing (32bit for lp64 or 64 bit for ilp64)
typedef MKL_INT64 li; // li is always 64 bit


#endif


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


#ifndef _SEAMASS_MATH_MKL_HPP_
#define _SEAMASS_MATH_MKL_HPP_


//#include "MatrixMKL.hpp" //Not sure why you have included this here... Recursive?!?!

#include <string>

//#define MKL_ILP64 // use 64 bit addressing (comment out for 32 bit)
#include <mkl.h>
#include <mkl_spblas.h>


typedef float fp; // fp is the selected floating point precision (float or double)
typedef MKL_INT ii; // ii is the selected addressing (32 or 64 bit)
typedef MKL_INT64 li; // li is always 64 bit


std::string getThreadInfo();
li getId();
void resetElapsedTime();
double getElapsedTime();
li getUsedMemory();
std::string getTimeStamp();
void setDebugLevel(int debugLevel);
int getDebugLevel();


#endif


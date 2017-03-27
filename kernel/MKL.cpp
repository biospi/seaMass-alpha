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


#include "MatrixSparseMKL.hpp"

#include <iomanip>
#include <iostream>
//#include <cassert>
//#include <cstring>
//#include <cmath>
#include <sstream>
//#include <algorithm>

#if defined(_OPENMP)
  #include <omp.h>
#endif

//#include <ippcore.h>
//#include <ipps.h>
//#include <ippi.h>

//#include "../kernel/NetcdfFile.hpp"


using namespace std;


std::string getThreadInfo()
{
    ostringstream out;
    out << "Config: " << 8 * sizeof(ii) << "bit MKL addressing, " << mkl_get_max_threads() << " MKL threads, ";
#if defined(_OPENMP)
    out << omp_get_max_threads() << " OpenMP threads";
#else
    out << "non-OpenMP build";
#endif
    return out.str();
}


static li id_ = -1;


li getId()
{
    return id_;
}


static double startTime_ = dsecnd();


void resetElapsedTime()
{
    startTime_ = dsecnd();
}


double getElapsedTime()
{
    return dsecnd() - startTime_;
}


li getUsedMemory()
{
    int allocatedBuffers;
    return mkl_mem_stat(&allocatedBuffers);
}


string getTimeStamp()
{
    ostringstream out;
    out << "[" << setw(9) << ++id_ << "," << fixed << internal << setw(9) << std::setprecision(3) << getElapsedTime() << "," << setw(9) << getUsedMemory()/1024.0/1024.0 << "] ";
    return out.str();
}


static int debugLevel_ = 0;


void setDebugLevel(int debugLevel)
{
    debugLevel_ = debugLevel;
}


int getDebugLevel()
{
    return debugLevel_;
}


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


#ifndef _SEAMASS_HPP_
#define _SEAMASS_HPP_

#if defined _WIN32 || defined __CYGWIN__
  #define SEAMASS_HELPER_DLL_IMPORT __declspec(dllimport)
  #define SEAMASS_HELPER_DLL_EXPORT __declspec(dllexport)
  #define SEAMASS_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define SEAMASS_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define SEAMASS_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define SEAMASS_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define SEAMASS_HELPER_DLL_IMPORT
    #define SEAMASS_HELPER_DLL_EXPORT
    #define SEAMASS_HELPER_DLL_LOCAL
  #endif
#endif

// Now we use the generic helper definitions above to define SEAMASS_API and SEAMASS_LOCAL.
// SEAMASS_API is used for the public API symbols. It either DLL imports or DLL exports (or does nothing for static build)
// SEAMASS_LOCAL is used for non-api symbols.

#ifdef SEAMASS_DLL // defined if SEAMASS is compiled as a DLL
  #ifdef SEAMASS_DLL_EXPORTS // defined if we are building the SEAMASS DLL (instead of using it)
    #define SEAMASS_API SEAMASS_HELPER_DLL_EXPORT
  #else
    #define SEAMASS_API SEAMASS_HELPER_DLL_IMPORT
  #endif // SEAMASS_DLL_EXPORTS
  #define SEAMASS_LOCAL SEAMASS_HELPER_DLL_LOCAL
#else // SEAMASS_DLL is not defined: this means SEAMASS is a static lib.
  #define SEAMASS_API
  #define SEAMASS_LOCAL
#endif // SEAMASS_DLL


#include <vector>

namespace seamass
{

SEAMASS_API
void notice();

SEAMASS_API
void process(const std::string& id,
             const std::string& config_id,
			 int instrument_type,
             std::vector<double>& rts,
             std::vector< std::vector<double> >& mzs,
             std::vector< std::vector<double> >& intensities,
             int mz0, int mz1,
             int rt0, int rt1,
             int shrinkage0, int shrinkage1,
             int tolerance0, int tolerance1,
		     int threads = 1, int debug = 0);

}

#endif // _SEAMASS_HPP_


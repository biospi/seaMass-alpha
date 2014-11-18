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


#ifndef _SEAMASS_RESTORATION_IMPL_HPP_
#define _SEAMASS_RESTORATION_IMPL_HPP_


#include "core.hpp"
#include "BasisFunctions.hpp"
#include "HDF5File.hpp"


class OptimiserASRL
{
protected:
    const vector<Basis*>& bases;
    vector<fp>& gs;
 
    vector< vector<fp> > cs;
    vector< vector<fp> > wcs;
    vector< vector<fp> > l2;
    
    ii accell;
    vector< vector<fp> > c0s;
    vector< vector<fp> > u0s;
    vector< vector<fp> > q0s;
    
    ii iteration;

public:    
	OptimiserASRL(const vector<Basis*>& _bases,
                  vector<fp>& gs,
                  ii accell = 2);
	~OptimiserASRL();
    
    double step(ii iteration, double strength);
    void threshold(double threshold);
    vector< vector<fp> >& get_cs() { return cs; }
    
    void write_h5(const HDF5File& file, const string& datafilename,
                  const vector<ii>& scale_bases, const vector<li>& is, const vector<ii>& js);
    void calc_error(const std::string& id);
};


#endif // _SPECTRUMLIST_SEAMASSRESTORATION_HPP_ 


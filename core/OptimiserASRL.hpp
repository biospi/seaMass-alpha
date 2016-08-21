//
// $Id$
//
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


#ifndef _SEAMASS_RESTORATION_IMPL_HPP_
#define _SEAMASS_RESTORATION_IMPL_HPP_


#include "BasisFunctions.hpp"


class OptimiserASRL
{
protected:
    const std::vector<Basis*>& bases;
	std::vector<fp>& gs;
 
	std::vector< std::vector<fp> > cs;
	std::vector< std::vector<fp> > wcs;
	std::vector< std::vector<fp> > l2;
    
    ii accell;
	std::vector< std::vector<fp> > c0s;
	std::vector< std::vector<fp> > u0s;
	std::vector< std::vector<fp> > q0s;
    
    ii iteration;

public:    
	OptimiserASRL(const std::vector<Basis*>& _bases, std::vector<fp>& gs, ii accell = 2);
	~OptimiserASRL();
    
	double step(ii iteration, fp lambda);

	void synthesis(std::vector<fp>& fs, ii return_cs = -1) const;
	void error(std::vector<fp>& fs) const;
	void analysis(std::vector< std::vector<fp> >& dcs, const std::vector<fp>& dfs) const;
	void shrinkage(std::vector< std::vector<fp> >& dcs, fp lambda);

    void threshold(double threshold);
	std::vector< std::vector<fp> >& get_cs() { return cs; }
	std::vector<fp>& get_gs() { return gs; }

	/*fp compute_norm_max_counts(ii n_core_bases = 0);
    
    void write_h5(const SMOWriter& file, const string& datafilename,
                  const vector<ii>& scale_bases, const vector<li>& is, const vector<ii>& js, const vector<fp>& gains);
    void calc_error(const std::string& id);*/
};


#endif // _SPECTRUMLIST_SEAMASSRESTORATION_HPP_ 


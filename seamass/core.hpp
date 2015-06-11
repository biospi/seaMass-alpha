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
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <float.h>
#include <cstdlib>
#include <omp.h>
#include <mkl.h>


using namespace std;


typedef float fp; // fp is the selected floating point precision
typedef MKL_INT ii; // ii is the selected indexing integer size
typedef long long li;


struct CoeffsMetadata
{
    CoeffsMetadata(ii d);
    ~CoeffsMetadata();
    
    ii size() const;
    void print(ostream& out) const;
    void operator=(const CoeffsMetadata& cm);
    
    ii d;        // dimension of the output coefficients
    vector<ii> l; // coefficient dyadic level for each dimension
    vector<ii> o; // coefficient offset
    vector<ii> n; // number of coefficients per set
};


void bin_mzs_intensities(vector< vector<double> >& mzs,
                         vector< vector<double> >& intensities,
                         ii instrument_type,
						 vector<fp>& gains);

void merge_bins(vector< vector<double> >& mzs,
                vector< vector<double> >& intensities,
				double width);

void create_gs(vector<fp>& gs,
               vector<li>& is,
               vector<ii>& js,
               const vector< vector<double> >& intensities);


namespace bspline
{
    double m(double x, int k, int i, vector<fp>& ks);
    double m(double x, int k, int i);
    double im(double x, int k);
    int factorial(int n);
}


#endif // _SEAMASSRESTORATION_AUX_HPP_


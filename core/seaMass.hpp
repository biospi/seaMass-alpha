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


#ifndef _SEAMASS_HPP_
#define _SEAMASS_HPP_


#include <vector>
#include <mkl.h>


typedef float fp; // fp is the selected floating point precision
typedef MKL_INT ii; // ii is the selected indexing integer size
typedef long long li;


/**
* seaMass fitting of a 1-d curve or 2-dimensional surface to the input spectr(um|a).
*/
class seaMass
{
public:
	static void notice();

	struct Input {
		std::vector<double> bin_counts;
		std::vector<double> bin_locations;
		std::vector<li> spectrum_index;
		std::vector<double> start_times;
		std::vector<double> finish_times;
		std::vector<double> exposures;
	};

	struct Output {
		std::vector<fp> coeffs;
		std::vector<short> levels;
		std::vector<long> offsets;
	};

	/**
	* Initialise seaMass fitting of a 1-d curve or 2-dimensional surface to the input spectr(um|a).
	* param mzs m/z values for the input spectra, where each value denotes the start m/z of the respective intensity bin; may be modified to save momo
	* param intensities ion counts for the input spectra, where each vector must have one less value than the respective vector in mzs (unless both are empty).
	* param stan_times Start scan times for the input spectra, must contain one more item than that of mzs and intensities (denoting the end scan time of the last spectrum)
	* param mz_res
	* param st_res
	* param shrinkage
	* param tolerance
	*/
	seaMass(Input& input, std::vector<int>& resolutions, double shrinkage, double tolerance);

	virtual bool iteration();

	void get_output(Output& output);

private:
	void create_gs(std::vector<fp>& gs, std::vector<li>& is, std::vector<ii>& js, const std::vector< std::vector<double> >& bin_locations);
};


#endif // _SEAMASS_HPP_


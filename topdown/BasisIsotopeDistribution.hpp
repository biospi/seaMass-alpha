//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  biospi Laboratory, EEE, University of Liverpool, UK
//
// This file is part of seaMass-TD.
//
// seaMass-TD is free software: you can redistribute it and/or modify
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
// along with seaMass-TD.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef _SEAMASS_BASISISOTOPEDISTRIBUTION_HPP_
#define _SEAMASS_BASISISOTOPEDISTRIBUTION_HPP_


#include "../core/BSpline.hpp"
#include "BasisChargeDistribution.hpp"


/*class BasisIsotopeDistribution : public Basis
{
private:
	// Input
	BasisChargeDistribution* parent;

	// Sparse basis matrix A (same for each input spectrum)
	SparseMatrix a;

	// Output (cs) metadata
	li ns, nc; // number of spectra and number of coefficients
	std::vector<ii> ois; //output indicies into coefficients
	std::vector<ii> gis; //group indicies into coefficients
	std::vector<fp> ave_masses;

public:
	BasisIsotopeDistribution(std::vector<Basis*>& bases, BasisChargeDistribution* parent,
		                     ii out_res, ii factor_res, ii max_z,
							 bool transient = false);

	~BasisIsotopeDistribution();

	void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true) const;
	void analysis(std::vector<fp>& es, const std::vector<fp>& fs) const;
	void l2norm(std::vector<fp>& es, const std::vector<fp>& fs) const;
	void shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage) const;
	li get_nc() const { return nc; }

	void write_cs(const std::vector<fp>& cs) const;
	void restrict_range(std::vector<fp>& cs, double mass0, double mass1) const;
};*/


#endif


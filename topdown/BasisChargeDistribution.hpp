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


#ifndef _SEAMASS_TOPDOWN_BASISCHARGEDISTRIBUTION_HPP_
#define _SEAMASS_TOPDOWN_BASISCHARGEDISTRIBUTION_HPP_


#include "../asrl/Basis.hpp"
#include "../core/Bspline.hpp"


/*class BasisChargeDistribution : public Basis
{
public:
	BasisChargeDistribution(std::vector<Basis*>& bases, const std::vector<fp>& counts, ii scale, ii offset, ii maxMass, ii binsPerDalton, bool isTransient = false);

	~BasisChargeDistribution();

	void synthesise(Matrix& f, const Matrix& x, bool accumulate) const;
	void analyse(Matrix& xE, const Matrix& fE, bool sqrA = false) const;

	ii getM() const;
	ii getN() const;

private:
	MatrixSparse a_;

	std::vector<li> groupXOffsets_;
	std::vector<short> groupZOffsets_;

	static double protonMass_; // mass of a proton (positive charge) in Daltons
};*/


#endif


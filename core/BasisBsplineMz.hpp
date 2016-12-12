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


#ifndef _SEAMASS_CORE_BASISBSPLINEMZ_HPP_
#define _SEAMASS_CORE_BASISBSPLINEMZ_HPP_


#include "BasisBspline.hpp"


class BasisBsplineMz : public BasisBspline
{
public:
	BasisBsplineMz(std::vector<Basis*>& bases, const std::vector<fp>& binCounts, const std::vector<li>& spectrumIndex,
		           const std::vector<double>& binEdges, short resolution, ii order = 3, bool isTransient = false);
	virtual ~BasisBsplineMz();

	void synthesis(MatrixSparse& f, const MatrixSparse& x, bool accumulate) const;
	void analysis(MatrixSparse& xE, const MatrixSparse& fE, bool sqrA = false) const;

private:
	std::vector<MatrixSparse> as_; // CSR sparse 'A' basis matrices
	std::vector<li> is_; // spectrum_index into 'g'
};


#endif


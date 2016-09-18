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


#ifndef _SEAMASS_CORE_BASISBSPLINE_HPP_
#define _SEAMASS_CORE_BASISBSPLINE_HPP_


#include "Basis.hpp"


class BasisBspline : public Basis
{
public:
	struct MeshInfo
	{
		ii dimensions;  	    // dimension of each b-spline coefficients mesh
		std::vector<ii> scale;  // dyadic scale for each dimension
		std::vector<ii> offset; // coefficient offset for each dimension
		std::vector<ii> extent; // number of coefficients for each dimension (make up the columns)
		ii n;                   // number of meshes (columns in resulting matrix)

		MeshInfo(ii dims);
		~MeshInfo();

		ii m() const;           // number of coefficients in a mesh (rows in resulting matrix)
		li size() const;        // number of coefficients across all meshes

		void operator=(const MeshInfo& meshInfo);
	};

	BasisBspline(std::vector<Basis*>& bases, ii dims, bool isTransient = false, ii parentIndex = -1);
	virtual ~BasisBspline();

	ii getM() const;
	ii getN() const;
	const MeshInfo& getMeshInfo() const;

protected:
	MeshInfo& meshInfo();

private:
	MeshInfo meshInfo_;
};

std::ostream& operator<<(std::ostream& os, const BasisBspline::MeshInfo& meshInfo);

#endif


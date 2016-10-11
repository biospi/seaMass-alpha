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


#include "BasisBsplineScale.hpp"

#include "Bspline.hpp"

#include <limits>
#include <iomanip>
#include <cmath>


using namespace std;


BasisBsplineScale::
BasisBsplineScale(vector<Basis*>& bases, ii parentIndex, ii _dimension, ii order, bool transient)
	: BasisBspline(bases, static_cast<BasisBspline*>(bases[parentIndex])->getMeshInfo().dimensions, transient, parentIndex), dimension(_dimension)
{
	const MeshInfo parentMeshInfo = static_cast<BasisBspline*>(bases[parentIndex])->getMeshInfo();
	meshInfo() = parentMeshInfo;
	meshInfo().scale[dimension] = parentMeshInfo.scale[dimension] - 1;
	meshInfo().offset[dimension] = parentMeshInfo.offset[dimension] / 2;
	meshInfo().extent[dimension] = (parentMeshInfo.offset[dimension] + parentMeshInfo.extent[dimension] - 1 - order) / 2 + order + 1 - meshInfo().offset[dimension];
	ii m = parentMeshInfo.extent[dimension];
	ii n = meshInfo().extent[dimension];

	ii stride = 1;
	for (ii j = 0; j < dimension; j++) stride *= meshInfo().extent[j];

	// create our kernel
	ii nh = order + 2;
	vector<fp> hs(nh);
	double sum = 0.0;
	for (ii i = 0; i < nh; i++)
	{
		hs[i] = (fp) (1.0 / pow(2.0, (double)order) * Bspline::factorial(order + 1) / (double)(Bspline::factorial(i)*Bspline::factorial(order + 1 - i)));
		sum += hs[i];
	}
	for (ii i = 0; i < nh; i++)
	{
		hs[i] /= (fp) sum;
	}

	// create A as a temporary COO matrix
	vector<fp> acoo(nh * n);
	vector<ii> rowind(nh * n);
	vector<ii> colind(nh * n);

	ii nnz = 0;
	ii offset = order + ((parentMeshInfo.offset[dimension] + 1) % 2);
	for (ii j = 0; j < n; j++)
	{
		for (ii i = 0; i < nh; i++)
		{
			rowind[nnz] = 2 * j + i - offset;
			if (rowind[nnz] < 0 || rowind[nnz] >= m) continue;
			acoo[nnz] = hs[i];
			colind[nnz] = j;

			nnz++;
		}
	}

	a.init(m, n, nnz, acoo.data(), rowind.data(), colind.data());
	aT.init(n, m, nnz, acoo.data(), colind.data(), rowind.data());

#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineScale";
	if (isTransient()) cout << " (t)";
	cout << " parent=" << getParentIndex() << " dimension=" << dimension << " " << meshInfo() << endl;
	cout << "  A" << a << " (";
    cout.unsetf(std::ios::floatfield);
    cout << setprecision(2) << (a.mem() + aT.mem()) / 1024.0 / 1024.0 << "Mb)" << endl;
    //cout << "  A" << a << " (" << defaultfloat << setprecision(2) << (a.mem() + aT.mem()) / 1024.0 / 1024.0 << "Mb)" << endl;
#endif
}


BasisBsplineScale::~BasisBsplineScale()
{
}


void
BasisBsplineScale::
synthesis(Matrix& f, const Matrix& c, bool accumulate) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineScale::synthesis" << endl;
#endif

	f.mul(a, c, accumulate, false);
}


void
BasisBsplineScale::
analysis(Matrix& cE, const Matrix& fE, bool sqrA) const
{
#ifndef NDEBUG
	cout << " " << getIndex() << " BasisBsplineScale::analysis" << endl;
#endif

	if (sqrA)
	{
		MatrixSparse aT_sqrd; aT_sqrd.elementwiseSqr(aT);
		cE.mul(aT_sqrd, fE, false, false);
	}
	else
	{
		cE.mul(aT, fE, false, false);
	}
}

//
// $Id$
//
//
// Original author: Ranjeet Bhamber <ranjeet <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Liverpool, UK
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

#ifndef SMPEAK_SMDATA_HPP_
#define SMPEAK_SMDATA_HPP_

#include "peakcore.hpp"


template<typename T=float, typename R=double>
struct DataAxis
{
	vector<R> rt;
	vector<double> mz;
	VecMat<T>* alpha;
protected:
	~DataAxis(){};
};


template
<
	template<class Operator> class MathOp,
	typename T=float,
	typename R=pair<int,double>
>
struct SMData1D : public DataAxis<T,R>, public MathOp<T>
{
	SMData1D(hsize_t dims[], int offset[], double mz_res, R rt,
			vector<T> &vec);
	~SMData1D(){delete this->alpha;};
};


template
<
	template<class Operator> class MathOp,
	typename T=float,
	typename R=double
>
struct SMData2D : public DataAxis<T,R>, public MathOp<T>
{
	SMData2D(vector<R> &_rt, vector<double> &_mz, vector<T> &vec);
	SMData2D(hsize_t dims[], int offset[], double mz_res, double rt_res,
			vector<T> &vec);
	~SMData2D(){delete this->alpha;};
};


template<template<class Operator> class MathOp, typename T, typename R>
SMData1D<MathOp,T,R>::SMData1D(hsize_t dims[], int offset[],
		double mz_res, R rt, vector<T> &vec)
{
	DataAxis<T,R>::alpha = new VecMat<T>(dims[0],dims[1],vec);
	apply(dims[0],dims[1],DataAxis<T,R>::alpha->m);
	axisMZ(dims[1], offset[0], mz_res, DataAxis<T,R>::mz);
	DataAxis<T,R>::rt.push_back(rt);
}


template<template<class Operator> class MathOp, typename T,typename R>
SMData2D<MathOp,T,R>::SMData2D(vector<R> &_rt, vector<double> &_mz,
		vector<T> &vec) : rt(_rt),mz(_mz)
{
	hsize_t row = hsize_t(_rt.size());
	hsize_t col = hsize_t(_mz.size());
	DataAxis<T,R>::alpha = new VecMat<T>(row,col,vec);
}

template<template<class Operator> class MathOp, typename T,typename R>
SMData2D<MathOp,T,R>::SMData2D(hsize_t dims[], int offset[],
		double mz_res, double rt_res, vector<T> &vec)
{
	DataAxis<T,R>::alpha = new VecMat<T>(dims[0],dims[1],vec);
	apply(dims[0],dims[1],DataAxis<T,R>::alpha->m);
	axisRT(dims[0], offset[1], rt_res, DataAxis<T,R>::rt);
	axisMZ(dims[1], offset[0], mz_res, DataAxis<T,R>::mz);
}

#endif /* SMPEAK_SMDATA_HPP_ */

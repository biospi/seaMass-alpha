//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Bristol, UK
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

#include "../io/VecMat.hpp"


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
	SMData1D(li dims[], int offset[], double mz_res, R rt,
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
	SMData2D(li dims[], int offset[], double mz_res, double rt_res,
			vector<T> &vec);
	~SMData2D(){delete this->alpha;};
};


#include "SMData.tpp"

#endif /* SMPEAK_SMDATA_HPP_ */

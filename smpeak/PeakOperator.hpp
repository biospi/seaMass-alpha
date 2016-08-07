//
// $Id$
//
//
// Author: Ranjeet Bhamber <ranjeet <a.t> liverpool.ac.uk>
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

#ifndef SMPEAK_PEAKOPERATOR_HPP_
#define SMPEAK_PEAKOPERATOR_HPP_

#include "../io/peakcore.hpp"


template<typename T, typename R>
class MathOp
{
protected:
	void calMidPoint(lli rtIdx, lli mzIdx, T** alpha,
					const vector<double> &mza, double &mz1, double &a1);
	double calT(const double y0, const double y1, const double y2);
	double calX(double t, double x0, double x1, double x2);
	vector<T> cal3rdMidPoint(lli rtIdx, lli mzIdx, T **P);
	T calPeakCount(vector<T> &ry, double t);
	void calPeakWidth(lli rtIdx,lli mzIdx, T** alpha, const vector<double> d2mz,
					double &mzlhs, double &mzrhs);
	void calPeakMZ(DataAxis<T,R> const *bs,DataAxis<T,R> const *dbs,DataAxis<T,R> const *d2bs,
					lli i, lli j, double &mzPeak, T &countMax,
					double &mzlhs, double &mzrhs, double &t0, lli &falsePeak);
	~MathOp(){};
};


template
<
	template <class Peak> class pPeak,
	template <class Data1, class Data2> class pData,
	typename T = float,
	typename R = double

>
class Centroid1D : public MathOp<T,R>
{
public:
	void calculate(pPeak<T> *peak, pData<R,T> *data, T threshold);
protected:
	~Centroid1D(){};
};


template
<
	template <class Peak> class pPeak,
	template <class Data1, class Data2> class pData,
	typename T = float,
	typename R = double
>
class Centroid2D : public MathOp<T,R>
{
public:
	void calculate(pPeak<T> *peak, pData<R,T> *data, T threshold);
protected:
	~Centroid2D(){};
};


template
<
	template <class Peak> class pPeak,
	template <class Data1, class Data2> class pData,
	typename T = float,
	typename R = double
>
class ExtractPeak : public MathOp<T,R>
{
public:
	void calculate(pPeak<T> *peak, pData<R,T> *data, T threshold);
protected:
	~ExtractPeak(){};
private:
	class DataPoint
	{
	public:
		double mz; // x
		double rt; // y
		T pk;      // z
		double t0;
	};
	void mulVecMat(const VecMat<T> &cs, const VecMat<T> &bp, VecMat<T> &pcs, lli i, lli j);
	int calInnerPeakMZ(DataAxis<T,R> const *bs,DataAxis<T,R> const *dhbs,
					lli i, lli j, double &mzPeak, T &countMax,
					double &t0, lli &falsePeak);
	void bezierP(vector<DataPoint> const &spt, vector<DataPoint> &p);
	void bezierT(vector<DataPoint> &p);
	DataPoint bezierO3(vector<DataPoint> &p);
	void calMidPointRT(lli rtIdx, lli mzIdx, T** alpha,
					const vector<double> &rta, double &rt1, double &a1);
	vector<T> cal3rdMidPointRT(lli rtIdx, lli mzIdx, T **P);
	void calPeakWidthRT(lli rtIdx,lli mzIdx, T** alpha, const vector<double> drt,
					double &rtlhs, double &rtrhs);
	void calPeakRT(DataAxis<T,R> const *bs,DataAxis<T,R> const *dbs,DataAxis<T,R> const *d2bs,
					lli i, lli j, double &mzPeak, T &countMax,
					double &mzlhs, double &mzrhs, double &t0, lli &falsePeak);
	void calPeakLenRT(lli rtIdx,lli mzIdx, T** alpha, const vector<double> d2rt,
					double &rtlhs, double &rtrhs);
};

#include "PeakOperator.tpp"

#endif /* SMPEAK_PEAKOPERATOR_HPP_ */

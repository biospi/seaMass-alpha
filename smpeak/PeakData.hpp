//
// $Id$
//
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

#ifndef SMPEAK_PEAKDATA_HPP_
#define SMPEAK_PEAKDATA_HPP_

#include "../io/iomath.hpp"
#include "../io/NetCDFile.hpp"

template<typename T = float>
struct Peak
{
public:
	Peak(double _mz, double _rt,T _pkcnt,
		pair<double,double> _mzW,
		pair<double,double> _rtW,
		double _t,
		lli _mz_idx, lli _rt_idx);
	double mz; // MZ value of Peak
	double rt; // RT value of Peak
	T pkcnt; // Peak count value
	pair<double,double> mzW; // MZ value of Peak Width [LHS,RHS]
	pair<double,double> rtW; // RT value of Peak Width [LHS,RHS]
	double t;  // Value of t at second derivative
	lli mz_idx; // MZ index value
	lli rt_idx; // RT index value, also scan number
};

template<typename T = float>
class PeakData
{
protected:
	typedef vector<Peak<T> > pdata;
	pdata peakData;
	lli falsePeak;
	lli falseWidth;
public:
	PeakData():falsePeak(0),falseWidth(0){}
	pdata* getPeakData(void);
	void addPeak(double _mz, double _rt, T _pkcnt,
		pair<double,double> _mzW, pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx);
	void addPeakArray(pdata* peakArray);
	void updateFalseData(lli fPeak, lli fWidth);
	size_t numOfPeaks(void);
	vector<double> getMZ(void);
	vector<double> getRT(void);
	vector<T> getPKcount(void);
	vector<double> getMZwidth(void);
	vector<double> getRTwidth(void);
	vector<double> getT(void);
	vector<lli> getMZIdx(void);
	vector<lli> getRTIdx(void);
	void getPeakMat(VecMat<double> &mz, VecMat<T> &pk, size_t maxRT, vector<size_t> &vecSize);
	void getPeakMatT(VecMat<double> &mz, VecMat<T> &pk, size_t maxRT, vector<size_t> &vecSize);
	void dumpPeakData(string filename, nc_type data_type_id=NC_FLOAT);
	//void dumpPeakData(string filename, const H5::DataType &data_type_id=H5::PredType::NATIVE_FLOAT);
	lli getFalsePeaks(void);
	lli getFalseWidths(void);
	void clear(void);
};

#include"PeakData.tpp"

#endif /* SMPEAK_PEAKDATA_HPP_ */

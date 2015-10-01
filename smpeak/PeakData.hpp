#ifndef SMPEAK_PEAKDATA_HPP_
#define SMPEAK_PEAKDATA_HPP_

#include"peakcore.hpp"

template<typename T = float>
struct Peak
{
	Peak(double _mz, double _rt, T _pkcnt,
		pair<double,double> _mzW, pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx);
	double mz; // MZ value of Peak
	double rt; // RT value of Peak
	T pkcnt; // Peak count value
	pair<double,double> mzW; // MZ value of Peak Width [LHS,RHS]
	pair<double,double> rtW; // RT value of Peak Width [LHS,RHS]
	double t;  // Value of t at second derivative
	lli mz_idx; // MZ index value
	lli rt_idx; // RT index value, also scan number
};


template<typename T>
Peak<T>::Peak(double _mz, double _rt, T _pkcnt,
		pair<double,double> _mzW, pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx):
		mz(_mz), rt(_rt), pkcnt(_pkcnt), mzW(_mzW), rtW(_rtW), t(_t),
		mz_idx(_mz_idx), rt_idx(_rt_idx){}

#endif /* SMPEAK_PEAKDATA_HPP_ */

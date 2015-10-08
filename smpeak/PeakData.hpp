#ifndef SMPEAK_PEAKDATA_HPP_
#define SMPEAK_PEAKDATA_HPP_

#include "peakcore.hpp"

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
public:
	pdata* getPeakData(void);
	void addPeak(double _mz, double _rt, T _pkcnt,
		pair<double,double> _mzW, pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx);
	size_t numOfPeaks(void);
	vector<double> getMZ(void);
	vector<double> getRT(void);
	vector<T> getPKcount(void);
	vector<double> getMZwidth(void);
	vector<double> getRTwidth(void);
	vector<double> getT(void);
	vector<lli> getMZIdx(void);
	vector<lli> getRTIdx(void);
};

#include"PeakData.tpp"

#endif /* SMPEAK_PEAKDATA_HPP_ */

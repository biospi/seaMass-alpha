#ifndef SMPEAK_PEAKDATA_HPP_
#define SMPEAK_PEAKDATA_HPP_

#include"peakcore.hpp"

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

template<typename T>
Peak<T>::Peak(double _mz, double _rt,T _pkcnt,
		pair<double,double> _mzW,pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx) :
		mz(_mz), rt(_rt), pkcnt(_pkcnt), mzW(_mzW),rtW(_rtW),t(_t),
		mz_idx(_mz_idx),rt_idx(_rt_idx){}


template<typename T>
vector<Peak<T> >* PeakData<T>::getPeakData(void)
{
	return *peakData;
}

template<typename T>
void PeakData<T>::addPeak(double _mz, double _rt,
		T _pkcnt,
		pair<double,double> _mzW,
		pair<double,double> _rtW,
		double _t,
		lli _mz_idx, lli _rt_idx)
{
	peakData.push_back(Peak<T>::Peak(_mz,_rt, _pkcnt, _mzW, _rtW, _t, _mz_idx, _rt_idx));
}

template<typename T>
size_t PeakData<T>::numOfPeaks(void)
{
	return peakData.size();
}

template<typename T>
vector<double> PeakData<T>::getMZ(void)
{
	vector<double> val(peakData.size(),0.0);
	for(lli i = 0; i < peakData.size(); ++i)
		val[i]=peakData[i].mz;
	return val;
}

template<typename T>
vector<double> PeakData<T>::getRT(void)
{
	vector<double> val(peakData.size(),0.0);
	for(lli i = 0; i < peakData.size(); ++i)
		val[i]=peakData[i].rt;
	return val;
}

template<typename T>
vector<T> PeakData<T>::getPKcount(void)
{
	vector<T> val(peakData.size(),0.0);
	for(lli i = 0; i < peakData.size(); ++i)
		val[i]=peakData[i].pkcnt;
	return val;
}

template<typename T>
vector<double> PeakData<T>::getMZwidth(void)
{
	pair<double,double> width(0.0,0.0);
	vector<double> val(peakData.size()*2,0.0);
	lli idx=0;
	for(lli i = 0; i < peakData.size(); ++i)
	{
		width=peakData[i].mzW;
		val[idx]=width.first;
		++idx;
		val[idx]=width.second;
		++idx;
	}
	return val;
}

template<typename T>
vector<double> PeakData<T>::getRTwidth(void)
{
	pair<double,double> width(0.0,0.0);
	vector<double> val(peakData.size()*2,0.0);
	lli idx=0;
	for(lli i = 0; i < peakData.size(); ++i)
	{
		width=peakData[i].rtW;
		val[idx]=width.first;
		++idx;
		val[idx]=width.second;
		++idx;
	}
	return val;
}

template<typename T>
vector<double> PeakData<T>::getT(void)
{
	vector<double> val(peakData.size(),0.0);
	for(lli i = 0; i < peakData.size(); ++i)
		val[i]=peakData[i].t;
	return val;
}

template<typename T>
vector<lli> PeakData<T>::getMZIdx(void)
{
	vector<lli> val(peakData.size(),0.0);
	for(lli i = 0; i < peakData.size(); ++i)
		val[i]=peakData[i].mz_idx;
	return val;
}

template<typename T>
vector<lli> PeakData<T>::getRTIdx(void)
{
	vector<lli> val(peakData.size(),0.0);
	for(lli i = 0; i < peakData.size(); ++i)
		val[i]=peakData[i].rt_idx;
	return val;
}

/*
class PeakData
{
private:
	typedef vector<Peak<> > pdata;
public:
	pdata* getPeakData(void);
	void addPeak(double _mz, double _rt, float _pkcnt,
		pair<double,double> _mzW, pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx);
private:
	pdata peakData;

};


template<typename T>
Peak<T>::Peak(double _mz, double _rt, T _pkcnt,
		pair<double,double> _mzW, pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx):
		mz(_mz), rt(_rt), pkcnt(_pkcnt), mzW(_mzW), rtW(_rtW), t(_t),
		mz_idx(_mz_idx), rt_idx(_rt_idx){}


vector<Peak<> >* PeakData::getPeakData(void)
{
	return *peakData;
}

void PeakData::addPeak(double _mz, double _rt, float _pkcnt,
		pair<double,double> _mzW, pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx)
{
	peakData.push_back(Peak<>(_mz,_rt, _pkcnt, _mzW, _rtW, _t, _mz_idx, _rt_idx));
}
*/

#endif /* SMPEAK_PEAKDATA_HPP_ */

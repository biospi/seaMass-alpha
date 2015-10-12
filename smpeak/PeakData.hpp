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
	void getPeakMat(VecMat<double> &mz, VecMat<T> &pk, size_t maxRT, vector<size_t> &vecSize);
	void getPeakMatT(VecMat<double> &mz, VecMat<T> &pk, size_t maxRT, vector<size_t> &vecSize);
};

#include"PeakData.tpp"


template<typename T>
void PeakData<T>::getPeakMat(VecMat<double> &mz, VecMat<T> &pk, size_t maxRT, vector<size_t> &vecSize)
{
	size_t N = this->peakData.size();
	vector<vector<double> > mzbuff;
	vector<vector<T> > pkbuff;
	//size_t maxRT=0;
	size_t maxMZ=0;

	//for(size_t i = 0; i < N; ++i)
		//if(peakData[i].rt_idx > maxRT) maxRT = this->peakData[i].rt_idx;

	mzbuff.resize(maxRT);
	pkbuff.resize(maxRT);

	for(size_t i = 0; i < N; ++i)
	{
		mzbuff[peakData[i].rt_idx].push_back(peakData[i].mz);
		pkbuff[peakData[i].rt_idx].push_back(peakData[i].pkcnt);
	}

	for(size_t i = 0;i < maxRT ; ++i)
		if(mzbuff[i].size() > maxMZ) maxMZ = mzbuff[i].size();

	mz.set(hsize_t(maxRT),hsize_t(maxMZ));
	pk.set(hsize_t(maxRT),hsize_t(maxMZ));

	for(size_t i = 0; i < maxRT; ++i)
	{
		vecSize.push_back(mzbuff[i].size());
		for(size_t j = 0; j < mzbuff[i].size(); ++j)
		{
			mz.m[i][j]=mzbuff[i][j];
			pk.m[i][j]=pkbuff[i][j];
		}
	}
}

template<typename T>
void PeakData<T>::getPeakMatT(VecMat<double> &mz, VecMat<T> &pk, size_t maxRT, vector<size_t> &vecSize)
{
	size_t N = peakData.size();
	vector<vector<double> > mzbuff;
	vector<vector<T> > pkbuff;
	//size_t maxRT=0;
	size_t maxMZ=0;

	//for(size_t i = 0; i < N; ++i)
		//if(peakData[i].rt_idx > maxRT) maxRT = peakData[i].rt_idx;

	mzbuff.resize(maxRT);
	pkbuff.resize(maxRT);

	for(size_t i = 0; i < N; ++i)
	{
		mzbuff[peakData[i].rt_idx].push_back(peakData[i].mz);
		pkbuff[peakData[i].rt_idx].push_back(peakData[i].pkcnt);
	}

	for(size_t i = 0; i < maxRT; ++i)
		if(mzbuff[i].size() > maxMZ) maxMZ = mzbuff[i].size();

	mz.set(hsize_t(maxMZ),hsize_t(maxRT));
	pk.set(hsize_t(maxMZ),hsize_t(maxRT));

	for(size_t i = 0; i < maxRT; ++i)
	{
		vecSize.push_back(mzbuff[i].size());
		for(size_t j = 0; j < mzbuff[i].size(); ++j)
		{
			mz.m[j][i]=mzbuff[i][j];
			pk.m[j][i]=pkbuff[i][j];
		}
	}
}

#endif /* SMPEAK_PEAKDATA_HPP_ */

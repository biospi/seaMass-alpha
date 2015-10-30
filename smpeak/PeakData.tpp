#ifndef SMPEAK_PEAKDATA_TPP_
#define SMPEAK_PEAKDATA_TPP_


template<typename T>
Peak<T>::Peak(double _mz, double _rt,T _pkcnt,
		pair<double,double> _mzW,pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx) :
		mz(_mz), rt(_rt), pkcnt(_pkcnt), mzW(_mzW),rtW(_rtW),t(_t),
		mz_idx(_mz_idx),rt_idx(_rt_idx){}


template<typename T>
vector<Peak<T> >* PeakData<T>::getPeakData(void)
{
	return &peakData;
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
void PeakData<T>::addPeakArray(pdata* peakArray)
{
	peakData.insert(peakData.end(),peakArray->begin(),peakArray->end());
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

#endif /* SMPEAK_PEAKDATA_TPP_ */

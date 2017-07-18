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

#ifndef SMPEAK_PEAKDATA_TPP_
#define SMPEAK_PEAKDATA_TPP_

#include "PeakData.hpp"

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
	peakData.push_back(typename Peak<T>::Peak(_mz,_rt, _pkcnt, _mzW, _rtW, _t, _mz_idx, _rt_idx));
}


template<typename T>
void PeakData<T>::addPeakArray(pdata* peakArray)
{
	peakData.insert(peakData.end(),peakArray->begin(),peakArray->end());
}

template<typename T>
void PeakData<T>::updateFalseData(lli fPeaks, lli fWidths)
{
	this->falsePeak+=fPeaks;
	this->falseWidth+=fWidths;
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

	mz.set(uli(maxRT),uli(maxMZ));
	pk.set(uli(maxRT),uli(maxMZ));

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

	mz.set(uli(maxMZ),uli(maxRT));
	pk.set(uli(maxMZ),uli(maxRT));

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

template<typename T>
void PeakData<T>::dumpPeakData(string filename, nc_type data_type_id)
{
	// Write data to SMP file.
	string outFileName=filename.substr(0,filename.size()-4)+".smp";
	FileNetcdf smpDataFile(string(outFileName),NC_NETCDF4);

	cout<<"\nSaving Peak Data to File: "<<outFileName<<endl;

	vector<T> tbuff;
	vector<double> dbuff;
	vector<lli> ibuff;

	dbuff=this->getMZ();
	smpDataFile.write_VecNC("Peak_mz",dbuff,NC_DOUBLE);
	dbuff=this->getMZwidth();
	smpDataFile.write_VecNC("Peak_mz_width",dbuff,NC_DOUBLE);
	dbuff=this->getRT();
	smpDataFile.write_VecNC("Peak_rt",dbuff,NC_DOUBLE);
	dbuff=this->getRTwidth();
	smpDataFile.write_VecNC("Peak_rt_width",dbuff,NC_DOUBLE);
	tbuff=this->getPKcount();
	smpDataFile.write_VecNC("Peak_Count",tbuff,data_type_id);
	ibuff=this->getMZIdx();
	smpDataFile.write_VecNC("Peak_mz_idx",ibuff,NC_INT64);
	ibuff=this->getRTIdx();
	smpDataFile.write_VecNC("Peak_rt_idx",ibuff,NC_INT64);
}

template<typename T>
void PeakData<T>::writePeakWidth(string filename, nc_type data_type_id)
{
    // Write data to seaMass Width file SMW.
    //string outFileName = filename.substr(0, filename.size() - 4) + ".smw";
    string outFileName = filename + ".smw";
    FileNetcdf smwDataFile(string(outFileName), NC_NETCDF4);

    cout << "\nSaving Peak Data to File: " << outFileName << endl;

    vector<double> peak;
    vector<double> absWidth;
    vector<T> count;

    peak = this->getMZ();
    smwDataFile.write_VecNC("Peak_mz", peak, data_type_id);
    peak = this->getMZwidth();
    smwDataFile.write_VecNC("Peak_mz_width_locations", peak, data_type_id);
    count = this->getPKcount();
    smwDataFile.write_VecNC("Peak_Count",count,NC_FLOAT);

    for(li i = 0; i < peak.size()-1; i=i+2)
    {
        double width=peak[i+1]-peak[i];
        absWidth.push_back(width);
    };

    smwDataFile.write_VecNC("Peak_mz_width", absWidth, data_type_id);
    peak = this->getRT();
    smwDataFile.write_VecNC("Peak_rt",peak,data_type_id);
    absWidth.clear();
    for(li i = 0; i < peak.size(); ++i)
    {
        absWidth.push_back(peak[i]);
        absWidth.push_back(peak[i]);
    }
    smwDataFile.write_VecNC("Peak_rt_width_locations",absWidth,data_type_id);
}

template<typename T>
lli PeakData<T>::getFalsePeaks(void)
{
	return falsePeak;
}

template<typename T>
lli PeakData<T>::getFalseWidths(void)
{
	return falseWidth;
}

template<typename T>
void PeakData<T>::clear(void)
{
	peakData.clear();
	vector<Peak<T> >(peakData).swap(peakData);
}

#endif /* SMPEAK_PEAKDATA_TPP_ */

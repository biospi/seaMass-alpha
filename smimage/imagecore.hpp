#ifndef _IMAGECORE_HPP
#define _IMAGECORE_HPP

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <cmath>

using namespace std;

//typedef unsigned long long int lli;
typedef long long lli;
typedef unsigned int ui;

struct MassSpecData
{
	MassSpecData(void);
	MassSpecData(int mz, int rt);

	void calMZi(void);
	void calRange(void);

	vector<double> rt; // StartTime from smj file
	vector<double> mz; // SpectrumMZ from rebinned smo file
	vector<float> sc; // SpectrumCount from rebinned smo file
	vector<lli> sci; // SpectrumCountIndex from rebinned smo file
	vector<lli> mzi; // SpectrumMZ Index calculated from sci
	vector<double> precursorMZ; // PrecursorMZ
	vector<lli> rti; // StartTime index for valid MZ.
	vector<lli> rtp; // StartTime index for Precursor MZ.
	lli N;

	pair<double,double> mzRange; // Viewable limits
	pair<double,double> rtRange; // Viewable limits
};

struct MassSpecImage
{
	MassSpecImage(void);
	MassSpecImage(pair<ui, ui> _xyview);

	vector<double> rt; // Start Time / RT axis
	vector<double> mz; // Spectrum MZ / MZ axis
	vector<float> sc; // Spectrum Count = I(MZ,RT)
	double drt;
	double dmz;

	pair<double,double> mzRange; // Viewable limits
	pair<double,double> rtRange; // Viewable limits
	pair<ui,ui> xypxl; // Physical Pixel dimensions of smimage
};

void scanMZ(MassSpecData &raw, MassSpecImage &imgbox, lli const idx,
			lli const row, double const scaleRT);

template<typename T1, typename T2>
pair<T1,T2> numStrBraket(string numket);

template<typename T>
T genAxis(vector<T> &x, T xmin, T xmax);

template<typename T>
T genAxisN(vector<T> &x, T xmin, T xmax);


template<typename T1, typename T2>
pair<T1,T2> numStrBraket(string numket)
{
	pair<T1,T2> range;
	size_t blhs=numket.find("[");
	size_t brhs=numket.find("]");
	size_t coma=numket.find(",");

	if(blhs == coma-1)
		range.first=T1(-1);
	else
		istringstream(numket.substr(1,coma-1))>>range.first;

	if(brhs == coma+1)
		range.second=T2(-1);
	else
		istringstream(numket.substr(coma+1,brhs-coma-1))>>range.second;

	return range;
}


template<typename T>
T genAxis(vector<T> &x, T xmin, T xmax)
{
	size_t N=x.size();
	T dx=fabs(xmax-xmin)/T(N-1);

	for(size_t i = 0; i < N; ++i)
	{
		x[i]=T(i)*dx+xmin;
	}
	return dx;
}


template<typename T>
T genAxisN(vector<T> &x, T xmin, T xmax)
{
	size_t N=x.size();
	T dx=fabs(xmax-xmin)/T(N);

	for(size_t i = 0; i < N; ++i)
	{
		x[i]=T(i)*dx+xmin;
	}
	return dx;
}


#endif //_IMAGECORE_HPP

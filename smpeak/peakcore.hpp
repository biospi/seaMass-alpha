#ifndef SMPEAK_PEAKCORE_HPP_
#define SMPEAK_PEAKCORE_HPP_

#include <iostream>
#include <vector>
#include <cmath>
#include <H5Cpp.h>

using namespace std;

typedef long long int lli;

template<typename T = float>
struct VecMat
{
	VecMat(void);
	VecMat(hsize_t _r, hsize_t _c, vector<T> &_vec);
	VecMat(hsize_t _r, hsize_t _c);
	vector<T> v; // Vector of Matrix data.
	T** m; // Data Matrix
	void set(hsize_t _r, hsize_t _c, vector<T> &_vec);
	void set(hsize_t _r, hsize_t _c);
	void getDims(hsize_t dims[]);
	void clear(void);
private:
	vector<T*> matIdx;
	hsize_t row;
	hsize_t col;
};

template<typename T>
void findVecString(vector<char> &vecStr,vector<T> &vec,
		const string subStr = "<spectrum index",
		const string endSubStr = "</spectrum>");

template<typename T>
vector<size_t> findSize(VecMat<T> data);

template<typename T>
vector<size_t> findSizeT(VecMat<T> data);

template<typename T>
void repackPeakDataT(VecMat<T> &peak, VecMat<T> &paw, vector<int> msl,
				vector<size_t> psize, vector<size_t> rsize);


template<typename T>
vector<size_t> findSize(VecMat<T> data)
{
	hsize_t dims[2];
	data.getDims(dims);
	vector<size_t> dsize;

	for(size_t idx=0; idx < dims[0]; ++idx)
	{
		size_t len=0;
		while(len < dims[1] && data.m[idx][len] > 0) ++len;
		dsize.push_back(len);
	}
	return dsize;
}

template<typename T>
vector<size_t> findSizeT(VecMat<T> data)
{
	hsize_t dims[2];
	data.getDims(dims);
	vector<size_t> dsize;

	for(size_t idx=0; idx < dims[1]; ++idx)
	{
		size_t len=0;
		while(len < dims[0] && data.m[len][idx] > 0) ++len;
		dsize.push_back(len);
	}
	return dsize;
}


template<typename T>
void repackPeakDataT(VecMat<T> &peak, VecMat<T> &raw, vector<int> msl,
				vector<size_t> psize, vector<size_t> rsize)
{
	vector<vector<T> > dataBuff;

	hsize_t pdim[2];
	hsize_t rdim[2];
	hsize_t r;
	hsize_t c;

	peak.getDims(pdim);
	raw.getDims(rdim);

	(pdim[0] > rdim[0]) ?  r = pdim[0]: r = rdim[0];
	(pdim[1] > rdim[1]) ?  c = pdim[1]: c = rdim[1];

	dataBuff.resize(r);

	for(hsize_t i = 0; i < r; ++i)
		dataBuff[i].resize(c,0);

	for(size_t idx=0, rdx=0; idx < msl.size(); ++idx)
	{
		if((msl[idx] == 1) && (rdx < psize.size() ))
		{
			for(int i = 0; i < psize[rdx]; ++i)
			{
				dataBuff[i][idx]=peak.m[i][rdx];
			}
			++rdx;
		}
		else if(msl[idx] == 2)
		{
			for(int i = 0; i < rsize[idx]; ++i)
			{
				dataBuff[i][idx]=raw.m[i][idx];
			}
		}
		else
		{
			cout << "Error in Repack!!! Reput Old data in missing ms-level!!!"<<endl;
			for(int i = 0; i < rsize[idx]; ++i)
			{
				dataBuff[i][idx]=raw.m[i][idx];
			}
		}
	}

	raw.clear();
	raw.set(r,c);

	for(int i=0; i < r; ++i)
		for(int j=0; j < c; ++j)
			raw.m[i][j]=dataBuff[i][j];
}

#include"peakcore.tpp"

#endif /* SMPEAK_PEAKCORE_HPP_ */

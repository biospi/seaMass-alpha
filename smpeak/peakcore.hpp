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
void repackPeakData(VecMat<T> &peak, VecMat<T> &paw, vector<int> msl,
				vector<size_t> psize, vector<size_t> rsize);

template<typename T>
void repackPeakDataT(VecMat<T> &peak, VecMat<T> &paw, vector<int> msl,
				vector<size_t> psize, vector<size_t> rsize);

#include"peakcore.tpp"

#endif /* SMPEAK_PEAKCORE_HPP_ */

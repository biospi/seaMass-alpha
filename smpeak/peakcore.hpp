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
	vector<T> v; // Vector of Matrix data.
	T** m; // Data Matrix
	void getDims(hsize_t dims[]);
	void set(hsize_t _r, hsize_t _c, vector<T> &_vec);
private:
	vector<T*> matIdx;
	hsize_t row;
	hsize_t col;
};

#include"peakcore.tpp"

#endif /* SMPEAK_PEAKCORE_HPP_ */

#ifndef SMPEAK_PEAKDATA_HPP_
#define SMPEAK_PEAKDATA_HPP_

#include"peakcore.hpp"
#include<H5Cpp.h>

template<typename T = float>
struct VecMat
{
	VecMat(void);
	VecMat(hsize_t _r, hsize_t _c, vector<T> &_vec);
	vector<T> v; // Vector of Matrix data.
	T** m; // Data Matrix
	void getDims(hsize_t &dims);
	void set(hsize_t _r, hsize_t _c, vector<T> &_vec);
private:
	vector<T*> matIdx;
	hsize_t row;
	hsize_t col;
};

#include"peakdata.tpp"

struct SMData
{
	SMData(vector<double> &_rt, vector<double> &_mz, vector<float> &vec);
	vector<double> rt;
	vector<double> mz;
	VecMat<float> alpha;
};

struct Peak
{
	double mz; // MZ value of Peak
	double rt; // RT value of Peak
	float pkcnt; // Peak count value
	pair<double,double> mzW; // MZ value of Peak Width [LHS,RHS]
	pair<double,double> rtW; // RT value of Peak Width [LHS,RHS]
	double t;  // Value of t at second derivative
	lli mz_idx; // MZ index value
	lli rt_idx; // RT index value, also scan number
};

#endif /* SMPEAK_PEAKDATA_HPP_ */

#ifndef SMPEAK_SMDATA_HPP_
#define SMPEAK_SMDATA_HPP_

#include"peakcore.hpp"

template<template<class Operator> class MathOp, typename T=float>
struct SMData : public MathOp<T>
{
	SMData(vector<double> &_rt, vector<double> &_mz, vector<T> &vec);
	SMData(hsize_t dims[], int offset[], double mz_res, double rt_res,
			vector<T> &vec);
	~SMData(){delete alpha;};
	vector<double> rt;
	vector<double> mz;
	VecMat<T>* alpha;
};

template<template<class Operator> class MathOp, typename T>
SMData<MathOp,T>::SMData(vector<double> &_rt, vector<double> &_mz,
		vector<T> &vec) : rt(_rt),mz(_mz)
{
	hsize_t row = hsize_t(_rt.size());
	hsize_t col = hsize_t(_mz.size());
	alpha = new VecMat<T>(row,col,vec);
}

template<template<class Operator> class MathOp, typename T>
SMData<MathOp,T>::SMData(hsize_t dims[], int offset[],
		double mz_res, double rt_res, vector<T> &vec)
{
	alpha = new VecMat<T>(dims[0],dims[1],vec);
	apply(dims[0],dims[1],this->alpha->m);
	axisRT(dims[0], offset[1], rt_res, this->rt);
	axisMZ(dims[1], offset[0], mz_res, this->mz);
}

template<typename T = float>
struct Peak
{
	Peak(double _mz, double _rt, T _pkcnt,
		pair<double,double> _mzW, pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx);
	double mz; // MZ value of Peak
	double rt; // RT value of Peak
	T pkcnt; // Peak count value
	pair<double,double> mzW; // MZ value of Peak Width [LHS,RHS]
	pair<double,double> rtW; // RT value of Peak Width [LHS,RHS]
	double t;  // Value of t at second derivative
	lli mz_idx; // MZ index value
	lli rt_idx; // RT index value, also scan number
};


template<typename T>
Peak<T>::Peak(double _mz, double _rt, T _pkcnt,
		pair<double,double> _mzW, pair<double,double> _rtW,
		double _t, lli _mz_idx, lli _rt_idx):
		mz(_mz), rt(_rt), pkcnt(_pkcnt), mzW(_mzW), rtW(_rtW), t(_t),
		mz_idx(_mz_idx), rt_idx(_rt_idx){}

#endif /* SMPEAK_SMDATA_HPP_ */

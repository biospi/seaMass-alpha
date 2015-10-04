#ifndef SMPEAK_SMDATA_HPP_
#define SMPEAK_SMDATA_HPP_

#include"peakcore.hpp"

template<typename T=float>
struct DataAxis
{
	vector<double> rt;
	vector<double> mz;
	VecMat<T>* alpha;
protected:
	~DataAxis(){};
};

template
<
	template<class Operator> class MathOp,
	typename T=float
>
struct SMData : public DataAxis<T>, public MathOp<T>
{
	SMData(vector<double> &_rt, vector<double> &_mz, vector<T> &vec);
	SMData(hsize_t dims[], int offset[], double mz_res, double rt_res,
			vector<T> &vec);
	~SMData(){delete this->alpha;};
};


template<template<class Operator> class MathOp, typename T>
SMData<MathOp,T>::SMData(vector<double> &_rt, vector<double> &_mz,
		vector<T> &vec) : rt(_rt),mz(_mz)
{
	hsize_t row = hsize_t(_rt.size());
	hsize_t col = hsize_t(_mz.size());
	DataAxis<T>::alpha = new VecMat<T>(row,col,vec);
}

template<template<class Operator> class MathOp, typename T>
SMData<MathOp,T>::SMData(hsize_t dims[], int offset[],
		double mz_res, double rt_res, vector<T> &vec)
{
	DataAxis<T>::alpha = new VecMat<T>(dims[0],dims[1],vec);
	apply(dims[0],dims[1],DataAxis<T>::alpha->m);
	axisRT(dims[0], offset[1], rt_res, DataAxis<T>::rt);
	axisMZ(dims[1], offset[0], mz_res, DataAxis<T>::mz);
}

#endif /* SMPEAK_SMDATA_HPP_ */

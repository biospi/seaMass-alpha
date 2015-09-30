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
	//alpha.set(hsize_t(_rt.size()),hsize_t(_mz.size()),vec);
	alpha = new VecMat<T>(row,col,vec);
	apply(row,col,alpha->m);
}

template<template<class Operator> class MathOp, typename T>
SMData<MathOp,T>::SMData(hsize_t dims[], int offset[],
		double mz_res, double rt_res, vector<T> &vec)
{
	alpha = new VecMat<T>(dims[0],dims[1],vec);
	apply(dims[0],dims[1],alpha->m);
	axisRT(dims[0], offset[1], rt_res, this->rt);
	axisMZ(dims[1], offset[0], mz_res, this->mz);
}

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

#endif /* SMPEAK_SMDATA_HPP_ */

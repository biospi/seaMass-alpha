#ifndef SMPEAK_SMDATA_HPP_
#define SMPEAK_SMDATA_HPP_

#include"peakcore.hpp"

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

#endif /* SMPEAK_SMDATA_HPP_ */

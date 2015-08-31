#ifndef SMPEAK_PEAKCORE_HPP_
#define SMPEAK_PEAKCORE_HPP_

#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

typedef long long int lli;

struct PeakData
{
	void add_peak(double _mz, double _rt, float _pVal, lli mzidx, lli rtidx);

	vector<double> mz; // MZ value of Peak
	vector<double> rt; // RT value of Peak
	vector<float> pVal; // Peak Value
	vector<lli> mz_idx; // MZ index value
	vector<lli> rt_idx; // RT index value, also scan number
};

vector<vector<float> > nabla(float **alpha, lli row, lli col);

vector<vector<float> > nabla2(float **alpha, lli row, lli col);

void calMZalpha(vector<double>& _mza, int offset);

void calRTalpha(vector<double>& _rt, int offset);

void calMidPoint(lli rtIdx, lli mzIdx, vector<vector<float> > &alpha,
		vector<double> &mza, double &mz1, double &a1);

double calT(const double y0, const double y1, const double y2);

double calX(double t, double x0, double x1, double x2);

#endif /* SMPEAK_PEAKCORE_HPP_ */

#ifndef SMPEAK_PEAKCORE_HPP_
#define SMPEAK_PEAKCORE_HPP_

#include<iostream>
#include<vector>
#include<cmath>
#include<H5Cpp.h>

using namespace std;

typedef long long int lli;

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

#include"peakcore.tpp"

//////////////////////////////////////////////////////////////////////////
struct PeakData
{
	void add_peak(double _mz, double _rt, float _pVal, double _t, lli mzidx, lli rtidx,
			double mzlhs, double mzrhs);

	vector<double> mz; // MZ value of Peak
	vector<double> mzW; // MZ value of Peak Width
	vector<double> rt; // RT value of Peak
	vector<double> rtW; // RT value of Peak Width
	vector<double> t;  // Value of t at second derivative
	vector<float> count; // Peak count value
	vector<lli> mz_idx; // MZ index value
	vector<lli> rt_idx; // RT index value, also scan number
};

//////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------
vector<vector<float> > nabla(float **alpha, lli row, lli col);
vector<vector<float> > nabla2(float **alpha, lli row, lli col);

void calDMZalpha(vector<double>& _mza, int _offset, double mz_rez);
void calD2MZalpha(vector<double>& _mza, int _offset, double mz_rez);

void calRTalpha(vector<double>& _rt, int _offset, double rt_rez);
//-----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void calMidPoint(lli rtIdx, lli mzIdx, vector<vector<float> > &alpha,
		vector<double> &mza, double &mz1, double &a1);

double calT(const double y0, const double y1, const double y2);

double calX(double t, double x0, double x1, double x2);

vector<float> cal3rdMidPoint(lli rtIdx, lli mzIdx, float **P);

float calPeakCount(vector<float> &ry, double t);

void calPeakWidth(lli rtIdx,lli mzIdx, vector<vector<float> > &alpha, vector<double> d2mz,
		double &mzlhs, double &mzrhs);
//////////////////////////////////////////////////////////////////////////

#endif /* SMPEAK_PEAKCORE_HPP_ */

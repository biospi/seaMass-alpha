#include"peakcore.hpp"

void PeakData::add_peak(double _mz, double _rt, float _pVal, lli mzidx, lli rtidx)
{
	mz.push_back(_mz);
	rt.push_back(_rt);
	pVal.push_back(_pVal);
	mz_idx.push_back(mzidx);
	rt_idx.push_back(rtidx);
}


vector<vector<float> >nabla(float **alpha, lli row, lli col)
{
	vector<vector<float> > dAlpha;

	dAlpha.resize(row);
	for(lli i=0; i < dAlpha.size(); ++i)
		dAlpha[i].resize(col,0.0);

	#pragma omp parallel for
	for(lli i=0; i < row; ++i)
		for(lli j=1; j < col; ++j)
			dAlpha[i][j]=alpha[i][j]-alpha[i][j-1];

	return dAlpha;
}

vector<vector<float> > nabla2(float **alpha, lli row, lli col)
{
	vector<vector<float> > d2Alpha;

	d2Alpha.resize(row);
	for(lli i=0; i < d2Alpha.size(); ++i)
		d2Alpha[i].resize(col,0.0);

	#pragma omp parallel for
	for(lli i=0; i < row; ++i)
		for(lli j=2; j < col; ++j)
			d2Alpha[i][j]=alpha[i][j]-2.0*alpha[i][j-1]+alpha[i][j-2];

	return d2Alpha;
}


void calMZalpha(vector<double>& _mza, int offset) {
	for (lli i = 0; i < _mza.size(); ++i) {
		double mz_rez = 1.0;
		double ppbmz = 1.0033548378 / (pow(2.0, mz_rez) * 60.0);
		_mza[i] = double((offset - 1 + i)) * ppbmz;
	}
}

void calRTalpha(vector<double>& _rt, int offset) {
	for (lli i = 0; i < _rt.size(); ++i) {
		double rt_rez = 4.0;
		double ppbrt = 1.0 / (pow(2.0, rt_rez));
		_rt[i] = double(offset + i) * ppbrt;
	}
}

void calMidPoint(lli rtIdx, lli mzIdx, vector<vector<float> > &alpha,
		vector<double> &mza, double &mz1, double &a1)
{
	mz1=0.5*(mza[mzIdx]+mza[mzIdx+1]);
	a1=double(0.5*(alpha[rtIdx][mzIdx]+alpha[rtIdx][mzIdx+1]));
}

double calT(const double y0, const double y1, const double y2)
{
	double a=y1-y0;
	double b=y0-2.0*y1+y2;
	double g=a/b;
	double r=sqrt(g*g-y0/b);
	double tp=-g+r;
	double tn=-g-r;

	if(tp >= 0 && tp <= 1) return tp;
	else if(tn >= 0 && tn <= 1)	return tn;
	else
	{
		cout<<"ERROR calculating T is out of range !!!"<<endl;
		return -1.0;
	}
}

double calX(double t, double x0, double x1, double x2)
{
	double c=(1.0-t);
	return c*c*x0+2.0*c*t*x1+t*t*x2;
}

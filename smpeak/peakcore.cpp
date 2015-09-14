#include"peakcore.hpp"

void PeakData::add_peak(double _mz, double _rt, float _pVal, double _t, lli mzidx, lli rtidx,
		double mzlhs, double mzrhs)
{
	mz.push_back(_mz);
	rt.push_back(_rt);
	count.push_back(_pVal);
	t.push_back(_t);
	mz_idx.push_back(mzidx);
	rt_idx.push_back(rtidx);
	mzW.push_back(mzlhs);
	mzW.push_back(mzrhs);
	rtW.push_back(_rt);
	rtW.push_back(_rt);
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


void calDMZalpha(vector<double>& _mza, int _offset, double mz_rez)
{
	double offset = double(_offset) - 1.5; // Factor due to Derivative along MZ
	double ppbmz = 1.0033548378 / (pow(2.0, mz_rez) * 60.0);
	for (lli i = 0; i < _mza.size(); ++i) {
		_mza[i] = (offset + i) * ppbmz;
	}
}

void calD2MZalpha(vector<double>& _mza, int _offset, double mz_rez)
{
	double offset = double(_offset) - 2; // Factor due to double Derivative along MZ
	double ppbmz = 1.0033548378 / (pow(2.0, mz_rez) * 60.0);
	for (lli i = 0; i < _mza.size(); ++i) {
		_mza[i] = (offset + i) * ppbmz;
	}
}

void calRTalpha(vector<double>& _rt, int _offset, double rt_rez)
{
	double offset = double(_offset) - 1; // Factor due to non-Derivative
	double ppbrt = 1.0 / (pow(2.0, rt_rez));
	for (lli i = 0; i < _rt.size(); ++i) {
		_rt[i] = (offset + i) * ppbrt;
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
		//cout<<"Insignificant False Peak Detected - Peak Ignored"<<endl;
		return -1.0;
	}
}

double calX(double t, double x0, double x1, double x2)
{
	double c=(1.0-t);
	return c*c*x0+2.0*c*t*x1+t*t*x2;
}

vector<float> cal3rdMidPoint(lli row, lli col, float **P)
{
	vector<float> ry(4,0.0);
	float p0=P[row][col-1];
	float p1=P[row][col];
	float p2=P[row][col+1];
	float p3=P[row][col+2];
	float u=1.0/3.0;
	float v=2.0/3.0;
	float w=1.0/2.0;
	float up=(1-u);
	float vp=(1-v);
	float wp=(1-w);

	ry[0]=wp*(vp*p0+v*p1)+w*(up*p1+u*p2);
	ry[1]=up*p1+u*p2;
	ry[2]=vp*p1+v*p2;
	ry[3]=wp*(vp*p1+v*p2)+w*(up*p2+u*p3);

	return ry;
}


float calPeakCount(vector<float> &ry, double t)
{
	double tp=1-t;
	float p0=tp*ry[0]+t*ry[1];
	float p1=tp*ry[1]+t*ry[2];
	float p2=tp*ry[2]+t*ry[3];
	float q0=tp*p0+t*p1;
	float q1=tp*p1+t*p2;
	return tp*q0+t*q1;
}

void calPeakWidth(lli rtIdx,lli mzIdx, vector<vector<float> > &alpha, vector<double> d2mz,
		double &mzlhs, double &mzrhs)
{
	lli lhsIdx=mzIdx;
	lli rhsIdx=mzIdx;

	do{--lhsIdx;} while( alpha[rtIdx][lhsIdx] < 0);
	do{++rhsIdx;} while( alpha[rtIdx][rhsIdx] < 0);

	double x1,x2,y1,y2,m;

	x1=d2mz[lhsIdx];
	y1=alpha[rtIdx][lhsIdx];
	x2=d2mz[lhsIdx+1];
	y2=alpha[rtIdx][lhsIdx+1];
	m=(y2-y1)/(x2-x1);
	mzlhs=(-y1/m)+x1;

	x1=d2mz[rhsIdx];
	y1=alpha[rtIdx][rhsIdx];
	x2=d2mz[rhsIdx-1];
	y2=alpha[rtIdx][rhsIdx-1];
	m=(y2-y1)/(x2-x1);
	mzrhs=(-y1/m)+x1;
}

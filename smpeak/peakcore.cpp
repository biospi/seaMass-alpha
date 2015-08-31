#include"peakcore.hpp"

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

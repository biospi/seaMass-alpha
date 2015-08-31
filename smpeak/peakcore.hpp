#ifndef SMPEAK_PEAKCORE_HPP_
#define SMPEAK_PEAKCORE_HPP_

#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

typedef long long int lli;


vector<vector<float> > nabla(float **alpha, lli row, lli col);

vector<vector<float> > nabla2(float **alpha, lli row, lli col);

void calMZalpha(vector<double>& _mza, int offset);

void calRTalpha(vector<double>& _rt, int offset);


#endif /* SMPEAK_PEAKCORE_HPP_ */

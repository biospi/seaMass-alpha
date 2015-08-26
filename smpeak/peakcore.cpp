#include "peakcore.hpp"
#include <iostream>

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

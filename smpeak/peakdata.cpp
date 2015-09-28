#include "peakdata.hpp"

SMData::SMData(vector<double> &_rt,	vector<double> &_mz,
		vector<float> &vec) : rt(_rt),mz(_mz)
{
	alpha.set(hsize_t(_rt.size()),hsize_t(_mz.size()),vec);
}

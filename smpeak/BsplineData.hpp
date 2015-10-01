#ifndef SMPEAK_BSPLINEDATA_HPP_
#define SMPEAK_BSPLINEDATA_HPP_

#include"peakcore.hpp"
#include"SMData.hpp"

class BsplineData
{
public:
	BsplineData(DataAxis<> &bs, DataAxis<> &dbs, DataAxis<> &d2bs);
	BsplineData(DataAxis<> &bs, DataAxis<> &dhbs, DataAxis<> &d2hbs,
				DataAxis<> &dvbs, DataAxis<> &d2vbs);
	vector<DataAxis<>* > get(void);
private:
	vector<DataAxis<>* > bspObjP;
};

#endif /* SMPEAK_BSPLINEDATA_HPP_ */

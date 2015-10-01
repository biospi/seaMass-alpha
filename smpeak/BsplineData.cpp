#include "BsplineData.hpp"

BsplineData::BsplineData(DataAxis<> &bs, DataAxis<> &dbs, DataAxis<> &d2bs)
{
	bspObjP.push_back(&bs);
	bspObjP.push_back(&dbs);
	bspObjP.push_back(&d2bs);
}


BsplineData::BsplineData(DataAxis<> &bs, DataAxis<> &dhbs, DataAxis<> &d2hbs,
		DataAxis<> &dvbs, DataAxis<> &d2vbs)
{
	bspObjP.push_back(&bs);
	bspObjP.push_back(&dhbs);
	bspObjP.push_back(&d2hbs);
	bspObjP.push_back(&dvbs);
	bspObjP.push_back(&d2vbs);
}

vector<DataAxis<>* > BsplineData::get(void)
 {
	 return bspObjP;
 }

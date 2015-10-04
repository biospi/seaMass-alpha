#ifndef SMPEAK_BSPLINEDATA_HPP_
#define SMPEAK_BSPLINEDATA_HPP_

#include"peakcore.hpp"
#include"SMData.hpp"

/*
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
*/

template<typename T = float>
class BsplineData
{
public:
	BsplineData(DataAxis<T> &bs, DataAxis<T> &dbs, DataAxis<T> &d2bs);
	BsplineData(DataAxis<T> &bs, DataAxis<T> &dhbs, DataAxis<T> &d2hbs,
				DataAxis<T> &dvbs, DataAxis<T> &d2vbs);
	vector<DataAxis<T>* > get(void);
private:
	vector<DataAxis<T>* > bspObjP;
};


template<typename T>
BsplineData<T>::BsplineData(DataAxis<T> &bs, DataAxis<T> &dbs, DataAxis<T> &d2bs)
{
	bspObjP.push_back(&bs);
	bspObjP.push_back(&dbs);
	bspObjP.push_back(&d2bs);
}

template<typename T>
BsplineData<T>::BsplineData(DataAxis<T> &bs, DataAxis<T> &dhbs, DataAxis<T> &d2hbs,
		DataAxis<T> &dvbs, DataAxis<T> &d2vbs)
{
	bspObjP.push_back(&bs);
	bspObjP.push_back(&dhbs);
	bspObjP.push_back(&d2hbs);
	bspObjP.push_back(&dvbs);
	bspObjP.push_back(&d2vbs);
}

template<typename T>
vector<DataAxis<T>* > BsplineData<T>::get(void)
{
	return bspObjP;
}


#endif /* SMPEAK_BSPLINEDATA_HPP_ */

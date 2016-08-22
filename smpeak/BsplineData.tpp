//
// $Id$
//
//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Bristol, UK
//
// This file is part of seaMass.
//
// seaMass is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// seaMass is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with seaMass.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef SMPEAK_BSPLINEDATA_TPP_
#define SMPEAK_BSPLINEDATA_TPP_


template<typename T>
BasisPatch<T>::BasisPatch()
{
	//T x[] = {-1.0, -2.0/3.0, -1.0/3.0, 0.0, 1.0/3.0};
	T x[] = {-1.0, -2.0/3.0, -1.0/3.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0};
	int t[] = {-4, -3, -2, -1, 0};
	int k=4;
	b.set(5,7);
	for(int i = 0; i < 5; ++i)
		for(int j = 0; j < 7; ++j)
			b.m[i][j] =  BSpline::m(x[j],k,t[i]);
}

template<typename T>
void BasisPatch<T>::set(vector<T> &x, vector<int> &t)
{
	int k =4;
	uli r = t.size();
	uli c = x.size();
	b.set(r,c);
	for(uli i = 0; i < r; ++i)
		for(uli j = 0; j < c; ++j)
			b.m[i][j] = m(x[j],k,t[i]);
}


template<typename R,typename T>
BsplineData<R,T>::BsplineData(DataAxis<T,R> &bs, DataAxis<T,R> &dbs, DataAxis<T,R> &d2bs)
{
	bspObjP.push_back(&bs);
	bspObjP.push_back(&dbs);
	bspObjP.push_back(&d2bs);
}


template<typename R,typename T>
BsplineBasisData<R,T>::BsplineBasisData(DataAxis<T,R> &bs,
			DataAxis<T,R> &dhbs, DataAxis<T,R> &d2hbs,
			DataAxis<T,R> &dvbs, DataAxis<T,R> &d2vbs,
			BasisPatch<T> &bp)
{
	bspObjP.push_back(&bs);
	bspObjP.push_back(&dhbs);
	bspObjP.push_back(&d2hbs);
	bspObjP.push_back(&dvbs);
	bspObjP.push_back(&d2vbs);
	bPat = &bp;
}

template<typename R,typename T>
BsplineBasisData<R,T>::BsplineBasisData(DataAxis<T,R> &bs,
			DataAxis<T,R> &dhbs, DataAxis<T,R> &d2hbs,
			DataAxis<T,R> &dvbs, BasisPatch<T> &bp)
{
	bspObjP.push_back(&bs);
	bspObjP.push_back(&dhbs);
	bspObjP.push_back(&d2hbs);
	bspObjP.push_back(&dvbs);
	bPat = &bp;
}

template<typename R,typename T>
vector<DataAxis<T,R>* > BsplineData<R,T>::get(void)
{
	return bspObjP;
}

template<typename R, typename T>
void BsplineBasisData<R,T>::get(vector<DataAxis<T,R>* > &bsDat, BasisPatch<T> *&bp)
{
	bsDat=bspObjP;
	bp=bPat;
}

template<typename R, typename T>
void BsplineBasisData<R,T>::dumpData(string filename, nc_type data_type_id)
{

	// Write data to SMD (debug) file.
	string outFileName=filename.substr(0,filename.size()-4)+".smd";
	NetCDFile smpDataFile(outFileName,NC_NETCDF4);

	DataAxis<T,R> const *bs=bspObjP[0];
	DataAxis<T,R> const *dhbs=bspObjP[1];
	DataAxis<T,R> const *dh2bs=bspObjP[2];
	DataAxis<T,R> const *dvbs=bspObjP[3];
	DataAxis<T,R> const *dv2bs=bspObjP[4];

	cout<<"\nSaving Peak Debugging Data to File:"<<endl;
	smpDataFile.write_VecNC("csOrig",bs->alpha->v,data_type_id);
	smpDataFile.write_VecNC("dhcs",dhbs->alpha->v,data_type_id);
	smpDataFile.write_VecNC("dh2cs",dh2bs->alpha->v,data_type_id);
	smpDataFile.write_VecNC("dvcs",dvbs->alpha->v,data_type_id);
	smpDataFile.write_VecNC("dv2cs",dv2bs->alpha->v,data_type_id);
}

#endif /* SMPEAK_BSPLINEDATA_TPP_ */

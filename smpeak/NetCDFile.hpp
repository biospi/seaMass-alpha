#ifndef SMPEAK_NETCDFILE_HPP_
#define SMPEAK_NETCDFILE_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <netcdf.h>

#include "peakcore.hpp"

#define ERRCODE 2
#define ERR(e) {cout<<"Error: "<<nc_strerror(e)<<endl; exit(ERRCODE);}

using namespace std;

class NetCDFile
{
public:
	NetCDFile(void) : fileStatus(false), ncid(0),retval(0){};
	NetCDFile(const string _fileName, int omode = NC_NOWRITE);
	void open(const string _fileName, int omode = NC_NOWRITE);
	void close(void);
	~NetCDFile();

	template<typename T>
	void read_VecNC(const string dataSet, vector<T> &vecData);
	template<typename T>
	void read_MatNC(const string dataSet, VecMat<T> &vm);
	template<typename T>
	void read_MatNCT(const string dataSet, VecMat<T> &vm);
private:
	string fileName;
	int ncid;
	int retval;
	bool fileStatus;
	void err(int e);
};


template<typename T>
void NetCDFile::read_VecNC(const string dataSet, vector<T> &vecData)
{

	int varid;
	int ndim;
	vector<int> dimid;
	vector<size_t> dimSize;
	size_t N=1;
	nc_type typId;

	if((retval = nc_inq_varid (ncid, dataSet.c_str(), &varid) ))
		ERR(retval);

	if((retval = nc_inq_vartype(ncid, varid, &typId) ))
		ERR(retval);

	if((retval = nc_inq_varndims(ncid,varid,&ndim) ))
		ERR(retval);

	dimid.resize(ndim);
	dimSize.resize(ndim);

	if((retval = nc_inq_vardimid(ncid, varid, &dimid[0]) ))
		ERR(retval);

	for(int i = 0; i < ndim; ++i)
	{
		if ((retval = nc_inq_dimlen(ncid, dimid[i], &dimSize[i]) ))
			ERR(retval);
	}

	for(int i = 0; i < ndim; ++i)
		N*=dimSize[i];

	vecData.resize(N);

	if (( retval = nc_get_var(ncid, varid, &vecData[0]) ))
		ERR(retval)
}


template<typename T>
void NetCDFile::read_MatNC(const string dataSet, VecMat<T> &vm)
{
	int varid;
	int ndim;
	vector<int> dimid;
	vector<size_t> dimSize;
	size_t N=1;
	nc_type typId;

	if((retval = nc_inq_varid (ncid, dataSet.c_str(), &varid) ))
		ERR(retval);

	if((retval = nc_inq_vartype(ncid, varid, &typId) ))
		ERR(retval);

	if((retval = nc_inq_varndims(ncid,varid,&ndim) ))
		ERR(retval);

	dimid.resize(ndim);
	dimSize.resize(ndim);

	if((retval = nc_inq_vardimid(ncid, varid, &dimid[0]) ))
		ERR(retval);

	for(int i = 0; i < ndim; ++i)
	{
		if ((retval = nc_inq_dimlen(ncid, dimid[i], &dimSize[i]) ))
			ERR(retval);
	}

	for(int i = 0; i < ndim; ++i)
		N*=dimSize[i];

	vm.set(dimSize[0],dimSize[1]);

	if (( retval = nc_get_var(ncid, varid, &vm.v[0]) ))
		ERR(retval)
}


template<typename T>
void NetCDFile::read_MatNCT(const string dataSet, VecMat<T> &vm)
{
	int varid;
	int ndim;
	vector<int> dimid;
	vector<size_t> dimSize;
	size_t N=1;
	nc_type typId;
	vector<T> vecbuff;

	if((retval = nc_inq_varid (ncid, dataSet.c_str(), &varid) ))
		ERR(retval);

	if((retval = nc_inq_vartype(ncid, varid, &typId) ))
		ERR(retval);

	if((retval = nc_inq_varndims(ncid,varid,&ndim) ))
		ERR(retval);

	dimid.resize(ndim);
	dimSize.resize(ndim);

	if((retval = nc_inq_vardimid(ncid, varid, &dimid[0]) ))
		ERR(retval);

	for(int i = 0; i < ndim; ++i)
	{
		if ((retval = nc_inq_dimlen(ncid, dimid[i], &dimSize[i]) ))
			ERR(retval);
	}

	for(int i = 0; i < ndim; ++i)
		N*=dimSize[i];

	vecbuff.resize(N);

	vm.set(dimSize[1],dimSize[0]);

	if (( retval = nc_get_var(ncid, varid, &vecbuff[0]) ))
		ERR(retval)

	for(int i=0; i < dimSize[0]; ++i)
	{
		for(int j=0; j < dimSize[1]; ++j)
		{
			vm.m[j][i]=vecbuff[i*dimSize[1]+j];
		}
	}
}

#endif /* SMPEAK_NETCDFILE_HPP_ */

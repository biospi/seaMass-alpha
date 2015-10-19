#ifndef SMPEAK_NETCDFILE_HPP_
#define SMPEAK_NETCDFILE_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <netcdf.h>
#include <fstream>

#include "peakcore.hpp"

#define ERRCODE 2
#define ERR(e) {cout<<"Error: "<<nc_strerror(e)<<endl; exit(ERRCODE);}

using namespace std;


struct InfoGrpVar
{
	InfoGrpVar(int _g, int _v, string _strg, string _strv):
		grpid(_g), varid(_v), grpName(_strg), varName(_strv){};
	int grpid;
	int varid;
	string grpName;
	string varName;
};

class NetCDFile
{
public:
	NetCDFile(void) : fileStatus(false), ncid(0),retval(0){};
	NetCDFile(const string _fileName, int omode = NC_NOWRITE);
	void open(const string _fileName, int omode = NC_NOWRITE);
	void close(void);
	~NetCDFile();

	template<typename T>
	void read_VecNC(const string dataSet, vector<T> &vecData, int grpid = 0);
	template<typename T>
	void read_VecNC(const string dataSet, T *vecData, int grpid = 0);

	template<typename T>
	void read_MatNC(const string dataSet, VecMat<T> &vm, int grpid = 0);
	template<typename T>
	void read_MatNCT(const string dataSet, VecMat<T> &vm, int grpid = 0);
	template<typename T>
	void read_MatNC(const string dataSet, vector<vector<T> > &vm, int grpid = 0);
	template<typename T>
	void read_MatNCT(const string dataSet, vector<vector<T> > &vm, int grpid = 0);

	int search_Group(const string dataSet, int grpid = 0);
	template<typename T>
	T search_Group(size_t level, int grpid = 0);

	//template<typename T>
	//void read_AttNC(const string dataSet, VecMat<T> &vm);

	template<typename T>
	void write_VecNC(const string dataSet, vector<T> &vec, nc_type xtype,
			size_t chunks = 4096, int deflate_level = 5, int shuffle = NC_SHUFFLE);
	template<typename T>
	void write_MatNC(const string dataSet, VecMat<T> &vm, nc_type xtype,
			size_t chunk = 64, int deflate_level = 5, int shuffle = NC_SHUFFLE);

	vector<InfoGrpVar> get_Info(void) {return dataSetList;};
	//template<typename T>
	//void write_AttNC();
private:
	string fileName;
	int ncid;
	int retval;
	bool fileStatus;
	vector<InfoGrpVar> dataSetList;
	void err(int e);
};

void mzMLdump(const string fileName, string data);

template<typename T>
void NetCDFile::read_VecNC(const string dataSet, vector<T> &vecData, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int varid;
	int ndim;
	vector<int> dimid;
	vector<size_t> dimSize;
	size_t N=1;
	nc_type typId;

	if((retval = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
		ERR(retval);

	if((retval = nc_inq_vartype(grpid, varid, &typId) ))
		ERR(retval);

	if((retval = nc_inq_varndims(grpid,varid,&ndim) ))
		ERR(retval);

	dimid.resize(ndim);
	dimSize.resize(ndim);

	if((retval = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
		ERR(retval);

	for(int i = 0; i < ndim; ++i)
	{
		if ((retval = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
			ERR(retval);
	}

	for(int i = 0; i < ndim; ++i)
		N*=dimSize[i];

	vecData.resize(N);

	if (( retval = nc_get_var(grpid, varid, &vecData[0]) ))
		ERR(retval)
}


template<typename T>
void NetCDFile::read_VecNC(const string dataSet, T *vecData, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int varid;
	int ndim;
	vector<int> dimid;
	vector<size_t> dimSize;
	size_t N=1;
	nc_type typId;

	if((retval = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
		ERR(retval);

	if((retval = nc_inq_vartype(grpid, varid, &typId) ))
		ERR(retval);

	if((retval = nc_inq_varndims(grpid,varid,&ndim) ))
		ERR(retval);

	dimid.resize(ndim);
	dimSize.resize(ndim);

	if((retval = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
		ERR(retval);

	for(int i = 0; i < ndim; ++i)
	{
		if ((retval = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
			ERR(retval);
	}

	for(int i = 0; i < ndim; ++i)
		N*=dimSize[i];

	vecData = new T[N];

	if (( retval = nc_get_var(grpid, varid, &vecData[0]) ))
		ERR(retval)
}


template<typename T>
void NetCDFile::read_MatNC(const string dataSet, VecMat<T> &vm, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int varid;
	int ndim;
	vector<int> dimid;
	vector<size_t> dimSize;
	nc_type typId;

	if((retval = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
		ERR(retval);

	if((retval = nc_inq_vartype(grpid, varid, &typId) ))
		ERR(retval);

	if((retval = nc_inq_varndims(grpid,varid,&ndim) ))
		ERR(retval);

	dimid.resize(ndim);
	dimSize.resize(ndim);

	if((retval = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
		ERR(retval);

	for(int i = 0; i < ndim; ++i)
	{
		if ((retval = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
			ERR(retval);
	}

	vm.set(dimSize[0],dimSize[1]);

	if (( retval = nc_get_var(grpid, varid, &vm.v[0]) ))
		ERR(retval)
}


template<typename T>
void NetCDFile::read_MatNCT(const string dataSet, VecMat<T> &vm, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int varid;
	int ndim;
	vector<int> dimid;
	vector<size_t> dimSize;
	size_t N=1;
	nc_type typId;
	vector<T> vecbuff;

	if((retval = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
		ERR(retval);

	if((retval = nc_inq_vartype(grpid, varid, &typId) ))
		ERR(retval);

	if((retval = nc_inq_varndims(grpid,varid,&ndim) ))
		ERR(retval);

	dimid.resize(ndim);
	dimSize.resize(ndim);

	if((retval = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
		ERR(retval);

	for(int i = 0; i < ndim; ++i)
	{
		if ((retval = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
			ERR(retval);
	}

	for(int i = 0; i < ndim; ++i)
		N*=dimSize[i];

	vecbuff.resize(N);

	vm.set(dimSize[1],dimSize[0]);

	if (( retval = nc_get_var(grpid, varid, &vecbuff[0]) ))
		ERR(retval)

	for(int i=0; i < dimSize[0]; ++i)
	{
		for(int j=0; j < dimSize[1]; ++j)
		{
			vm.m[j][i]=vecbuff[i*dimSize[1]+j];
		}
	}
}

template<typename T>
void NetCDFile::read_MatNC(const string dataSet, vector<vector<T> > &vm, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int varid;
	int ndim;
	vector<int> dimid;
	vector<size_t> dimSize;
	size_t N=1;
	nc_type typId;
	vector<T> vecbuff;

	if((retval = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
		ERR(retval);

	if((retval = nc_inq_vartype(grpid, varid, &typId) ))
		ERR(retval);

	if((retval = nc_inq_varndims(grpid,varid,&ndim) ))
		ERR(retval);

	dimid.resize(ndim);
	dimSize.resize(ndim);

	if((retval = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
		ERR(retval);

	for(int i = 0; i < ndim; ++i)
	{
		if ((retval = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
			ERR(retval);
	}

	for(int i = 0; i < ndim; ++i)
		N*=dimSize[i];

	vecbuff.resize(N);

	vm.resize(dimSize[0]);
	for(size_t i = 0; i < dimSize[0] ; ++i)
		vm[i].resize(dimSize[1]);

	if (( retval = nc_get_var(grpid, varid, &vecbuff[0]) ))
		ERR(retval)

	for(int i=0; i < dimSize[0]; ++i)
	{
		for(int j=0; j < dimSize[1]; ++j)
		{
			vm[i][j]=vecbuff[i*dimSize[1]+j];
		}
	}
}

template<typename T>
void NetCDFile::read_MatNCT(const string dataSet, vector<vector<T> > &vm, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int varid;
	int ndim;
	vector<int> dimid;
	vector<size_t> dimSize;
	size_t N=1;
	nc_type typId;
	vector<T> vecbuff;

	if((retval = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
		ERR(retval);

	if((retval = nc_inq_vartype(grpid, varid, &typId) ))
		ERR(retval);

	if((retval = nc_inq_varndims(grpid,varid,&ndim) ))
		ERR(retval);

	dimid.resize(ndim);
	dimSize.resize(ndim);

	if((retval = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
		ERR(retval);

	for(int i = 0; i < ndim; ++i)
	{
		if ((retval = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
			ERR(retval);
	}

	for(int i = 0; i < ndim; ++i)
		N*=dimSize[i];

	vecbuff.resize(N);

	vm.resize(dimSize[1]);
	for(size_t i = 0; i < dimSize[1] ; ++i)
		vm[i].resize(dimSize[0]);

	if (( retval = nc_get_var(grpid, varid, &vecbuff[0]) ))
		ERR(retval)

	for(int i=0; i < dimSize[0]; ++i)
	{
		for(int j=0; j < dimSize[1]; ++j)
		{
			vm[j][i]=vecbuff[i*dimSize[1]+j];
		}
	}
}


template<typename T>
T NetCDFile::search_Group(size_t level, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int nGrps;
	vector<int> ngrpids;
	string strGrp;
	size_t strGrpL=0;

	int nVars;
	int varid;

	T val;
	size_t ds=0;

	for(hsize_t i = 0; i <= level; ++i)
	{
		size_t cidx=strGrpL;
		// How many Groups
		if(( retval = nc_inq_grps(grpid, &nGrps, NULL) ))
			err(retval);

		if(nGrps > 0 )
		{
			ngrpids.resize(nGrps);
			if(( retval = nc_inq_grps(grpid, NULL, &ngrpids[0]) ))
				err(retval);
		}

		if(( retval = nc_inq_grpname_len(grpid, &strGrpL)))
			err(retval);
		strGrp.resize(strGrpL);
		if(( retval = nc_inq_grpname_full(grpid,NULL,&strGrp[0])))
			err(retval);

		ds=strGrpL-cidx;

		grpid = ngrpids[0];
	}

	istringstream(strGrp.substr(strGrpL-ds+1,ds-1))>>val;
	return val;
}


template<typename T>
void NetCDFile::write_VecNC(const string dataSet, vector<T> &vec, nc_type xtype,
			size_t chunks, int deflate_level, int shuffle)
{
	int varid;
	int dimid;
	int ndim = 1;
	int deflate = 0;
	size_t N = vec.size();
	string dimName = "dim_" + dataSet;

	/* Set chunking, shuffle, and deflate. */
	shuffle = NC_SHUFFLE;
	if(deflate_level > 0 && deflate_level < 10)
		deflate = 1;

	/* Define the dimensions. */
	if ((retval = nc_def_dim(ncid, dimName.c_str(), N, &dimid)))
	   ERR(retval);

	/* Define the variable. */
	if ((retval = nc_def_var(ncid, dataSet.c_str(), xtype, ndim,
	                         &dimid, &varid)))
	   ERR(retval);
	if ((retval = nc_def_var_chunking(ncid, varid, NC_CHUNKED, &chunks)))
	   ERR(retval);
	if ((retval = nc_def_var_deflate(ncid, varid, shuffle, deflate,
	                                 deflate_level)))
	   ERR(retval);
	/* No need to explicitly end define mode for netCDF-4 files. Write
	 * the pretend data to the file. */
	if ((retval = nc_put_var(ncid, varid, &vec[0])))
	   ERR(retval);
}

template<typename T>
void NetCDFile::write_MatNC(const string dataSet, VecMat<T> &vm, nc_type xtype,
			size_t chunk, int deflate_level, int shuffle)
{
	int varid;
	int dimid[2];
	size_t chunks[2];
	int ndim = 2;
	int deflate = 0;
	size_t N[2];
	hsize_t buffN[2];
	string dimName1 = "row_" + dataSet;
	string dimName2 = "col_" + dataSet;

	/* Set chunking, shuffle, and deflate. */
	shuffle = NC_SHUFFLE;
	if(deflate_level > 0 && deflate_level < 10)
		deflate = 1;

	vm.getDims(buffN);

	N[0]=size_t(buffN[0]);
	N[1]=size_t(buffN[1]);

	/* Define the dimensions. */
	if ((retval = nc_def_dim(ncid,dimName1.c_str(),N[0],&dimid[0])))
	   ERR(retval);
	if ((retval = nc_def_dim(ncid,dimName2.c_str(),N[1],&dimid[1])))
	   ERR(retval);

	chunks[0] = chunk;
	chunks[1] = chunk;

	/* Define the variable. */
	if ((retval = nc_def_var(ncid,dataSet.c_str(),xtype,ndim,&dimid[0],&varid)))
	   ERR(retval);

	if ((retval = nc_def_var_chunking(ncid,varid,NC_CHUNKED,&chunks[0])))
	   ERR(retval);

	if ((retval = nc_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level)))
	   ERR(retval);

	/* No need to explicitly end define mode for netCDF-4 files. Write
	 * the pretend data to the file. */
	if ((retval = nc_put_var(ncid, varid, &vm.v[0])))
	   ERR(retval);
}




#endif /* SMPEAK_NETCDFILE_HPP_ */

/*
 * NetCDFile.tpp
 *
 *  Created on: 20 Oct 2015
 *      Author: ranjeet
 */

#ifndef SMPEAK_NETCDFILE_TPP_
#define SMPEAK_NETCDFILE_TPP_


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
void NetCDFile::read_AttNC(const string attName, int varid, vector<T> &attVal, int grpid)
{
	if(grpid == 0) grpid = ncid;

	nc_type xtype;
	size_t len;

	if(( retval = nc_inq_att(grpid, varid, attName.c_str(), &xtype, &len) ))
		err(retval);

	attVal.resize(len);

	if(( retval =  nc_get_att(grpid, varid, attName.c_str(), &attVal[0]) ))
		err(retval);
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
			size_t chunks, int deflate_level, int shuffle, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int varid;
	int dimid;
	int ndim = 1;
	int deflate = 0;
	size_t N = vec.size();
	string dimName = "dim_" + dataSet;

	// Set chunking, shuffle, and deflate.
	shuffle = NC_SHUFFLE;
	if(deflate_level > 0 && deflate_level < 10)
		deflate = 1;

	// Define the dimensions.
	if((retval = nc_def_dim(grpid, dimName.c_str(), N, &dimid)))
		ERR(retval);

	// Define the variable.
	if((retval = nc_def_var(grpid, dataSet.c_str(), xtype, ndim,
	                         &dimid, &varid)))
		ERR(retval);

	if(N < chunks) chunks = N;

	if((retval = nc_def_var_chunking(grpid, varid, NC_CHUNKED, &chunks)))
		ERR(retval);
	if((retval = nc_def_var_deflate(grpid, varid, shuffle, deflate,
	                                 deflate_level)))
		ERR(retval);

	if((retval = nc_enddef(grpid)))
		ERR(retval);

	// No need to explicitly end define mode for netCDF-4 files. Write
	// the pretend data to the file.
	if ((retval = nc_put_var(grpid, varid, &vec[0])))
		ERR(retval);
}


template<typename T>
void NetCDFile::write_MatNC(const string dataSet, VecMat<T> &vm, nc_type xtype,
			size_t chunk, int deflate_level, int shuffle, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int varid;
	int dimid[2];
	size_t chunks[2];
	int ndim = 2;
	int deflate = 0;
	size_t N[2];
	hsize_t buffN[2];
	string dimName1 = "row_" + dataSet;
	string dimName2 = "col_" + dataSet;

	// Set chunking, shuffle, and deflate.
	shuffle = NC_SHUFFLE;
	if(deflate_level > 0 && deflate_level < 10)
		deflate = 1;

	vm.getDims(buffN);

	N[0]=size_t(buffN[0]);
	N[1]=size_t(buffN[1]);

	// Define the dimensions.
	if((retval = nc_def_dim(grpid,dimName1.c_str(),N[0],&dimid[0])))
	   ERR(retval);
	if((retval = nc_def_dim(grpid,dimName2.c_str(),N[1],&dimid[1])))
	   ERR(retval);

	//if(N[0]*N[1] < chunk) chunk = N[0]*N[1];
	if(N[1] < chunk) chunk = N[1];
	chunks[0] = 1;
	chunks[1] = chunk;

	// Define the variable.
	if((retval = nc_def_var(grpid,dataSet.c_str(),xtype,ndim,&dimid[0],&varid)))
	   ERR(retval);

	if((retval = nc_def_var_chunking(grpid,varid,NC_CHUNKED,&chunks[0])))
	   ERR(retval);

	if((retval = nc_def_var_deflate(grpid,varid,shuffle,deflate,deflate_level)))
	   ERR(retval);

	if((retval = nc_enddef(grpid)))
		ERR(retval);

	// No need to explicitly end define mode for netCDF-4 files. Write
	// the pretend data to the file.
	if((retval = nc_put_var(grpid, varid, &vm.v[0])))
	   ERR(retval);
}

template<typename T>
void NetCDFile::write_AttNC(const string dataSet, const string attName,
			vector<T> &attVal, nc_type xtype, int grpid)
{
	if(grpid == 0) grpid = ncid;

	int varid;
	nc_type typatt;
	size_t len;

	if((retval = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
		ERR(retval);

	len = attVal.size();

	if((retval = nc_put_att(grpid, varid, attName.c_str(), xtype, len, &attVal[0]) ))
		ERR(retval);
}

#endif /* SMPEAK_NETCDFILE_TPP_ */

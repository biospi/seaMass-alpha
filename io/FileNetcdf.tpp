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

#ifndef SEAMASS_IO_FILENETCDF_TPP
#define SEAMASS_IO_FILENETCDF_TPP


#include "FileNetcdf.hpp"
#include <sstream>


template<typename T>
void FileNetcdf::read_VecNC(const string dataSet, vector<T> &vecData, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    vector<size_t> dimSize;
    size_t N=1;
    nc_type typId;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if((retval_ = nc_inq_vartype(grpid, varid, &typId) ))
        err(retval_);

    if((retval_ = nc_inq_varndims(grpid,varid,&ndim) ))
        err(retval_);

    dimid.resize(ndim);
    dimSize.resize(ndim);

    if((retval_ = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
        err(retval_);

    for(int i = 0; i < ndim; ++i)
    {
        if ((retval_ = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
            err(retval_);
    }

    for(int i = 0; i < ndim; ++i)
        N*=dimSize[i];

    vecData.resize(N);

    if (( retval_ = nc_get_var(grpid, varid, &vecData[0]) ))
        err(retval_);
}

template<typename T>
void FileNetcdf::read_VecNC(const string dataSet, T *vecData, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    vector<size_t> dimSize;
    size_t N=1;
    nc_type typId;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if((retval_ = nc_inq_vartype(grpid, varid, &typId) ))
        err(retval_);

    if((retval_ = nc_inq_varndims(grpid,varid,&ndim) ))
        err(retval_);

    dimid.resize(ndim);
    dimSize.resize(ndim);

    if((retval_ = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
        err(retval_);

    for(int i = 0; i < ndim; ++i)
    {
        if ((retval_ = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
            err(retval_);
    }

    for(int i = 0; i < ndim; ++i)
        N*=dimSize[i];

    vecData = new T[N];

    if (( retval_ = nc_get_var(grpid, varid, &vecData[0]) ))
        err(retval_);
}


template<typename T>
void FileNetcdf::read_MatNC(const string dataSet, VecMat<T> &vm, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    vector<size_t> dimSize;
    nc_type typId;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if((retval_ = nc_inq_vartype(grpid, varid, &typId) ))
        err(retval_);

    if((retval_ = nc_inq_varndims(grpid,varid,&ndim) ))
        err(retval_);

    dimid.resize(ndim);
    dimSize.resize(ndim);

    if((retval_ = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
        err(retval_);

    for(int i = 0; i < ndim; ++i)
    {
        if ((retval_ = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
            err(retval_);
    }

    vm.set(dimSize[0],dimSize[1]);

    if (( retval_ = nc_get_var(grpid, varid, &vm.v[0]) ))
        err(retval_);
}

template<typename T>
void FileNetcdf::read_MatNCT(const string dataSet, VecMat<T> &vm, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    vector<size_t> dimSize;
    size_t N=1;
    nc_type typId;
    vector<T> vecbuff;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if((retval_ = nc_inq_vartype(grpid, varid, &typId) ))
        err(retval_);

    if((retval_ = nc_inq_varndims(grpid,varid,&ndim) ))
        err(retval_);

    dimid.resize(ndim);
    dimSize.resize(ndim);

    if((retval_ = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
        err(retval_);

    for(int i = 0; i < ndim; ++i)
    {
        if ((retval_ = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
            err(retval_);
    }

    for(int i = 0; i < ndim; ++i)
        N*=dimSize[i];

    vecbuff.resize(N);

    vm.set(dimSize[1],dimSize[0]);

    if (( retval_ = nc_get_var(grpid, varid, &vecbuff[0]) ))
        err(retval_);

    for(int i=0; i < dimSize[0]; ++i)
    {
        for(int j=0; j < dimSize[1]; ++j)
        {
            vm.m[j][i]=vecbuff[i*dimSize[1]+j];
        }
    }
}

template<typename T>
void FileNetcdf::read_MatNC(const string dataSet, vector<vector<T> > &vm, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    vector<size_t> dimSize;
    size_t N=1;
    nc_type typId;
    vector<T> vecbuff;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if((retval_ = nc_inq_vartype(grpid, varid, &typId) ))
        err(retval_);

    if((retval_ = nc_inq_varndims(grpid,varid,&ndim) ))
        err(retval_);

    dimid.resize(ndim);
    dimSize.resize(ndim);

    if((retval_ = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
        err(retval_);

    for(int i = 0; i < ndim; ++i)
    {
        if ((retval_ = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
            err(retval_);
    }

    for(int i = 0; i < ndim; ++i)
        N*=dimSize[i];

    vecbuff.resize(N);

    vm.resize(dimSize[0]);
    for(size_t i = 0; i < dimSize[0] ; ++i)
        vm[i].resize(dimSize[1]);

    if (( retval_ = nc_get_var(grpid, varid, &vecbuff[0]) ))
        err(retval_);

    for(int i=0; i < dimSize[0]; ++i)
    {
        for(int j=0; j < dimSize[1]; ++j)
        {
            vm[i][j]=vecbuff[i*dimSize[1]+j];
        }
    }
}

template<typename T>
void FileNetcdf::read_MatNCT(const string dataSet, vector<vector<T> > &vm, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    vector<size_t> dimSize;
    size_t N=1;
    nc_type typId;
    vector<T> vecbuff;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if((retval_ = nc_inq_vartype(grpid, varid, &typId) ))
        err(retval_);

    if((retval_ = nc_inq_varndims(grpid,varid,&ndim) ))
        err(retval_);

    dimid.resize(ndim);
    dimSize.resize(ndim);

    if((retval_ = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
        err(retval_);

    for(int i = 0; i < ndim; ++i)
    {
        if ((retval_ = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
            err(retval_);
    }

    for(int i = 0; i < ndim; ++i)
        N*=dimSize[i];

    vecbuff.resize(N);

    vm.resize(dimSize[1]);
    for(size_t i = 0; i < dimSize[1] ; ++i)
        vm[i].resize(dimSize[0]);

    if (( retval_ = nc_get_var(grpid, varid, &vecbuff[0]) ))
        err(retval_);

    for(int i=0; i < dimSize[0]; ++i)
    {
        for(int j=0; j < dimSize[1]; ++j)
        {
            vm[j][i]=vecbuff[i*dimSize[1]+j];
        }
    }
}

template<typename T>
vector<size_t> FileNetcdf::read_DimNC(const string dataSet, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    vector<size_t> dimSize;
    nc_type typId;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if((retval_ = nc_inq_vartype(grpid, varid, &typId) ))
        err(retval_);

    if((retval_ = nc_inq_varndims(grpid,varid,&ndim) ))
        err(retval_);

    dimid.resize(ndim);
    dimSize.resize(ndim);

    if((retval_ = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
        err(retval_);

    for(int i = 0; i < ndim; ++i)
    {
        if ((retval_ = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
            err(retval_);
    }

    return dimSize;
}


template<typename T>
void FileNetcdf::read_AttNC(const string attName, int varid, vector<T> &attVal, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    nc_type xtype;
    size_t len;

    if(( retval_ = nc_inq_att(grpid, varid, attName.c_str(), &xtype, &len) ))
        err(retval_);

    attVal.resize(len);

    if(( retval_ =  nc_get_att(grpid, varid, attName.c_str(), &attVal[0]) ))
        err(retval_);
}


template<typename T>
void FileNetcdf::read_HypVecNC(const string dataSet, vector<T> &vm,
        size_t *rcIdx, size_t *len, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    int dimid;
    size_t dimSize;
    nc_type typId;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if((retval_ = nc_inq_vartype(grpid, varid, &typId) ))
        err(retval_);

    if((retval_ = nc_inq_varndims(grpid,varid,&ndim) ))
        err(retval_);

    if((retval_ = nc_inq_vardimid(grpid, varid, &dimid) ))
        err(retval_);

    if ((retval_ = nc_inq_dimlen(grpid, dimid, &dimSize) ))
        err(retval_);

    vm.resize(*len);

    if(typeid(vector<float>) == typeid(vm))
    {
        if (( retval_ = nc_get_vara_float(grpid, varid, rcIdx, len, reinterpret_cast<float*>(&vm[0])) ))
            err(retval_);
    }
    else if(typeid(vector<double>) == typeid(vm))
    {
        if (( retval_ = nc_get_vara_double(grpid, varid, rcIdx, len, reinterpret_cast<double*>(&vm[0])) ))
            err(retval_);
    }
    else
    {
        if (( retval_ = nc_get_vara(grpid, varid, rcIdx, len, &vm[0]) ))
            err(retval_);
    }
}

template<typename T>
void FileNetcdf::read_HypMatNC(const string dataSet, VecMat<T> &vm,
        size_t *rcIdx, size_t *len, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    vector<size_t> dimSize;
    size_t N=1;
    nc_type typId;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if((retval_ = nc_inq_vartype(grpid, varid, &typId) ))
        err(retval_);

    if((retval_ = nc_inq_varndims(grpid,varid,&ndim) ))
        err(retval_);

    dimid.resize(ndim);
    dimSize.resize(ndim);

    if((retval_ = nc_inq_vardimid(grpid, varid, &dimid[0]) ))
        err(retval_);

    for(int i = 0; i < ndim; ++i)
    {
        if ((retval_ = nc_inq_dimlen(grpid, dimid[i], &dimSize[i]) ))
            err(retval_);
    }

    for(int i = 0; i < ndim; ++i)
    {
        if(len[i] == 0) len[i]=dimSize[i]-rcIdx[i];
        N *=len[i];
    }

    vm.set(uli(len[0]), uli(len[1]));


    if(typeid(vector<float>) == typeid(vm.v))
    {
        if (( retval_ = nc_get_vara_float(grpid, varid, rcIdx, len, reinterpret_cast<float*>(&vm.v[0])) ))
            err(retval_);
    }
    else if(typeid(vector<double>) == typeid(vm.v))
    {
        if (( retval_ = nc_get_vara_double(grpid, varid, rcIdx, len, reinterpret_cast<double*>(&vm.v[0])) ))
            err(retval_);
    }
    else
    {
        if (( retval_ = nc_get_vara(grpid, varid, rcIdx, len, &vm.v[0]) ))
            err(retval_);
    }

}


template<typename T>
T FileNetcdf::search_Group(size_t level, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int nGrps;
    vector<int> ngrpids;
    string strGrp;
    size_t strGrpL=0;

    int nVars;
    int varid;

    T val;
    size_t ds=0;

    for(uli i = 0; i <= level; ++i)
    {
        size_t cidx=strGrpL;
        // How many Groups
        if(( retval_ = nc_inq_grps(grpid, &nGrps, NULL) ))
            err(retval_);

        if(nGrps > 0 )
        {
            ngrpids.resize(nGrps);
            if(( retval_ = nc_inq_grps(grpid, NULL, &ngrpids[0]) ))
                err(retval_);
        }

        if(( retval_ = nc_inq_grpname_len(grpid, &strGrpL)))
            err(retval_);
        strGrp.resize(strGrpL);
        if(( retval_ = nc_inq_grpname_full(grpid,NULL,&strGrp[0])))
            err(retval_);

        ds=strGrpL-cidx;

        grpid = ngrpids[0];
    }

    istringstream(strGrp.substr(strGrpL-ds+1,ds-1))>>val;
    return val;
}


template<typename T>
int FileNetcdf::write_VecNC(const string dataSet, const vector<T> &vec, nc_type xtype,
            int grpid, bool unlim,
            size_t chunks, int deflate_level, int shuffle)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int dimid;
    int ndim = 1;
    int deflate = 0;
    size_t N;
    //string dimName = "dim_" + dataSet;
    if(unlim == false)
    {
        N = vec.size();
        if(N < chunks) chunks = N;
    }
    else
    {
        N = NC_UNLIMITED;
    }

    // Set chunking, shuffle, and deflate.
    shuffle = NC_SHUFFLE;
    if(deflate_level > 0 && deflate_level < 10)
        deflate = 1;

    ostringstream oss;
    oss << dataSet << "_col";
    // Define the dimensions.
    if((retval_ = nc_def_dim(grpid, oss.str().c_str(), N, &dimid)))
        err(retval_);

    // Define the variable.
    if((retval_ = nc_def_var(grpid, dataSet.c_str(), xtype, ndim,
                             &dimid, &varid)))
        err(retval_);

    if((retval_ = nc_def_var_chunking(grpid, varid, NC_CHUNKED, &chunks)))
        err(retval_);
    if((retval_ = nc_def_var_deflate(grpid, varid, shuffle, deflate,
                                     deflate_level)))
        err(retval_);

    if(unlim == true)
    {
        T attVal[1] = {0};
        if((retval_ = nc_put_att(grpid, varid,"_FillValue", xtype, 1, attVal) ))
            err(retval_);
    }

    if((retval_ = nc_enddef(grpid)))
        err(retval_);

    // No need to explicitly end define mode for netCDF-4 files. Write
    // the data to the file.

    if(unlim == true)
    {
        N=vec.size();
        size_t cIdx=0;
        if((retval_ = nc_put_vara(grpid, varid,&cIdx, &N, &vec[0])))
            err(retval_);
    }
    else
    {
        if((retval_ = nc_put_var(grpid, varid, &vec[0])))
            err(retval_);
    }

    return varid;
}


template<typename T>
int FileNetcdf::write_VecNC(const string dataSet, const T *vec, size_t len, nc_type xtype,
            int grpid, bool unlim,
            size_t chunks, int deflate_level, int shuffle)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int dimid;
    int ndim = 1;
    int deflate = 0;
    size_t N;
    //string dimName = "dim_" + dataSet;
    if(unlim == false)
    {
        N = len;
        if(N < chunks) chunks = N;
    }
    else
    {
        N = NC_UNLIMITED;
    }

    // Set chunking, shuffle, and deflate.
    shuffle = NC_SHUFFLE;
    if(deflate_level > 0 && deflate_level < 10)
        deflate = 1;

    // Define the dimensions.
    if((retval_ = nc_def_dim(grpid, dataSet.c_str(), N, &dimid)))
        err(retval_);

    // Define the variable.
    if((retval_ = nc_def_var(grpid, dataSet.c_str(), xtype, ndim,
                             &dimid, &varid)))
        err(retval_);

    if((retval_ = nc_def_var_chunking(grpid, varid, NC_CHUNKED, &chunks)))
        err(retval_);
    if((retval_ = nc_def_var_deflate(grpid, varid, shuffle, deflate,
                                     deflate_level)))
        err(retval_);

    if(unlim == true)
    {
        T attVal[1] = {0};
        if((retval_ = nc_put_att(grpid, varid,"_FillValue", xtype, 1, attVal) ))
            err(retval_);
    }

    if((retval_ = nc_enddef(grpid)))
        err(retval_);

    // No need to explicitly end define mode for netCDF-4 files. Write
    // the data to the file.

    if(unlim == true)
    {
        N=len;
        size_t cIdx=0;
        if((retval_ = nc_put_vara(grpid, varid,&cIdx, &N, &vec[0])))
            err(retval_);
    }
    else
    {
        if((retval_ = nc_put_var(grpid, varid, &vec[0])))
            err(retval_);
    }

    return varid;
}


template<typename T>
int FileNetcdf::write_MatNC(const string dataSet, const VecMat<T> &vm, nc_type xtype,
            int grpid, const string rowY, const string colX,
            size_t chunk, int deflate_level, int shuffle)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int dimid[2];
    size_t chunks[2];
    int ndim = 2;
    int deflate = 0;
    size_t N[2];
    uli buffN[2];

    // Set chunking, shuffle, and deflate.
    shuffle = NC_SHUFFLE;
    if(deflate_level > 0 && deflate_level < 10)
        deflate = 1;

    vm.getDims(buffN);

    N[0]=size_t(buffN[0]);
    N[1]=size_t(buffN[1]);

    // Define the dimensions.
    if(rowY.size() != 0)
    {
        int tmpvar;
        int tmpdims[1];
        if((retval_ = nc_inq_varid(grpid,rowY.c_str(),&tmpvar) ))
            err(retval_);
        if((retval_ = nc_inq_vardimid(grpid,tmpvar,tmpdims) ))
            err(retval_);
        dimid[0]=tmpdims[0];

    }
    else
    {
        string dimName1 = dataSet+"_row";
        if((retval_ = nc_def_dim(grpid,dimName1.c_str(),N[0],&dimid[0])))
            err(retval_);
    }

    if(colX.size() != 0)
    {
        int tmpvar;
        int tmpdims[1];
        if((retval_ = nc_inq_varid(grpid,colX.c_str(),&tmpvar) ))
            err(retval_);
        if((retval_ = nc_inq_vardimid(grpid,tmpvar,tmpdims) ))
            err(retval_);
        dimid[1]=tmpdims[0];

    }
    else{
        string dimName2 = dataSet+"_col";
        if((retval_ = nc_def_dim(grpid,dimName2.c_str(),N[1],&dimid[1])))
            err(retval_);
    }

    //if(N[0]*N[1] < chunk) chunk = N[0]*N[1];
    if(N[0] < chunk) chunk = N[0];
    chunks[0] = 1;
    chunks[1] = chunk;

    // Define the variable.
    if((retval_ = nc_def_var(grpid,dataSet.c_str(),xtype,ndim,&dimid[0],&varid)))
       err(retval_);

    if((retval_ = nc_def_var_chunking(grpid,varid,NC_CHUNKED,&chunks[0])))
       err(retval_);

    if((retval_ = nc_def_var_deflate(grpid,varid,shuffle,deflate,deflate_level)))
       err(retval_);

    if((retval_ = nc_enddef(grpid)))
        err(retval_);

    // No need to explicitly end define mode for netCDF-4 files. Write
    // the data to the file.
    if((retval_ = nc_put_var(grpid, varid, &vm.v[0])))
       err(retval_);

    return varid;
}


template<typename T,typename X, typename Y>
int FileNetcdf::write_MatAxisNC(const string dataSet, const VecMat<T> &vm, nc_type ztype,
                   vector<X> colAxisX, nc_type xtype,
                   vector<Y> rowAxisY, nc_type ytype,
                   const string colX, const string rowY,
                   int grpid,
                   size_t chunk, int deflate_level,
                   int shuffle)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int dimid[2];
    int axisVarid[2];
    size_t chunks[2];
    size_t xchunk, ychunk;
    int ndim = 2;
    int vecDim=1;
    int deflate = 0;
    size_t N[2];
    uli buffN[2];

    string dimName1;
    string dimName2;

    // Set chunking, shuffle, and deflate.
    shuffle = NC_SHUFFLE;
    if(deflate_level > 0 && deflate_level < 10)
        deflate = 1;

    vm.getDims(buffN);

    N[0]=size_t(buffN[0]);
    N[1]=size_t(buffN[1]);

    // Define the dimensions.
    if(rowY.size() == 0)
    {
        dimName1 = dataSet+"_row";
    }
    else
    {
        dimName1=rowY;
    }
    if(colX.size() == 0)
    {
        dimName2 = dataSet+"_col";
    }
    else
    {
        dimName2=colX;
    }
    if((retval_ = nc_def_dim(grpid,dimName1.c_str(),N[0],&dimid[0])))
        err(retval_);
    if((retval_ = nc_def_dim(grpid,dimName2.c_str(),N[1],&dimid[1])))
        err(retval_);

    // Chunking for Axis Data
    if(colAxisX.size() < chunk) xchunk=colAxisX.size()-1;
    if(rowAxisY.size() < chunk) ychunk=rowAxisY.size()-1;

    // Chunking for Matrix
    //if(N[0]*N[1] < chunk) chunk = N[0]*N[1];
    if(N[1] < chunk) chunk = N[1];
    chunks[0] = 1;
    chunks[1] = chunk;

    // Define Axises variables.
    if((retval_ = nc_def_var(grpid,dimName1.c_str(),ytype,vecDim,&dimid[0],&axisVarid[0])))
       err(retval_);
    if((retval_ = nc_def_var_chunking(grpid,axisVarid[0],NC_CHUNKED,&ychunk)))
       err(retval_);
    if((retval_ = nc_def_var_deflate(grpid,axisVarid[0],shuffle,deflate,deflate_level)))
       err(retval_);

    if((retval_ = nc_def_var(grpid,dimName2.c_str(),xtype,vecDim,&dimid[1],&axisVarid[1])))
       err(retval_);
    if((retval_ = nc_def_var_chunking(grpid,axisVarid[1],NC_CHUNKED,&xchunk)))
       err(retval_);
    if((retval_ = nc_def_var_deflate(grpid,axisVarid[1],shuffle,deflate,deflate_level)))
       err(retval_);


    // Define the Matrix variable.
    if((retval_ = nc_def_var(grpid,dataSet.c_str(),ztype,ndim,&dimid[0],&varid)))
       err(retval_);

    if((retval_ = nc_def_var_chunking(grpid,varid,NC_CHUNKED,&chunks[0])))
       err(retval_);

    if((retval_ = nc_def_var_deflate(grpid,varid,shuffle,deflate,deflate_level)))
       err(retval_);

    if((retval_ = nc_enddef(grpid)))
        err(retval_);

    // Write the data to the file.
    if((retval_ = nc_put_var(grpid, axisVarid[0], &rowAxisY[0])))
        err(retval_);
    if((retval_ = nc_put_var(grpid, axisVarid[1], &colAxisX[0])))
        err(retval_);
    if((retval_ = nc_put_var(grpid, varid, &vm.v[0])))
       err(retval_);

    return varid;
}


template<typename T>
void FileNetcdf::write_AttNC(const string dataSet, const string attName,
                             const vector<T> &attVal, nc_type xtype, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    size_t len;

    if (dataSet == "")
    {
        varid = NC_GLOBAL;
    }
    else
    {
        if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
            err(retval_);
    }

    len = attVal.size();

    if((retval_ = nc_put_att(grpid, varid, attName.c_str(), xtype, len, &attVal[0]) ))
        err(retval_);
}


template<typename T>
void FileNetcdf::write_DefHypVecNC(const string dataSet, nc_type xtype, int grpid,
                       size_t chunk, int deflate_level, int shuffle)
{
    if(grpid == 0) grpid = ncid_;

    int dimid;
    int varid;
    int ndim = 1;
    int deflate = 0;
    size_t N = NC_UNLIMITED;;

    // Set chunking, shuffle, and deflate.
    shuffle = NC_SHUFFLE;
    if(deflate_level > 0 && deflate_level < 10)
        deflate = 1;

    // Define the dimensions.
    if((retval_ = nc_def_dim(grpid, dataSet.c_str(), N, &dimid)))
        err(retval_);

    // Define the variable.
    if((retval_ = nc_def_var(grpid, dataSet.c_str(), xtype, ndim,
                             &dimid, &varid)))
        err(retval_);

    if((retval_ = nc_def_var_chunking(grpid, varid, NC_CHUNKED, &chunk)))
        err(retval_);
    if((retval_ = nc_def_var_deflate(grpid, varid, shuffle, deflate,
                                     deflate_level)))
        err(retval_);

    T attVal[1] = {0};
    if((retval_ = nc_put_att(grpid, varid,"_FillValue", xtype, 1, attVal) ))
        err(retval_);

    if((retval_ = nc_enddef(grpid)))
        err(retval_);
}

template<typename T>
void FileNetcdf::write_PutHypVecNC(const string dataSet, const vector<T> &vec,
                       size_t idx, size_t len, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if ((retval_ = nc_put_vara(grpid, varid, idx, len, &vec[0]) ))
        err(retval_);
}

template<typename T>
void FileNetcdf::write_PutHypVecNC(const string dataSet, const T* vec,
                       size_t idx, size_t len, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if ((retval_ = nc_put_vara(grpid, varid, idx, len, vec) ))
        err(retval_);
}

template<typename T>
void FileNetcdf::write_CatHypVecNC(const string dataSet, const vector<T> &vec, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    size_t len = vec.size();
    int dimid;
    int varid;
    size_t idx;

    if((retval_ = nc_inq_dimid(grpid, dataSet.c_str(), &dimid) ))
        err(retval_);

    if((retval_ = nc_inq_dim(grpid, dimid, NULL, &idx) ))
        err(retval_);

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if ((retval_ = nc_put_vara(grpid, varid, &idx, &len, &vec[0]) ))
        err(retval_);
}

template<typename T>
void FileNetcdf::write_CatHypVecNC(const string dataSet, const T* vec,size_t len,int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int dimid;
    int varid;
    size_t idx;

    if((retval_ = nc_inq_dimid(grpid, dataSet.c_str(), &dimid) ))
        err(retval_);

    if((retval_ = nc_inq_dim(grpid, dimid, NULL, &idx) ))
        err(retval_);

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if ((retval_ = nc_put_vara(grpid, varid, &idx, &len, vec) ))
        err(retval_);
}

template<typename T>
void FileNetcdf::write_DefHypMatNC(const string dataSet, size_t dims[], nc_type xtype,
            int grpid,
            size_t chunk, int deflate_level, int shuffle)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int dimid[2];
    size_t chunks[2];
    int ndim = 2;
    int deflate = 0;
    string dimName1 = dataSet + "_row";
    string dimName2 = dataSet + "_col";

    // Set chunking, shuffle, and deflate.
    shuffle = NC_SHUFFLE;
    if(deflate_level > 0 && deflate_level < 10)
        deflate = 1;

    // Define the dimensions.
    if((retval_ = nc_def_dim(grpid,dimName1.c_str(),dims[0],&dimid[0])))
        err(retval_);
    if((retval_ = nc_def_dim(grpid,dimName2.c_str(),dims[1],&dimid[1])))
        err(retval_);

    //if(N[0]*N[1] < chunk) chunk = N[0]*N[1];
    if(dims[1] < chunk)
    {
        if(dims[1] != NC_UNLIMITED) chunk = dims[1];
    }

    chunks[0] = 1;
    chunks[1] = chunk;

    // Define the variable.
    if((retval_ = nc_def_var(grpid,dataSet.c_str(),xtype,ndim,&dimid[0],&varid)))
       err(retval_);

    if((retval_ = nc_def_var_chunking(grpid,varid,NC_CHUNKED,&chunks[0])))
       err(retval_);

    if((retval_ = nc_def_var_deflate(grpid,varid,shuffle,deflate,deflate_level)))
       err(retval_);

    T attVal[1] = {0};
    if((retval_ = nc_put_att(grpid, varid,"_FillValue", xtype, 1, attVal) ))
        err(retval_);

    if((retval_ = nc_enddef(grpid)))
        err(retval_);
}

template<typename T>
void FileNetcdf::write_DefHypMatNC(const string dataSet, const string rowY, const string colX,
        nc_type xtype, int grpid,
        size_t chunk, int deflate_level, int shuffle)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int dimid[2];
    size_t chunks[2];
    size_t dims[2];
    int ndim = 2;
    int deflate = 0;

    // Set chunking, shuffle, and deflate.
    shuffle = NC_SHUFFLE;
    if(deflate_level > 0 && deflate_level < 10)
        deflate = 1;

    // Define the dimensions.
    int cordvar;
    int corddims[1];

    if ((retval_ = nc_inq_varid(grpid,rowY.c_str(),&cordvar) ))
        err(retval_);
    if ((retval_ = nc_inq_vardimid(grpid,cordvar,corddims) ))
        err(retval_);
    dimid[0]=corddims[0];
    if ((retval_ = nc_inq_dimlen(grpid, dimid[0], &dims[0]) ))
        err(retval_);

    if ((retval_ = nc_inq_varid(grpid,colX.c_str(),&cordvar) ))
        err(retval_);
    if ((retval_ = nc_inq_vardimid(grpid,cordvar,corddims) ))
        err(retval_);
    dimid[1]=corddims[0];
    if ((retval_ = nc_inq_dimlen(grpid, dimid[1], &dims[1]) ))
        err(retval_);

    //if(N[0]*N[1] < chunk) chunk = N[0]*N[1];
    if(dims[1] < chunk)
    {
        if(dims[1] != NC_UNLIMITED) chunk = dims[1];
    }

    chunks[0] = 1;
    chunks[1] = chunk;

    // Define the variable.
    if((retval_ = nc_def_var(grpid,dataSet.c_str(),xtype,ndim,&dimid[0],&varid)))
       err(retval_);

    if((retval_ = nc_def_var_chunking(grpid,varid,NC_CHUNKED,&chunks[0])))
       err(retval_);

    if((retval_ = nc_def_var_deflate(grpid,varid,shuffle,deflate,deflate_level)))
       err(retval_);


    T attVal[1] = {0};
    if((retval_ = nc_put_att(grpid, varid,"_FillValue", xtype, 1, attVal) ))
        err(retval_);

    if((retval_ = nc_enddef(grpid)))
        err(retval_);
}

template<typename T>
void FileNetcdf::write_PutHypMatNC(const string dataSet, const VecMat<T> &vm,
        size_t rcIdx[2], size_t len[2], int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if ((retval_ = nc_put_vara(grpid, varid, rcIdx, len, &vm.v[0]) ))
        err(retval_);
}

template<typename T>
void FileNetcdf::write_PutHypMatNC(const string dataSet, const T *vm,
        size_t rcIdx[2], size_t len[2], int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
        err(retval_);

    if ((retval_ = nc_put_vara(grpid, varid, rcIdx, len, vm) ))
        err(retval_);
}


#endif

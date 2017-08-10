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


#include "FileNetcdf.hpp"
#include <cstring>


FileNetcdf::FileNetcdf(const string _fileName, int omode) : fileName_(_fileName)
{
    switch(omode)
    {
    case NC_NOWRITE:
        if ((retval_ = nc_open(fileName_.c_str(), omode, &ncid_)))
            err(retval_);
        fileStatus_ = true;
        break;
    case NC_WRITE:
        if ((retval_ = nc_open(fileName_.c_str(), omode, &ncid_)))
            err(retval_);
        fileStatus_ = true;
        break;
    case NC_NETCDF4:
       if ((retval_ = nc_create(fileName_.c_str(), omode|NC_CLOBBER, &ncid_)))
          err(retval_);
       fileStatus_ = true;
       break;
    }
}

void FileNetcdf::open(const string _fileName, int omode)
{
    if (fileStatus_ == false)
    {
        fileName_=_fileName;
        switch(omode)
        {
        case NC_NOWRITE:
            if ((retval_ = nc_open(fileName_.c_str(), omode, &ncid_)))
                err(retval_);
            fileStatus_ = true;
            break;
        case NC_WRITE:
            if ((retval_ = nc_open(fileName_.c_str(), omode, &ncid_)))
                err(retval_);
            fileStatus_ = true;
            break;
        case NC_NETCDF4:
            if ((retval_ = nc_create(fileName_.c_str(), omode|NC_CLOBBER, &ncid_)))
                err(retval_);
            fileStatus_ = true;
            break;
        }
    }
    else
    {
        throw runtime_error("Error: File already opened");
    }
}

void FileNetcdf::close(void)
{
    if (fileStatus_ == true)
    {
        if ((retval_ = nc_close(ncid_)))
            err(retval_);
        fileStatus_ = false;
    }
    else
    {
        throw runtime_error("Error: No file to close");
    }
}

FileNetcdf::~FileNetcdf()
{
    if(fileStatus_ == true)
    {
        if((retval_ = nc_close(ncid_)))
            err(retval_);
    }
}

int FileNetcdf::read_VarIDNC(const string dataSet, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;

    retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid);

    if(retval_ != 0) varid = -1;

    return varid;
}


vector<size_t> FileNetcdf::read_DimNC(const string dataSet, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    vector<size_t> dimSize;

    if((retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid) ))
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

int FileNetcdf::create_Group(const string name, int grpid)
{
    if (grpid == 0) grpid = ncid_;

    int varid;
    if(( retval_ = nc_def_grp(grpid, name.c_str(), &varid) ))
        err(retval_);

    return varid;
}

int FileNetcdf::open_Group(const string name, int grpid)
{
    if (grpid == 0) grpid = ncid_;

    int ncid;
    retval_ = nc_inq_ncid(grpid, name.c_str(), &ncid);

    if(retval_ != 0) ncid = -1;

    return ncid;
}

int FileNetcdf::search_Group(const string dataSet, int grpid)
{
    if(grpid == 0) grpid = ncid_;

    int nGrps;
    vector<int> ngrpids;
    string strGrp;
    size_t strGrpL;

    int nVars;
    int varid;

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

    // Scan Groups
    for(int i = 0; i < nGrps; ++i)
    {
        search_Group(dataSet, ngrpids[i]);
    }
    // Scan Variables
    if(( retval_ = nc_inq_varids(grpid, &nVars, NULL) ))
        err(retval_);
    if(nVars > 0)
    {
        retval_ = nc_inq_varid(grpid, dataSet.c_str(), &varid);
        if(retval_ == NC_NOERR)
        {
            this->dataSetList_.push_back(InfoGrpVar(grpid,varid,strGrp,dataSet));
            return 1;
        }
    }
    return 0;
}


void FileNetcdf::read(Matrix& a, const string name, int grpid)
{
    vector<size_t> dims = read_DimNC(name, grpid);
    if (dims.size() == 1)
    {
        // HDF5 1D dataset is treated as a row vector
        vector<fp> v;
        read_VecNC(name, v, grpid);
        a.importFromArray(1, v.size(), v.data());
    }
    else
    {
        VecMat<fp> vm;
        read_MatNC(name, vm, grpid);
        uli dims[2];
        vm.getDims(dims);
        a.importFromArray(dims[0], dims[1], vm.v.data());
    }
}


void FileNetcdf::write(const Matrix& a, const string name, int grpid)
{
    if (a.m() == 1)
    {
        // HDF5 1D dataset is treated as a row vector
        write_VecNC(name, a.vs(), a.n(), sizeof(fp) == 4 ? NC_FLOAT : NC_DOUBLE, grpid);
    }
    else
    {
        VecMat<fp> vm(a.m(), a.n());
        memcpy(vm.v.data(), a.vs(), sizeof(fp) * vm.v.size());
        write_MatNC(name, vm, sizeof(fp) == 4 ? NC_FLOAT : NC_DOUBLE, grpid);
    }
}


void FileNetcdf::read(MatrixSparse& a, const string name, int grpid)
{
    int grpidMat = open_Group(name, grpid);

    vector<ii> m(1);
    read_AttNC("m", NC_GLOBAL, m, grpidMat);

    vector<ii> n(1);
    read_AttNC("n", NC_GLOBAL, n, grpidMat);

    if (read_VarIDNC("v", grpidMat) != -1)
    {
        vector<ii> rowind;
        read_VecNC("i", rowind, grpidMat);

        vector<ii> colind;
        read_VecNC("j", colind, grpidMat);

        vector<fp> acoo;
        read_VecNC("v", acoo, grpidMat);

        a.importFromCoo(m[0], n[0], acoo.size(), rowind.data(), colind.data(), acoo.data());
    }
    else
    {
        a.init(m[0], n[0]);
    }
}


int FileNetcdf::write(const MatrixSparse& a, const string name, int grpid)
{
    int grpidMat = create_Group(name, grpid);

    vector<ii> m(1); m[0] = a.m();
    write_AttNC("", "m", m, sizeof(ii) == 4 ? NC_INT : NC_INT64, grpidMat);

    vector<ii> n(1); n[0] = a.n();
    write_AttNC("", "n", n, sizeof(ii) == 4 ? NC_INT : NC_INT64, grpidMat);

    vector<ii> rowind(a.nnz());
    vector<ii> colind(a.nnz());
    vector<fp> acoo(a.nnz());
    if (a.nnz() > 0)
        a.exportToCoo(rowind.data(), colind.data(), acoo.data());

    write_VecNC("i", rowind, sizeof(ii) == 4 ? NC_INT : NC_INT64, grpidMat, true);
    write_VecNC("j", colind, sizeof(ii) == 4 ? NC_INT : NC_INT64, grpidMat, true);
    write_VecNC("v", acoo, sizeof(fp) == 4 ? NC_FLOAT : NC_DOUBLE, grpidMat, true);

    return grpidMat;
}


void FileNetcdf::err(int e)
{
    throw runtime_error("ERROR: '" + string(nc_strerror(e)) + "' processing " + fileName_);
}



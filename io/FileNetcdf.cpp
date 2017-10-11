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


FileNetcdf::FileNetcdf() : fileStatus_(false)
{
};


FileNetcdf::FileNetcdf(const string& filename, int omode) : filename_(filename), fileStatus_(false)
{
    open(filename, omode);
}


FileNetcdf::~FileNetcdf()
{
    close();
}


void FileNetcdf::open(const string& filename, int omode)
{
    if (!fileStatus_)
    {
        filename_ = filename;

        switch(omode)
        {
        case NC_NOWRITE:
            if ((retval_ = nc_open(filename_.c_str(), omode, &ncid_)) != NC_NOERR)
                err(retval_);
            break;
        case NC_WRITE:
            if ((retval_ = nc_open(filename_.c_str(), omode, &ncid_)) != NC_NOERR)
                err(retval_);
            break;
        case NC_NETCDF4:
            if ((retval_ = nc_create(filename_.c_str(), omode|NC_CLOBBER, &ncid_)) != NC_NOERR)
                err(retval_);
            break;
        }

        fileStatus_ = true;
    }
    else
    {
        throw runtime_error("Error: File already opened");
    }
}


void FileNetcdf::close()
{
    if (fileStatus_)
    {
        if ((retval_ = nc_close(ncid_)) != NC_NOERR)
            err(retval_);

        fileStatus_ = false;
    }
    else
    {
        throw runtime_error("Error: No file to close");
    }
}


bool FileNetcdf::exists(const string& variable, int parentId)
{
    if (parentId == 0)
        parentId = ncid_;

    int varid;
    return nc_inq_varid(parentId, variable.c_str(), &varid) == NC_NOERR;
}


void FileNetcdf::readExtent(vector<size_t>& extent, const string& name, int parentId)
{
    if (parentId == 0)
        parentId = ncid_;

    int varid;
    int ndim;
    vector<int> dimid;
    nc_type typId;

    if ((retval_ = nc_inq_varid(parentId, name.c_str(), &varid) ) != NC_NOERR)
        err(retval_);

    if ((retval_ = nc_inq_vartype(parentId, varid, &typId)) != NC_NOERR)
        err(retval_);

    if ((retval_ = nc_inq_varndims(parentId,varid,&ndim)) != NC_NOERR)
        err(retval_);

    dimid.resize(ndim);
    extent.resize(ndim);

    if ((retval_ = nc_inq_vardimid(parentId, varid, &dimid[0])) != NC_NOERR)
        err(retval_);

    for(int i = 0; i < ndim; ++i)
    {
        if ((retval_ = nc_inq_dimlen(parentId, dimid[i], &extent[i])) != NC_NOERR)
            err(retval_);
    }
}

void FileNetcdf::readAttribute(string& data, const string& name, const string& dataset, int parentId)
{
    if(parentId == 0)
        parentId = ncid_;

    int varid;
    if (dataset == "")
    {
        varid = NC_GLOBAL;
    }
    else
    {
        if((retval_ = nc_inq_varid(parentId, dataset.c_str(), &varid)) != NC_NOERR)
            err(retval_);
    }

    nc_type xtype;
    size_t len;
    if((retval_ = nc_inq_att(parentId, varid, name.c_str(), &xtype, &len)) != NC_NOERR)
        err(retval_);

    data.resize(len);

    if((retval_ =  nc_get_att_text(parentId, varid, name.c_str(), &data[0])) != NC_NOERR)
        err(retval_);
}


void FileNetcdf::writeAttribute(const string& data, const string& attribute, const string& dataset, int parentId)
{
    if(parentId == 0)
        parentId = ncid_;

    int varid;
    if (dataset == "")
    {
        varid = NC_GLOBAL;
    }
    else
    {
        if((retval_ = nc_inq_varid(parentId, dataset.c_str(), &varid)) != NC_NOERR)
            err(retval_);
    }

    if((retval_ = nc_put_att_text(parentId, varid, attribute.c_str(), data.size(), data.data())) != NC_NOERR)
        err(retval_);
}



size_t FileNetcdf::readSize(const string& name, int parentId)
{
    vector<size_t> extent;
    readExtent(extent, name, parentId);

    size_t size = 1;
    for(size_t i = 0; i < extent.size(); ++i)
        size *= extent[i];

    return size;
}


int FileNetcdf::createGroup(const string& name, int parentId)
{
    if (parentId == 0)
        parentId = ncid_;

    int groupId;
    if((retval_ = nc_def_grp(parentId, name.c_str(), &groupId)) != NC_NOERR)
        err(retval_);

    return groupId;
}


int FileNetcdf::openGroup(const string& name, int parentId)
{
    if (parentId == 0)
        parentId = ncid_;

    int groupId;
    retval_ = nc_inq_ncid(parentId, name.c_str(), &groupId);

    if(retval_ != NC_NOERR)
        groupId = -1;

    return groupId;
}


int FileNetcdf::searchGroup(const string dataSet, int parentId)
{
    if(parentId == 0) parentId = ncid_;

    int nGrps;
    vector<int> ngrpids;
    string strGrp;
    size_t strGrpL;

    int nVars;
    int varid;

    // How many Groups
    if(( retval_ = nc_inq_grps(parentId, &nGrps, NULL) ))
        err(retval_);
    if(nGrps > 0 )
    {
        ngrpids.resize(nGrps);
        if(( retval_ = nc_inq_grps(parentId, NULL, &ngrpids[0]) ))
            err(retval_);
    }

    if(( retval_ = nc_inq_grpname_len(parentId, &strGrpL)))
        err(retval_);
    strGrp.resize(strGrpL);
    if(( retval_ = nc_inq_grpname_full(parentId,NULL,&strGrp[0])))
        err(retval_);

    // Scan Groups
    for(int i = 0; i < nGrps; ++i)
    {
        searchGroup(dataSet, ngrpids[i]);
    }
    // Scan Variables
    if(( retval_ = nc_inq_varids(parentId, &nVars, NULL) ))
        err(retval_);
    if(nVars > 0)
    {
        retval_ = nc_inq_varid(parentId, dataSet.c_str(), &varid);
        if(retval_ == NC_NOERR)
        {
            this->dataSetList_.push_back(InfoGrpVar(parentId,varid,strGrp,dataSet));
            return 1;
        }
    }
    return 0;
}


void FileNetcdf::readMatrix(Matrix &a, const string& name, int grpid)
{
    vector<size_t> dims;
    readExtent(dims, name, grpid);
    if (dims.size() == 1)
    {
        // HDF5 1D dataset is treated as a row vector
        vector<fp> v;
        readVector(v, name, grpid);
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


void FileNetcdf::writeMatrix(const Matrix& a, const string& name, int parentId)
{
    if (a.m() == 1)
    {
        // HDF5 1D dataset is treated as a row vector
        writeVector(a.vs(), a.n(), name, parentId);
    }
    else
    {
        VecMat<fp> vm(a.m(), a.n());
        memcpy(vm.v.data(), a.vs(), sizeof(fp) * vm.v.size());
        write_MatNC(name, vm, parentId);
    }
}


int FileNetcdf::readMatrixSparseCoo(MatrixSparse& a, const string& name, int parentId)
{
    int matrixId = openGroup(name, parentId);

    ii m = readAttribute<ii>("m", "", matrixId);
    ii n = readAttribute<ii>("n", "", matrixId);

    vector<ii> is;
    readVector(is, "i", matrixId);

    vector<ii> js;
    readVector(js, "j", matrixId);

    vector<fp> vs;
    readVector(vs, "v", matrixId);

    a.importFromCoo(m, n, vs.size(), is.data(), js.data(), vs.data());

    return matrixId;
}


int FileNetcdf::writeMatrixSparseCoo(const MatrixSparse& a, const string& name, int parentId)
{
    int matrixId = createGroup(name, parentId);

    writeAttribute(a.m(), "m", "", matrixId);
    writeAttribute(a.n(), "n", "", matrixId);

    vector<ii> is(a.nnz());
    vector<ii> js(a.nnz());
    vector<fp> vs(a.nnz());
    if (a.nnz() > 0)
        a.exportToCoo(is.data(), js.data(), vs.data());

    writeVector(is, "i", matrixId);
    writeVector(js, "j", matrixId);
    writeVector(vs, "v", matrixId);

    return matrixId;
}


int FileNetcdf::readMatrixSparseCsr(MatrixSparse& a, const string& name, int parentId)
{
    int matrixId = openGroup(name, parentId);

    ii m = ii(readSize("ij", matrixId)) - 1;
    ii n = readAttribute<ii>("n", "", matrixId);
    ii nnz = ii(readSize("j", matrixId));

    a.createCsr(m, n, nnz);

    readVector(a.ijs(), "ij", matrixId);
    readVector(a.js(), "j", matrixId);
    readVector(a.vs(), "v", matrixId);

    a.commitCsr(true);

    return matrixId;
}


int FileNetcdf::writeMatrixSparseCsr(const MatrixSparse& a, const string& name, int parentId)
{
    int matrixId = createGroup(name, parentId);

    writeAttribute(a.n(), "n", "", matrixId);

    if (a.initCsr(a.nnz() > 0))
    {
        writeVector(a.ijs(), a.m() + 1, "ij", matrixId);
    }
    else
    {
        vector<ii> ijs(a.m() + 1, 0);
        writeVector(ijs, "ij", matrixId);
    }

    writeVector(a.js(), a.nnz(), "j", matrixId);
    writeVector(a.vs(), a.nnz(), "v", matrixId);

    return matrixId;
}


void FileNetcdf::err(int e)
{
    throw runtime_error("ERROR: '" + string(nc_strerror(e)) + "' processing " + filename_);
}




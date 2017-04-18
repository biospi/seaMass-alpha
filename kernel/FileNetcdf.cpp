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


void FileNetcdf::write(const MatrixSparse& a, const string name, int grpid)
{
	int grpidMat = create_Group(name, grpid);

    vector<ii> m(1); m[0] = a.m();
    write_AttNC("", "m", m, sizeof(ii) == 4 ? NC_INT : NC_INT64, grpidMat);

    vector<ii> n(1); n[0] = a.n();
    write_AttNC("", "n", n, sizeof(ii) == 4 ? NC_INT : NC_INT64, grpidMat);

    if (a.nnz() > 0)
    {
        vector<ii> is;
        vector<ii> js;
        vector<fp> vs;
        a.exportTo(is, js, vs);

        write_VecNC("i", is, sizeof(ii) == 4 ? NC_INT : NC_INT64, grpidMat);
        write_VecNC("j", js, sizeof(ii) == 4 ? NC_INT : NC_INT64, grpidMat);
        write_VecNC("v", vs, sizeof(fp) == 4 ? NC_FLOAT : NC_DOUBLE, grpidMat);
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
        vector<ii> is;
        read_VecNC("i", is, grpidMat);

        vector<ii> js;
        read_VecNC("j", js, grpidMat);

        vector<fp> vs;
        read_VecNC("v", vs, grpidMat);

        a.copy(m[0], n[0], is, js, vs);
    }
    else
    {
        a.init(m[0], n[0]);
    }
}



void FileNetcdf::err(int e)
{
    throw runtime_error("ERROR: " + string(nc_strerror(e)));
}



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

#ifndef SMPEAK_NETCDFILE_HPP_
#define SMPEAK_NETCDFILE_HPP_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <netcdf.h>
#include <fstream>
#include <typeinfo>

#include "iomath.hpp"

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

	int read_VarIDNC(const string dataSet, int grpid = 0);

	template<typename T>
	void read_MatNC(const string dataSet, VecMat<T> &vm, int grpid = 0);
	template<typename T>
	void read_MatNCT(const string dataSet, VecMat<T> &vm, int grpid = 0);
	template<typename T>
	void read_MatNC(const string dataSet, vector<vector<T> > &vm, int grpid = 0);
	template<typename T>
	void read_MatNCT(const string dataSet, vector<vector<T> > &vm, int grpid = 0);
	template<typename T>
	vector<size_t> read_DimNC(const string dataSet, int grpid = 0);
	template<typename T>
	void read_AttNC(const string attName, int varid, vector<T> &attVal, int grpid = 0);
	vector<size_t> read_DimNC(const string dataSet, int grpid = 0);

	template<typename T>
	void read_HypVecNC(const string dataSet, vector<T> &vm,
			size_t *rcIdx, size_t *len, int grpid = 0);
	template<typename T>
	void read_HypMatNC(const string dataSet, VecMat<T> &vm,
			size_t *rcIdx, size_t *len, int grpid = 0);

	int search_Group(const string dataSet, int grpid = 0);
	template<typename T>
	T search_Group(size_t level, int grpid = 0);

	template<typename T>
	int write_VecNC(const string dataSet, vector<T> &vec, nc_type xtype,
			int grpid = 0, bool unlim = false,
			size_t chunks = 4096, int deflate_level = 1,
			int shuffle = NC_SHUFFLE);
	template<typename T>
	int write_MatNC(const string dataSet, VecMat<T> &vm, nc_type xtype,
			int grpid = 0, const string rowY="", const string colX="",
			size_t chunk = 4096, int deflate_level = 1,
			int shuffle = NC_SHUFFLE);
	template<typename T>
	void write_AttNC(const string dataSet, const string attName,
			vector<T> &attVal, nc_type xtype, int grpid = 0);

	template<typename T>
	void write_DefHypMatNC(const string dataSet, size_t dims[], nc_type xtype,
			int grpid = 0,
			size_t chunk = 4096, int deflate_level = 1, int shuffle = NC_SHUFFLE);
	template<typename T>
	void write_DefHypMatNC(const string dataSet, const string rowY, const string colX, nc_type xtype,
			int grpid = 0,
			size_t chunk = 4096, int deflate_level = 1, int shuffle = NC_SHUFFLE);
	template<typename T>
	void write_PutHypMatNC(const string dataSet, VecMat<T> &vm,
		size_t rcIdx[2], size_t len[2], int grpid = 0);
	template<typename T>
	void write_PutHypMatNC(const string dataSet, T *vm,
		size_t rcIdx[2], size_t len[2], int grpid = 0);

	vector<InfoGrpVar> get_Info(void) {return dataSetList;};
private:
	string fileName;
	bool fileStatus;
	int ncid;
	int retval;
	vector<InfoGrpVar> dataSetList;
	void err(int e);
};

#include "NetCDFile.tpp"

#endif /* SMPEAK_NETCDFILE_HPP_ */

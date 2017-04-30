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

#ifndef SEAMASS_IO_FILENETCDF_HPP
#define SEAMASS_IO_FILENETCDF_HPP


#include "../kernel/VecMat.hpp"
#include "../kernel/kernel.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <netcdf.h>
#include <fstream>
#include <typeinfo>


struct InfoGrpVar
{
    InfoGrpVar(int _g, int _v, string _strg, string _strv):
        grpid(_g), varid(_v), grpName(_strg), varName(_strv){};
    int grpid;
    int varid;
    string grpName;
    string varName;
};


class FileNetcdf
{
public:
    FileNetcdf(void) : fileStatus_(false), ncid_(0),retval_(0){};
    FileNetcdf(const string _fileName, int omode = NC_NOWRITE);
    void open(const string _fileName, int omode = NC_NOWRITE);
    void close(void);
    ~FileNetcdf();

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

    int create_Group(const string name, int grpid = 0);
    int open_Group(const string name, int grpid = 0);

    int search_Group(const string dataSet, int grpid = 0);
    template<typename T>
    T search_Group(size_t level, int grpid = 0);

    template<typename T>
    int write_VecNC(const string dataSet, const vector<T> &vec, nc_type xtype,
            int grpid = 0, bool unlim = false,
            size_t chunks = 1048576, int deflate_level = 4,
            int shuffle = NC_SHUFFLE);
    template<typename T>
    int write_VecNC(const string dataSet, const T *vec, size_t len, nc_type xtype,
            int grpid = 0, bool unlim = false,
            size_t chunks = 1048576, int deflate_level = 4,
            int shuffle = NC_SHUFFLE);
    template<typename T>
    int write_MatNC(const string dataSet, const VecMat<T> &vm, nc_type xtype,
            int grpid = 0, const string rowY="", const string colX="",
            size_t chunk = 1048576, int deflate_level = 4,
            int shuffle = NC_SHUFFLE);
    template<typename T,typename X, typename Y>
    int write_MatAxisNC(const string dataSet, const VecMat<T> &vm, nc_type ztype,
                    vector<X> colAxisX, nc_type xtype,
                    vector<Y> rowAxisY, nc_type ytype,
                    const string colX="", const string rowY="",
                    int grpid = 0,
                    size_t chunk = 1048576, int deflate_level = 4,
                    int shuffle = NC_SHUFFLE);
    template<typename T>
    void write_AttNC(const string dataSet, const string attName,
                     const vector<T> &attVal, nc_type xtype, int grpid = 0);

    template<typename T>
    void write_DefHypVecNC(const string dataSet, nc_type xtype, int grpid = 0,
            size_t chunk = 1048576, int deflate_level = 4, int shuffle = NC_SHUFFLE);
    template<typename T>
    void write_PutHypVecNC(const string dataSet, const vector<T> &vec,
        size_t idx, size_t len, int grpid = 0);
    template<typename T>
    void write_PutHypVecNC(const string dataSet, const T* vec,
        size_t idx, size_t len, int grpid = 0);
    template<typename T>
    void write_CatHypVecNC(const string dataSet, const vector<T> &vec, int grpid = 0);
    template<typename T>
    void write_CatHypVecNC(const string dataSet, const T* vec,size_t len,int grpid = 0);

    template<typename T>
    void write_DefHypMatNC(const string dataSet, size_t dims[], nc_type xtype,
            int grpid = 0,
            size_t chunk = 1048576, int deflate_level = 4, int shuffle = NC_SHUFFLE);
    template<typename T>
    void write_DefHypMatNC(const string dataSet, const string rowY, const string colX, nc_type xtype,
            int grpid = 0,
            size_t chunk = 1048576, int deflate_level = 4, int shuffle = NC_SHUFFLE);
    template<typename T>
    void write_PutHypMatNC(const string dataSet, const VecMat<T> &vm,
        size_t rcIdx[2], size_t len[2], int grpid = 0);
    template<typename T>
    void write_PutHypMatNC(const string dataSet, const T *vm,
        size_t rcIdx[2], size_t len[2], int grpid = 0);

    void read(Matrix& a, const string name, int grpid = 0);
    void write(const Matrix& a, const string name, int grpid = 0);

    void read(MatrixSparse& a, const string name, int grpid = 0);
    void write(const MatrixSparse& a, const string name, int grpid = 0);

    vector<InfoGrpVar> get_Info(void) {return dataSetList_;};
private:
    string fileName_;
    bool fileStatus_;
    int ncid_;
    int retval_;
    vector<InfoGrpVar> dataSetList_;
    void err(int e);
};


#include "FileNetcdf.tpp"


#endif

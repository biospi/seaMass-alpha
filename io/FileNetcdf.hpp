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


#include "VecMat.hpp"
#include <Matrix.hpp>
#include <MatrixSparse.hpp>
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
    FileNetcdf();
    FileNetcdf(const string& filename, int omode = NC_NOWRITE);
    ~FileNetcdf();

    void open(const string& filename, int omode = NC_NOWRITE);
    void close(void);

    bool exists(const string& variable, int parentId = 0);

    int createGroup(const string& name, int parentId = 0);
    int openGroup(const string& name, int parentId = 0);

    void readExtent(vector<size_t>& extent, const string& name, int parentId = 0);
    size_t readSize(const string& name, int parentId = 0);

    template<typename T>
    void readAttribute(vector<T>& data, const string& name, const string& dataset, int parentId = 0);
    template<typename T>
    void writeAttribute(const vector<T> &data, const string& attribute, const string& dataset,
                        int parentId = 0);

    template<typename T>
    T readAttribute(const string& name, const string& dataset, int parentId = 0);
    template<typename T>
    void writeAttribute(T data, const string& attribute, const string& dataset, int parentId = 0);

    template<typename T>
    void readVector(T* data, const string& dataset, int parentId = 0);
    template<typename T>
    void writeVector(const T* data, size_t length, const string& dataset,
                     int parentId = 0, size_t chunkExtent = 1048576, int deflateLevel = 4);

    template<typename T>
    void readVector(vector<T> &data, const string& dataset, int parentId = 0);
    template<typename T>
    void writeVector(const vector<T>& data, const string& dataset,
                     int parentId = 0, size_t chunkExtent = 1048576, int deflateLevel = 4);

    void readMatrix(Matrix &a, const string& name, int grpid = 0);
    void writeMatrix(const Matrix &a, const string& name, int parentId = 0);

    int readMatrixSparseCoo(MatrixSparse &a, const string& name, int parentId = 0);
    int writeMatrixSparseCoo(const MatrixSparse &a, const string& name, int parentId = 0);

    int readMatrixSparseCsr(MatrixSparse &a, const string& dataset, int parentId = 0);
    int writeMatrixSparseCsr(const MatrixSparse &a, const string& name, int parentId = 0);


    // make clearer?

    int searchGroup(const string dataSet, int parentId = 0);


    // deprecate?

    template<typename T>
    T searchGroup(size_t level, int parentId = 0);


    template<typename T>
    void read_MatNC(const string dataSet, VecMat<T> &vm, int grpid = 0);
    template<typename T>
    void read_MatNCT(const string dataSet, VecMat<T> &vm, int grpid = 0);
    template<typename T>
    void read_MatNC(const string dataSet, vector<vector<T> > &vm, int grpid = 0);
    template<typename T>
    void read_MatNCT(const string dataSet, vector<vector<T> > &vm, int grpid = 0);


    template<typename T>
    void read_HypVecNC(const string dataSet, vector<T> &vm,
            size_t *rcIdx, size_t *len, int grpid = 0);
    template<typename T>
    void read_HypMatNC(const string dataSet, VecMat<T> &vm,
            size_t *rcIdx, size_t *len, int grpid = 0);




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

    vector<InfoGrpVar> get_Info(void) {return dataSetList_;};

private:
    template<typename T>
    int getType() const;

    void err(int e);

    string filename_;
    bool fileStatus_;
    int ncid_;
    int retval_;
    vector<InfoGrpVar> dataSetList_;
};


#include "FileNetcdf.tpp"


#endif

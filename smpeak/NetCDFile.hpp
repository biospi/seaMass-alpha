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
	void read_VecNC(const string dataSet, vector<T> &vecData, int grpid = NULL);
	template<typename T>
	void read_VecNC(const string dataSet, T *vecData, int grpid = NULL);

	template<typename T>
	void read_MatNC(const string dataSet, VecMat<T> &vm, int grpid = NULL);
	template<typename T>
	void read_MatNCT(const string dataSet, VecMat<T> &vm, int grpid = NULL);
	template<typename T>
	void read_MatNC(const string dataSet, vector<vector<T> > &vm, int grpid = NULL);
	template<typename T>
	void read_MatNCT(const string dataSet, vector<vector<T> > &vm, int grpid = NULL);
	template<typename T>
	void read_AttNC(const string attName, int varid, vector<T> &attVal, int grpid = NULL);
	vector<size_t> read_DimNC(const string dataSet, int grpid = NULL);

	template<typename T>
	void read_HypVecNC(const string dataSet, vector<T> &vm,
			size_t *rcIdx, size_t *len, int grpid);
	template<typename T>
	void read_HypMatNC(const string dataSet, VecMat<T> &vm,
			size_t *rcIdx, size_t *len, int grpid);

	int search_Group(const string dataSet, int grpid = NULL);

	template<typename T>
	T search_Group(size_t level, int grpid = NULL);

	template<typename T>
	int write_VecNC(const string dataSet, vector<T> &vec, nc_type xtype,
			int grpid = NULL, bool unlim = false,
			size_t chunks = 4096, int deflate_level = 4,
			int shuffle = NC_SHUFFLE);
	template<typename T>
	int write_MatNC(const string dataSet, VecMat<T> &vm, nc_type xtype,
			int grpid = NULL, const string rowY="", const string colX="",
			size_t chunk = 4096, int deflate_level = 4,
			int shuffle = NC_SHUFFLE);
	template<typename T>
	void write_AttNC(const string dataSet, const string attName,
			vector<T> &attVal, nc_type xtype, int grpid = NULL);

	template<typename T>
	void write_DefUMatNC(const string dataSet, size_t dims[], nc_type xtype,
			int grpid = NULL,
			size_t chunk = 4096, int deflate_level = 4, int shuffle = NC_SHUFFLE);	template<typename T>
	void write_DefUMatNC(const string dataSet, const string rowY, const string colX, nc_type xtype,
			int grpid = NULL,
			size_t chunk = 4096, int deflate_level = 4, int shuffle = NC_SHUFFLE);
	template<typename T>
	void write_PutUMatNC(const string dataSet, VecMat<T> &vm,
		size_t rcIdx[2], size_t len[2], int grpid = NULL);

	vector<InfoGrpVar> get_Info(void) {return dataSetList;};
private:
	string fileName;
	int ncid;
	int retval;
	bool fileStatus;
	vector<InfoGrpVar> dataSetList;
	void err(int e);
};

void mzMLdump(const string fileName, string data);

#include "NetCDFile.tpp"

#endif /* SMPEAK_NETCDFILE_HPP_ */

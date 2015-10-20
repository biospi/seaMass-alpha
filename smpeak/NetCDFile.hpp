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

	template<typename T>
	void read_AttNC(const string attName, int varid, vector<T> &attVal, int grpid = 0);

	template<typename T>
	void write_VecNC(const string dataSet, vector<T> &vec, nc_type xtype,
			size_t chunks = 4096, int deflate_level = 5, int shuffle = NC_SHUFFLE, int grpid = 0);
	template<typename T>
	void write_MatNC(const string dataSet, VecMat<T> &vm, nc_type xtype,
			size_t chunk = 64, int deflate_level = 5, int shuffle = NC_SHUFFLE, int grpid = 0);
	template<typename T>
	void write_AttNC(const string dataSet, const string attName,
			vector<T> &attVal, nc_type xtype, int grpid = 0);

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

//
// Created by ranjeet.
//

#ifndef SEAMASS_MSFILEDATA_HPP
#define SEAMASS_MSFILEDATA_HPP

#include <vector>
#include <string>
#include "NetCDFile.hpp"

struct spectrum
{
    size_t index;
    unsigned short preset_config;
    double precursor_mz;
    double scan_start_time;
    size_t scan_start_time_index;
    size_t count;
};

class MassSpecFile;

class FileFactory
{
public:
    //FileFactory(string fName);
    //virtual ~FileFactory(void){};
    static MassSpecFile* createFileObj(string fName);
//protected:
//    string fileName;
//    string id;
};


class MassSpecFile
{
protected:
    virtual void extractData(void) = 0;
    vector<spectrum> msSpec;
    string fileName;
public:
    virtual vector<spectrum> getSpectrum() = 0;
    virtual void getScanMZ(vector<double> &mz, size_t index, size_t count) = 0;
    virtual void getScanIntensities(vector<double> &intensities, size_t index, size_t count) = 0;
    virtual ~MassSpecFile() {};
};


class MSmzMLb3: public MassSpecFile
{
public:
    MSmzMLb3(string fileName);
    vector<spectrum> getSpectrum();
    void getScanMZ(vector<double> &mz, size_t index, size_t count);
    void getScanIntensities(vector<double> &intensities, size_t index, size_t count);
    //vector<double> getScanMZ();
    //vector<double> getScanIntensities();
    virtual ~MSmzMLb3(){};
protected:
    void extractData(void);
    NetCDFile mzMLb3File;
    vector<char> mzMLBuff;
    vector<InfoGrpVar> dataSetList;
    vector<spectrum> spectra;
	vector<size_t> hypIdx;
	vector<size_t> rdLen;
};


/*
class MScsv: public SMFile
{
public:
    virtual vector<spectrum> getSpectrum() = 0;
    virtual vector<double> getScanMZ() = 0;
    virtual vector<double> getScanIntensities() = 0;
    virtual ~MScsv(){};
protected:
};
*/

#endif //SEAMASS_MSFILEDATA_HPP

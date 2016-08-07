//
// Created by ranjeet.
//

#ifndef SEAMASS_MSFILEDATA_HPP
#define SEAMASS_MSFILEDATA_HPP

#include <vector>
#include <string>
#include "NetCDFile.hpp"
#include "../core/seaMass.hpp"

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
    static MassSpecFile* createFileObj(string fName);
};


class MassSpecFile
{
protected:
    virtual void extractData(void) = 0;
    vector<spectrum> msSpec;
    string fileName;
    unsigned long instrument_type;
public:
    virtual vector<spectrum> getSpectrum() = 0;
    virtual void getScanMZ(vector<double> &mz, size_t index, size_t count) = 0;
    virtual void getScanIntensities(vector<double> &intensities, size_t index, size_t count) = 0;
    virtual unsigned long getInstrument(void) = 0;
    virtual ~MassSpecFile() {};
};


class MSmzMLb3: public MassSpecFile
{
public:
    MSmzMLb3(string fileName);
    vector<spectrum> getSpectrum();
    void getScanMZ(vector<double> &mz, size_t index, size_t count);
    void getScanIntensities(vector<double> &intensities, size_t index, size_t count);
    unsigned long getInstrument(void);
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


class MSFile
{
public:
	MSFile(string fileName);
	~MSFile();

	bool next(seaMass::Input& output);

protected:
	MassSpecFile* msFile;
	vector<spectrum> spectra;
	vector<double> scan_start_times;
	unsigned long instrument_type;
	size_t i;

	void bin_mzs_intensities();
	static bool scan_start_time_order(const spectrum& lhs, const spectrum& rhs);
	static bool seamass_order(const spectrum& lhs, const spectrum& rhs);
};


/*
class MScsv: public SMFile
{
public:
    vector<spectrum> getSpectrum() = 0;
    getScanMZ(vector<double> &mz, size_t index, size_t count) = 0;
    getScanIntensities(vector<double> &intensities, size_t index, size_t count) = 0;
    ~MScsv(){};
protected:
};
*/

#endif //SEAMASS_MSFILEDATA_HPP

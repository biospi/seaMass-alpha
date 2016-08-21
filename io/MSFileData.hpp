//
// Created by ranjeet.
//

#ifndef SEAMASS_MSFILEDATA_HPP
#define SEAMASS_MSFILEDATA_HPP

#include <vector>
#include <string>
#include "NetCDFile.hpp"
#include "../core/seaMass.hpp"


struct spectrumMetaData
{
    size_t index;
	size_t count;
	double start_time;
	double finish_time;
	size_t preset_config;
	bool positive_polarity;
    double precursor_mz;
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
    vector<spectrumMetaData> msSpec;
    string fileName;
    unsigned long instrument_type;
public:
    virtual vector<spectrumMetaData>& getSpectraMetaData() = 0;
    virtual void getScanMZs(vector<double> &mz, size_t index, size_t count) = 0;
    virtual void getScanIntensities(vector<double> &intensities, size_t index, size_t count) = 0;
    virtual unsigned long getInstrument(void) = 0;
    virtual ~MassSpecFile() {};
};


class MSmzMLb3: public MassSpecFile
{
public:
    MSmzMLb3(string fileName);
	vector<spectrumMetaData>& getSpectraMetaData();
    void getScanMZs(vector<double> &mz, size_t index, size_t count);
    void getScanIntensities(vector<double> &intensities, size_t index, size_t count);
    unsigned long getInstrument(void);
    virtual ~MSmzMLb3(){};
protected:
    void extractData(void);
    NetCDFile mzMLb3File;
    vector<char> mzMLBuff;
    vector<InfoGrpVar> dataSetList;
    vector<spectrumMetaData> spectraMetaData;
	vector<size_t> hypIdx;
	vector<size_t> rdLen;
};


class MSmzMLb: public MassSpecFile
{
public:
    MSmzMLb(string fileName);
    vector<spectrumMetaData>& getSpectraMetaData();
    void getScanMZs(vector<double> &mz, size_t index, size_t count);
    void getScanIntensities(vector<double> &intensities, size_t index, size_t count);
    unsigned long getInstrument(void);
    //vector<double> getScanStartTimes(void);
    virtual ~MSmzMLb(){};
protected:
    void extractData(void);
    NetCDFile mzMLbFile;
    vector<char> mzMLBuff;
    vector<InfoGrpVar> dataSetList;
    vector<spectrumMetaData> spectraMetaData;
	vector<size_t> hypIdx;
	vector<size_t> rdLen;
	//vector<double> scan_start_times;
	vector<size_t> mzIdx;
	vector<size_t> intensitiesIdx;
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////


class InputFile
{
public:
	virtual bool next(seaMass::Input& output, std::string& id) = 0;
	virtual ~InputFile(){};
};


class mzMLbInputFile : public InputFile
{
public:
	mzMLbInputFile(string fileName);
	~mzMLbInputFile();

	virtual bool next(seaMass::Input& output, std::string& id);

protected:
	MassSpecFile* msFile;
	vector<spectrumMetaData> spectraMetaData;
	unsigned long instrument_type;
	size_t i;

	static bool scan_start_time_order(const spectrumMetaData& lhs, const spectrumMetaData& rhs);
	static bool seamass_order(const spectrumMetaData& lhs, const spectrumMetaData& rhs);
};


class SMIInputFile : public InputFile
{
public:
	SMIInputFile(string fileName) {}

	virtual bool next(seaMass::Input& output, std::string& id) { return false; }
};


#endif //SEAMASS_MSFILEDATA_HPP

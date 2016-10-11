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

#ifndef SEAMASS_MSFILEDATA_HPP
#define SEAMASS_MSFILEDATA_HPP

#include <vector>
#include <string>
#include "NetCDFile.hpp"
#include "../core/SeaMass.hpp"


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
	virtual bool next(SeaMass::Input& output, std::string& id) = 0;
	virtual ~InputFile(){};
};


class mzMLbInputFile : public InputFile
{
public:
	mzMLbInputFile(string fileName);
	~mzMLbInputFile();

	virtual bool next(SeaMass::Input& output, std::string& id);

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

	virtual bool next(SeaMass::Input& output, std::string& id) { return false; }
};


#endif //SEAMASS_MSFILEDATA_HPP

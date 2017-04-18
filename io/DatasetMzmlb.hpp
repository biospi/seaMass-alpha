//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
//
// Copyright (C) 2016  biospi Laboratory, University of Bristol, UK
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

#ifndef SEAMASS_DATASETMZMLB_HPP
#define SEAMASS_DATASETMZMLB_HPP


#include <vector>
#include <string>
#include "Dataset.hpp"
#include "../kernel/FileNetcdf.hpp"
#include <pugixml.hpp>

namespace xml = pugi;


class DatasetMzmlb: public Dataset
{
public:
	struct SpectrumMetadata
	{
		size_t mzmlSpectrumIndex; // index of spectrum in original mzML <SpectrumList> tag
    	std::string id; // id differentiates which set of spectra this spectrum is in for seaMass

		bool isProfileMode;

		double startTime;
		double finishTime;
		string startTimeString;

		std::string config;

		enum DataType { Unknown, IonCount, IonCurrent } dataType;

		size_t defaultArrayLength;
		std::string mzsDataset;
		size_t mzsOffset;
		std::string intensitiesDataset;
		size_t intensitiesOffset;
	};

	DatasetMzmlb(std::string &filename);
    virtual ~DatasetMzmlb();

    virtual bool read(Seamass::Input &input, std::string &id);
    virtual bool read(Seamass::Input &input, Seamass::Output &output, std::string &id) { return false; }

	virtual void write(const Seamass::Input& input, const std::string& id) {}
    virtual void write(const Seamass::Input& input, const Seamass::Output& output, const std::string& id);

	virtual void writeData(Seamass &sm_, Seamass::Input &input_, std::string id, bool centriod_, double threshold_);

private:
    static bool startTimeOrder(const SpectrumMetadata &lhs, const SpectrumMetadata &rhs);
    static bool seamassOrder(const SpectrumMetadata &lhs, const SpectrumMetadata &rhs);

    template<typename T>
    T getXmlValue(xml::xml_document &scan, string xpath, string attrib);
    template<typename T>
    void setXmlValue(xml::xml_document &scan, string xpath, string attrib,T value);

	void writeVecData(vector<fp>& data_);
	void writeXmlData();
	size_t idxDataArrayOffSet_;
	vector<uli> specIdx_;
	vector<uli> newSpecIdx_;

	void writePeakData(VecMat<double>& mzPeak_, VecMat<float>& pkPeak_,
					   vector<size_t>& mzpkVecSize_);
	void writePeakXmlData(vector<size_t>& mzpkVecSize_);
	void writeChromatogramXmlEnd();

    FileNetcdf file_;
	FileNetcdf fileOut_;

    vector<SpectrumMetadata> metadata_; // this will be sorted for 'next()'
    li spectrumIndex_;
    li lastSpectrumIndex_;
	li extent_;

};

#include "DatasetMzmlb.tpp"

#endif

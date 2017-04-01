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

#ifndef SEAMASS_MZMLBINPUTFILE_HPP
#define SEAMASS_MZMLBINPUTFILE_HPP


#include <vector>
#include <string>
#include "Dataset.hpp"
#include "../kernel/FileNetcdf.hpp"


struct MzmlbSpectrumMetadata
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


class DatasetMzmlb: public Dataset
{
public:
    DatasetMzmlb(std::string &filename);
    virtual ~DatasetMzmlb();

    virtual bool next(SeamassCore::Input& output, std::string& id);

private:
    static bool startTimeOrder(const MzmlbSpectrumMetadata &lhs, const MzmlbSpectrumMetadata &rhs);
    static bool seamassOrder(const MzmlbSpectrumMetadata &lhs, const MzmlbSpectrumMetadata &rhs);

    FileNetcdf file_;

    vector<MzmlbSpectrumMetadata> metadata_; // this will be sorted for 'next()'
    li spectrumIndex_;
    li lastSpectrumIndex_;
};


#endif

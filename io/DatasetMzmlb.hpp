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


struct MetadataMzmlbSpectrum
{
    size_t mzmlSpectrumIndex;
    string id;

	bool isProfileMode;
	bool isPositivePolarity;

	double startTime;
	double finishTime;

    enum DataType { Unknown, IonCount, IonCurrent } dataType;

	size_t defaultArrayLength;

	//double mz0;
	//double mz1;

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
    static bool startTimeOrder(const MetadataMzmlbSpectrum &lhs, const MetadataMzmlbSpectrum &rhs);
    static bool seamassOrder(const MetadataMzmlbSpectrum &lhs, const MetadataMzmlbSpectrum &rhs);

    FileNetcdf file_;

    vector<MetadataMzmlbSpectrum> metadata_; // this will be sorted for 'next()'
    size_t spectrumIndex_;
};


#endif

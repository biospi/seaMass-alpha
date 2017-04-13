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

#include "DatasetSmi.hpp"

DatasetSmi::DatasetSmi(std::string &filename_):processed(false)
{
	size_t lastdot = filename_.find_last_of(".");
	string outFileName=filename_.substr(0,lastdot)+".smv";

	file_.open(filename_);
	fileOut_.open(outFileName,NC_NETCDF4);
}


bool DatasetSmi::next(SeamassCore::Input &output_, std::string &id_)
{

	if(processed == true) return false;

	// clear out
	if(output_.binCounts.size() > 0)
	{
		vector<li>().swap(output_.binCountsIndex);
		vector<double>().swap(output_.startTimes);
		vector<double>().swap(output_.finishTimes);
		vector<fp>().swap(output_.exposures);
		vector<fp>().swap(output_.binCounts);
		vector<double>().swap(output_.binEdges);
	}
	file_.read_VecNC("binCounts",output_.binCounts);
	file_.read_VecNC("binCountsIndex",output_.binCountsIndex);
	file_.read_VecNC("binEdges",output_.binEdges);
	file_.read_VecNC("startTimes",output_.startTimes);
	file_.read_VecNC("finishTimes",output_.finishTimes);
	file_.read_VecNC("exposures",output_.exposures);

	processed = true;

	return true;
}

void DatasetSmi::writeData(SeamassCore &sm_, SeamassCore::Input &input_, bool centriod_,double threshold_)
{
	// New smv file creation will go here...
}

DatasetSmi::~DatasetSmi(){}

//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> manchester.ac.uk>
//
// Copyright (C) 2013  CADET Bioinformatics Laboratory, University of Manchester, UK
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


#ifndef _SEAMASSRESTORATION_SMVWRITER_HPP_
#define _SEAMASSRESTORATION_SMVWRITER_HPP_


#include "../core/seaMass.hpp"
#include <SpatialIndex.h>


/*class SMVFile
{
public:
	class Iterator
	{
	public:
		Iterator(SMVFile& parent, const vector<ii>& low, const vector<ii>& high); // Constructor 
		Iterator(SMVFile& parent, const vector<double>& low, const vector<double>& high);

		bool next(size_t n);

	protected:
		const vector<double>& low;
		const vector<double>& high;
	};

	SMVFile(const string& filename, const vector<ii> baseline_size, const vector<ii> baseline_scale, const vector<ii> baseline_offset);
	~SMVFile();

protected:
	std::string basename;
	SpatialIndex::IStorageManager* diskfile;
	SpatialIndex::StorageManager::IBuffer* file;
	SpatialIndex::ISpatialIndex* tree;

	const vector<ii> baseline_size;
	const vector<ii> baseline_scale;
	const vector<ii> baseline_offset;
};*/


#include <string>
#include <vector>
#include <SpatialIndex.h>


struct Coef
{
	int lx, ly, x, y;
	fp v;

	bool operator()(Coef& a, Coef& b);
};

class RTreeReader
{
protected:
	SpatialIndex::IStorageManager* diskfile;
	SpatialIndex::StorageManager::IBuffer* file;
	SpatialIndex::ISpatialIndex* tree;

public:
	RTreeReader(const std::string& filename);
	~RTreeReader();

	void read(SeaMass::Output& output);
};


#endif // _SEAMASSRESTORATION_SMVWRITER_HPP_


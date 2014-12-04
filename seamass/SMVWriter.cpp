//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> manchester.ac.uk>
//
// Copyright (C) 2013  CADET Laboratory for Medical Bioinformatics, University of Manchester, UK
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

#include "SMVWriter.hpp"

#include <boost/filesystem.hpp>
#include <SpatialIndex.h>


using namespace SpatialIndex;
using namespace boost;

class MyDataStream : public IDataStream
{
public:
	vector<Basis*>& bases;
	ii n_core_bases;
	vector< vector<fp> >& cs;
	ii j;
	ii y;
	ii x;
	ii id;
	IData* next;

	MyDataStream(vector<Basis*>& _bases, ii _n_core_bases, vector< vector<fp> >& _cs) : 	
		bases(_bases),
		n_core_bases(_n_core_bases),
		cs(_cs)
	{
		rewind();
		getNext();
	}

	virtual ~MyDataStream() {}

	virtual IData* getNext()
	{
		IData* ret = next;

		if (x == bases[j]->get_cm().n[0]) { x = 0; y++; }
		if (y == bases[j]->get_cm().n[1]) { y = 0; j--; }

		for (; j >= n_core_bases; j--)
        if (!bases[j]->is_transient())
		{
			for (; y < bases[j]->get_cm().n[1]; ++y)
			{
				for (; x < bases[j]->get_cm().n[0]; ++x)
				{
					double i = cs[j][x + y*bases[j]->get_cm().n[0]];
					if (i > 0.0)
					{
						double low[3] = {
							pow(2.0, -bases[j]->get_cm().l[0]) * (bases[j]->get_cm().o[0]+x-3),
							pow(2.0, -bases[j]->get_cm().l[1]) * (bases[j]->get_cm().o[1]+y-3),
							-i
						};
						double high[3] = {
							pow(2.0, -bases[j]->get_cm().l[0]) * (bases[j]->get_cm().o[0]+x+1),
							pow(2.0, -bases[j]->get_cm().l[1]) * (bases[j]->get_cm().o[1]+y+1),
							-i
						};				
						Region r(low, high, 3);

						x++;						
						next = new RTree::Data(0, 0, r, id++);
						return ret;
					}
				}
				x = 0;
			}
			y = 0;		
		}
		
		next = 0;
		return ret;
	}

	virtual bool hasNext()
	{
		return (next != 0);
	}

	virtual uint32_t size()
	{
		throw Tools::NotSupportedException("Operation not supported.");
	}

	virtual void rewind()
	{
		j = bases.size()-1;
		x = 0;
		y = 0;
		id = 0;
	}
};


SMVWriter::
SMVWriter(const string& _directory) :
	directory(_directory)
{
	filesystem::path viz_path(directory);
	if (filesystem::exists(viz_path))
	{
		for (filesystem::directory_iterator end_dir_it, it(viz_path); it!=end_dir_it; ++it)
		{
			filesystem::remove_all(it->path());
		}
	}
	else
	{
		filesystem::create_directories(viz_path);
	}
}


SMVWriter::
~SMVWriter()
{
}


// only supports cm with dimension 2 at present
void
SMVWriter::
write_cs(const string& basename, vector<Basis*>& bases, ii n_core_bases, vector< vector<fp> >& cs,
         double mz_min, double mz_max, double rt_min, double rt_max, double max_intensity) const
{
    cout << "Writing " << directory << "/" << basename << ".txt" << endl;

	boost::filesystem::path txt_path(directory);
	ostringstream oss; oss << basename << ".txt";
	txt_path /= oss.str();
	ofstream ofs(txt_path.string().c_str());
	ofs << mz_min << " " << mz_max << " " << rt_min << " " << rt_max << " " << max_intensity;

    cout << "Writing " << directory << "/" << basename << ".idx" << endl;
    cout << "Writing " << directory << "/" << basename << ".dat" << endl;

	// Create a new storage manager with the provided base name and a 4K page size.
	boost::filesystem::path viz_path(directory);
	viz_path /= basename;
	string vp = viz_path.string();
	SpatialIndex::IStorageManager* diskfile = StorageManager::createNewDiskStorageManager(vp, 4096);

	SpatialIndex::StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
	// applies a main memory random buffer on top of the persistent storage manager
	// (LRU buffer, etc can be created the same way).
 
	MyDataStream stream(bases, n_core_bases, cs);

	// Create and bulk load a new RTree with dimensionality 2, using "file" as
	// the StorageManager and the RSTAR splitting policy.
	id_type indexIdentifier;
	ISpatialIndex* tree = RTree::createAndBulkLoadNewRTree(RTree::BLM_STR, stream, *file, 0.7, 100, 100, 3, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);

	cout << "Max Intensity: " << max_intensity << endl;
	cout << *tree;
	cout << "Buffer hits: " << file->getHits() << endl;
	cout << "Index ID: " << indexIdentifier << endl;

	bool ret = tree->isIndexValid();
	if (ret == false) cout << "ERROR: Structure is invalid!" << endl;
	else cout << "The stucture seems O.K." << endl;

	delete tree;
	delete file;
	delete diskfile;
}

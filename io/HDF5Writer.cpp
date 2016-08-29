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

#include "HDF5Writer.hpp"
#include <iostream>
#include <boost/filesystem.hpp>
#include <SpatialIndex.h>

using namespace std;
using namespace SpatialIndex;
using namespace boost;


class MyDataStream : public IDataStream
{
public:
	const seaMass::Output& output;
	ii index;
	ii dimensions;
	vector<double> low;
	vector<double> high;

	MyDataStream(const seaMass::Output& _output) :
		output(_output),
		index(0),
		dimensions(output.baseline_size.size()),
		low(dimensions + 1),
		high(dimensions + 1)
	{
	}

	virtual ~MyDataStream() {}

	virtual IData* getNext()
	{
		double i = output.weights[index];
		for (ii j = 0; j < dimensions; j++)
		{
			i *= pow(2.0, -output.scales[j][index]);
			low[j] = pow(2.0, -output.scales[j][index]) * (output.offsets[j][index] - 3);
			high[j] = pow(2.0, -output.scales[j][index]) * (output.offsets[j][index] + 1);
		}
		low[dimensions] = -i;
		high[dimensions] = -i;

		return new RTree::Data(0, 0, Region(low.data(), high.data(), dimensions + 1), index++);
	}

	virtual bool hasNext()
	{
		return index <= output.weights.size();
	}

	virtual uint32_t size()
	{
		return output.weights.size();
	}

	virtual void rewind()
	{
		index = 0;
	}
};


HDF5Writer::
HDF5Writer(const string& _filename) :
	filename(_filename)
{
	file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file < 0)
	{
		// throw exception
		cerr << "problem creating smo" << endl;
		throw "problem creating smo";
	}
}


HDF5Writer::
~HDF5Writer()
{
	if (H5Fclose(file) < 0)
    {
        // throw exception
        cerr << "problem closing smo" << endl;
		throw "problem closing smo";
	}
}


void
HDF5Writer::
write_input(const seaMass::Input& input) const
{
	write("bin_edges", input.bin_edges);
	write("bin_counts", input.bin_counts);
	if (input.spectrum_index.size() > 0) write("spectrum_index", input.spectrum_index);
	if (input.start_times.size() > 0) write("start_times", input.start_times);
	if (input.finish_times.size() > 0) write("finish_times", input.finish_times);
	if (input.exposures.size() > 0) write("exposures", input.exposures);
}

void
HDF5Writer::
write_output(const seaMass::Output& output, ii shrinkage, ii tolerance, ii page_size) const
{
	// placeholder idx and dat datasets
	hsize_t dims = 1;
	hid_t fspace = H5Screate_simple(1, &dims, &dims);
	hid_t type_id = H5Tcreate(H5T_OPAQUE, 1);
	hid_t dataset = H5Dcreate(file, "seamass_index", type_id, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	{
		hsize_t dims = output.baseline_size.size();
		hid_t scale_fspace = H5Screate_simple(1, &dims, &dims);
		hid_t scale_attr = H5Acreate(dataset, "size", H5T_NATIVE_INT, scale_fspace, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(scale_attr, H5T_NATIVE_INT, &(output.baseline_size[0]));
		H5Sclose(scale_fspace);
		H5Aclose(scale_attr);
	}
	{
		hsize_t dims = output.baseline_size.size();
		hid_t scale_fspace = H5Screate_simple(1, &dims, &dims);
		hid_t scale_attr = H5Acreate(dataset, "scale", H5T_NATIVE_INT, scale_fspace, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(scale_attr, H5T_NATIVE_INT, &(output.baseline_scale[0]));
		H5Sclose(scale_fspace);
		H5Aclose(scale_attr);
	}
	{
		hsize_t dims = output.baseline_size.size();
		hid_t scale_fspace = H5Screate_simple(1, &dims, &dims);
		hid_t scale_attr = H5Acreate(dataset, "offset", H5T_NATIVE_INT, scale_fspace, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(scale_attr, H5T_NATIVE_INT, &(output.baseline_offset[0]));
		H5Sclose(scale_fspace);
		H5Aclose(scale_attr);
	}
	{
		hsize_t dims = 1;
		hid_t scale_fspace = H5Screate_simple(1, &dims, &dims);
		hid_t scale_attr = H5Acreate(dataset, "shrinkage", H5T_NATIVE_INT, scale_fspace, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(scale_attr, H5T_NATIVE_INT, &shrinkage);
		H5Sclose(scale_fspace);
		H5Aclose(scale_attr);
	}
	{
		hsize_t dims = 1;
		hid_t scale_fspace = H5Screate_simple(1, &dims, &dims);
		hid_t scale_attr = H5Acreate(dataset, "tolerance", H5T_NATIVE_INT, scale_fspace, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(scale_attr, H5T_NATIVE_INT, &tolerance);
		H5Sclose(scale_fspace);
		H5Aclose(scale_attr);
	}
	H5Dclose(dataset);

	hid_t dataset2 = H5Dcreate(file, "seamass_data", type_id, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(fspace);
	H5Dclose(dataset2);

	// Create a new storage manager with the provided base name and page size.
	string _filename = filename;
	SpatialIndex::IStorageManager* diskfile = StorageManager::createNewDiskStorageManager(_filename, page_size);

	// applies a main memory random buffer on top of the persistent storage manager (LRU buffer, etc can be created the same way).
	SpatialIndex::StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);

	// Create and bulk load a new RTree, using "file" as the StorageManager and the RSTAR splitting policy.
	MyDataStream stream(output);
	id_type indexIdentifier;
	ISpatialIndex* tree = RTree::createAndBulkLoadNewRTree(RTree::BLM_STR, stream, *file, 0.7, 100, 100, output.baseline_size.size() + 1, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);

	cout << "RTREE OUTPUT" << endl;
	cout << *tree;
	cout << "Buffer hits: " << file->getHits() << endl;
	cout << "Index ID: " << indexIdentifier << endl;

	bool ret = tree->isIndexValid();
	if (ret == false) cout << "ERROR: Structure is invalid!" << endl;

	delete tree;
	delete file;
	delete diskfile;
}


// only supports cm with dimension 2 at present
void
HDF5Writer::
write_output_control_points(const seaMass::ControlPoints& control_points) const
{
    hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl_id, 1);
    
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    float fillval = -3.0/8.0;
    H5Pset_fill_value(dcpl_id, H5T_NATIVE_FLOAT, &fillval);

	hsize_t cdims[2] = { control_points.size[1] < 128 ? static_cast<hsize_t>(control_points.size[1]) : 128, control_points.size[0] < 128 ? static_cast<hsize_t>(control_points.size[0]) : 128 };
    H5Pset_chunk(dcpl_id, 2, cdims);
    H5Pset_shuffle(dcpl_id);
    H5Pset_deflate(dcpl_id, 1);

	hsize_t dims[2] = { static_cast<hsize_t>(control_points.size[1]), static_cast<hsize_t>(control_points.size[0]) };
    hid_t fspace = H5Screate_simple(2, dims, dims);
    hid_t dataset = H5Dcreate(file, "control_points", H5T_NATIVE_FLOAT, fspace, lcpl_id,  dcpl_id, H5P_DEFAULT);

	// write
	hsize_t mdims[2] = { static_cast<hsize_t>(control_points.size[1]), static_cast<hsize_t>(control_points.size[0]) };
    hid_t mspace = H5Screate_simple(2, mdims, mdims);
	H5Dwrite(dataset, H5T_NATIVE_FLOAT, mspace, fspace, H5P_DEFAULT, control_points.coeffs.data());
    
    H5Sclose(fspace);
    H5Sclose(mspace);

    hsize_t two = 2;
    hid_t scale_fspace = H5Screate_simple(1, &two, &two);
	hid_t scale_attr = H5Acreate(dataset, "scale", H5T_NATIVE_INT, scale_fspace, H5P_DEFAULT, H5P_DEFAULT);
	int scale_val[2] = { control_points.scale[0], control_points.scale[1] };
	H5Awrite(scale_attr, H5T_NATIVE_INT, &scale_val);
	H5Sclose(scale_fspace);
	H5Aclose(scale_attr);

	hid_t offset_fspace = H5Screate_simple(1, &two, &two);
	hid_t offset_attr = H5Acreate(dataset, "offset", H5T_NATIVE_INT, offset_fspace, H5P_DEFAULT, H5P_DEFAULT);
	int offset_val[2] = { control_points.offset[0], control_points.offset[1] };
	H5Awrite(offset_attr, H5T_NATIVE_INT, &offset_val);
	H5Sclose(offset_fspace);
	H5Aclose(offset_attr);
    
    if(H5Dclose(dataset)<0) cout << "ARGH" << endl;
}


void
HDF5Writer::
write(const string& objectname, const vector<unsigned char>& cdata) const
{
	write(objectname, cdata, H5T_NATIVE_UCHAR);
}


void
HDF5Writer::
write(const string& objectname, const vector<float>& cdata) const
{
    write(objectname, cdata, H5T_NATIVE_FLOAT);
}


void
HDF5Writer::
write(const string& objectname, const vector<double>& cdata) const
{
	write(objectname, cdata, H5T_NATIVE_DOUBLE);
}


void
HDF5Writer::
write(const string& objectname, const vector<long long>& cdata) const
{
    write(objectname, cdata, H5T_NATIVE_LLONG);
}


void
HDF5Writer::
write(const string& objectname, const vector<long>& cdata) const
{
	write(objectname, cdata, H5T_NATIVE_INT);
}

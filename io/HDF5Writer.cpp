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

using namespace std;


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
write_output(const seaMass::Output& output) const
{
	write("weights", output.weights);
	// write scales as matrix
	// write offsets as matrix
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

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


#include "HDF5File.hpp"


HDF5File::
HDF5File(const string& _filename)
{
    ostringstream oss; oss << _filename << ".h5";
    filename = oss.str();

    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0)
    {
        // throw exception
    }
}


HDF5File::
~HDF5File()
{
	if (H5Fclose(file) < 0)
    {
        // throw exception
    }
}


// only supports cm with dimension 2 at present
void
HDF5File::
write_cs(const string& objectname, const CoeffsMetadata& cm, const vector<fp>& cs) const
{
    cout << "Writing " << filename << "/" << objectname << endl;
 
    hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl_id, 1);
    
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    float fillval = -3.0/8.0;
    H5Pset_fill_value(dcpl_id, H5T_NATIVE_FLOAT, &fillval);

    hsize_t cdims[2] = {cm.n[1] < 128 ? cm.n[1] : 128, cm.n[0] < 128 ? cm.n[0] : 128};
    H5Pset_chunk(dcpl_id, 2, cdims);
    H5Pset_shuffle(dcpl_id);
    H5Pset_deflate(dcpl_id, 1);

    hsize_t dims[2] = {cm.n[1], cm.n[0]};
    hid_t fspace = H5Screate_simple(2, dims, dims);
    hid_t dataset = H5Dcreate(file, objectname.c_str(), H5T_NATIVE_FLOAT, fspace, lcpl_id,  dcpl_id, H5P_DEFAULT);

	// write
    hsize_t mdims[2] = {cm.n[1], cm.n[0]};
    hid_t mspace = H5Screate_simple(2, mdims, mdims);
    H5Dwrite(dataset, H5T_NATIVE_FLOAT, mspace, fspace, H5P_DEFAULT, cs.data());
    
    H5Sclose(fspace);
    H5Sclose(mspace);

    hsize_t two = 2;
    hid_t offset_fspace = H5Screate_simple(1, &two, &two);
    hid_t offset_attr = H5Acreate(dataset, "Offset", H5T_NATIVE_INT, offset_fspace, H5P_DEFAULT, H5P_DEFAULT);
    int offset_val[2] = {cm.o[0],cm.o[1]};
    H5Awrite(offset_attr, H5T_NATIVE_INT, &offset_val);
    H5Sclose(offset_fspace);
    H5Aclose(offset_attr);
    
    H5Dclose(dataset);
}


void
HDF5File::
write_fs(const string& objectname,
         const vector<fp>& fs,
         const vector<li>& is,
         const vector<ii>& js) const
{
    cout << "Writing " << filename << "/" << objectname << "/" << 0 << ":" << js.size()-1 << endl;
    
    for (ii j = 0; j < js.size(); j++)
    {
        ii n = is[j+1] - is[j];
        
        hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(lcpl_id, 1);
        
        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        float fillval = -3.0/8.0;
        H5Pset_fill_value(dcpl_id, H5T_NATIVE_FLOAT, &fillval);
        
        hsize_t cdims[1] = {n < 16384 ? n : 16384};
        H5Pset_chunk(dcpl_id, 1, cdims);
        H5Pset_shuffle(dcpl_id);
        H5Pset_deflate(dcpl_id, 1);
        
        hsize_t dims[1] = {n};
        hid_t fspace = H5Screate_simple(1, dims, dims);
        ostringstream oss; oss << objectname << "/" << j;
        hid_t dataset = H5Dcreate(file, oss.str().c_str(), H5T_NATIVE_FLOAT, fspace, lcpl_id,  dcpl_id, H5P_DEFAULT);
        
        // write
        hsize_t mdims[1] = {n};
        hid_t mspace = H5Screate_simple(1, mdims, mdims);
        H5Dwrite(dataset, H5T_NATIVE_FLOAT, mspace, fspace, H5P_DEFAULT, &fs.data()[is[j]]);
        
        H5Sclose(fspace);
        H5Sclose(mspace);
        
        H5Dclose(dataset);
    }
}




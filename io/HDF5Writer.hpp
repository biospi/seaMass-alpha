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


#ifndef _SEAMASSRESTORATION_SMOWRITER_HPP_
#define _SEAMASSRESTORATION_SMOWRITER_HPP_


#include <hdf5.h>
#include <iostream>
#include "../core/seaMass.hpp"

class HDF5Writer
{
protected:
    std::string filename;
    int file;
    template<typename T>
    void write(const std::string& _objectname, const std::vector<T>& _cdata, hid_t& _data_type_id) const;
    
public:
	HDF5Writer(const std::string& filename);
	~HDF5Writer();
    
	void write_input(const SeaMass::Input& input) const;
	void write_output(const SeaMass::Output& output, ii shrinkage, ii tolerance, ii page_size) const;
	void write_output_control_points(const SeaMass::ControlPoints& control_points) const;
   
	void write(const std::string& objectname, const std::vector<unsigned char>& cdata) const;
	void write(const std::string& objectname, const std::vector<float>& cdata) const;
	void write(const std::string& objectname, const std::vector<double>& cdata) const;
	void write(const std::string& objectname, const std::vector<long long>& cdata) const;
	void write(const std::string& objectname, const std::vector<long>& cdata) const;
};


template<typename T>
void
HDF5Writer::
write(const std::string& _objectname, const std::vector<T>& _cdata, hid_t& _data_type_id) const
{
    hsize_t n = _cdata.size();

    hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl_id, 1);

    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    float fillval = -3.0/8.0;
    H5Pset_fill_value(dcpl_id, _data_type_id, &fillval);

    hsize_t cdims[1] ={n < 16384 ? n : 16384};
    H5Pset_chunk(dcpl_id, 1, cdims);
    H5Pset_shuffle(dcpl_id);
    H5Pset_deflate(dcpl_id, 1);

    hsize_t dims[1] = {n};
    hid_t fspace = H5Screate_simple(1, dims, dims);
    hid_t dataset = H5Dcreate(file, _objectname.c_str(), _data_type_id, fspace, lcpl_id,  dcpl_id, H5P_DEFAULT);

    // write
    hsize_t mdims[1] = {n};
    hid_t mspace = H5Screate_simple(1, mdims, mdims);
    H5Dwrite(dataset, _data_type_id, mspace, fspace, H5P_DEFAULT, &_cdata[0]);

    H5Sclose(fspace);
    H5Sclose(mspace);

   if(H5Dclose(dataset)<0) std::cout << "ARGH" << std::endl;
}


#endif // _SEAMASSRESTORATION_SMOWRITER_HPP_


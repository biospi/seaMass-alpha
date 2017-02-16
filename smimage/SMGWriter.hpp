//
// $Id$
//
//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Bristol, UK
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

#ifndef SMIMAGE_SMGWRITER_HPP_
#define SMIMAGE_SMGWRITER_HPP_

#include <iostream>
#include <vector>
#include <H5Cpp.h>

using namespace std;

template<typename T>
bool read_VecH5(string filename, string datagroup, vector<T>& vec_data,
				const H5::DataType &data_type_id);

template<typename T>
bool read_AttH5(string filename, string dataSetName, string dataAttName,T& attribute,
				const H5::DataType &data_type_id);

class SMGWriter
{
protected:
	string filename;
	H5::H5File *h5file;

public:
	SMGWriter(string _filename);
	~SMGWriter();

	template<typename T>
	void write_VecMatH5(string group, vector<T> const &data_set,
					const vector<hsize_t> dims,
					const H5::DataType &data_type_id);

	template<typename T>
	void write_VecMatH5(string group, vector<vector<T> > const &data_set,
					const H5::DataType &data_type_id);
};


template<typename T>
bool read_VecH5(string filename, string datagroup, vector<T>& vec_data,
				const H5::DataType &data_type_id)
{
	cout<<"Loading DATA: "<< datagroup <<" from file: "<<filename <<endl;
	// Try block to detect exceptions raised by any of the calls inside it
	try
	{
		// Turn off the auto-printing when failure occurs so that we can
		// handle the errors appropriately
		H5::Exception::dontPrint();

		// Open an existing file and dataset.
		H5::H5File file(filename, H5F_ACC_RDONLY);
		H5::DataSet dataset = file.openDataSet(datagroup);

		// Get dataspace of the dataset.
		H5::DataSpace dataspace = dataset.getSpace();

		// Get the number of dimensions in the dataspace.
		int rank = dataspace.getSimpleExtentNdims();

		// Get the dimension size of each dimension in the dataspace and
		// display them.
		hsize_t dims_out[rank];
		int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
		cout << "rank " << rank << ", dimensions " << (unsigned long)(dims_out[0]) << endl;

		vec_data.resize(dims_out[0],0);
		// Read the data to the dataset using default memory space, file
		// space, and transfer properties.
		dataset.read(&vec_data[0], data_type_id);

		file.close();

	}  // end of try block

	// catch failure caused by the H5File operations
	catch(H5::FileIException& error)
	{
		cout<<"ERROR HDF5 FILE"<<endl;
		error.printError();
		return false;
	}

	// catch failure caused by the DataSet operations
	catch(H5::DataSetIException& error)
	{
		cout<<"ERROR HDF5 DATA"<<endl;
		error.printError();
		return false;
	}

	// catch failure caused by the DataSpace operations
	catch(H5::DataSpaceIException& error)
	{
		error.printError();
		return false;
	}

	// catch failure caused by the Attribute operations
	catch(H5::AttributeIException& error)
	{
		error.printError();
		return false;
	}
	return true;
}

template<typename T>
bool read_AttH5(string filename, string dataSetName, string dataAttName, T &attribute,
				const H5::DataType &data_type_id)
{
	cout<<"Loading Attribute from Dataset: "<< dataSetName <<" from file: "<<filename <<endl;
	// Try block to detect exceptions raised by any of the calls inside it
	try
	{
		// Turn off the auto-printing when failure occurs so that we can
		// handle the errors appropriately
		H5::Exception::dontPrint();

		// Open an existing file and dataset.
		H5::H5File file(filename, H5F_ACC_RDONLY);
		H5::DataSet dataset = file.openDataSet(dataSetName);

	    H5::Attribute att = dataset.openAttribute("instrumentType");
	    att.read(data_type_id, &attribute);

		file.close();

	}  // end of try block

	// catch failure caused by the H5File operations
	catch(H5::FileIException& error)
	{
		cout<<"ERROR HDF5 FILE"<<endl;
		error.printError();
		return false;
	}

	// catch failure caused by the DataSet operations
	catch(H5::DataSetIException& error)
	{
		cout<<"ERROR HDF5 DATA"<<endl;
		error.printError();
		return false;
	}

	// catch failure caused by the DataSpace operations
	catch(H5::DataSpaceIException& error)
	{
		error.printError();
		return false;
	}

	// catch failure caused by the Attribute operations
	catch(H5::AttributeIException& error)
	{
		error.printError();
		return false;
	}
	return true;
}

template<typename T>
void SMGWriter::write_VecMatH5(string group, vector<T> const &data_set,
					const vector<hsize_t> dims,
					const H5::DataType &data_type_id)
{
	// Write out both Vector or Matrix data in HDF5 file, with an packed vector
	// data input.

	int rank=dims.size();

	try{
		H5::Group h5group;
		H5::DataSet h5dataset;

		vector<hsize_t> cdims;

		for (int i=0; i < rank; ++i)
		{
			cdims.push_back({dims[i] < 16384 ? dims[i] : 16384});
		}

		H5::DSetCreatPropList dataplist;
	    dataplist.setChunk(rank, &cdims[0]);
	    dataplist.setDeflate(7);

		//string data_group = "/RootGroup/"+group;
		//h5file->createGroup("/RootGroup");

		H5::DataSpace dataspace(rank, &dims[0]);
		H5::DataSet *dataset = new H5::DataSet(
				h5file->createDataSet(group.c_str(),data_type_id,dataspace,dataplist));
		dataset->write(&data_set[0], data_type_id);

		dataset->close();
		delete dataset;
	}
	catch(const H5::FileIException & error){
		error.printError();
	}
	catch(const H5::GroupIException & error)
	{
		error.printError();
	}
	catch(const H5::DataSetIException & error)
	{
		error.printError();
	}
}

template<typename T>
void SMGWriter::write_VecMatH5(string group, vector<vector<T> > const &data_set,
					const H5::DataType &data_type_id)
{
	const int rank=2;
	hsize_t nrow = data_set.size();
	//hsize_t min_ncol=data_set[0].size();
	hsize_t max_ncol=0;
	hsize_t dim[2];
	// Seamass HDF5 file containing data for Image
	typename vector<vector<T> >::const_iterator itrow;

	for(itrow = data_set.begin(); itrow != data_set.end(); ++itrow )
	{
		hsize_t vec_size = itrow->size();
		if(vec_size > max_ncol) max_ncol=vec_size;
		//if(vec_size < min_ncol) min_ncol=vec_size;
	}
	//cout<<"Length of row vector [Min, Max]       : ["<<min_ncol<<","<<max_ncol<<"]"<<endl;
	cout<<"Dimensions of output Matrix [row,col]: ["<<nrow<<","<<max_ncol<<"]"<<endl;

	dim[0]=nrow;
	dim[1]=max_ncol;

	vector<T> matvec(nrow*max_ncol, -1.0);
	for(hsize_t i=0; i < data_set.size(); ++i)
		for(hsize_t j=0; j < data_set[i].size(); ++j)
			matvec[max_ncol*i+j]=data_set[i][j];

	try{
		H5::Group h5group;
		H5::DataSet h5dataset;
		string data_group = "/RootGroup/"+group;

		h5file->createGroup("/RootGroup");

		H5::DataSpace dataspace(rank, dim);
		H5::DataSet dataset = h5file->createDataSet(data_group.c_str(),data_type_id,dataspace);
		dataset.write(&matvec[0], data_type_id);

		dataset.close();
	}
	catch(const H5::FileIException & error){
		error.printError();
	}
	catch(const H5::GroupIException & error)
	{
		error.printError();
	}
	catch(const H5::DataSetIException & error)
	{
		error.printError();
	}
}

#endif /* SMIMAGE_SMGWRITER_HPP_ */

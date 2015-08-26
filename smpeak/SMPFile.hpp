#ifndef SMPEAK_SMPFILE_HPP_
#define SMPEAK_SMPFILE_HPP_

#include<iostream>
#include<vector>
#include<H5Cpp.h>
#include"peakcore.hpp"

using namespace std;

class ReadSMFile
{
private:
	string filename;
	H5::H5File *h5file;
	vector<string> dataSetList;
public:
	ReadSMFile(string _filename);
	~ReadSMFile();
	void searchGroup(const string group, const string dataSet);
	vector<string> getDataSetName(void);
	void open(string _file);
	void close(void);

	template<typename T>
	void read_VecH5(const string dataSetName, vector<T> &vec_data,
				const H5::DataType &data_type_id);
	template<typename T>
	void read_MatH5(const string dataSetName, vector<T> &mat_data,
				lli &row, lli &col, const H5::DataType &data_type_id);
	template<typename T>
	void read_AttH5(const string dataSetName, const string dataAttName, T& attribute,
					const H5::DataType &data_type_id);
};

class SMPFile
{
protected:
	string filename;
	H5::H5File *h5file;

public:
	SMPFile(string _filename);
	~SMPFile();

	template<typename T>
	void write_VecMatH5(string group, vector<T> const &data_set,
					const vector<hsize_t> dims,
					const H5::DataType &data_type_id);

	template<typename T>
	void write_VecMatH5(string group, vector<vector<T> > const &data_set,
					const H5::DataType &data_type_id);
};


template<typename T>
void ReadSMFile::read_VecH5(const string dataSetName, vector<T> &vec_data,
				const H5::DataType &data_type_id)
{
	cout<<"Loading DATA: "<< dataSetName <<" from file: "<<filename <<endl;
	// Try block to detect exceptions raised by any of the calls inside it
	try
	{
		// Open an existing dataset.
		H5::DataSet dataset = h5file->openDataSet(dataSetName);

		// Get dataspace of the dataset.
		H5::DataSpace dataspace = dataset.getSpace();

		// Get the number of dimensions in the dataspace.
		int rank = dataspace.getSimpleExtentNdims();

		// Get the dimension size of each dimension in the dataspace and
		// display them.
		hsize_t dims_out[rank];
		int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
		cout << "Rank " << rank << ", Dimensions: " << (unsigned long)(dims_out[0]) << endl;

		vec_data.resize(dims_out[0],0);
		// Read the data to the dataset using default memory space, file
		// space, and transfer properties.
		dataset.read(&vec_data[0], data_type_id);
	}
	// catch failure caused by the H5File operations
	catch(H5::FileIException& error)
	{
		cout<<"ERROR HDF5 FILE"<<endl;
		error.printError();
	}
	// catch failure caused by the DataSet operations
	catch(H5::DataSetIException& error)
	{
		cout<<"ERROR HDF5 DATA"<<endl;
		error.printError();
	}
	// catch failure caused by the DataSpace operations
	catch(H5::DataSpaceIException& error)
	{
		error.printError();
	}
	// catch failure caused by the Attribute operations
	catch(H5::AttributeIException& error)
	{
		error.printError();
	}
}


template<typename T>
void ReadSMFile::read_MatH5(const string dataSetName, vector<T> &mat_data,
				lli &row, lli &col, const H5::DataType &data_type_id)
//void ReadSMFile::read_MatH5(const string dataGroup, vector<vector<T> > &mat_data,
//				const H5::DataType &data_type_id)
{
	cout<<"Loading DATA: "<< dataSetName <<" from file: "<<filename <<endl;
	// Try block to detect exceptions raised by any of the calls inside it
	try
	{
		// Open an existing dataset.
		H5::DataSet dataset = h5file->openDataSet(dataSetName);

		// Get dataspace of the dataset.
		H5::DataSpace dataspace = dataset.getSpace();

		// Get the number of dimensions in the dataspace.
		int rank = dataspace.getSimpleExtentNdims();

		// Get the dimension size of each dimension in the dataspace and
		// display them.
		hsize_t dims_out[rank];
		//vector<hsize_t> dims_out(rank);
		int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
		cout << "Rank " << rank << ", Dimensions Row: " << (dims_out[0]) << endl;
		cout << "Rank " << rank << ", Dimensions Col: " << (dims_out[1]) << endl;

		//H5::DataSpace memspace(rank,dims_out, NULL);

		row=dims_out[0];
		col=dims_out[1];
		mat_data.resize(row*col);

		cout<<"Matrix Row: "<<row<<endl;
		cout<<"Matrix Col: "<<col<<endl;

		// Read the data to the dataset using default memory space, file
		// space, and transfer properties.
		dataset.read(&mat_data[0], data_type_id);
	}
	// catch failure caused by the H5File operations
	catch(H5::FileIException& error)
	{
		cout<<"ERROR HDF5 FILE"<<endl;
		error.printError();
	}
	// catch failure caused by the DataSet operations
	catch(H5::DataSetIException& error)
	{
		cout<<"ERROR HDF5 DATA"<<endl;
		error.printError();
	}
	// catch failure caused by the DataSpace operations
	catch(H5::DataSpaceIException& error)
	{
		error.printError();
	}
	// catch failure caused by the Attribute operations
	catch(H5::AttributeIException& error)
	{
		error.printError();
	}
}


template<typename T>
void ReadSMFile::read_AttH5(const string dataSetName, const string dataAttName, T &attribute,
				const H5::DataType &data_type_id)
{
	cout<<"Loading Attribute from DataSet: "<< dataSetName <<" from file: "<<filename <<endl;
	// Try block to detect exceptions raised by any of the calls inside it
	try
	{
		H5::DataSet dataset = h5file->openDataSet(dataSetName);

	    H5::Attribute att = dataset.openAttribute(dataAttName);
	    att.read(data_type_id, &attribute);

	}  // end of try block

	// catch failure caused by the H5File operations
	catch(H5::FileIException& error)
	{
		cout<<"ERROR HDF5 FILE"<<endl;
		error.printError();
	}

	// catch failure caused by the DataSet operations
	catch(H5::DataSetIException& error)
	{
		cout<<"ERROR HDF5 DATA"<<endl;
		error.printError();
	}

	// catch failure caused by the DataSpace operations
	catch(H5::DataSpaceIException& error)
	{
		error.printError();
	}

	// catch failure caused by the Attribute operations
	catch(H5::AttributeIException& error)
	{
		error.printError();
	}
}


template<typename T>
void SMPFile::write_VecMatH5(string group, vector<T> const &data_set,
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
void SMPFile::write_VecMatH5(string group, vector<vector<T> > const &data_set,
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

#endif /* SMPEAK_SMPFILE_HPP_ */

/*
 * SMGWriter.cpp
 *
 *      Author: Ranjeet Singh Bhamber
 */

#include "SMGWriter.hpp"

SMGWriter::SMGWriter(string _filename)
{
	filename=_filename+".smg";
	try{
		H5::Exception::dontPrint();
		h5file = new H5::H5File(filename.c_str(), H5F_ACC_TRUNC);
	}
	catch(const H5::FileIException & error){
		error.printError();
	}
}

SMGWriter::~SMGWriter()
{
	h5file->close();
	delete h5file;
}

/*
 * SMPFile.cpp
 *
 *      Author: Ranjeet Bhamber
 */


#include "SMPFile.hpp"

SMPFile::SMPFile(string _filename)
{
	filename=_filename+".smp";
	try{
		H5::Exception::dontPrint();
		h5file = new H5::H5File(filename.c_str(), H5F_ACC_TRUNC);
	}
	catch(const H5::FileIException & error){
		error.printError();
	}
}

SMPFile::~SMPFile()
{
	h5file->close();
	delete h5file;
}

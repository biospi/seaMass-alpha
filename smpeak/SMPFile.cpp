//
// $Id$
//
//
// Author: Ranjeet Bhamber <ranjeet <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Liverpool, UK
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

#include "SMPFile.hpp"

ReadSMFile::ReadSMFile(string _filename):filename(_filename)
{
	try{
		H5::Exception::dontPrint();
		h5file = new H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	}
	catch(const H5::FileIException & error){
		error.printError();
	}

}

void ReadSMFile::open(string _file)
{
	try{
		h5file->openFile(filename.c_str(), H5F_ACC_RDONLY);
	}
	catch(const H5::FileIException & error){
		error.printError();
	}
}

ReadSMFile::~ReadSMFile()
{
	h5file->close();
}

vector<string> ReadSMFile::getDataSetName(void)
{
	return dataSetList;
}

void ReadSMFile::close(void)
{
	h5file->close();
}

void ReadSMFile::searchGroup(const string group, const string setName)
{
	H5::Group dataGroup = h5file->openGroup(group);
	hsize_t nObj = dataGroup.getNumObjs();
	for (hsize_t i = 0; i < nObj; ++i) {
		string obj = dataGroup.getObjnameByIdx(i);
		int oType = dataGroup.getObjTypeByIdx(i);
		//cout << "Object Name: " << obj << endl;
		//cout << "Object Type: " << oType << endl;
		string newObj;
		newObj=group+obj+"/";
		switch(oType)
		{
			case H5G_GROUP:
				//cout<<"Group NAME: "<<newObj<<endl;
				searchGroup(newObj, setName);
				break;
			case H5G_DATASET:
				//cout<<"Group NAME: "<<newObj<<endl;
				//cout<<"DataSet NAME: "<<obj<<endl;
				int n = newObj.find(setName);
				if(n > 0)
					dataSetList.push_back(group+obj);
				break;
		}
	}
}

double ReadSMFile::searchGroup(const string group, int level)
{
	double val;
	string newObj(group);
	string obj;
	H5::Group dataGroup = h5file->openGroup(group);
	hsize_t nObj = dataGroup.getNumObjs();
	if(nObj > 1)
	{
		cout<<"Error Specify Single Group!"<<endl;
		return NULL;
	}
	for (hsize_t i = 0; i <= level; ++i) {
		obj = dataGroup.getObjnameByIdx(0);
		int oType = dataGroup.getObjTypeByIdx(0);
		//cout << "Object Name: " << obj << endl;
		//cout << "Object Type: " << oType << endl;
		switch(oType)
		{
			case H5G_GROUP:
				//cout<<"Group NAME: "<<newObj<<endl;
				newObj=newObj+obj+"/";
				dataGroup = h5file->openGroup(newObj);
				nObj = dataGroup.getNumObjs();
				break;
			case H5G_DATASET:
				//cout<<"Group NAME: "<<newObj<<endl;
				cout<<"Error! Level: "<<level<<" too deep!"<<endl;
				return NULL;
				break;
		}
	}
	istringstream(obj)>>val;
	return val;
}


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

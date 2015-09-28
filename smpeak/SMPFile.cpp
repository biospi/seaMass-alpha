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

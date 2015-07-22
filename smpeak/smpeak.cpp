/*
 * smpeak.cpp
 *
 *      Author: Ranjeet Bhamber
 */
#include<iostream>
#include"peakcore.hpp"
#include"SMPFile.hpp"

void scanGroup(const H5::H5File& file,
				const string group, const string dataSet,
				vector<string> &setMatch)
{
	H5::Group dataGroup = file.openGroup(group);
	hsize_t nodj = dataGroup.getNumObjs();
	for (hsize_t i = 0; i < nodj; ++i) {
		string obj = dataGroup.getObjnameByIdx(i);
		int otype = dataGroup.getObjTypeByIdx(i);
		cout << "Object Name: " << obj << endl;
		cout << "Object Type: " << otype << endl;

		string newObj;
		newObj=group+obj+"/";
		switch(otype)
		{
			case H5G_GROUP:
				cout<<"Group NAME: "<<newObj<<endl;
				scanGroup(file, newObj, dataSet, setMatch);
				break;
			case H5G_DATASET:
				cout<<"Group NAME: "<<newObj<<endl;
				cout<<"DataSet NAME: "<<obj<<endl;
				int n = newObj.find(dataSet);
				if(n > 0)
					setMatch.push_back(newObj+obj);
				break;
		}
	}
}

int main(int argc, char **argv)
{

	if(argc != 2)
	{
		cout<<"Usage:"<<endl;
		cout<<"   smpeak <smo File>"<<endl;
		return -1;
	}

	int i=1;
	string filename=argv[i];

	cout<<"Peak Detection"<<endl;
	try
	{
		H5::Exception::dontPrint();

		// Open an existing file and Group.
		H5::H5File file(filename.c_str(), H5F_ACC_RDONLY);

		cout << "List all groups within file: " << filename << endl;

		vector<string> dataSetMatch;
		scanGroup(file,"/","fs",dataSetMatch);

		cout<<"Matching Dataset Found: "<<dataSetMatch[0]<<endl;
		file.close();

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

	return 0;
}

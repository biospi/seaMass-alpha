/*
 * smpeak.cpp
 *
 *      Author: Ranjeet Bhamber
 */
#include<iostream>
#include"peakcore.hpp"
#include"SMPFile.hpp"


int main(int argc, char **argv)
{

	if(argc != 2)
	{
		cout<<"Usage:"<<endl;
		cout<<"   smpeak <smo File>"<<endl;
		return -1;
	}

	cout<<"Peak Detection"<<endl;
	int i=1;
	string filename=argv[i];
	ReadSMFile dataFile(filename);
	vector<double> vecData;

	cout << "List all groups within file: " << filename << endl;
	vector<string> groupList;

	dataFile.searchGroup("/","/cs");
	dataFile.searchGroup("/","/fs");
	dataFile.searchGroup("/","/fcs");
	dataFile.searchGroup("/","/SpectrumMZ");
	groupList = dataFile.getDataSetName();

	for(int i=0; i < groupList.size(); ++i)
	{
		cout<<"DataSets found ["<< i<<"]: "<<groupList[i]<<endl;
	}

	vector<float> vecMat;
	lli row,col;
	dataFile.read_VecH5(groupList[3], vecData, H5::PredType::NATIVE_DOUBLE);
	dataFile.read_MatH5(groupList[0], vecMat, row, col, H5::PredType::NATIVE_FLOAT);

	cout<<"Size of Vector: "<<vecMat.size()<<endl;

	vector<float*> idx;
	idx.reserve(row);
	for(lli i=0; i < row; ++i)
		idx.push_back(&vecMat[i*col]);

	float **matData = &idx[0];

	for(lli i=0; i<row; ++i){
		for(lli j=0; j<col; ++j)
			cout<<matData[i][j]<<"  ";
		cout<<endl;
	}

	//cout<<"Size of Matrix: ["<<matData.size()<<","<<matData[0].size()<<"]"<<endl;

	return 0;
}

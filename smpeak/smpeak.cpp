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
	string fileName=argv[i];
	ReadSMFile dataFile(fileName);
	string outFileName=fileName.substr(0,fileName.size()-4);

	cout << "List all groups within file: " << fileName << endl;
	vector<string> dataSetList;

	dataFile.searchGroup("/","/cs");
	dataFile.searchGroup("/","/fs");
	dataFile.searchGroup("/","/fcs");
	dataFile.searchGroup("/","/SpectrumMZ");
	dataSetList = dataFile.getDataSetName();

	for(int i=0; i < dataSetList.size(); ++i)
		cout<<"DataSets found ["<< i<<"]: "<<dataSetList[i]<<endl;

	vector<float> csVecMat;
	int offset[2];
	hsize_t row,col;
	vector<hsize_t> dim;
	dataFile.read_MatH5(dataSetList[0], csVecMat, row, col, H5::PredType::NATIVE_FLOAT);
	dataFile.read_AttH5(dataSetList[0],"Offset", offset, H5::PredType::NATIVE_INT);
	dim.push_back(row);
	dim.push_back(col);

	// Calculate the MZs of control points.
	vector<double> mza(col,0.0);
	calMZalpha(mza, offset[0]);

	// Calculate the RTs of control points.
	vector<double> rt(row,0.0);
	calRTalpha(rt, offset[1]);

	cout<<"Size of Vector: "<<csVecMat.size()<<endl;
	vector<float*> csIdx;
	csIdx.reserve(row);
	for(lli i=0; i < row; ++i)
		csIdx.push_back(&csVecMat[i*col]);
	float **csMat = &csIdx[0];

	for(lli i=0; i<6; ++i){
		for(lli j=0; j<8; ++j)
			cout<<csMat[i][j]<<"  ";
		cout<<endl;
	}

	vector<vector<float> > dcs=nabla(csMat, row, col);
	cout<<"\nFirst Derivative of Coeffs"<<endl;
	for(lli i=0; i<6; ++i){
		for(lli j=0; j<8; ++j)
			cout<<dcs[i][j]<<"  ";
		cout<<endl;
	}

	vector<vector<float> > d2cs=nabla2(csMat, row, col);
	cout<<"\nSecond Derivative of Coeffs"<<endl;
	for(lli i=0; i<6; ++i){
		for(lli j=0; j<8; ++j)
			cout<<d2cs[i][j]<<"  ";
		cout<<endl;
	}


	// Write data to SMP file.
	SMPFile smpDataFile(outFileName);

	smpDataFile.write_VecMatH5("csOrig",csVecMat,dim,H5::PredType::NATIVE_FLOAT);
	smpDataFile.write_MatH5("dcs",dcs,H5::PredType::NATIVE_FLOAT);
	smpDataFile.write_MatH5("d2cs",d2cs,H5::PredType::NATIVE_FLOAT);

	return 0;
}

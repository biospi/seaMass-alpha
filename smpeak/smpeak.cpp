#include<iostream>
#include<cmath>
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
	vector<string> dataSetList;

	dataFile.searchGroup("/","/cs");
	dataFile.searchGroup("/","/fs");
	dataFile.searchGroup("/","/fcs");
	dataFile.searchGroup("/","/SpectrumMZ");
	dataSetList = dataFile.getDataSetName();

	for(int i=0; i < dataSetList.size(); ++i)
	{
		cout<<"DataSets found ["<< i<<"]: "<<dataSetList[i]<<endl;
	}

	vector<float> csVecMat;
	int offset[2];
	lli row,col;
	dataFile.read_VecH5(dataSetList[3], vecData, H5::PredType::NATIVE_DOUBLE);
	dataFile.read_MatH5(dataSetList[1], csVecMat, row, col, H5::PredType::NATIVE_FLOAT);
	dataFile.read_AttH5(dataSetList[1],"Offset", offset, H5::PredType::NATIVE_INT);

	// Calculate the mz of control points.
	vector<double> mza(col,0.0);
	for(lli i=0; i < mza.size(); ++i)
	{
		double mz_rez=1.0;
		double ppbmz=1.0033548378/(pow(2.0,mz_rez)*60.0);
		mza[i]=double((offset[0]-1+i))*ppbmz;
	}

	// Calculate the rt of control points.
	vector<double> rt(row,0.0);
	for(lli i = 0; i < rt.size(); ++i)
	{
		double rt_rez=4.0;
		double ppbrt=1.0/(pow(2.0,rt_rez));
		rt[i]=double(offset[1]+i)*ppbrt;
	}

	cout<<"Size of Vector: "<<csVecMat.size()<<endl;

	vector<float*> csIdx;
	csIdx.reserve(row);
	for(lli i=0; i < row; ++i)
		csIdx.push_back(&csVecMat[i*col]);

	float **csMat = &csIdx[0];

	for(lli i=0; i<3; ++i){
		for(lli j=0; j<col; ++j)
			cout<<csMat[i][j]<<"  ";
		cout<<endl;
	}

	vector<vector <float> > dcs=nabla(csMat, row, col);

	cout<<"\nDerivative of Coeffs"<<endl;

	for(lli i=0; i<3; ++i){
		for(lli j=0; j<col; ++j)
			cout<<dcs[i][j]<<"  ";
		cout<<endl;
	}


	return 0;
}

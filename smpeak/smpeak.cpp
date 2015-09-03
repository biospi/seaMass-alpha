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
	calMZalpha(mza, offset[0],1.0);

	// Calculate the RTs of control points.
	vector<double> rta(row,0.0);
	calRTalpha(rta, offset[1],4.0);

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


	PeakData centriodPeak;

	int falsePeak=0;
	// Find Peaks and exact MZ values.
	cout<<"Find Peaks along MZ axis"<<endl;
	for(lli i = 0; i < row; ++i)
	{
		for(lli j = 2; j < col-2; ++j)
		{
			if((dcs[i][j] > 0) && (dcs[i][j+1] < 0) )
			{
				double pa1=0.0;
				double pmz1=0.0;
				calMidPoint(i,j,dcs,mza,pmz1,pa1);
				if(pa1 < 0){
					double pa0=0.0;
					double pmz0=0.0;
					calMidPoint(i,j-1,dcs,mza,pmz0,pa0);
					double t0 = calT(pa0,double(dcs[i][j]),pa1);
					double mzPeak=calX(t0,pmz0,mza[j],pmz1);
					centriodPeak.add_peak(mzPeak,rta[i],csMat[i][j],i,j);
				}
				else
				{
					double pa2=0.0;
					double pmz2=0.0;
					calMidPoint(i,j+1,dcs,mza,pmz2,pa2);
					double t0 = calT(pa1,double(dcs[i][j+1]),pa2);
					if (t0>=0)
					{
						double mzPeak=calX(t0,pmz1,mza[j+1],pmz2);
						centriodPeak.add_peak(mzPeak,rta[i],csMat[i][j],i,j);
					}
					else
					{
						++falsePeak;
					}
				}
			}
		}
	}

	cout<<"Found ["<<centriodPeak.mz.size()<<"] Peaks."<<endl;
	if(falsePeak > 0)
		cout<<"Warning !!! Found ["<<falsePeak<<"] Insignificant False Peak Detected - Peaks Ignored"<<endl;
	vector<hsize_t> vecN;
	vecN.push_back(0.0);

	cout<<"\nSaving Data to File:"<<endl;
	vecN[0]=centriodPeak.mz.size();
	smpDataFile.write_VecMatH5("Peak_mz",centriodPeak.mz,vecN,H5::PredType::NATIVE_DOUBLE);
	vecN[0]=centriodPeak.rt.size();
	smpDataFile.write_VecMatH5("Peak_rt",centriodPeak.rt,vecN,H5::PredType::NATIVE_DOUBLE);
	vecN[0]=centriodPeak.pVal.size();
	smpDataFile.write_VecMatH5("Peak_Count",centriodPeak.pVal,vecN,H5::PredType::NATIVE_FLOAT);
	vecN[0]=centriodPeak.mz_idx.size();
	smpDataFile.write_VecMatH5("Peak_mz_idx",centriodPeak.mz_idx,vecN,H5::PredType::NATIVE_LLONG);
	vecN[0]=centriodPeak.rt_idx.size();
	smpDataFile.write_VecMatH5("Peak_rt_idx",centriodPeak.rt_idx,vecN,H5::PredType::NATIVE_LLONG);

	return 0;
}

#include <iostream>
#include <iterator>
#include <boost/program_options.hpp>
#include "peakcore.hpp"
#include "SMPFile.hpp"
#include "SMData.hpp"
#include "MathOperator.hpp"
#include "BsplineData.hpp"
#include "PeakOperator.hpp"
#include "PeakData.hpp"
#include "PeakManager.hpp"
#include "NetCDFile.hpp"


namespace po = boost::program_options;

int main(int argc, char **argv)
{
	string fileName;
	string mzMLFileName;
	string outMZFileName;
	bool centroid;

	po::options_description desc("Usage: smpeak [OPTION...] [SMO FILE] [mzMLb3 FILE]");

	desc.add_options()
		("help,h", "Produce help message")
		("smo,s", po::value<string>(&fileName), "Filename of input smo data file")
		("mzMLb3,z", po::value<string>(&mzMLFileName), "Filename of input mzMLb3 data file")
		("output,o", po::value<string>(&outMZFileName), "Filename of output peak data file")
		("centroid,c", "Centroid mode Peak transformation along M/Z");

	try
	{
		po::positional_options_description pod;
		pod.add("smo", 1);
		pod.add("mzMLb3", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
		po::notify(vm);

		if(vm.count("help"))
		{
			cout<<desc<<endl;
			return 0;
		}
		if(vm.count("centroid"))
		{
			cout<<"Transforming data to peak detection Centroid mode."<<endl;
			centroid=true;
		}
		else
		{
			centroid=false;
		}
		if(vm.count("smo"))
		{
			cout<<"Opening SMO file: "<<vm["smo"].as<string>()<<endl;
		}
		else
		{
			cout<<desc<<endl;
			throw "SMO file was not give...";
		}
		if(vm.count("mzMLb3"))
		{
			cout<<"Opening mzMLb3 file: "<<vm["mzMLb3"].as<string>()<<endl;
		}
		else
		{
			cout<<desc<<endl;
			throw "mzMLb3 file was not give...";
		}
		if(vm.count("output"))
		{
			cout<<"Output Peak file: "<<vm["output"].as<string>()<<endl;
		}
		else
		{
			outMZFileName=fileName.substr(0,fileName.size()-4)+"_Peak.mzMLb3";
		}
	}
	catch(exception& e)
	{
		cerr<<"error: " << e.what() <<endl;
		return 1;
	}
	catch(const char* msg)
	{
		cerr<<"error: "<<msg<<endl;
		return 1;
	}
	catch(...)
		{
		cerr<<"Exception of unknown type!\n";
	}

	// Test Read data
	NetCDFile mzMLFile;

	mzMLFile.open(mzMLFileName);

	VecMat<double> mzData;
	VecMat<double> mzDataT;
	vector<char> mzMLbuff;

	mzMLFile.read_MatNC("spectrum_MS_1000514",mzData);
	mzMLFile.read_MatNCT("spectrum_MS_1000514",mzDataT);
	mzMLFile.read_VecNC("mzML",mzMLbuff);

	// Test Write data
	NetCDFile hammerNCDF4(outMZFileName,NC_NETCDF4);

	hammerNCDF4.write_VecNC("mzML",mzMLbuff,NC_CHAR);
	hammerNCDF4.write_MatNC("spectrum_MS_1000514",mzData,NC_DOUBLE);

	string buff(&mzMLbuff[0]);
	size_t loc=0;
	loc = buff.find("<spectrum index");
	loc = buff.find("<spectrum index", 586900);
	loc = buff.find("<spectrum index", 590943);

	cout<<"Spectrum Index: "<<loc<<endl;

	vector<unsigned int> specIdx;

	findVecString(mzMLbuff, specIdx);

	//cout<<endl<<buff<<endl;
	/*
	for(int i = 0; i < mzMLbuff.size(); ++i)
	{
		cout<<mzMLbuff[i];
	}
	*/

	cout<<"Centroid Peak Data Set..."<<endl;
	ReadSMFile dataFile(fileName);
	string outFileName=fileName.substr(0,fileName.size()-4);

	cout << "List all groups within file: " << fileName << endl;
	vector<string> dataSetList;

	dataFile.searchGroup("/","/cs");
	dataSetList = dataFile.getDataSetName();

	for(int i=0; i < dataSetList.size(); ++i)
		cout<<"DataSets found ["<< i<<"]: "<<dataSetList[i]<<endl;

	cout<<"List Parameters from SMO:"<<endl;
	double mzRes=dataFile.searchGroup("/",1);
	cout<<"MZ_resolution: "<<mzRes<<endl;
	double rtRes=dataFile.searchGroup("/",2);
	cout<<"RT resolution: "<<rtRes<<endl;

	vector<float> rawCoeff;
	int offset[2];
	hsize_t row,col;
	dataFile.read_MatH5(dataSetList[0], rawCoeff, row, col, H5::PredType::NATIVE_FLOAT);
	dataFile.read_AttH5(dataSetList[0],"Offset", offset, H5::PredType::NATIVE_INT);

	// Write data to SMP file.
	SMPFile smpDataFile(outFileName);

	hsize_t dims[2];
	dims[0]=row;
	dims[1]=col;

	SMData<OpUnit> A(dims,offset,mzRes,rtRes,rawCoeff);
	SMData<OpNablaH> dhA(dims,offset,mzRes,rtRes,rawCoeff);
	SMData<OpNabla2H> d2hA(dims,offset,mzRes,rtRes,rawCoeff);
	SMData<OpNablaV> dvA(dims,offset,mzRes,rtRes,rawCoeff);
	SMData<OpNabla2V> d2vA(dims,offset,mzRes,rtRes,rawCoeff);

	BsplineData<> bsData(A,dhA,d2hA);
	BsplineData<> bsPeakData(A,dhA,d2hA,dvA,d2vA);

	PeakManager<PeakData,BsplineData,Centroid> centriodPeak(bsData);
	centriodPeak.execute();

	PeakManager<PeakData,BsplineData,ExtractPeak> extractPeak(bsPeakData);
	extractPeak.execute();

	cout<<"\nSaving Data to File:"<<endl;

	vector<hsize_t> vecN;
	vecN.push_back(0.0);
	vecN[0]=centriodPeak.peak->getMZ().size();
	smpDataFile.write_VecMatH5("Peak_mz",centriodPeak.peak->getMZ(),vecN,H5::PredType::NATIVE_DOUBLE);
	vecN[0]=centriodPeak.peak->getMZwidth().size();
	smpDataFile.write_VecMatH5("Peak_mz_width",centriodPeak.peak->getMZwidth(),vecN,H5::PredType::NATIVE_DOUBLE);
	vecN[0]=centriodPeak.peak->getRT().size();
	smpDataFile.write_VecMatH5("Peak_rt",centriodPeak.peak->getRT(),vecN,H5::PredType::NATIVE_DOUBLE);
	vecN[0]=centriodPeak.peak->getRTwidth().size();
	smpDataFile.write_VecMatH5("Peak_rt_width",centriodPeak.peak->getRTwidth(),vecN,H5::PredType::NATIVE_DOUBLE);
	vecN[0]=centriodPeak.peak->getPKcount().size();
	smpDataFile.write_VecMatH5("Peak_Count",centriodPeak.peak->getPKcount(),vecN,H5::PredType::NATIVE_FLOAT);
	vecN[0]=centriodPeak.peak->getMZIdx().size();
	smpDataFile.write_VecMatH5("Peak_mz_idx",centriodPeak.peak->getMZIdx(),vecN,H5::PredType::NATIVE_LLONG);
	vecN[0]=centriodPeak.peak->getRTIdx().size();
	smpDataFile.write_VecMatH5("Peak_rt_idx",centriodPeak.peak->getRTIdx(),vecN,H5::PredType::NATIVE_LLONG);

	smpDataFile.write_VecMatH5("csOrig",A.alpha->v,dims,H5::PredType::NATIVE_FLOAT);
	smpDataFile.write_VecMatH5("dcs",dhA.alpha->v,dims,H5::PredType::NATIVE_FLOAT);
	smpDataFile.write_VecMatH5("d2cs",d2hA.alpha->v,dims,H5::PredType::NATIVE_FLOAT);
	smpDataFile.write_VecMatH5("dvcs",dvA.alpha->v,dims,H5::PredType::NATIVE_FLOAT);
	smpDataFile.write_VecMatH5("d2vcs",d2vA.alpha->v,dims,H5::PredType::NATIVE_FLOAT);

	return 0;
}

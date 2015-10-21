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
#include <pugixml.hpp>
#include <sstream>

namespace po = boost::program_options;
namespace xml = pugi;

int main(int argc, char **argv)
{
	string smoFileName;
	string mzMLFileName;
	string outMZFileName;
	bool centroid, debug;

	po::options_description desc("Usage: smpeak [OPTION...] [SMO FILE] [mzMLb3 FILE]");

	desc.add_options()
		("help,h", "Produce help message")
		("centroid,c", "Centroid mode Peak transformation along M/Z")
		("debug,g", "Output debugging Peak data")
		("smo,s", po::value<string>(&smoFileName), "Filename of input smo data file")
		("mzMLb3,z", po::value<string>(&mzMLFileName), "Filename of input mzMLb3 data file")
		("output,o", po::value<string>(&outMZFileName), "Filename of output peak data file");

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
			outMZFileName=smoFileName.substr(0,smoFileName.size()-4)+"_Peak.mzMLb3";
		}
		if(vm.count("debug"))
		{
			debug = true;
		}
		else
		{
			debug = false;
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


	//---------------------------------------------------------------------
	// NetCDF4 Read data
	//---------------------------------------------------------------------
	cout<<"Centroid Peak Data Set..."<<endl;

	//----------------------------
	// Load SMO Data that we need.
	//----------------------------
	NetCDFile dataFile(smoFileName);

	VecMat<float> rawCoeff;
	vector<int> offset;
	hsize_t dims[2];
	vector<InfoGrpVar> dataSetList;
	string outFileName=smoFileName.substr(0,smoFileName.size()-4);

	cout << "List all groups within file: " << smoFileName << endl;

	dataFile.search_Group("fcs");
	dataSetList = dataFile.get_Info();

	for(int i=0; i < dataSetList.size(); ++i)
		cout<<"DataSets found ["<< i<<"]: "<<dataSetList[i].grpName<<"/"
		<<dataSetList[0].varName<<endl;

	cout<<"List Parameters from SMO:"<<endl;
	double mzRes=dataFile.search_Group<double>(2);
	cout<<"MZ_resolution: "<<mzRes<<endl;
	double rtRes=dataFile.search_Group<double>(3);
	cout<<"RT resolution: "<<rtRes<<endl;

	dataFile.read_MatNC(dataSetList[0].varName,rawCoeff,dataSetList[0].grpid);
	dataFile.read_AttNC("Offset",dataSetList[0].varid,offset,dataSetList[0].grpid);
	rawCoeff.getDims(dims);

	dataFile.close();

	//-------------------------------
	// Load mzMLb3 Data that we need.
	//-------------------------------

	vector<char> mzMLbuff;
	VecMat<double> rawMZ;
	VecMat<float> rawPK;

	dataFile.open(mzMLFileName);

	dataFile.read_VecNC("mzML",mzMLbuff);
	dataFile.read_MatNC("spectrum_MS_1000514",rawMZ);
	dataFile.read_MatNC("spectrum_MS_1000515",rawPK);

	vector<size_t> rawSize = findSizeT(rawMZ);

	//-------------------------------------------
	// XML scan and read Start Time and ms-level.
	//-------------------------------------------
	vector<int> msLevel;
	vector<pair<int,double> > rtRaw;

	xml::xml_document docmzML;
	size_t xmlSize=sizeof(char)*mzMLbuff.size();
	xml::xml_parse_result result = docmzML.load_buffer_inplace(&mzMLbuff[0],xmlSize);


	xml::xpath_node_set tools =
		docmzML.select_nodes("mzML/run/spectrumList/spectrum/cvParam[@name='ms level']");

	cout<<"XML: "<<tools.size()<<endl;;
	int idx=0;
	for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr, ++idx)
	{
		int ms;
		double rt;
		//xml::xpath_node n = *i;
		//cout<<"Spectrum: "<<
			//n.node().parent().attribute("index").value()<<"  "<<
			//n.node().attribute("value").value()<<endl;
		//cout<<"Spectrum: "<<
			//i->node().parent().attribute("index").value()<<"  "<<
			//i->node().attribute("value").value()<<"\t";
		//cout<<"Start Time: "<<
			//i->node().parent().child("scanList").child("scan").child("cvParam").attribute("value").value()<<endl;

		istringstream(itr->node().attribute("value").value()) >> ms;
		msLevel.push_back(ms);
		if(ms == 1)
		{
			istringstream(itr->node().parent().child("scanList").child("scan").child("cvParam").attribute("value").value()) >> rt;
			rtRaw.push_back(make_pair(idx,rt));
		}
	}


	//---------------------------------------------------------------------
	// Centroid Peak Data Set...
	//---------------------------------------------------------------------
	SMData<OpUnit> A(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
	SMData<OpNablaH> dhA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
	SMData<OpNabla2H> d2hA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
	//SMData<OpNablaV> dvA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
	//SMData<OpNabla2V> d2vA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
	rawCoeff.clear();

	for(size_t i = 0; i < A.rt.size(); ++i)
	{
		A.rt[i] = rtRaw[i].second;
		dhA.rt[i] = rtRaw[i].second;
		d2hA.rt[i] = rtRaw[i].second;
	}

	BsplineData<> bsData(A,dhA,d2hA);
	//BsplineData<> bsPeakData(A,dhA,d2hA,dvA,d2vA);

	PeakManager<PeakData,BsplineData,Centroid> centriodPeak(bsData);
	centriodPeak.execute();

	//PeakManager<PeakData,BsplineData,ExtractPeak> extractPeak(bsPeakData);
	//extractPeak.execute();

	//---------------------------------------------------------------------
	// XML Modify
	//---------------------------------------------------------------------
	VecMat<double> mzPeak;
	VecMat<float> pkPeak;
	vector<size_t> vecSize;

	centriodPeak.peak->getPeakMatT(mzPeak, pkPeak, A.rt.size(), vecSize);

	repackPeakDataT(mzPeak, rawMZ, msLevel, vecSize, rawSize);
	repackPeakDataT(pkPeak, rawPK, msLevel, vecSize, rawSize);

	tools =	docmzML.select_nodes("mzML/run/spectrumList/spectrum");

	cout<<"XML: "<<tools.size()<<endl;;
	for(int i = 0, j = 0; i < tools.size(); ++i)
	{
		if( (msLevel[i] == 1) && (j < vecSize.size()) )
		{
			ostringstream buff;
			buff << vecSize[j];
			string newVal(buff.str());
			tools[i].node().attribute("defaultArrayLength").set_value(newVal.c_str());
			++j;
		}
	}

	stringstream newmzML;
	docmzML.save(newmzML);
	string output = newmzML.str();

	// Free mzMLbuff
	vector<char>().swap(mzMLbuff);
	mzMLbuff.assign(output.begin(),output.end());

	// free up memory
	newmzML.clear();
	newmzML.str(std::string());

	//vector<char> vs(output.begin(),output.end());

	cout<<"It works!"<<endl;
	//---------------------------------------------------------------------
	//---------------------------------------------------------------------

	//---------------------------------------------------------------------
	// NetCDF4 Test Write data
	//---------------------------------------------------------------------
	/*
	NetCDFile hammerNCDF4(outMZFileName,NC_NETCDF4);

	//hammerNCDF4.write_VecNC("mzML",mzMLbuff,NC_BYTE);
	hammerNCDF4.write_VecNC("mzML",mzMLbuff,NC_BYTE,4096,0);
	hammerNCDF4.write_MatNC("spectrum_MS_1000514",mzData,NC_DOUBLE);

	vector<double> attVal;

	attVal.push_back(123.422);
	attVal.push_back(576.452);
	attVal.push_back(34.4542);
	attVal.push_back(983.232);

	double *p;
	p = &attVal[0];

	cout<<p[1]<<endl;

	hammerNCDF4.write_AttNC("mzML","MyTest",attVal, NC_DOUBLE);

	vector<unsigned int> specIdx;
	findVecString(mzMLbuff, specIdx);
	*/


	if(debug)
	{
		// Write data to SMP file.
		SMPFile smpDataFile(outFileName);

		mzMLdump(string(outFileName+"_mzML.txt"),output);

		cout<<"\nSaving Peak Debugging Data to File:"<<endl;

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
		//smpDataFile.write_VecMatH5("dvcs",dvA.alpha->v,dims,H5::PredType::NATIVE_FLOAT);
		//smpDataFile.write_VecMatH5("d2vcs",d2vA.alpha->v,dims,H5::PredType::NATIVE_FLOAT);

		mzPeak.getDims(dims);
		smpDataFile.write_VecMatH5("mat_Peak_MZ",mzPeak.v,dims,H5::PredType::NATIVE_DOUBLE);
		smpDataFile.write_VecMatH5("mat_Peak_PK",pkPeak.v,dims,H5::PredType::NATIVE_FLOAT);

		rawMZ.getDims(dims);
		smpDataFile.write_VecMatH5("mat_raw_Peak_MZ",rawMZ.v,dims,H5::PredType::NATIVE_DOUBLE);
		smpDataFile.write_VecMatH5("mat_raw_Peak_PK",rawPK.v,dims,H5::PredType::NATIVE_FLOAT);
	}

	return 0;
}




	//---------------------------------------------------------------------
	// HDF5 Load Peak Data Set...
	//---------------------------------------------------------------------
	/*cout<<"Centroid Peak Data Set..."<<endl;
	ReadSMFile dataFile(fileName);
	string outFileName=fileName.substr(0,fileName.size()-4);

	cout << "List all groups within file: " << fileName << endl;
	vector<string> dataSetList;

	dataFile.searchGroup("/","/fcs");
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

	hsize_t dims[2];
	dims[0]=row;
	dims[1]=col;
	*/


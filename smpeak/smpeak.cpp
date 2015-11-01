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
#include <omp.h>

namespace po = boost::program_options;
namespace xml = pugi;

int main(int argc, char **argv)
{
	string smoFileName;
	string mzMLFileName;
	string outMZFileName;
	bool centroid, debug;
	typedef pair<int,double> rtIdxData;

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

	//-------------------------------
	// Load mzMLb3 Data that we need.
	//-------------------------------
	vector<char> mzMLbuff;
	NetCDFile mzMLb3DF(mzMLFileName);
	mzMLb3DF.read_VecNC("mzML",mzMLbuff);

	//-------------------------------------------
	// XML scan and read Start Time and ms-level.
	//-------------------------------------------
	vector<int> msLevel;
	vector<rtIdxData> rtRaw;
	vector<size_t> rawSize;
	size_t maxSize=0;

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
		size_t len;

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

		istringstream(itr->node().parent().attribute("defaultArrayLength").value()) >> len;
		rawSize.push_back(len);
		if(rawSize.back() > maxSize) maxSize = rawSize.back();

		if(ms == 1)
		{
			istringstream(itr->node().parent().child("scanList").child("scan").child("cvParam").attribute("value").value()) >> rt;
			rtRaw.push_back(make_pair(idx,rt));
		}
	}

	//----------------------------
	// Load SMO Data that we need.
	//----------------------------
	NetCDFile smoDF(smoFileName);

	VecMat<float> rawCoeffDEL;

	vector<int> offset;
	hsize_t dimsDel[2];
	vector<hsize_t> dims(2,0);
	vector<size_t> fcsLen;
	vector<InfoGrpVar> dataSetList;

	cout << "List all groups within file: " << smoFileName << endl;

	smoDF.search_Group("fcs");
	dataSetList = smoDF.get_Info();

	for(int i=0; i < dataSetList.size(); ++i)
		cout<<"DataSets found ["<< i<<"]: "<<dataSetList[i].grpName<<"/"
		<<dataSetList[0].varName<<endl;

	cout<<"List Parameters from SMO:"<<endl;
	double mzRes=smoDF.search_Group<double>(2);
	cout<<"MZ_resolution: "<<mzRes<<endl;
	double rtRes=smoDF.search_Group<double>(3);
	cout<<"RT resolution: "<<rtRes<<endl;

	smoDF.read_AttNC("Offset",dataSetList[0].varid,offset,dataSetList[0].grpid);
	fcsLen = smoDF.read_DimNC(dataSetList[0].varName,dataSetList[0].grpid);

	//=====================================================================
	smoDF.read_MatNC(dataSetList[0].varName,rawCoeffDEL,dataSetList[0].grpid);
	rawCoeffDEL.getDims(dimsDel);

	//---------------------------------------------------------------------
	// NetCDF4 Read data TEST Hyperslab
	VecMat<float> hypTest;
	size_t rc[2]={7,4};
	size_t lenA[2]={6,1};

	smoDF.read_HypVecNC(dataSetList[0].varName,hypTest.v,rc,lenA,dataSetList[0].grpid);
	hypTest.clear();

	lenA[0]=1;
	lenA[1]=0;
	smoDF.read_HypVecNC(dataSetList[0].varName,hypTest.v,rc,lenA,dataSetList[0].grpid);
	hypTest.clear();

	lenA[0]=0;
	lenA[1]=1;
	smoDF.read_HypVecNC(dataSetList[0].varName,hypTest.v,rc,lenA,dataSetList[0].grpid);
	hypTest.clear();

	lenA[0]=5;
	lenA[1]=3;
	smoDF.read_HypMatNC(dataSetList[0].varName,hypTest,rc,lenA,dataSetList[0].grpid);
	hypTest.clear();

	rc[0]=3;
	rc[1]=1192;
	lenA[0]=3;
	lenA[1]=0;
	smoDF.read_HypMatNC(dataSetList[0].varName,hypTest,rc,lenA,dataSetList[0].grpid);
	hypTest.clear();
	//=====================================================================

	//---------------------------------------------------------------------
	// Centroid Peak Data Set...
	// Process Scan by Scan
	//---------------------------------------------------------------------
	vector<size_t> hypIdx(2);
	vector<size_t> rdLen(2);

	rdLen[0]=1;
	rdLen[1]=fcsLen[1];
	dims[0]=hsize_t(rdLen[0]);
	dims[1]=hsize_t(rdLen[1]);
	hypIdx[1]=0; // = Always read from first Column;

	PeakData<> totalPeaks;
	vector<PeakData<> *> peakThreads(omp_get_max_threads());

	cout<<"Extract Peaks from Mass Spec Data"<<endl;
	// Manual reduction of STL container as STLs are not thread safe...
	#pragma omp parallel
	{
		int nthrd = omp_get_num_threads();
		PeakData<> localPeaks;
		int thrdid=omp_get_thread_num();
		peakThreads[thrdid] = &localPeaks;

		#pragma omp for firstprivate(hypIdx,rdLen)
		for(size_t rt_idx = 0; rt_idx < fcsLen[0]; ++rt_idx)
		{
			vector<float> rawCoeff;
			hypIdx[0]=rt_idx;

			smoDF.read_HypVecNC(dataSetList[0].varName,rawCoeff,&hypIdx[0],&rdLen[0],
					dataSetList[0].grpid);

			SMData1D<OpUnit> A(&dims[0],&offset[0],mzRes,rtRaw[rt_idx],rawCoeff);
			SMData1D<OpNablaH> dhA(&dims[0],&offset[0],mzRes,rtRaw[rt_idx],rawCoeff);
			SMData1D<OpNabla2H> d2hA(&dims[0],&offset[0],mzRes,rtRaw[rt_idx],rawCoeff);

			BsplineData<rtIdxData> bsData(A,dhA,d2hA);

			PeakManager<PeakData,BsplineData,Centroid1D,rtIdxData> centriodPeak(bsData);
			centriodPeak.execute();

			localPeaks.addPeakArray(centriodPeak.peak->getPeakData());
		}

		#pragma omp single
		{
			for(int i = 0; i < nthrd; ++i)
			{
				cout<<"Thread ["<<i<<"] peaks found: "<<peakThreads[i]->numOfPeaks()<<endl;
				totalPeaks.addPeakArray(peakThreads[i]->getPeakData());
			}
		}
	}
	cout<<"Total ["<<totalPeaks.numOfPeaks()<<"] peaks found."<<endl;

	//=====================================================================
	SMData2D<OpUnit> Adel(dimsDel,&offset[0],mzRes,rtRes,rawCoeffDEL.v);
	SMData2D<OpNablaH> dhAdel(dimsDel,&offset[0],mzRes,rtRes,rawCoeffDEL.v);
	SMData2D<OpNabla2H> d2hAdel(dimsDel,&offset[0],mzRes,rtRes,rawCoeffDEL.v);
	rawCoeffDEL.clear();

	for(size_t i = 0; i < Adel.rt.size(); ++i)
	{
		Adel.rt[i] = rtRaw[i].second;
		dhAdel.rt[i] = rtRaw[i].second;
		d2hAdel.rt[i] = rtRaw[i].second;
	}

	BsplineData<> bsDatadel(Adel,dhAdel,d2hAdel);

	PeakManager<PeakData,BsplineData,Centroid2D> centriodPeakdel(bsDatadel);
	centriodPeakdel.execute();
	//=====================================================================

	//SMData<OpNablaV> dvA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
	//SMData<OpNabla2V> d2vA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
	//BsplineData<> bsPeakData(A,dhAdel,d2hAdel,dvA,d2vA);
	//PeakManager<PeakData,BsplineData,ExtractPeak> extractPeak(bsPeakData);
	//extractPeak.execute();

	//---------------------------------------------------------------------
	// Raw MZ and Peak Data
	//---------------------------------------------------------------------
	VecMat<double> rawMZ;
	VecMat<float> rawPK;
	mzMLb3DF.read_MatNC("spectrum_MS_1000514",rawMZ);
	mzMLb3DF.read_MatNC("spectrum_MS_1000515",rawPK);

	VecMat<double> mzPeak;
	VecMat<float> pkPeak;
	vector<size_t> vecSize;

	centriodPeakdel.peak->getPeakMatT(mzPeak, pkPeak, Adel.rt.size(), vecSize);

	repackPeakDataT(mzPeak, rawMZ, msLevel, vecSize, rawSize);
	repackPeakDataT(pkPeak, rawPK, msLevel, vecSize, rawSize);

	//---------------------------------------------------------------------
	// XML Modify
	//---------------------------------------------------------------------
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

	cout<<"It works!"<<endl;

	vector<unsigned int> mzpkSpecIdx;
	vector<unsigned int> chromatSpecIdx;
	findVecString(mzMLbuff, mzpkSpecIdx);
	findVecString(mzMLbuff, chromatSpecIdx,"<chromatogram index","</chromatogram>");

	//---------------------------------------------------------------------
	// Close Files...
	smoDF.close();

	//---------------------------------------------------------------------
	// NetCDF4 Write data to Peak mzMLb3 File.
	//---------------------------------------------------------------------
	vector<double> chromatTime;
	vector<float> chromatInten;
	hsize_t rdims[2];
	rawMZ.getDims(rdims);
	vector<double> specMZAxisRow(rdims[0],0.0);
	vector<double> specMZAxisCol(rdims[1],0.0);
	vector<float> specPKAxisRow(rdims[0],0.0);
	vector<float> specPKAxisCol(rdims[1],0.0);

	// Load rest of data from mzMLb3 file...
	mzMLb3DF.read_VecNC("chromatogram_MS_1000595",chromatTime);
	mzMLb3DF.read_VecNC("chromatogram_MS_1000515",chromatInten);

	NetCDFile mzMLb3NCDF4(outMZFileName,NC_NETCDF4);

	mzMLb3NCDF4.write_VecNC("mzML",mzMLbuff,NC_BYTE);
	vector<unsigned int> mzMLver;
	mzMLver.push_back(3);
	mzMLb3NCDF4.write_AttNC("mzML","mzMLb_version",mzMLver, NC_UINT);

	mzMLb3NCDF4.write_VecNC("mzML_chromatogramIndex",chromatSpecIdx,NC_UINT);
	mzMLb3NCDF4.write_VecNC("chromatogram_MS_1000595",chromatTime,NC_DOUBLE);
	mzMLb3NCDF4.write_VecNC("chromatogram_MS_1000515",chromatInten,NC_FLOAT);

	mzMLb3NCDF4.write_VecNC("mzML_spectrumIndex",mzpkSpecIdx,NC_UINT);
	mzMLb3NCDF4.write_VecNC("spectrum_MS_1000514_row",specMZAxisRow,NC_DOUBLE);
	mzMLb3NCDF4.write_VecNC("spectrum_MS_1000514_col",specMZAxisCol,NC_DOUBLE);
	mzMLb3NCDF4.write_MatNC("spectrum_MS_1000514",rawMZ,NC_DOUBLE,NULL,
			"spectrum_MS_1000514_row","spectrum_MS_1000514_col");
	mzMLb3NCDF4.write_VecNC("spectrum_MS_1000515_row",specPKAxisRow,NC_FLOAT);
	mzMLb3NCDF4.write_VecNC("spectrum_MS_1000515_col",specPKAxisCol,NC_FLOAT);
	mzMLb3NCDF4.write_MatNC("spectrum_MS_1000515",rawPK,NC_FLOAT,NULL,
			"spectrum_MS_1000515_row","spectrum_MS_1000515_col");

	//---------------------------------------------------------------------
	// Test Unlimited write of Vector.

	vector<double> mzAxisMZ(0,0.0);
	vector<double> mzAxisRT(rdims[1],0.0);
	mzMLb3NCDF4.write_VecNC("mzAxisMZ",mzAxisMZ,NC_DOUBLE,NULL,true);
	mzMLb3NCDF4.write_VecNC("mzAxisRT",mzAxisRT,NC_DOUBLE);

	// Test UnLimited write of Matrix.

	VecMat<double> testUMat1;
	VecMat<float> testUMat2;
	size_t udims[2];


	testUMat1.set(1,30);
	testUMat2.set(1,30);
	for(int i=0; i < 30; ++i)
	{
		testUMat1.m[0][i]=rawMZ.m[i][0];
		testUMat2.m[0][i]=rawPK.m[i][0];
	}

	udims[0]=NC_UNLIMITED;
	udims[1]=size_t(rdims[1]);
	mzMLb3NCDF4.write_DefUMatNC<double>("UMZTest","mzAxisMZ","mzAxisRT",NC_DOUBLE);
	mzMLb3NCDF4.write_DefUMatNC<float>("UPKTest",udims,NC_FLOAT);

	size_t wrcIdx[2];
	size_t wlen[2];

	wrcIdx[0]=0;
	wrcIdx[1]=0;
	wlen[0]=30;
	wlen[1]=1;

	mzMLb3NCDF4.write_PutUMatNC("UMZTest",testUMat1,wrcIdx,wlen);
	mzMLb3NCDF4.write_PutUMatNC("UPKTest",testUMat2,wrcIdx,wlen);

	testUMat1.clear();
	testUMat1.set(1,33);
	testUMat2.clear();
	testUMat2.set(1,33);
	for(int i=0; i < 33; ++i)
	{
		testUMat1.m[0][i]=rawMZ.m[i][1];
		testUMat2.m[0][i]=rawPK.m[i][1];
	}
	wrcIdx[0]=3;
	wrcIdx[1]=2;
	wlen[0]=33;
	wlen[1]=1;

	mzMLb3NCDF4.write_PutUMatNC("UMZTest",testUMat1,wrcIdx,wlen);
	mzMLb3NCDF4.write_PutUMatNC("UPKTest",testUMat2,wrcIdx,wlen);

	//---------------------------------------------------------------------

	if(debug)
	{
		string outFileName=smoFileName.substr(0,smoFileName.size()-4);
		// Write data to SMP file.
		SMPFile smpDataFile(outFileName);

		mzMLdump(string(outFileName+"_mzML.txt"),output);

		cout<<"\nSaving Peak Debugging Data to File:"<<endl;

		vector<hsize_t> vecN;
		vecN.push_back(0.0);
		vecN[0]=centriodPeakdel.peak->getMZ().size();
		smpDataFile.write_VecMatH5("Peak_mz",centriodPeakdel.peak->getMZ(),vecN,H5::PredType::NATIVE_DOUBLE);
		vecN[0]=centriodPeakdel.peak->getMZwidth().size();
		smpDataFile.write_VecMatH5("Peak_mz_width",centriodPeakdel.peak->getMZwidth(),vecN,H5::PredType::NATIVE_DOUBLE);
		vecN[0]=centriodPeakdel.peak->getRT().size();
		smpDataFile.write_VecMatH5("Peak_rt",centriodPeakdel.peak->getRT(),vecN,H5::PredType::NATIVE_DOUBLE);
		vecN[0]=centriodPeakdel.peak->getRTwidth().size();
		smpDataFile.write_VecMatH5("Peak_rt_width",centriodPeakdel.peak->getRTwidth(),vecN,H5::PredType::NATIVE_DOUBLE);
		vecN[0]=centriodPeakdel.peak->getPKcount().size();
		smpDataFile.write_VecMatH5("Peak_Count",centriodPeakdel.peak->getPKcount(),vecN,H5::PredType::NATIVE_FLOAT);
		vecN[0]=centriodPeakdel.peak->getMZIdx().size();
		smpDataFile.write_VecMatH5("Peak_mz_idx",centriodPeakdel.peak->getMZIdx(),vecN,H5::PredType::NATIVE_LLONG);
		vecN[0]=centriodPeakdel.peak->getRTIdx().size();
		smpDataFile.write_VecMatH5("Peak_rt_idx",centriodPeakdel.peak->getRTIdx(),vecN,H5::PredType::NATIVE_LLONG);

		smpDataFile.write_VecMatH5("csOrig",Adel.alpha->v,dimsDel,H5::PredType::NATIVE_FLOAT);
		smpDataFile.write_VecMatH5("dcs",dhAdel.alpha->v,dimsDel,H5::PredType::NATIVE_FLOAT);
		smpDataFile.write_VecMatH5("d2cs",d2hAdel.alpha->v,dimsDel,H5::PredType::NATIVE_FLOAT);
		//smpDataFile.write_VecMatH5("dvcs",dvA.alpha->v,dims,H5::PredType::NATIVE_FLOAT);
		//smpDataFile.write_VecMatH5("d2vcs",d2vA.alpha->v,dims,H5::PredType::NATIVE_FLOAT);

		mzPeak.getDims(dimsDel);
		smpDataFile.write_VecMatH5("mat_Peak_MZ",mzPeak.v,dims,H5::PredType::NATIVE_DOUBLE);
		smpDataFile.write_VecMatH5("mat_Peak_PK",pkPeak.v,dims,H5::PredType::NATIVE_FLOAT);

		rawMZ.getDims(dimsDel);
		smpDataFile.write_VecMatH5("mat_raw_Peak_MZ",rawMZ.v,dims,H5::PredType::NATIVE_DOUBLE);
		smpDataFile.write_VecMatH5("mat_raw_Peak_PK",rawPK.v,dims,H5::PredType::NATIVE_FLOAT);
	}

	return 0;
}

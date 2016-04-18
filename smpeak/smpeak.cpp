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

#include <iostream>
#include <iterator>
#include <boost/program_options.hpp>
#include <pugixml.hpp>
#include <sstream>
#include <omp.h>

#include "peakcore.hpp"
#include "NetCDFile.hpp"
#include "SMPFile.hpp"

#include "SMData.hpp"
#include "MathOperator.hpp"
#include "BsplineData.hpp"
#include "PeakOperator.hpp"
#include "PeakData.hpp"
#include "PeakManager.hpp"

namespace po = boost::program_options;
namespace xml = pugi;

int main(int argc, char **argv)
{
	enum Algorithm{NOPROC,CENT1DIM,CENT2DIM,PEAKPICK};
	string smoFileName;
	string mzMLFileName;
	string outMZFileName;
	string optfcscs;
	Algorithm process;
	double threshold;
	int threads;
	bool debug;
	typedef pair<int,double> rtIdxData;

	po::options_description general("Usage: smpeak [OPTION...] [SMO FILE]\n"
			"Options");

	general.add_options()
		("help,h", "Produce help message")
		("centroid,c", "Centroid mode Peak transformation along M/Z")
		("peak,p", "Peak Pick mode, across 2D RT vs M/Z data set")
		("threshold,l",po::value<double>(&threshold)->default_value(0.0), "Threshold Level for Peak Counts. Default = 0")
		("smo,s", po::value<string>(&smoFileName), "Filename of input smo data file")
		("mzMLb3,z", po::value<string>(&mzMLFileName), "Filename of input mzMLb3 data file")
		("output,o", po::value<string>(&outMZFileName), "Filename of output peak data file")
		("threads,t",po::value<int>(&threads), "Number of CPU threads - default system max")
		("debug,g", "Output debugging Peak data");

	po::options_description hidden("Hidden options");
	hidden.add_options()
		("hidden","Hidden help message")
		("centroid-2d","Centroid mode Peak transformation along M/Z, using 2D algorithm");

	po::options_description examples("Eamples:\n"
			"Centroid mode:\n"
			"\tsmpeak -c -s [SMO FILE] -z [mzMLb3 FILE]\n"
			"\tsmpeak -s [SMO FILE] -z [mzMLb3 FILE]\n"
			"\tsmpeak -c [SMO FILE] -z [mzMLb3 FILE]\n"
			"\tsmpeak [SMO FILE] -z [mzMLb3 FILE]\n"
			"2D Peak Pick mode:\n"
			"\tsmpeak -p -s [SMO FILE]\n"
			"\tsmpeak -p [SMO FILE]\n"
			"Default mode if -c or -p not given - Centroid");

	po::options_description cmdline;
	cmdline.add(general).add(hidden).add(examples);

	po::options_description desc;
	desc.add(general).add(examples);

	try
	{
		po::positional_options_description pod;
		pod.add("smo", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(cmdline).positional(pod).run(), vm);
		po::notify(vm);

		if(vm.count("help"))
		{
			cout<<desc<<endl;
			return 0;
		}
		if(vm.count("hidden"))
		{
			cout<<cmdline<<endl;
			return 0;
		}
		if(vm.count("centroid"))
		{
			cout<<"Transforming data to peak detection Centroid, 1D scan mode."<<endl;
			process=CENT1DIM;
			optfcscs="fcs";
		}
		else if(vm.count("centroid-2d"))
		{
			cout<<"Transforming data to peak detection Centroid, 2D scan mode."<<endl;
			process=CENT2DIM;
			optfcscs="fcs";
		}
		else if(vm.count("peak"))
		{
			cout<<"Extract 2D Peak from RT vs M/Z data set."<<endl;
			process=PEAKPICK;
			optfcscs="cs";
		}
		else if(vm.count("smo"))
		{
			cout<<"Transforming data to peak detection Centroid, 1D scan mode."<<endl;
			process=CENT1DIM;
			optfcscs="fcs";
		}
		else
		{
			process=NOPROC;
			optfcscs="fcs";
		}
		if(vm.count("threads"))
		{
			threads=vm["threads"].as<int>();
		}
		else
		{
			threads=omp_get_max_threads();
		}
		if(vm.count("smo"))
		{
			cout<<"Opening SMO file: "<<vm["smo"].as<string>()<<endl;
		}
		else
		{
			throw "SMO file was not give...";
		}
		if((vm.count("mzMLb3") && process == CENT1DIM) || (vm.count("mzMLb3") && process == CENT2DIM))
		{
			cout<<"Opening mzMLb3 file: "<<vm["mzMLb3"].as<string>()<<endl;
		}
		else if(process == CENT1DIM || process == CENT2DIM)
		{
			throw "mzMLb3 file was not give...";
		}
		else if(vm.count("mzMLb3") && process == PEAKPICK)
		{
			throw "Invalid program options for given operation...";
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
		cout<<desc<<endl;
		return 1;
	}
	catch(const char* msg)
	{
		cerr<<"error: "<<msg<<endl;
		cout<<desc<<endl;
		return 1;
	}
	catch(...)
		{
		cerr<<"Exception of unknown type!\n";
	}

	omp_set_num_threads(threads);
	cout<<"Number of threads: "<<threads<<endl;

	//---------------------------------------------------------------------
	// NetCDF4 Read data
	//---------------------------------------------------------------------
	cout<<"Processing SMO Data Set..."<<endl;

	vector<char> mzMLbuff;
	NetCDFile mzMLb3DF;

	vector<int> msLevel;
	vector<rtIdxData> rtRaw;
	vector<size_t> rawSize;
	bool constPrecursor=true;
	xml::xml_document docmzML;
	xml::xpath_node_set tools;

	if(process == CENT1DIM || process == CENT2DIM)
	{
		//-------------------------------
		// Load mzMLb3 Data that we need.
		//-------------------------------
		mzMLb3DF.open(mzMLFileName);
		mzMLb3DF.read_VecNC("mzML",mzMLbuff);

		//-------------------------------------------
		// XML scan and read Start Time and ms-level.
		//-------------------------------------------
		size_t xmlSize=sizeof(char)*mzMLbuff.size();
		xml::xml_parse_result result = docmzML.load_buffer_inplace(&mzMLbuff[0],xmlSize);

		tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/precursorList/precursor/selectedIonList/selectedIon/cvParam[@accession='MS:1000744']");

		if(tools.empty())
		{
			constPrecursor = false;
		}
		else
		{
			for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end()-1; ++itr)
			{
				double prc0;
				double prc1;
				istringstream(itr->node().attribute("value").value()) >> prc0;
				istringstream((itr+1)->node().attribute("value").value()) >> prc1;
				if(prc0 != prc1)
				{
					constPrecursor = false;
					break;
				}
			}
		}

		tools =
			docmzML.select_nodes("mzML/run/spectrumList/spectrum/cvParam[@name='ms level']");

		if(tools.empty())
		{
			cout<<"Error reading XML data from mzML file..."<<endl;
		}
		else
		{
			int idx=0;
			for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
			{
				int ms;
				double rt;
				size_t len;

				istringstream(itr->node().attribute("value").value()) >> ms;
				if(constPrecursor)
					msLevel.push_back(3);
				else
					msLevel.push_back(ms);

				istringstream(itr->node().parent().attribute("defaultArrayLength").value()) >> len;
				rawSize.push_back(len);

				if(msLevel.back() == 1 || msLevel.back() == 3 )
				{
					istringstream(itr->node().parent().child("scanList").child("scan").find_child_by_attribute("accession","MS:1000016").attribute("value").value()) >> rt;
					rtRaw.push_back(make_pair(idx,rt));
					++idx;
				}
			}
		}
	}

	//----------------------------
	// Load SMO Data that we need.
	//----------------------------
	NetCDFile smoDF(smoFileName);

	vector<int> offset;
	vector<size_t> dataMatLen;
	vector<InfoGrpVar> dataSetList;

	cout << "List all groups within file: " << smoFileName << endl;

	smoDF.search_Group(optfcscs);
	dataSetList = smoDF.get_Info();

	for(int i=0; i < dataSetList.size(); ++i)
		cout<<"DataSets found ["<< i<<"]: "<<dataSetList[i].grpName<<"/"
		<<dataSetList[i].varName<<endl;

	cout<<"List Parameters from SMO:"<<endl;
	double mzRes=smoDF.search_Group<double>(2);
	cout<<"MZ_resolution: "<<mzRes<<endl;
	double rtRes=smoDF.search_Group<double>(3);
	cout<<"RT resolution: "<<rtRes<<endl;

	smoDF.read_AttNC("Offset",dataSetList[0].varid,offset,dataSetList[0].grpid);
	dataMatLen = smoDF.read_DimNC(dataSetList[0].varName,dataSetList[0].grpid);

	VecMat<double> mzPeak;
	VecMat<float> pkPeak;
	vector<size_t> mzpkVecSize;

	//---------------------------------------------------------------------
	// Centroid Peak Data Set...
	// Process Scan by Scan
	//---------------------------------------------------------------------
	if(process == CENT1DIM)
	{
		vector<hsize_t> dims(2,0);
		vector<size_t> hypIdx(2);
		vector<size_t> rdLen(2);

		rdLen[0]=1;
		rdLen[1]=dataMatLen[1];
		dims[0]=hsize_t(rdLen[0]);
		dims[1]=hsize_t(rdLen[1]);
		hypIdx[1]=0; // = Always read from first Column;

		PeakData<> totalPeaks;
		vector<PeakData<> *> peakThreads(omp_get_max_threads());

		cout<<"Extract Peaks from Mass Spec Data"<<endl;
		// Manual reduction of STL container as STLs are not thread safe...
		int run=0;
		#pragma omp parallel
		{
			int nthrd=omp_get_num_threads();
			PeakData<> localPeaks;
			int thrdid=omp_get_thread_num();
			peakThreads[thrdid] = &localPeaks;
			run=0;

			#pragma omp for firstprivate(hypIdx,rdLen) schedule(dynamic)
			for(ii rt_idx = 0; rt_idx < dataMatLen[0]; ++rt_idx)
			{
				vector<float> rawCoeff;
				hypIdx[0]=rt_idx;

				if(thrdid == 0)
					cout<<"\r"<<"Processing Scan: "<<dataMatLen[0]<<"/"<<run<<flush;

				smoDF.read_HypVecNC(dataSetList[0].varName,rawCoeff,&hypIdx[0],&rdLen[0],
						dataSetList[0].grpid);

				SMData1D<OpUnitS> A(&dims[0],&offset[0],mzRes,rtRaw[rt_idx],rawCoeff);
				SMData1D<OpNablaHS> dhA(&dims[0],&offset[0],mzRes,rtRaw[rt_idx],rawCoeff);
				SMData1D<OpNabla2HS> d2hA(&dims[0],&offset[0],mzRes,rtRaw[rt_idx],rawCoeff);

				BsplineData<rtIdxData> bsData(A,dhA,d2hA);

				PeakManager<PeakData,BsplineData,Centroid1D,rtIdxData> centriodPeak(bsData,threshold);
				centriodPeak.execute();

				localPeaks.addPeakArray(centriodPeak.peak->getPeakData());
				localPeaks.updateFalseData(centriodPeak.peak->getFalsePeaks(),centriodPeak.peak->getFalseWidths());
				#pragma omp atomic
					++run;
			}

			#pragma omp single
			{
				cout<<"\r"<<"Processing Scan: "<<dataMatLen[0]<<"/"<<dataMatLen[0]<<endl;
				cout<<"Gathering Peaks from all threads..."<<endl;
				for(int i = 0; i < nthrd; ++i)
				{
					//if(debug)
					cout<<"Thread ["<<i<<"] peaks found: "<<peakThreads[i]->numOfPeaks()<<endl;
					totalPeaks.addPeakArray(peakThreads[i]->getPeakData());
					totalPeaks.updateFalseData(peakThreads[i]->getFalsePeaks(),peakThreads[i]->getFalseWidths());
					peakThreads[i]->clear();
				}
			}
		}

		totalPeaks.getPeakMat(mzPeak, pkPeak, dataMatLen[0], mzpkVecSize);
		cout<<"Total Peaks found - "<<"["<<totalPeaks.numOfPeaks()<<"]"<<endl;
		if(totalPeaks.getFalsePeaks() > 0 || debug == true)
			cout<<"Total insignificant false Peaks detected - "
				<<"["<<totalPeaks.getFalsePeaks()<<"] Peaks Ignored"<<endl;
		if(totalPeaks.getFalsePeaks() > 0 || debug == true)
			cout<<"Total false Peak detected with incorrect Peak Widths - "
				<<"["<<totalPeaks.getFalseWidths()<<"] Peaks Ignored"<<endl;

		if(debug) totalPeaks.dumpPeakData(smoFileName,H5::PredType::NATIVE_FLOAT);
	}
	else if(process == CENT2DIM)
	{
		hsize_t dims[2];
		VecMat<float> rawCoeff;
		smoDF.read_MatNC(dataSetList[0].varName,rawCoeff,dataSetList[0].grpid);
		rawCoeff.getDims(dims);

		SMData2D<OpUnit> A(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
		SMData2D<OpNablaH> dhA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
		SMData2D<OpNabla2H> d2hA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);

		for(size_t i = 0; i < A.rt.size(); ++i)
		{
			A.rt[i] = rtRaw[i].second;
			dhA.rt[i] = rtRaw[i].second;
			d2hA.rt[i] = rtRaw[i].second;
		}

		BsplineData<> bsData(A,dhA,d2hA);

		PeakManager<PeakData,BsplineData,Centroid2D> centriodPeak(bsData,threshold);
		centriodPeak.execute();
		centriodPeak.peak->getPeakMat(mzPeak, pkPeak, dataMatLen[0], mzpkVecSize);

		if(debug) centriodPeak.peak->dumpPeakData(smoFileName);
	}
	else if(process == PEAKPICK)
	{
		hsize_t dims[2];
		VecMat<float> rawCoeff;
		smoDF.read_MatNC(dataSetList[0].varName,rawCoeff,dataSetList[0].grpid);
		rawCoeff.getDims(dims);

		BasisPatch<> bpatch;

		SMData2D<OpUnit> A(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
		SMData2D<OpNablaH> dhA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
		SMData2D<OpNabla2H> d2hA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
		SMData2D<OpNablaV> dvA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
		SMData2D<OpNabla2V> d2vA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);

		BsplineBasisData<> bsPeakData(A,dhA,d2hA,dvA,d2vA,bpatch);
		PeakManager<PeakData,BsplineBasisData,ExtractPeak> extractPeak(bsPeakData,threshold);
		extractPeak.execute();

		cout<<"Total Peaks found - "<<"["<<extractPeak.peak->numOfPeaks()<<"]"<<endl;
		if(extractPeak.peak->getFalsePeaks() > 0 || debug == true)
			cout<<"Total insignificant false Peaks detected - "
				<<"["<<extractPeak.peak->getFalsePeaks()<<"] Peaks Ignored"<<endl;
		if(extractPeak.peak->getFalseWidths() > 0 || debug == true)
			cout<<"Total false Peak detected with incorrect Peak Widths - "
				<<"["<<extractPeak.peak->getFalseWidths()<<"] Peaks Ignored"<<endl;

		extractPeak.peak->dumpPeakData(smoFileName,H5::PredType::NATIVE_FLOAT);
	}
	else if(process == NOPROC)
	{
		cout<<"No process performed on Data Set..."<<endl;
		return 0;
	}

	//---------------------------------------------------------------------
	// Close Files...
	smoDF.close();

	if(process == CENT1DIM || process == CENT2DIM)
	{
		//---------------------------------------------------------------------
		// XML Modify
		//---------------------------------------------------------------------
		tools =	docmzML.select_nodes("mzML/run/spectrumList/spectrum");

		for(int i = 0, j = 0; i < tools.size(); ++i)
		{
			if( (msLevel[i] == 1) && (j < mzpkVecSize.size()) )
			{
				ostringstream buff;
				buff << mzpkVecSize[j];
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

		if(debug)
		{
			cout<<"\nSaving Peak Debugging Data to File:"<<endl;
			string outFileName=smoFileName.substr(0,smoFileName.size()-4);
			mzMLdump(string(outFileName+"_mzML.txt"),output);
		}

		// Free up memory
		newmzML.clear();
		newmzML.str(std::string());
		string(output).swap(output);

		//---------------------------------------------------------------------
		// NetCDF4 Write data to Peak mzMLb3 File.
		//---------------------------------------------------------------------

		cout<<"Writing Data to new mzMLb3 file:"<<endl;

		vector<double> chromatTime;
		vector<float> chromatInten;
		vector<size_t> rdims = mzMLb3DF.read_DimNC("spectrum_MS_1000514");

		vector<double> specMZAxisRow(rdims[0],0.0);
		vector<float> specPKAxisRow(rdims[0],0.0);
		vector<double> specMZAxisCol;
		vector<float> specPKAxisCol;

		vector<unsigned int> mzpkSpecIdx;
		vector<unsigned int> chromatSpecIdx(1,0);
		findVecString(mzMLbuff, mzpkSpecIdx);

		// Load rest of data from mzMLb3 file...
		NetCDFile mzMLb3NCDF4(outMZFileName,NC_NETCDF4);

		mzMLb3NCDF4.write_VecNC("mzML",mzMLbuff,NC_BYTE);
		vector<unsigned int> mzMLver;
		mzMLver.push_back(3);
		mzMLb3NCDF4.write_AttNC("mzML","mzMLb_version",mzMLver, NC_UINT);

		if(mzMLb3DF.read_VarIDNC("chromatogram_MS_1000595") >= 0)
		{
			findVecString(mzMLbuff, chromatSpecIdx,"<chromatogram index","</chromatogram>");
			mzMLb3DF.read_VecNC("chromatogram_MS_1000595",chromatTime);
			mzMLb3DF.read_VecNC("chromatogram_MS_1000515",chromatInten);
			mzMLb3NCDF4.write_VecNC("chromatogram_MS_1000595",chromatTime,NC_DOUBLE);
			mzMLb3NCDF4.write_VecNC("chromatogram_MS_1000515",chromatInten,NC_FLOAT);
		}
		mzMLb3NCDF4.write_VecNC("mzML_chromatogramIndex",chromatSpecIdx,NC_UINT);

		mzMLb3NCDF4.write_VecNC("mzML_spectrumIndex",mzpkSpecIdx,NC_UINT);
		// Unlimited dimensions
		mzMLb3NCDF4.write_VecNC("spectrum_MS_1000514_row",specMZAxisRow,NC_DOUBLE);
		mzMLb3NCDF4.write_VecNC("spectrum_MS_1000515_row",specPKAxisRow,NC_FLOAT);
		mzMLb3NCDF4.write_VecNC("spectrum_MS_1000514_col",specMZAxisCol,NC_DOUBLE,NULL,true);
		mzMLb3NCDF4.write_VecNC("spectrum_MS_1000515_col",specPKAxisCol,NC_FLOAT,NULL,true);

		// Define the unlimited ragged matrix dimensions
		mzMLb3NCDF4.write_DefHypMatNC<double>("spectrum_MS_1000514",
				"spectrum_MS_1000514_row","spectrum_MS_1000514_col",NC_DOUBLE);
		mzMLb3NCDF4.write_DefHypMatNC<float>("spectrum_MS_1000515",
				"spectrum_MS_1000515_row","spectrum_MS_1000515_col",NC_FLOAT);

		// Write ragged unlimited matrix scan by scan
		hsize_t centrc[2];
		mzPeak.getDims(centrc);
		for(size_t rt_idx = 0, ms1 = 0; rt_idx < rdims[0]; ++rt_idx)
		{
			size_t wrcIdx[2]={rt_idx,0};
			size_t wlen[2]={1,0};
			vector<double> writeMZCoeff;
			vector<float> writePKCoeff;
			double *scanMZ;
			float *scanPK;

			cout<<"\r"<<"Writing Scan: "<<rdims[0]<<"/"<<rt_idx+1<<flush;
			if((msLevel[rt_idx] == 1  || msLevel[rt_idx] == 3) && ms1 < centrc[0])
			{
				scanMZ = &mzPeak.m[ms1][0];
				scanPK = &pkPeak.m[ms1][0];

				wlen[1]=mzpkVecSize[ms1];
				mzMLb3NCDF4.write_PutHypMatNC("spectrum_MS_1000514",scanMZ,wrcIdx,wlen);
				mzMLb3NCDF4.write_PutHypMatNC("spectrum_MS_1000515",scanPK,wrcIdx,wlen);
				++ms1;
			}
			else if(msLevel[rt_idx] == 2)
			{
				size_t rawLen[2]={1,rawSize[rt_idx]};

				mzMLb3DF.read_HypVecNC("spectrum_MS_1000514",writeMZCoeff,wrcIdx,rawLen);
				mzMLb3DF.read_HypVecNC("spectrum_MS_1000515",writePKCoeff,wrcIdx,rawLen);

				mzMLb3NCDF4.write_PutHypMatNC("spectrum_MS_1000514",&writeMZCoeff[0],wrcIdx,rawLen);
				mzMLb3NCDF4.write_PutHypMatNC("spectrum_MS_1000515",&writePKCoeff[0],wrcIdx,rawLen);
			}
			else
			{
				size_t rawLen[2]={1,rawSize[rt_idx]};

				mzMLb3DF.read_HypVecNC("spectrum_MS_1000514",writeMZCoeff,wrcIdx,rawLen);
				mzMLb3DF.read_HypVecNC("spectrum_MS_1000515",writePKCoeff,wrcIdx,rawLen);

				mzMLb3NCDF4.write_PutHypMatNC("spectrum_MS_1000514",&writeMZCoeff[0],wrcIdx,rawLen);
				mzMLb3NCDF4.write_PutHypMatNC("spectrum_MS_1000515",&writePKCoeff[0],wrcIdx,rawLen);
			}
		}
		cout<<"\nCentroid complete."<<endl;
	}

	return 0;
}

#include "imagecore.hpp"
#include "SMGWriter.hpp"
#include "core.hpp"
#include <utility>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char** argv){

	string smjFileName, channel;
	int mz_res=0;
	bool reBin, ms2;
	pair<double,double> mzInputRange;
	pair<double,double> rtInputRange;
	pair<ui,ui> xyview;

	po::options_description usage("Usage: smimage [OPTION...] [SMJ FILE]\n"
			"Example: smimage -m [200,300] -r [20,140] -o [3000,1000] --ms2 datafile.smj\n"
			"         smimage -o [,500] -m[,300] -r [12,] -f datafile.smj\n\n"
			"Options");

	usage.add_options()
		("help,h", "Produce help message")
		("file,f",
			po::value<string>(&smjFileName),
				"Filename - Load data from seaMass input file format (smj).")
		("output-image,o",
			po::value<string>(),
				"Output range of generated SMG file (default values [1280,1024]). "
				"The argument takes the following format: [MZ Pixels,RT Pixels], "
				"omitting a value will result in the appropriate default values "
				".i.e. [,Val] -> will give a range [1280,Val], "
				"[Val,] -> will give a range [Val,1024].")
		("mz-range,m",
			po::value<string>(),
				"MZ range units (m/z Th) argument takes the following format: [Min Value,Max Value], "
				"omitting a value will result in the appropriate [Min,Max] value "
				".i.e. [,Val] -> will give a range [Min,Val], "
				"[Val,] -> will give a range [Val,Max].")
		("rt-range,r",
			po::value<string>(),
				"RT range units (s) argument takes the following format: [Min Value,Max Value], "
				"omitting a value will result in the appropriate [Min,Max] value "
				".i.e. [,Val] -> will give a range [Min,Val], "
				"[Val,] -> will give a range [Val,Max].")
		("ms2", "Include MS2 data.")
		("mz-rebin,b",
				po::value<int>(&mz_res),
				"Re-Bin MZ for speed. MZ resolution given as: "
				"\"b-splines per Th = 2^mz_res * 60 / 1.0033548378\" "
				"guidelines: between 0 to 2 for ToF, 3 for Orbitrap.");
	try
	{
		po::positional_options_description pod;
		pod.add("file", -1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(usage).positional(pod).run(), vm);
		po::notify(vm);

		if(vm.count("help"))
		{
			cout<<usage<<endl;
			return 0;
		}
		if(vm.count("output-image"))
		{
			xyview=numStrBraket<ui,ui>(vm["output-image"].as<string>());
			if(xyview.first == -1) xyview.first = 1280;
			if(xyview.second == -1) xyview.second = 1024;
		}
		else
		{
			xyview = make_pair(1280,1024);
		}
		if(vm.count("mz-range"))
		{
			mzInputRange=numStrBraket<double,double>(vm["mz-range"].as<string>());
		}
		else
		{
			mzInputRange = make_pair(-1.0,-1.0);
		}
		if(vm.count("rt-range"))
		{
			rtInputRange=numStrBraket<double,double>(vm["rt-range"].as<string>());
			if(rtInputRange.first != -1) rtInputRange.first = rtInputRange.first/60.0;
			if(rtInputRange.second != -1) rtInputRange.second = rtInputRange.second/60.0;
		}
		else
		{
			rtInputRange = make_pair(-1.0,-1.0);
		}
		if(vm.count("ms2"))
			ms2=true;
		else
			ms2=false;
		if(vm.count("mz-rebin"))
			reBin=true;
		else
			reBin=false;
		if(vm.count("file"))
		{
			cout<<"Opening SMJ file: "<<vm["file"].as<string>()<<endl;
		}
		else
		{
			cout<<usage<<endl;
			throw "Input SMJ file was not give...";
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

	string outFileName=smjFileName.substr(0,smjFileName.size()-4);
	SMGWriter smgFile(outFileName);

	cout<<"Load raw data from data filies:"<<endl;
	ostringstream datah5;
	MassSpecData raw;
	MassSpecImage imgBox(xyview);

	// Load Data from smj file
	datah5 << "/" << "StartTime";
	cout<<"Load Start Time from SMJ:"<<endl;
	read_VecH5(smjFileName, datah5.str(), raw.rt, H5::PredType::NATIVE_DOUBLE);
	datah5.str("");
	datah5.clear();
	raw.N = raw.rt.size();

	// Load Data from smj file
	datah5 << "/" << "PrecursorMZ";
	cout<<"Load Precursor MZ from SMJ:"<<endl;
	read_VecH5(smjFileName, datah5.str(), raw.precursorMZ, H5::PredType::NATIVE_DOUBLE);
	datah5.str("");
	datah5.clear();

	vector<double> vmz;
	// Load Data from smj file
	datah5 << "/" << "SpectrumMZ";
	cout<<"Load Spectrim MZ from SMJ:"<<endl;
	read_VecH5(smjFileName, datah5.str(), vmz, H5::PredType::NATIVE_DOUBLE);
	datah5.str("");
	datah5.clear();

	vector<double> vsc;
	// Load Data from smj file
	datah5 << "/" << "SpectrumIntensity";
	cout<<"Load Spectrum Count from SMJ:"<<endl;
	read_VecH5(smjFileName, datah5.str(), vsc, H5::PredType::NATIVE_DOUBLE);
	datah5.str("");
	datah5.clear();

	vector<lli> vsci;
	// Load Data from smj file
	datah5 << "/" <<"SpectrumIndex";
	cout<<"Load Spectrum Count Index from SMJ:"<<endl;
	read_VecH5(smjFileName, datah5.str(), vsci, H5::PredType::NATIVE_ULONG);
	datah5.str("");
	datah5.clear();

	vector<vector <double> > mzs;
	vector<vector <double> > intensities;

	mzs.resize(raw.N);
	intensities.resize(raw.N);

	for(lli i = 0; i < raw.N; ++i)
	{
		lli elements=vsci[i+1]-vsci[i]-1;
		if(elements > 0)
		{
			mzs[i].reserve(elements);
			intensities[i].reserve(elements);
		}
		for(lli j = vsci[i]; j < vsci[i+1]; ++j)
		{
			mzs[i].push_back(vmz[j]);
			intensities[i].push_back(vsc[j]);
		}
	}

	// difference between carbon12 and carbon13
	double rc_mz = pow(2.0, (double) mz_res) * 60 / 1.0033548378;
	unsigned long instrument_type = 1;
	read_AttH5(smjFileName,"/SpectrumIntensity","instrumentType",
			instrument_type, H5::IntType(H5::PredType::NATIVE_USHORT));

	// Ensure the raw data is in binned format and compute exposures
	vector<fp> exposures;
	bin_mzs_intensities(mzs, intensities, instrument_type, exposures);

	if(reBin)
	{
		// for speed only, merge bins if rc_mz is set higher than twice bin width
		cout<<"Re-scaling MZ for speed..."<<endl;
		merge_bins(mzs, intensities, 0.5 / rc_mz);
	}

	// Convert intensities into format used in algorithm for gs
	vector<ii> nonZeroI;
	create_gs(raw.sc, raw.sci, nonZeroI, intensities);

	for (ii i = 0; i < (ii) intensities.size(); i++) vector<double>().swap(intensities[i]);

	for(lli i = 0; i < mzs.size(); ++i)
	{
		for(lli j = 0; j < mzs[i].size(); ++j)
		{
			raw.mz.push_back(mzs[i][j]);
		}
	}

	for (ii i = 0; i < (ii) mzs.size(); i++) vector<double>().swap(mzs[i]);

	cout<<"Processing data to create new SMG image file..."<<endl;
	// Calculate MZ index
	raw.calMZi();
	// Find data limits in mass spec data.
	raw.calRange();

	raw.rtp.resize(raw.precursorMZ.size());

	if(ms2)
	{
		for(lli i=0; i < raw.precursorMZ.size(); ++i)
		{
			raw.rtp[i]=i;
		}
	}
	else
	{
		for(lli i=0; i < raw.precursorMZ.size(); ++i)
		{
			if(raw.precursorMZ[i] == 0.0)
			{
				raw.rtp[i]=i;
			}
			else
			{
				raw.rtp[i]=-1;
			}
		}
	}

	raw.rti.resize(raw.N,-1);
	for(lli i=0; i < nonZeroI.size(); ++i)
	{
		if(raw.rtp[nonZeroI[i]] >= 0)
			raw.rti[nonZeroI[i]]=i;
	}

	// Clip Box if needed
	// Min MZ Value
	(mzInputRange.first < raw.mzRange.first || mzInputRange.first < 0 ) ?
			imgBox.mzRange.first=raw.mzRange.first:
			imgBox.mzRange.first=mzInputRange.first;
	// Max MZ Value
	(mzInputRange.second > raw.mzRange.second || mzInputRange.second < 0 ) ?
			imgBox.mzRange.second=raw.mzRange.second:
			imgBox.mzRange.second=mzInputRange.second;
	// Min RT Value
	(rtInputRange.first < raw.rtRange.first || rtInputRange.first < 0 ) ?
			imgBox.rtRange.first=raw.rtRange.first:
			imgBox.rtRange.first=rtInputRange.first;
	// Max RT Value
	(rtInputRange.second > raw.rtRange.second || rtInputRange.second < 0 ) ?
			imgBox.rtRange.second=raw.rtRange.second:
			imgBox.rtRange.second=rtInputRange.second;

	imgBox.dmz=genAxis(imgBox.mz, imgBox.mzRange.first, imgBox.mzRange.second);
	imgBox.drt=genAxis(imgBox.rt, imgBox.rtRange.first, imgBox.rtRange.second);

	// Find RT range index n raw data, as data has irregular DRT we have to scan manually.
	lli rtidxBegin=0, rtidxEnd=0;
	for(lli i = 0; i < raw.rt.size(); ++i)
	{
		if(imgBox.rt[0] < raw.rt[i])
		{
			rtidxBegin=i-1;
			break;
		}
	}
	for(lli i = rtidxBegin; i < raw.rt.size(); ++i)
	{
		if(imgBox.rt[imgBox.xypxl.second-1] < raw.rt[i])
		{
			rtidxEnd=i;
			break;
		}
		else if(i == raw.rt.size()-1)
		{
			rtidxEnd=raw.rt.size()-1;
		}
	}

	for(lli idx=rtidxBegin; idx <=rtidxEnd; ++idx)
	{
		double scaleRT=0.0;
		double drt=fabs(raw.rt[idx+1]-raw.rt[idx]);
		double yn=(raw.rt[idx]-imgBox.rt[0])/imgBox.drt;
		double yp=(raw.rt[idx+1]-imgBox.rt[0])/imgBox.drt;
		lli rowN=floor(yn);
		lli rowP=floor(yp);

		if(rowP > imgBox.xypxl.second-2) rowP=imgBox.xypxl.second-1;
		if(rowN >= imgBox.xypxl.second-1) break;
		// Overlapping image 1st box on boundary
		if(rowN < 0 && rowP >= 0)
		{
			rowN=0;
			for(lli row=rowN; row <= rowP; ++row)
			{
				double boxBegin=0.0;
				double boxEnd=0.0;
				// Cut top overlap else its contained
				(imgBox.rt[row] < raw.rt[idx]) ? boxBegin=raw.rt[idx] : boxBegin=imgBox.rt[row];
				// Cut bottom overlap else its contained
				(imgBox.rt[row+1] < raw.rt[idx+1]) ? boxEnd=imgBox.rt[row+1] : boxEnd=raw.rt[idx+1];
				scaleRT=(boxEnd-boxBegin)/drt;
				// Run MZ scan calculations
				if(raw.rti[idx] >= 0) scanMZ(raw, imgBox, raw.rti[idx], row, scaleRT);
			}
		}
		// Multiple smaller image boxes within single raw data scanline
		else if(rowN < rowP && rowN >=0)
		{
			for(lli row=rowN; row <= rowP; ++row)
			{
				double boxBegin=0.0;
				double boxEnd=0.0;
				// Cut top overlap else its contained
				(imgBox.rt[row] < raw.rt[idx]) ? boxBegin=raw.rt[idx] : boxBegin=imgBox.rt[row];
				// Cut bottom overlap else its contained
				(imgBox.rt[row+1] < raw.rt[idx+1]) ? boxEnd=imgBox.rt[row+1] : boxEnd=raw.rt[idx+1];
				scaleRT=(boxEnd-boxBegin)/drt;
				// Run MZ scan calculations
				if(raw.rti[idx] >= 0) scanMZ(raw, imgBox, raw.rti[idx], row, scaleRT);
			}
		}
		// Raw data scanline contained within a Image box hence no scale needed.
		else if(rowN == rowP)
		{
			scaleRT=1.0;
			// Run MZ scan calculations
			if(raw.rti[idx] >= 0) scanMZ(raw, imgBox, raw.rti[idx], rowN, scaleRT);
		}
	}

	cout<<"Writing out data to seaMass image data file: "<<endl;
	cout<<"    "<<outFileName+".smg"<<endl;

	vector<hsize_t> mzDims;
	mzDims.push_back(imgBox.mz.size());
	smgFile.write_VecMatH5("SpectrumMZ", imgBox.mz, mzDims, H5::PredType::NATIVE_DOUBLE);

	vector<hsize_t> scDims;
	scDims.push_back(imgBox.xypxl.second);
	scDims.push_back(imgBox.xypxl.first);
	smgFile.write_VecMatH5("SpectrumCount", imgBox.sc, scDims, H5::PredType::NATIVE_FLOAT);

	//Convert RT to seconds...
	for(int i=0;i < imgBox.rt.size(); ++i)
	{
		imgBox.rt[i]=imgBox.rt[i]*60.0;
	}

	vector<hsize_t> rtDims;
	rtDims.push_back(imgBox.rt.size()-1);
	smgFile.write_VecMatH5("StartTime", imgBox.rt, rtDims, H5::PredType::NATIVE_DOUBLE);

	return 0;
}

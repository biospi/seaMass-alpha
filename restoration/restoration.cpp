//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> manchester.ac.uk>
//
// Copyright (C) 2013  CADET Laboratory for Medical Bioinformatics, University of Manchester, UK
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
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <H5Cpp.h>

#include "seamass.hpp"

using namespace std;


struct spectrum
{
	size_t index;
	string preset_scan_configuration;
	int ms_level;
	double scan_start_time;
	size_t scan_start_time_index;
};


bool scan_start_time_order(const spectrum& lhs, const spectrum& rhs)
{
	return lhs.scan_start_time <= rhs.scan_start_time;
}


bool seamass_order(const spectrum& lhs, const spectrum& rhs)
{
	if (lhs.preset_scan_configuration == rhs.preset_scan_configuration)
	{
		if (lhs.ms_level == rhs.ms_level)
		{
			return lhs.scan_start_time <= rhs.scan_start_time;
		}
		else
		{
			return lhs.ms_level < rhs.ms_level;
		}
	}
	else
	{
		return lhs.preset_scan_configuration < rhs.preset_scan_configuration;
	}
}


int main(int argc, char *argv[])
{
	H5::Exception::dontPrint();
	seamass::notice();
	if (argc < 8)
	{
		cout << "Usage" << endl;
		cout << "-----" << endl;
		cout << "seamass <in_file> <mz_res> <rt_res> <shrinkage> <tol> <threads> <out_type>" << endl;
		cout << endl;
		cout << "<in_file>:  Raw input file in seaMass Input format (smi)" << endl;
		cout << "            guidelines: Use pwiz-seamass to convert from mzML or vendor format" << endl;
		cout << "<mz_res>:   MS resolution given as: \"b-splines per Th = 2^mz_res * 60 / 1.0033548378\"" << endl;
		cout << "            guidelines: between 0 to 2 for ToF, 3 for Orbitrap" << endl;
		cout << "<rt_res>:   LC resolution given as: \"b-splines per minute = 2^rt_res\"" << endl;
		cout << "            guidelines: around 4" << endl;
		cout << "<shrink>:   Amount of denoising given as: \"L1 shrinkage = 2^shrinkage\"" << endl;
		cout << "            guidelines: around -4" << endl;
		cout << "<tol>:      Convergence tolerance, given as: \"gradient <= 2^tol\"" << endl;
		cout << "            guidelines: around -9" << endl;
		cout << "<threads>:  Number of OpenMP threads to use" << endl;
		cout << "            guidelines: set to amount of CPU cores or 4, whichever is smaller" << endl;
		cout << "<out_type>: Type of output desired" << endl;
		cout << "            guidelines: 0 = just viz_client input; 1 = also smo; 2 = also debug" << endl;
		return 0;
	}
	string in_file = argv[1];
	int mz_res = atoi(argv[2]);
	int rt_res = atoi(argv[3]);
	int shrink = atoi(argv[4]);
	int tol = atoi(argv[5]);
	int threads = atoi(argv[6]);
	int out_type = atoi(argv[7]);

	int lastdot = in_file.find_last_of("."); 
	string id = (lastdot == string::npos) ? in_file : in_file.substr(0, lastdot); 

	// open H5 file
	H5::H5File file(in_file, H5F_ACC_RDONLY);

	// find instrument type - VERY ROUGH ATM
	int instrument_type = 0;
	try
	{
		H5::Group as = file.openGroup("/instrumentConfigurationList/instrumentConfiguration/0/componentList/analyzer");
		hsize_t na = as.getNumObjs();
		for (size_t i = 0; i < na; i++)
		{
			try
			{
				ostringstream oss; oss << as.getObjnameByIdx(i) << "/cvParam/MS:1000484";
				H5::Group a = as.openGroup(oss.str());
				instrument_type = 1;
				break;
			}
			catch  (const H5::GroupIException& e) {}
		}		
	}
	catch (const H5::FileIException& e) {}

	// open spectrumList
	H5::Group ss = file.openGroup("/mzML/run/spectrumList/spectrum");
	hsize_t ns = ss.getNumObjs();

	// query necessary metadata
	cout << "Querying metadata..." << endl;
	vector<spectrum> spectra(ns);
	for (size_t i = 0; i < ns; i++)
	{
		spectra[i].index = i;
		H5::Group s = ss.openGroup(ss.getObjnameByIdx(i));

		try
		{
			H5::Group cv = s.openGroup("scanList/scan/0/cvParam/MS:1000616");
			H5::Attribute sst_a = cv.openAttribute("value");
			sst_a.read(H5::StrType(0, H5T_VARIABLE), spectra[i].preset_scan_configuration);
		}
		catch (const H5::GroupIException& e)
		{
			spectra[i].preset_scan_configuration = "0";
		}

		try
		{
			H5::Group cv = s.openGroup("cvParam/MS:1000579");
			spectra[i].ms_level = 1;
		}
		catch (const H5::GroupIException& e)
		{
			spectra[i].ms_level = 2;
		}

		try
		{
			H5::Group cv = s.openGroup("scanList/scan/0/cvParam/MS:1000016");
			H5::Attribute sst_a = cv.openAttribute("value");
			H5std_string sst_s("");
			sst_a.read(H5::StrType(0, H5T_VARIABLE), sst_s);
			spectra[i].scan_start_time = atof(sst_s.c_str());	
			
			H5::Attribute ua_a = cv.openAttribute("unitAccession");
			H5std_string units("");
			ua_a.read(H5::StrType(0, H5T_VARIABLE), units);
			if (units == "UO:0000010") spectra[i].scan_start_time /= 60.0; // seconds
		}
		catch (const H5::GroupIException& e)
		{
			spectra[i].scan_start_time = -1.0;
		}
	}

	// determine start_scan_time order of spectra
	sort(spectra.begin(), spectra.end(), scan_start_time_order);
	for (size_t i = 0; i < spectra.size(); i++)
	{
		spectra[i].scan_start_time_index = i;
	}

	// save scan_start_times and sort into seamass processing order:
	// preset_scan_configuration -> ms_level -> scan_start_time
	vector<double> scan_start_times(ns);
	for (size_t i = 0; i < spectra.size(); i++)
	{
		scan_start_times[i] = spectra[i].scan_start_time;
	}
	sort(spectra.begin(), spectra.end(), seamass_order);

	// load spectra and process
	cout << "Loading spectra..." << endl;
	vector< std::vector<double> > mzs(ns);
	vector< std::vector<double> > intensities(ns);
	int loaded = 0;
	for (size_t i = 0; i < spectra.size(); i++)
	{
		H5::Group s = ss.openGroup(ss.getObjnameByIdx(spectra[i].index));	

		// load spectra into mzs and intensities vectors if profile & ms1
		try
		{
			s.openGroup("cvParam/MS:1000128"); // profile mode
					
			hsize_t n;

			H5::DataSet mzs_ds = s.openDataSet("binaryDataArrayList/binaryDataArray/MS:1000514/binary");
			mzs_ds.getSpace().getSimpleExtentDims(&n);
			mzs[spectra[i].scan_start_time_index].resize(n);			
			mzs_ds.read(mzs[spectra[i].scan_start_time_index].data(), H5::DataType(H5::PredType::NATIVE_DOUBLE));

			H5::DataSet intensities_ds = s.openDataSet("binaryDataArrayList/binaryDataArray/MS:1000515/binary");
			intensities_ds.getSpace().getSimpleExtentDims(&n);
			intensities[spectra[i].scan_start_time_index].resize(n);			
			intensities_ds.read(intensities[spectra[i].scan_start_time_index].data(), H5::DataType(H5::PredType::NATIVE_DOUBLE));
			
			loaded++;
			//cout << spectra[i].index << "," << spectra[i].preset_scan_configuration << "," << spectra[i].ms_level << "," << spectra[i].scan_start_time << "," << spectra[i].scan_start_time_index << ":" << n << endl;
		}
		catch (const H5::GroupIException& e) {}			
			
		if (loaded && (i == spectra.size()-1 ||
			spectra[i].preset_scan_configuration != spectra[i+1].preset_scan_configuration ||
			spectra[i].ms_level != spectra[i+1].ms_level))
		{				
			if (spectra[i].ms_level == 1)
			{
				// run seamass 2D
				seamass::process(id,
					spectra[i].preset_scan_configuration,
					instrument_type,
					scan_start_times, mzs, intensities,
					mz_res, mz_res,
					rt_res, rt_res,
					shrink, shrink,
					tol, tol,
					threads, out_type);
			}
			else
			{
				// run seamass 1D (todo)
			}
		}

	}


















	/*


	// this maps between our memory efficient spectrum ordering and spectra ordered by scan start time
	std::vector<size_t> sorted_indicies(ns);

	// data arrays passed on to seaMass













	// sort scan start times of all spectra in dataset
	vector< std::pair<double, size_t> > ssts(ns);
	vector<int> ms1s(ns);
	vector<double> precursors(ns);
	for (size_t i = 0; i < ns; i++)
	{
		H5::Group& s = ss.openGroup(ss.getObjnameByIdx(i));

		// read scan start times
		try
		{
			H5::Group& cv = s.openGroup("scanList/scan/0/cvParam/MS:1000016");
			Attribute& sst_a = cv.openAttribute("value");
			H5std_string sst_s("");
			sst_a.read(StrType(0, H5T_VARIABLE), sst_s);
			ssts[i] = pair<double, size_t>(atof(sst_s.c_str()),i);			
		}
		catch (const FileIException& e)
		{
			ssts[i] = pair<double, size_t>(-1.0,i);
		}

		// read if MS1 spectrum
		try
		{
			H5::Group& cv = s.openGroup("cvParam/MS:1000579");
			ms1s[i] = true;
		}
		catch (const FileIException& e)
		{
			ms1s[i] = false;
		} 
	}
	sort(ssts.begin(), ssts.end());
	// write scan_start_times
	for (size_t i = 0; i < ns; i++) scan_start_times[i] = ssts[i].first;
    
	// get indicies from our special sorted list to the scan start time order
	vector< std::pair<size_t, size_t> > sis(ns);
	for (size_t i = 0; i < ns; i++)
	{
		sis[i] = std::pair<size_t, size_t>(ssts[i].second, i);
	}
	std::sort(sis.begin(), sis.end());
	// write sorted_indicies
	for (size_t i = 0; i < ns; i++) sorted_indicies[i] = sis[i].second;
    


	for (size_t i = 0; i < scan_start_times.size(); i++)
	{
		cout << scan_start_times[i] << endl;
	}

	//seaMass = new seaMassRestoration(msd.run.id, instrument_type, debug);


	/*vector<hsize_t> ids(ns); // indices
	vector<double> ssts(ns); // start scan times
	vector<string> pscs(ns); // preset scan configurations
	vector<bool> ms1s(ns);   // is it an ms1 spectrum?
	vector<bool> pmss(ns);   // is it a profile mode spectrum?
	for (hsize_t i = 0; i < ns; i++)
	{
		ids[i] = i;



		cout << ids[i] << "," << ssts[i] << "," << pscs[i] << "," << ms1s[i] << "," << pmss[i] << endl;
	*/

	return 0;
}

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
#include <iomanip>
#include <H5Cpp.h>
#include <boost/program_options.hpp>
#include <pugixml.hpp>
#include <omp.h>

#include "seamass.hpp"
#include "NetCDFile.hpp"

#include "MSFileData.hpp"

using namespace std;
namespace po = boost::program_options;
namespace xml = pugi;

bool scan_start_time_order(const spectrum& lhs, const spectrum& rhs)
{
	return lhs.scan_start_time <= rhs.scan_start_time;
}


bool seamass_order(const spectrum& lhs, const spectrum& rhs)
{
	if (lhs.preset_config == rhs.preset_config)
	{
		if (lhs.precursor_mz == rhs.precursor_mz)
		{
			return lhs.scan_start_time <= rhs.scan_start_time;
		}
		else
		{
			return lhs.precursor_mz < rhs.precursor_mz;
		}
	}
	else
	{
		return lhs.preset_config < rhs.preset_config;
	}
}


int main(int argc, char *argv[])
{
	H5::Exception::dontPrint();
	seamass::notice();
	string in_file;
	int mz_res;
	int rt_res;
	int shrink;
	int tol;
	int out_type;
	int threads;
	double precoursorEvt;

	// *******************************************************************

	po::options_description general("Usage\n"
			"-----\n"
			"restoration [OPTIONS...] [MZMLB]\n"
			"restoration <-f in_file> <-m mz_res> <-r rt_res> <-s shrinkage> <-l tol> <-t threads> <-o out_type>\n"
			"restoration -m 1 -r 4 -s -4 -l -9 -t 4 -o 0");

	general.add_options()
		("help,h", "Produce help message")
		("file,f", po::value<string>(&in_file),
			"Raw input file in seaMass Input format (mzMLb, csv etc.) "
			"guidelines: Use pwiz-seamass to convert from mzML or vendor format")
		("mz_res,m",po::value<int>(&mz_res)->default_value(1),
			"MS resolution given as: \"b-splines per Th = 2^mz_res * 60 / 1.0033548378\" "
			"guidelines: between 0 to 1 for ToF (e.g. 1 is suitable for 30,000 resolution), 3 for Orbitrap, "
			"default: 1")
		("rt_res,r",po::value<int>(&rt_res)->default_value(4),
			"LC resolution given as: \"b-splines per minute = 2^rt_res\" "
			"guidelines: around 4, "
			"default: 4")
		("shrink,s",po::value<int>(&shrink)->default_value(-4),""
			"Amount of denoising given as: \"L1 shrinkage = 2^shrinkage\" "
			"guidelines: around -4, "
			"default: -4")
		("tol,l",po::value<int>(&tol)->default_value(-9),
			"Convergence tolerance, given as: \"gradient <= 2^tol\" "
			"guidelines: around -9, "
			"default: -9")
		("threads,t",po::value<int>(&threads)->default_value(4),
			"Number of OpenMP threads to use, "
			"guidelines: set to amount of CPU cores or 4, whichever is smaller, "
			"default: 4")
		("precursor,p",po::value<double>(&precoursorEvt),
			"Precursor value to single out Swath Window, "
			"guidelines: 0 = ms1 spectrum; > 0 = ms2 spectrum")
		("out_type,o",po::value<int>(&out_type)->default_value(0),
			"Type of output desired, "
			"guidelines: 0 = just viz_client input; 1 = also smo; 2 = also debug, "
			"default: 0");

	po::options_description desc;
	desc.add(general);

	try
	{
		po::positional_options_description pod;
		pod.add("file", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
		po::notify(vm);

		if(vm.count("help"))
		{
			cout<<desc<<endl;
			return 0;
		}
			if(vm.count("threads"))
		{
			threads=vm["threads"].as<int>();
		}
		else
		{
			threads=omp_get_max_threads();
		}
		if(vm.count("file"))
		{
			cout<<"Opening file: "<<vm["file"].as<string>()<<endl;
		}
		else
		{
			throw "Valid seamass input file was not give...";
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

	MassSpecFile* msFile = FileFactory::createFileObj(in_file);
    vector<spectrum> spectra = msFile->getSpectrum();
	unsigned long instrument_type = msFile->getInstrument();
	hsize_t ns = spectra.size();

    int lastdot = in_file.find_last_of(".");
    string id = (lastdot == string::npos) ? in_file : in_file.substr(0, lastdot);

	// determine start_scan_time order of spectra
	sort(spectra.begin(), spectra.end(), scan_start_time_order);
	for (size_t i = 0; i < spectra.size(); i++)
	{
		spectra[i].scan_start_time_index = i;
	}
	vector<double> scan_start_times(ns);
	for (size_t i = 0; i < spectra.size(); i++)
	{
		scan_start_times[i] = spectra[i].scan_start_time;
	}
	sort(spectra.begin(), spectra.end(), seamass_order);
	vector< std::vector<double> > mzs(ns-1);
	vector< std::vector<double> > intensities(ns-1);

	int loaded = 0;
    bool precursor_mz_is_constant = true;
	for (size_t i = 0; i < spectra.size(); i++)
	{
        if (loaded > 1 && spectra[i].precursor_mz != spectra[i-1].precursor_mz)
            precursor_mz_is_constant = false;

		if(spectra[i].count > 0 && spectra[i].index != spectra.size()-1)
		{
			msFile->getScanMZ(mzs[spectra[i].scan_start_time_index],spectra[i].index,spectra[i].count);
			msFile->getScanIntensities(intensities[spectra[i].scan_start_time_index],spectra[i].index,spectra[i].count);
		}

        loaded++;
		if ((i == spectra.size()-1 ||
			spectra[i].preset_config != spectra[i+1].preset_config ||
            (spectra[i].precursor_mz != spectra[i+1].precursor_mz && spectra[i].precursor_mz == 0.0)))
		{
            if (loaded > 1 && precursor_mz_is_constant)
			{
                ostringstream oss; oss << spectra[i].preset_config << "_" << spectra[i].precursor_mz;

				// run seamass 2D
				seamass::process(id,
					oss.str().c_str(),
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
            loaded = 0;
            precursor_mz_is_constant = true;
        }
	}

	delete msFile;
	return 0;
}

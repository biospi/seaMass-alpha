//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> bristol.ac.uk>
//
// Copyright (C) 2016  biospi Laboratory, University of Bristol, UK
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


#include <boost/program_options.hpp>
#include <omp.h>

#include "../io/MSFileData.hpp"
#include "../core/seaMass.hpp"


using namespace std;
namespace po = boost::program_options;


int main(int argc, char *argv[])
{
	seaMass::notice();

	string in_file;
	vector<int> res(2);
	int shrink;
	int tol;
	int threads;

	// *******************************************************************

	po::options_description general("Usage\n"
			"-----\n"
			"seamass [OPTIONS...] [MZMLB]\n"
			"seamass <-f in_file> <-m mz_res> <-r rt_res> <-s shrinkage> <-l tol> <-t threads> <-o out_type>\n"
			"seamass -m 1 -r 4 -s -4 -l -9 -t 4 -o 0");

	general.add_options()
		("help,h", "Produce help message")
		("file,f", po::value<string>(&in_file),
			"Raw input file in seaMass Input format (mzMLb, csv etc.) "
			"guidelines: Use pwiz-seamass to convert from mzML or vendor format")
		("mz_res,m",po::value<int>(&res[0])->default_value(1),
			"MS resolution given as: \"b-splines per Th = 2^mz_res * 60 / 1.0033548378\" "
			"guidelines: between 0 to 1 for ToF (e.g. 1 is suitable for 30,000 resolution), 3 for Orbitrap, "
			"default: 1")
		("rt_res,r",po::value<int>(&res[1])->default_value(4),
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
			"default: 4");

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
			throw "Valid seamass input file was not given...";
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

	mzMLbInputFile msFile(in_file);
	seaMass::Input input;
	while (msFile.next(input))
	{
		cout << input.spectrum_index.size() << ":" << input.bin_counts.size() << ":" << input.bin_locations.size() << endl;
		seaMass sm(input, res, shrink, tol);
		do
		{
			// if debug save 
		}
		while (sm.next());
		// save 
	}

	return 0;
}

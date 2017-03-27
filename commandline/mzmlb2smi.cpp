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


#include <limits>
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>
//#include <omp.h>

#include "../io/NetcdfWriter.hpp"
#include "../io/DatasetMzmlb.hpp"

using namespace std;
namespace po = boost::program_options;


int main(int argc, char *argv[])
{
#ifndef NDEBUG
	try
#endif
	{
		string fileName;
		int debugLevel;

		po::options_description general("Usage\n"
												"-----\n"
												"mzmlb2smi [OPTIONS...] [MZMLB]\n"
												"mzmlb2smi <-f inFile>\n");

		general.add_options()
				("help,h", "Produce help message")
				("file,f", po::value<string>(&fileName),
				 "input file in mzMLb format"
						 "guidelines: Use pwiz-mzmlb (https://github.com/biospi/mzmlb) to convert from mzML or vendor format")
				("debug_level,d", po::value<int>(&debugLevel)->default_value(0),
				 "debug level"
						 "guidelines: 1+ for convergence stats, 2+ for performance stats, 3+ to write intermediate iterations to disk, 4 for all math"
						 "default: 0");


		po::options_description desc;
		desc.add(general);

		po::positional_options_description pod;
		pod.add("file", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
		po::notify(vm);

		if(vm.count("help"))
		{
			cout << desc << endl;
			return 0;
		}
		if(!vm.count("file"))
		{
			throw runtime_error("Error: Valid mzMLb input file was not given");
		}

		setDebugLevel(debugLevel);

        if (getDebugLevel() % 10 == 0) cout << "Reading " << fileName << endl;
        DatasetMzmlb msFile(fileName);
		SeamassCore::Input input;
		string id;
		for (int i = 0; msFile.next(input, id); i++)
		{
			string smiFileName = boost::filesystem::change_extension(fileName, "").string() + "." + id + ".smi";
			if (getDebugLevel() % 10 == 0) cout << "Writing " << smiFileName << endl;
			NetcdfWriter netcdfWriter(smiFileName);
			netcdfWriter.writeSmi(input);
		}
	}
#ifndef NDEBUG
	catch(exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}
#endif

	return 0;
}

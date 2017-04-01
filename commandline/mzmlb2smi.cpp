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


int main(int argc, const char * const * argv)
{
#ifdef NDEBUG
	try
#endif
	{
		string filePath;
		int debugLevel;

        po::options_description general(
            "Usage\n"
            "-----\n"
            "mzmlb2smi [OPTIONS...] [MZMLB FILE]\n"
            "mzmlb2smi <file>"
        );

		general.add_options()
            ("help,h", "Produce help message")
			("file,f", po::value<string>(&filePath),
             "Input file in mzMLb format. Use pwiz-mzmlb (https://github.com/biospi/mzmlb) to convert from mzML or vendor format.")
			("debug_level,d", po::value<int>(&debugLevel)->default_value(0),
             "Debug level. Use 1+ for stats on DIA output, 2+ for all output, 3+ for stats on input spectra, ")
        ;

		po::options_description desc;
		desc.add(general);

		po::positional_options_description pod;
		pod.add("file", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
		po::notify(vm);

        if(vm.count("help") || !vm.count("file")) {
            cout << desc << endl;
            return 0;
        }

		setDebugLevel(debugLevel);

        DatasetMzmlb msFile(filePath);
        string fileName = boost::filesystem::path(filePath).stem().string();
        bool success = true;
		for (int i = 0; success; i++)
		{
            SeamassCore::Input input;
            string id;
            success = msFile.next(input, id);

            if (success)
            {
                string smiFileName = fileName + "." + id + ".smi";
                if (input.binCountsIndex.size() >= 1000 && getDebugLevel() % 10 >= 1)
                    cout << getTimeStamp() << "  Writing " << smiFileName << " ..." << endl;
                NetcdfWriter netcdfWriter(smiFileName);
                netcdfWriter.writeSmi(input);
            }
		}
	}
#ifdef NDEBUG
	catch(exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}
#endif

	return 0;
}

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
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>
#include <kernel.hpp>
#include "../core/DatasetSeamass.hpp"
using namespace std;
using namespace kernel;
namespace po = boost::program_options;


int main(int argc, const char * const * argv)
{
#ifdef NDEBUG
    try
#endif
    {
        string filePathIn;
        int debugLevel;

        // *******************************************************************

        po::options_description general(
            "Usage\n"
            "-----\n"
            "Restore an mzMLb/smv file to mzMLb/smb with optional resampling or centroiding.\n"
            "\n"
            "seamass-restore [OPTIONS...] <file>\n"
        );

        general.add_options()
            ("help,h",
             "Produce this help message")
            ("file,f", po::value<string>(&filePathIn),
             "Input file in mzMLv or smv format produced by 'seamass'.")
            ("debug,d", po::value<int>(&debugLevel)->default_value(0),
             "Debug level.")
        ;

        po::options_description desc;
        desc.add(general);

        po::positional_options_description pod;
        pod.add("file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
        po::notify(vm);

        cout << endl;
        cout << "seamass-restore : Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
        cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
        cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
        cout << endl;
        initKernel(debugLevel);

        if (vm.count("help") || !vm.count("file"))
        {
            cout << desc << endl;
            return 0;
        }

        string fileStemOut = boost::filesystem::path(filePathIn).stem().string();
        Dataset* dataset = FileFactory::createFileObj(filePathIn, fileStemOut, Dataset::WriteType::Input);
        if (!dataset)
            throw runtime_error("ERROR: Input file is missing or incorrect");

        Seamass::Input input;
        Seamass::Output output;
        string id;

        while (dataset->read(input, output, id))
        {
            if (debugLevel % 10 == 0)
                cout << "Processing " << id << endl;

            // load back into Seamass
            Seamass seamassCore(input, output);

            seamassCore.getOutputBinCounts(input.counts);

            dataset->write(input, id);

            if (debugLevel % 10 == 0)
                cout << endl;
        }

        delete dataset;
        cout << endl;
    }
#ifdef NDEBUG
    catch(exception& e)
    {
        cerr << e.what() << endl;
        cout << endl;
        return 1;
    }
#endif
    return 0;
}

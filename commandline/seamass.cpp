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
#include "../io/DatasetSeamass.hpp"
using namespace std;
using namespace kernel;
namespace po = boost::program_options;


int main(int argc, const char * const * argv)
{
#ifdef NDEBUG
	try
#endif
	{
        string fileName;
        vector<char> scale(2);
        int shrinkageExponent;
        int toleranceExponent;
        int debugLevel;

        // *******************************************************************

        po::options_description general(
            "Usage\n"
            "-----\n"
            "seamass [OPTIONS...] [MZMLB FILE]\n"
            "seamass <-m mz_scale> <-s st_scale> <-l lambda> <-t tol> <file>"
        );

        general.add_options()
            ("help,h",
             "Produce this help message")
            ("file,f", po::value<string>(&fileName),
             "Input file in mzMLb or binned smb format. Use pwiz-mzmlb (https://github.com/biospi/mzmlb) to convert from mzML/vendor format to mzMLb.")
            ("mz_scale,m", po::value<char>(&scale[0]),
             "Output m/z resolution given as \"b-splines per Th = 2^mz_scale * 60 / 1.0033548378\". "
             "Use 0 or 1 for ToF (e.g. 1 is suitable for 30,000 resolution), 3 for Orbitrap. "
             "Default is to autodetect.")
            ("st_scale,s", po::value<char>(&scale[1]),
             "output scan-time resolution given as \"b-splines per second = 2^st_scale\". Use around 4. "
             "Default is to autodetect.")
            ("lambda,l", po::value<int>(&shrinkageExponent)->default_value(0),
             "Amount of denoising given as \"L1 lambda = 2^shrinkage\". Use around 0.")
            ("tol,t", po::value<int>(&toleranceExponent)->default_value(-10),
             "Convergence tolerance, given as \"gradient <= 2^tol\". Use around -10.")
            ("debug,d", po::value<int>(&debugLevel)->default_value(0),
             "Debug level. Use 1+ for convergence stats, 2+ for performance stats, 3+ for sparsity info, "
             "4 to output all maths, +10 to write intermediate results to disk.")
        ;

        po::options_description desc;
        desc.add(general);

        po::positional_options_description pod;
		pod.add("file", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
		po::notify(vm);

        cout << endl;
        Seamass::notice();
        initKernel(debugLevel);

        if(vm.count("help") || !vm.count("file"))
		{
			cout << desc << endl;
			return 0;
		}
        cout << endl;

        if(!vm.count("mz_scale"))
            scale[0] = numeric_limits<char>::max();

        if(!vm.count("st_scale"))
            scale[1] = numeric_limits<char>::max();

       Dataset* dataset = FileFactory::createFileObj(fileName);
        if (!dataset)
            throw runtime_error("ERROR: Input file is missing or incorrect");

        Seamass::Input input;
        string id;
        double tolerance = pow(2.0, (double)toleranceExponent);
        double shrinkage = pow(2.0, (double)shrinkageExponent);

        while (dataset->read(input, id))
        {
            if (getDebugLevel() % 10 == 0)
                cout << "Processing " << id << endl;

            Seamass seamassCore(input, scale, shrinkage, tolerance);

            do
            {
                if (getDebugLevel() >= 10)
                {
                    Seamass::Output output;
                    seamassCore.getOutput(output);

                    // write intermediate output in seaMass format
                    DatasetSeamass datasetOut(fileName, true);
                    ostringstream oss; oss << id << "_" << setfill('0') << setw(4) << seamassCore.getIteration();
                    datasetOut.write(input, output, oss.str());
                }
            }
            while (seamassCore.step());

            // write output
            Seamass::Output output;
            seamassCore.getOutput(output);
            dataset->write(input, output, id);

            if (getDebugLevel() % 10 == 0)
                cout << endl;
        }

        delete dataset;
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

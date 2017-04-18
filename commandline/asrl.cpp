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


#include "../asrl/Asrl.hpp"
#include "../kernel/FileNetcdf.hpp"
#include <boost/filesystem/convenience.hpp>
#include <boost/program_options.hpp>
#include <iomanip>
using namespace std;
using namespace kernel;
namespace po = boost::program_options;


int main(int argc, const char * const * argv)
{
#ifdef NDEBUG
	try
#endif
	{
        po::options_description general(
                "Usage\n"
                "-----\n"
                "asrl [OPTIONS...] [SAI FILE]\n"
                "asrl <-l lambda> <-no_taper> <-t tol> <file>"
        );

        string filePath;
        int shrinkageExponent;
        int toleranceExponent;
        int debugLevel;
        bool noTaperLambda;

        general.add_options()
                ("help,h",
                 "Produce this help message")
                ("file,f", po::value<string>(&filePath),
                 "HDF5 or NetCDF4 input file in SAI format.")
                ("lambda,l", po::value<int>(&shrinkageExponent)->default_value(0),
                 "Amount of denoising given as \"L1 lambda = 2^shrinkage\". Use around 0.")
                ("no_taper", po::bool_switch(&noTaperLambda)->default_value(false),
                 "Use this to stop tapering of lambda to 0 before finishing.")
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
        Asrl::notice();
        initKernel(debugLevel);

		if(vm.count("help") || !vm.count("file"))
		{
			cout << desc << endl;
			return 0;
		}

        double tolerance = pow(2.0, (double)toleranceExponent);
        double shrinkage = pow(2.0, (double)shrinkageExponent);

        // read input file
        if (getDebugLevel() % 10 >= 1)
            cout << getTimeStamp() << "Reading " << filePath << " ..." << endl;
        Asrl::Input input;
        {
            FileNetcdf fileIn(filePath);

            vector<ii> aExtent(2);
            fileIn.read_AttNC("extent", fileIn.read_VarIDNC("A_v"), aExtent);
            input.aM = aExtent[0]; input.aN = aExtent[1];
            fileIn.read_VecNC("A_v", input.aVs);
            fileIn.read_VecNC("A_i", input.aIs);
            fileIn.read_VecNC("A_j", input.aJs);

            fileIn.read_VecNC("b", input.bs);

            try
            {
                fileIn.read_VecNC("x", input.xs);
            }
            catch(runtime_error& e) {}

            input.gN = 0;
            try
            {
                vector<ii> gExtent(2);
                fileIn.read_AttNC("extent", fileIn.read_VarIDNC("G_v"), gExtent);
                input.gM = gExtent[0]; input.gN = gExtent[1];
           }
            catch(runtime_error& e) {}
            if (input.gM > 0 && input.gN > 0)
            {
                fileIn.read_VecNC("G_v", input.gVs);
                fileIn.read_VecNC("G_i", input.gIs);
                fileIn.read_VecNC("G_j", input.gJs);
            }
        }

        // optimise!
        Asrl asrl(input, shrinkage, !noTaperLambda, tolerance);
        do
        {
            if (getDebugLevel() >= 10)
            {
                // write intermediate results
                Asrl::Output output;
                asrl.getOutput(output);
                if (output.xs.size() > 0)
                {
                    ostringstream oss;
                    oss << "." << setfill('0') << setw(4) << asrl.getIteration() << ".sao";
                    string fileName = boost::filesystem::path(filePath).stem().string() + oss.str();

                    if (getDebugLevel() % 10 >= 1)
                        cout << getTimeStamp() << "  Writing " << fileName << " ..." << endl;

                    FileNetcdf fileOut(fileName, NC_NETCDF4);
                    fileOut.write_VecNC("x", output.xs, sizeof(output.xs[0]) == 4 ? NC_FLOAT : NC_DOUBLE);
                    fileOut.write_VecNC("Ax", output.aXs, sizeof(output.aXs[0]) == 4 ? NC_FLOAT : NC_DOUBLE);
                    if (output.gXs.size() > 0)
                        fileOut.write_VecNC("Gx", output.gXs, sizeof(output.gXs[0]) == 4 ? NC_FLOAT : NC_DOUBLE);
                }
            }
        }
        while (asrl.step());

        Asrl::Output output;
        asrl.getOutput(output);
        if (output.xs.size() > 0)
        {
            // open and write output file
            string fileName = boost::filesystem::path(filePath).stem().string() + ".sao";

            if (getDebugLevel() % 10 >= 1)
                cout << getTimeStamp() << "Writing " << fileName << " ..." << endl;

            FileNetcdf fileOut(fileName, NC_NETCDF4);
            fileOut.write_VecNC("x", output.xs, sizeof(output.xs[0]) == 4 ? NC_FLOAT : NC_DOUBLE);
            fileOut.write_VecNC("Ax", output.aXs, sizeof(output.aXs[0]) == 4 ? NC_FLOAT : NC_DOUBLE);
            if (output.gXs.size() > 0)
                fileOut.write_VecNC("Gx", output.gXs, sizeof(output.gXs[0]) == 4 ? NC_FLOAT : NC_DOUBLE);
        }

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

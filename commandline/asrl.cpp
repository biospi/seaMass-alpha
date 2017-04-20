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
        int lambdaExponent;
        int lambdaGroupExponent;
        int toleranceExponent;
        int debugLevel;
        bool noTaperLambda;

        general.add_options()
                ("help,h",
                 "Produce this help message")
                ("file,f", po::value<string>(&filePath),
                 "HDF5 or NetCDF4 input file in SAI format.")
                ("lambda,l", po::value<int>(&lambdaExponent)->default_value(0),
                 "Amount of individual lambda given as \"L1_lambda = 2^lambda\". Use around 0.")
                ("group_lambda,g", po::value<int>(&lambdaGroupExponent)->default_value(0),
                 "Amount of group lambda given as \"L2_group_lambda = 2^lambda_group\". "
                 "Ignored if no groups are specified in the input. Use around 0.")
                ("no_taper", po::bool_switch(&noTaperLambda)->default_value(false),
                 "Use this to stop tapering of lambda to 0 before finishing.")
                ("tol,t", po::value<int>(&toleranceExponent)->default_value(-15),
                 "Convergence tolerance, given as \"gradient <= 2^tol\". Use around -15.")
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

        // read input file
        if (getDebugLevel() % 10 >= 1)
            cout << getTimeStamp() << "Reading " << filePath << " ..." << endl;

        FileNetcdf fileIn(filePath);

        Asrl::Input input;

        input.aT.resize(1);
        fileIn.read(input.aT[0], "At");

        input.bT.resize(1);
        fileIn.read(input.bT[0], "Bt");

        try
        {
            input.xT.resize(1);
            fileIn.read(input.xT[0], "Xt");
        }
        catch(runtime_error& e)
        {
            input.xT.resize(0);
        }

        try
        {
            input.gT.resize(1);
            fileIn.read(input.gT[0], "Gt");
        }
        catch(runtime_error& e)
        {
            input.gT.resize(0);
        }

        double tolerance = pow(2.0, (double)toleranceExponent);
        double lambda = pow(2.0, (double)lambdaExponent);
        double lambdaGroup = input.gT.size() > 0 ? lambdaGroup = pow(2.0, (double)lambdaGroupExponent) : 0.0;

        string fileStemOut = boost::filesystem::path(filePath).stem().string();

        // optimise!
        Asrl asrl(input, lambda, lambdaGroup, !noTaperLambda, tolerance);
        do
        {
            if (getDebugLevel() >= 10)
            {
                // write intermediate results
                Asrl::Output output;
                asrl.getOutput(output);

                ostringstream oss;
                oss << "." << setfill('0') << setw(4) << asrl.getIteration() << ".sao";
                string fileNameOut = fileStemOut + oss.str();

                if (getDebugLevel() % 10 >= 1)
                    cout << getTimeStamp() << "  Writing " << fileNameOut << " ..." << endl;

                FileNetcdf fileOut(fileNameOut, NC_NETCDF4);
                fileOut.write(output.xT[0], "Xt");
                fileOut.write(output.aTxT[0], "AtXt");
                if (output.gTxT.size() > 0)
                    fileOut.write(output.gTxT[0], "GtXt");
            }
        }
        while (asrl.step());

        Asrl::Output output;
        asrl.getOutput(output);
        ostringstream oss;

        string fileNameOut = fileStemOut + ".sao";

        if (getDebugLevel() % 10 >= 1)
            cout << getTimeStamp() << "  Writing " << fileNameOut << " ..." << endl;

        FileNetcdf fileOut(fileNameOut, NC_NETCDF4);
        fileOut.write(output.xT[0], "Xt");
        fileOut.write(output.aTxT[0], "AtXt");
        if (output.gTxT.size() > 0)
            fileOut.write(output.gTxT[0], "GtXt");

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

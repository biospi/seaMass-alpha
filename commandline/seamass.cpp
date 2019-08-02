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


#include "../kernel/Subject.hpp"
#include "../core/DatasetSeamass.hpp"
#include <kernel.hpp>
#include <limits>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>
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
        string libPathIn;
        int scaleX;
        int scaleY;
        int scaleZ;
        int lambdaExponent;
        int lambdaGroupExponent;
        bool noTaperLambda;
        int toleranceExponent;
        short chargeStates;
        double peakFwhm;
        int debugLevel;

        // *******************************************************************

        po::options_description general(
            "Usage\n"
            "-----\n"
            "seamass [OPTIONS...] [INPUT FILE]\n"
            "seamass <-l lambda> <file>"
        );

        general.add_options()
            ("help,h",
             "Produce this help message")
            ("file,f", po::value<string>(&filePathIn),
             "Input file in sml format.")
            ("lib,l", po::value<string>(&libPathIn)->default_value(""),
             "OPTIONAL: Spectral library in sml format to use for library deconvolution and identification. ")
            ("lambda,l", po::value<int>(&lambdaExponent)->default_value(0),
             "OPTIONAL: Amount of denoising given as \"L1 lambda = 2^shrinkage\", default 0.")
            ("group_lambda,g", po::value<int>(&lambdaGroupExponent)->default_value(0),
             "OPTIONAL: Amount of group lambda given as \"L2_group_lambda = 2^lambda_group\", default 0. "
             "Ignored if no groups are specified in the spectral library file.")
            ("OPTIONAL: no_taper", po::bool_switch(&noTaperLambda)->default_value(false),
             "Use this to stop tapering of lambda to 0 before finishing.")
            ("x", po::value<int>(&scaleX),
             "OPTIONAL: Output resolution for the x dimension (i.e. m/z for LC-MS data).")
            ("y", po::value<int>(&scaleY),
             "OPTIONAL: Output resolution for the y dimension (i.e. scantime for LC-MS data).")
            ("z", po::value<int>(&scaleZ),
             "OPTIONAL: Output resolution for the z dimension.")
            ("charges,c", po::value<short>(&chargeStates)->default_value(0),
             "OPTIONAL: Highest charge state to deconvolute. Default is no deconvolution.")
            ("fwhm,w", po::value<double>(&peakFwhm)->default_value(0.0),
             "OPTIONAL: FWHM for deconcolution of peak spread. Default is no deconvolution.")
            ("tol,t", po::value<int>(&toleranceExponent)->default_value(-10),
             "OPTIONAL: Convergence tolerance, given as \"gradient <= 2^tol\". Use around -10.")
            ("debug,d", po::value<int>(&debugLevel)->default_value(0),
             "OPTIONAL: Debug level. Use 1+ for convergence stats, 2+ for performance stats, 3+ for sparsity info, "
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
        cout << endl;
        initKernel(debugLevel);

        Subject::setDebugLevel(debugLevel);
        Observer* observer = 0;
        if (debugLevel % 10 >= 1)
            Subject::registerObserver(observer = new Observer());

        ObserverMatrix* observerMatrix = 0;
        ObserverMatrixSparse* observerMatrixSparse = 0;
        if (debugLevel / 10 >= 2)
        {
            SubjectMatrix::registerObserver(observerMatrix = new ObserverMatrix());
            SubjectMatrixSparse::registerObserver(observerMatrixSparse = new ObserverMatrixSparse());
        }

        if(vm.count("help") || !vm.count("file"))
        {
            cout << desc << endl;
            return 0;
        }

        vector<short> scale(3);

        if(vm.count("x"))
            scale[0] = short(scaleX);
        else
            scale[0] = numeric_limits<short>::max();

        if(vm.count("y"))
            scale[1] = short(scaleY);
        else
            scale[1] = numeric_limits<short>::max();

        if(vm.count("z"))
            scale[2] = short(scaleZ);
        else
            scale[2] = numeric_limits<short>::max();
        
       string fileStemOut = boost::filesystem::path(filePathIn).stem().string();
        Dataset* dataset = FileFactory::createFileObj(filePathIn, fileStemOut, Dataset::WriteType::InputOutput);
        if (!dataset)
            throw runtime_error("ERROR: Input file is missing or incorrect");

        Seamass::Input input;
        string id;
        fp tolerance = pow(2.0, fp(toleranceExponent));
        fp lambda = pow(2.0, fp(lambdaExponent));
        fp lambdaGroup = pow(2.0, fp(lambdaGroupExponent));

        while (dataset->read(input, id))
        {
            if (debugLevel % 10 == 0)
                cout << "Processing " << id << endl;

            Seamass seamass(input, libPathIn, scale, lambda, lambdaGroup, !noTaperLambda, tolerance,
                            peakFwhm, chargeStates);

            if (debugLevel / 10 >= 1)
            {
                /*Seamass::Input input2;
                input2.countsIndex = input.countsIndex;
                input2.startTimes = input.startTimes;
                input2.finishTimes = input.finishTimes;
                input2.exposures = input.exposures;
                input2.type = input.type;
                seamass.getInput(input2);

                // write input in seaMass format
                ostringstream oss; oss << fileStemOut << ".input";
                DatasetSeamass datasetOut("", oss.str(), Dataset::WriteType::Input);
                datasetOut.write(input2, id);*/
            }

            do
            {
                if (debugLevel / 10 >= 2)
                {
                    Seamass::Output output;
                    seamass.getOutput(output, true);

                    // write intermediate output in seaMass format
                    ostringstream oss; oss << fileStemOut << ".synthesized." << setfill('0') << setw(4) << seamass.getIteration();
                    DatasetSeamass datasetOut("", oss.str(), Dataset::WriteType::InputOutput);
                    datasetOut.write(input, output, id);
                }
            }
            while (seamass.step());

            // write output
            {
                Seamass::Output output;
                seamass.getOutput(output, false);
                dataset->write(input, output, id);
            }

            if (debugLevel / 10 >= 1)
            {
                Seamass::Output output;
                seamass.getOutput(output, true);

                ostringstream oss; oss << fileStemOut << ".synthesized";
                DatasetSeamass datasetOut("", oss.str(), Dataset::WriteType::InputOutput);
                datasetOut.write(input, output, id);
            }

            if (debugLevel % 10 == 0)
                cout << endl;
        }

        delete dataset;
        if (observer) delete observer;
        if (observerMatrix) delete observerMatrix;
        if (observerMatrixSparse) delete observerMatrixSparse;
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

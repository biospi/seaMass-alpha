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
        string dbFilename;
        int scaleMz;
        int scaleSt;
        double lambdaExp;
        double lambdaModExp;
        bool noTaperLambda;
        bool noUnknowns;
        double toleranceExp;
        double peakFwhm;
        short chargeStates;
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
            ("file,f", po::value<string>(&filePathIn),
                "Input file in mzMLb or binned smb format. Use pwiz-mzmlb (https://github.com/biospi/mzmlb) to convert "
                "from mzML/vendor format to mzMLb.")
            ("db,b", po::value<string>(&dbFilename),
                "Spectral library database in smd format. generate_unknowns will generate a generic db.")
            ("mz_scale,m", po::value<int>(&scaleMz),
                "Output mz resolution given as \"2^mz_scale * log2(mz - polarity*1.007276466879)\". "
                "Default is to autodetect.")
            ("st_scale,s", po::value<int>(&scaleSt),
                "output scantime resolution given as \"2^st_scale\"."
                "Default is to autodetect.")
            ("lambda,l", po::value<double>(&lambdaExp)->default_value(0),
                "Amount of shrinkage (denoising) given as \"L1 shrinkage = 2^lambda\"."
                "Use around 0.")
            ("lambda_mod,lm", po::value<double>(&lambdaModExp)->default_value(0),
                "modifier for amount of shrinkage of the unknowns given as \"L1 shrinkage_db = 2^(lambda+lambda_mod)\"."
                "For low resolution data, this should be >0, otherwise 0 is fine.")
            ("no_taper", po::bool_switch(&noTaperLambda)->default_value(false),
                "Use this to stop tapering of lambda to 0 before finishing.")
            ("no_unknowns", po::bool_switch(&noUnknowns)->default_value(false),
                "Use this to turn off fitting of unknowns.")
            ("tol,t", po::value<double>(&toleranceExp)->default_value(-10),
                "Convergence tolerance, given as \"gradient <= 2^tol\". Use around -10.")
            ("charge_states,c", po::value<short>(&chargeStates)->default_value(0),
                "Highest charge state to deconvolute. Default is no charge state deconvolution. Currently unimplemented.")
            ("fwhm,w", po::value<double>(&peakFwhm)->default_value(0.0),
                "Full width at half maximum of peak width. Currently unimplemented.")
            ("debug,d", po::value<int>(&debugLevel)->default_value(0),
                "Debug level. Use 1 for convergence stats, 2 for performance stats, 3 for sparsity info, "
                "4 to output all maths, +10 to write synthesised output, +20 to also write intermediate results to disk.")
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

        vector<short> scale(2);

        if(vm.count("mz_scale"))
            scale[0] = short(scaleMz);
        else
            scale[0] = numeric_limits<short>::max();

        if(vm.count("st_scale"))
            scale[1] = short(scaleSt);
        else
            scale[1] = numeric_limits<short>::max();

        string fileStemOut = boost::filesystem::path(filePathIn).stem().string();
        Dataset* dataset = FileFactory::createFileObj(filePathIn, fileStemOut, Dataset::WriteType::InputOutput);
        if (!dataset)
            throw runtime_error("ERROR: Input file is missing or incorrect");

        Seamass::Input input;
        string id;
        fp tolerance = pow(2.0, fp(toleranceExp));
        fp lambda = pow(2.0, fp(lambdaExp));
        fp lambdaScale = pow(2.0, fp(lambdaModExp));

        while (dataset->read(input, id))
        {
            if (debugLevel % 10 == 0)
                cout << "Processing " << id << endl;

            Seamass seamass(input, dbFilename, scale, lambda, lambdaScale, noTaperLambda, noUnknowns, peakFwhm, chargeStates, tolerance);

            do
            {
                if (debugLevel / 10 >= 2)
                {
                    {
                        Seamass::Output output;
                        seamass.getOutput(output, false);

                        // write intermediate output in seaMass format
                        ostringstream oss; oss << fileStemOut << "." << setfill('0') << setw(4) << seamass.getIteration();
                        DatasetSeamass datasetOut("", oss.str(), Dataset::WriteType::InputOutput);
                        datasetOut.write(input, output, id);
                    }
                    {
                        Seamass::Output output;
                        seamass.getOutput(output, true);

                        // write intermediate output in seaMass format
                        ostringstream oss; oss << fileStemOut << ".synthesized." << setfill('0') << setw(4) << seamass.getIteration();
                        DatasetSeamass datasetOut("", oss.str(), Dataset::WriteType::InputOutput);
                        datasetOut.write(input, output, id);
                    }
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
                {
                    Seamass::Output output;
                    seamass.getOutput(output, true, seamass.bases_mask_library_);

                    ostringstream oss; oss << fileStemOut << ".synthesized.library";
                    DatasetSeamass datasetOut("", oss.str(), Dataset::WriteType::InputOutput);
                    datasetOut.write(input, output, id);                
                }
                {
                    Seamass::Output output;
                    seamass.getOutput(output, true, seamass.bases_mask_unknowns_);

                    ostringstream oss; oss << fileStemOut << ".synthesized.unknowns";
                    DatasetSeamass datasetOut("", oss.str(), Dataset::WriteType::InputOutput);
                    datasetOut.write(input, output, id);
                }
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

//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> bristol.ac.uk>
//
// Copyright (C) 2018  biospi Laboratory, University of Bristol, UK
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
#include "../io/FileNetcdf.hpp"
#include "../core/DatasetSeamass.hpp"
#include <kernel.hpp>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <netcdf.h>
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
                        "driams_raw2smb [OPTIONS...] [txt FILE]\n"
                        "driams_raw2smb <file>"
        );

        general.add_options()
            ("help,h", "Produce help message")
            ("file,f", po::value<string>(&filePathIn),
                "Input file in driams format.")
            ("debug,d", po::value<int>(&debugLevel)->default_value(0),
                "Debug level. Use 1+ for stats on DIA output, 2+ for all output, 3+ for stats on input spectra.")
            ;

        po::options_description desc;
        desc.add(general);

        po::positional_options_description pod;
        pod.add("file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
        po::notify(vm);

        if(vm.count("help") || !vm.count("file"))
        {
            cout << endl;
            cout << "driams_raw2smb : Copyright (C) 2025 - biospi Laboratory, University of Bristol, UK" << endl;
            cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
            cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
            cout << endl;

            cout << desc << endl;
            return 0;
        }


        // read input spectrum
        typedef boost::tokenizer< boost::escaped_list_separator<char> > so_tokenizer;
        ifstream raw(filePathIn, ios_base::in);
        std::vector<double> mzs;
        std::vector<double> intensities;
        string line;
        do
        {
            boost::trim(line);
            if (line.size() == 0 || line[0] == '#')
                continue;

            so_tokenizer tok(line, boost::escaped_list_separator<char>("", " \t", "\"\'"));

            so_tokenizer::iterator toki = tok.begin();
            double mz = atof(toki->c_str());
            ++toki;
            double intensity = atof(toki->c_str());

            mzs.push_back(mz);
            intensities.push_back(intensity);            
        }
        while(getline(raw, line));


        // write output spectrum
        Subject::setDebugLevel(debugLevel);
        Observer* observer = 0;
        if (debugLevel % 10 >= 1)
            Subject::registerObserver(observer = new Observer());

        ObserverMatrix* observerMatrix = 0;
        ObserverMatrixSparse* observerMatrixSparse = 0;
        if (debugLevel / 10 >= 1)
        {
            SubjectMatrix::registerObserver(observerMatrix = new ObserverMatrix());
            SubjectMatrixSparse::registerObserver(observerMatrixSparse = new ObserverMatrixSparse());
        }

        if (vm.count("help") || !vm.count("file") || intensities.size() == 0)
        {
            cout << desc << endl;
            return 0;
        }

        Seamass::Input out;
        out.type = Seamass::Input::Type::Binned;
        out.polarity = 1;

        // use all the intensities except first and last
        for (ii k = 1; k < (ii)intensities.size() - 1; k++)
        {
            if (intensities[k] != 0.0 || intensities[k - 1] != 0.0) // merge zeros
            {
                out.locations.push_back(0.5 * (mzs[k - 1] + mzs[k]));
                out.counts.push_back(((fp)intensities[k]));
            }
        }
        out.locations.push_back(0.5 * (mzs[mzs.size() - 2] + mzs.back()));

        string fileStemOut = boost::filesystem::path(filePathIn).stem().string();
        DatasetSeamass datasetOut("", fileStemOut, Dataset::WriteType::Input);
        datasetOut.write(out, "");

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

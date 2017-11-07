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


#include "../core/DatasetMzmlb.hpp"
#include "../core/DatasetSeamass.hpp"
#include <kernel.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/program_options.hpp>
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
        string dirPathIn;
        int debugLevel;

        po::options_description general(
            "Usage\n"
            "-----\n"
            "Converts a set of smb files back to an mzMLb file, given the original mzMLb file.\n"
            "\n"
            "smb2mzmlb [OPTIONS...] <file>\n"
        );

        general.add_options()
            ("help,h", "Produce help message")
            ("file,f", po::value<string>(&filePathIn),
             "Input file in mzMLb format. Use pwiz-mzmlb (https://github.com/biospi/mzmlb) to convert from mzML or vendor format.")
            ("dir,i", po::value<string>(&dirPathIn)->default_value("."),
             "Input directory containing files in smb format [default=.]")
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

        cout << endl;
        cout << "smb2mzmlb : Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
        cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
        cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
        cout << endl;
        initKernel(debugLevel);

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

        if(vm.count("help") || !vm.count("file"))
        {
            cout << desc << endl;
            return 0;
        }

        boost::filesystem::path smbPathStem = boost::filesystem::path(dirPathIn) / boost::filesystem::path(filePathIn).stem();

        boost::filesystem::path fileNameOut = boost::filesystem::path(filePathIn).filename();
        // Need to fix this it seems that boost filesystem doesnt like some special charactors! .e.g. () {} %
        //if (equivalent(fileNameOut, filePathIn))
        //    throw runtime_error("ERROR: Make sure the input mzMLb file is not in the working directory.");

        string eek = fileNameOut.replace_extension("").string();
        DatasetMzmlb datasetMzmlb(filePathIn, eek, Dataset::WriteType::Input);

        Seamass::Input input;
        string id;
        int injected = 0;
        while(datasetMzmlb.read(input, id))
        {
            DatasetSeamass* datasetSeamass;
            string smbPathFile = smbPathStem.string() + "." + id + ".smb";
            try
            {
                datasetSeamass = new DatasetSeamass(smbPathFile, "");
            }
            catch (runtime_error r)
            {
                datasetSeamass = 0;
            }

            // replace input with smb file input if available
            if (datasetSeamass)
            {
                string nullId;
                datasetSeamass->read(input, nullId);

                 if (debugLevel % 10 >= 1)
                     cout << getTimeStamp() << " ";
                 else
                     cout << endl;
                cout << " Injected " << smbPathFile << endl;

                injected++;
           }

            datasetMzmlb.write(input, id);
        }

        if (injected == 0)
        {
            cerr << "WARNING: No smb files injected" << endl;
        }
        else
        {
            if (debugLevel % 10 >= 1)
                cout << getTimeStamp() << " ";
            else
                cout << endl;
            cout <<  "Injected " << injected << " smb file" << (injected > 0 ? "s" : "") << endl;
        }

        cout << endl;
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

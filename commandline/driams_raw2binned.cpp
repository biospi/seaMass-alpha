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
#include "../core/Bspline.hpp"
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


double PROTON_MASS = 1.007276466879;


int main(int argc, const char * const * argv)
{
#ifdef NDEBUG
    try
#endif
    {
        string filePathIn;
        int mzScale;
        double mzMin;
        double mzMax;
        float energy = -1.0;

        // *******************************************************************

        po::options_description general(
                "Usage\n"
                        "-----\n"
                        "driams_raw2binned [OPTIONS...] [txt FILE]\n"
                        "driams_raw2binned <-m mz_scale> <file>"
        );

        general.add_options()
                ("help,h",
                 "Produce this help message")
                ("file,f", po::value<string>(&filePathIn),
                 "Input spectum file in raw DRIAMS txt format.")
                ("mz_scale,m", po::value<int>(&mzScale)->default_value(12),
                 "Output mz resolution given as \"2^mz_scale * log2(mz - 1.007276466879)\". ")
                ("mz_min,0", po::value<double>(&mzMin)->default_value(2000.0),
                 "Minimum product ion mz. ")
                ("mz_max,1", po::value<double>(&mzMax)->default_value(20000.0),
                 "Maximum product ion mz. ")
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
            cout << "driams_raw2binned : Copyright (C) 2022 - biospi Laboratory, University of Bristol, UK" << endl;
            cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
            cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
            cout << endl;

            cout << desc << endl;
            return 0;
        }

        typedef boost::tokenizer< boost::escaped_list_separator<char> > so_tokenizer;

        ifstream raw(filePathIn, ios_base::in);
        ii m = 0;
        ii offset = ii(floor(log2(mzMin - PROTON_MASS) * (1L << mzScale))) - 1;
        ii n = (ii(ceil(log2(mzMax - PROTON_MASS) * (1L << mzScale))) + 1) - offset + 1;

        // output spectrum
        float* vs = new float[n];
        ii* js = new ii[n];
        for (ii j = 0; j < n; j++) vs[j] = 0.0f;

        Bspline bspline(3, 65536); // bspline basis function lookup table
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

            double bin = log2(mz - PROTON_MASS) * (1L << mzScale) - offset;

            fp b0 = 0.0f;
            fp b1 = ceil(bin) - bin;
            fp b2 = b1 + 1.0f;
            fp b3 = b2 + 1.0f;
            fp b4 = b3 + 1.0f;
            fp b5 = 4.0;

            ii ibin = ii(bin);
            if(ibin-2 >= 0 && ibin-2 < n) vs[ibin-2] += intensity * fp(bspline.ibasis(b1) - bspline.ibasis(b0));
            if(ibin-1 >= 0 && ibin-1 < n) vs[ibin-1] += intensity * fp(bspline.ibasis(b2) - bspline.ibasis(b1));
            if(ibin   >= 0 && ibin   < n) vs[ibin  ] += intensity * fp(bspline.ibasis(b3) - bspline.ibasis(b2));
            if(ibin+1 >= 0 && ibin+1 < n) vs[ibin+1] += intensity * fp(bspline.ibasis(b4) - bspline.ibasis(b3));
            if(ibin+2 >= 0 && ibin+2 < n) vs[ibin+2] += intensity * fp(bspline.ibasis(b5) - bspline.ibasis(b4));
        }
        while(getline(raw, line));

        cout << "bin counts" << endl;
        for (ii j = 0; j < n; j++) cout << offset + j << " " << vs[j] << endl;

        delete[] vs;
        delete[] js;
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

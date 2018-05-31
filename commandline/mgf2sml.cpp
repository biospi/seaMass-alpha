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
                        "mgf2sml [OPTIONS...] [mgf FILE]\n"
                        "mgf2sml <-m mz_scale> <file>"
        );

        general.add_options()
                ("help,h",
                 "Produce this help message")
                ("file,f", po::value<string>(&filePathIn),
                 "Input spectral library file in Mascot mgf format.")
                ("mz_scale,m", po::value<int>(&mzScale)->default_value(10),
                 "Output mz resolution given as \"2^mz_scale * log2(mz - 1.007276466879)\". ")
                ("mz_min,0", po::value<double>(&mzMin)->default_value(50.0),
                 "Minimum product ion mz. ")
                ("mz_max,1", po::value<double>(&mzMax)->default_value(3000.0),
                 "Maximum product ion mz. ")
                ("energy,e", po::value<float>(&energy),
                 "Override HCD collision energy written to sml.")
                ;

        po::options_description desc;
        desc.add(general);

        po::positional_options_description pod;
        pod.add("file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
        po::notify(vm);

        cout << endl;
        cout << "mgf2sml : Copyright (C) 2018 - biospi Laboratory, University of Bristol, UK" << endl;
        cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
        cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
        cout << endl;

        if(vm.count("help") || !vm.count("file"))
        {
            cout << desc << endl;
            return 0;
        }

        // output SML
        string fileOut = boost::filesystem::path(filePathIn).stem().replace_extension("sml").string();
        FileNetcdf sml(fileOut, NC_NETCDF4);

        li zero = 0;
        int peptideIdIndex_id = sml.write_VecNC("peptideIdIndex", &zero, 1, NC_LONG, 0, true);
        int peptideIds_id = sml.write_VecNC("peptideIds", vector<char>(), NC_CHAR, 0, true);
        int titleIndex_id = sml.write_VecNC("titleIndex", &zero, 1, NC_LONG, 0, true);
        int titles_id = sml.write_VecNC("titles", vector<char>(), NC_CHAR, 0, true);
        int charge_id = sml.write_VecNC("charge", vector<char>(), NC_BYTE, 0, true);
        int precursorMz_id = sml.write_VecNC("precursorMz", vector<double>(), NC_DOUBLE, 0, true);
        int collisionEnergy_id = sml.write_VecNC("collisionEnergy", vector<float>(), NC_FLOAT, 0, true);
        int startTime_id = sml.write_VecNC("startTime", vector<float>(), NC_FLOAT, 0, true);

        ostringstream oss; oss << "mzScale=" << mzScale;
        int spectraGroup_id = sml.create_Group(oss.str());
        int spectraGroup_i_id = sml.write_VecNC("spectrumIndex", &zero, 1, NC_LONG, spectraGroup_id, true);
        int spectraGroup_j_id = sml.write_VecNC("binLocations", vector<long>(), NC_LONG, spectraGroup_id, true);
        int spectraGroup_v_id = sml.write_VecNC("binIntensities", vector<float>(), NC_FLOAT, spectraGroup_id, true);

        typedef boost::tokenizer< boost::escaped_list_separator<char> > so_tokenizer;

        ifstream mgf(filePathIn, ios_base::in);
        ii m = 0;
        ii offset = ii(floor(log2(mzMin - PROTON_MASS) * (1L << mzScale))) - 1;
        ii n = (ii(ceil(log2(mzMax - PROTON_MASS) * (1L << mzScale))) + 1) - offset + 1;
        li nnz = 0;
        li nPeptideIdChar = 0;
        li nTitleChar = 0;

        // output spectrum
        float* vs = new float[n];
        ii* js = new ii[n];

        Bspline bspline(3, 65536); // bspline basis function lookup table
        for(string line; getline(mgf, line); )
        {
            boost::trim(line);
            if (line.size() == 0 || line[0] == '#')
                continue;

            if(line.compare(0, 10, "BEGIN IONS") != 0)
                throw runtime_error("ERROR: MGF file malformed");

             // metadata
            string peptideId;
            string title;
            float collisionEnergy = -1.0;
            float startTime = -1.0;
            double precursorMZ = 0.0;
            char charge = 0;

            // read header block
            for(; getline(mgf, line); )
            {
                boost::trim(line);
                if (line.size() == 0 || line[0] == '#')
                    continue;

                if(line.find_first_of("=") == string::npos)
                    break;

                if(line.compare(0, 4, "SEQ=") == 0)
                {
                    peptideId = line.substr(4);
                    boost::trim(peptideId);
                }
                else if(line.compare(0, 6, "TITLE=") == 0)
                {
                    title = line.substr(6);
                    boost::trim(title);
                }
                else if(line.compare(0, 12, "RTINSECONDS=") == 0)
                {
                    line = line.substr(12);
                    startTime = atof(line.c_str());
                }
                else if(line.compare(0, 7, "CHARGE=") == 0)
                {
                    line = line.substr(7);
                    charge = atoi(line.c_str());
                }
                else if(line.compare(0, 8, "PEPMASS=") == 0)
                {
                    line = line.substr(8);
                    precursorMZ = atof(line.c_str());
                }
                else if(energy < 0.0)
                {
                    if(line.compare(0, 17, "COLLISION_ENERGY=") == 0)
                    {
                        line = line.substr(17);
                        collisionEnergy = atof(line.c_str());
                    }
                }
                else
                {
                    collisionEnergy = energy;
                }
            }

            if(energy < 0.0 && collisionEnergy < 0.0)
                throw runtime_error("ERROR: MGF does not contain collision energies. You must supply mfg2sml with the '--energy' argument.");

            // read peaks
            for (ii j = 0; j < n; j++) vs[j] = 0.0f;
            do
            {
                boost::trim(line);
                if (line.size() == 0 || line[0] == '#')
                    continue;

                if(line.compare(0, 8, "END IONS") == 0)
                    break;

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
            while(getline(mgf, line));

            // make sparse
            ii spectrum_nnz = 0;
            ii jOut = 0;
            for(ii jIn = 0; jIn < n; jIn++)
            {
                if(vs[jIn] > 0.0f)
                {
                    vs[jOut] = vs[jIn];
                    js[jOut] = jIn;

                    spectrum_nnz++;
                    jOut++;
                }
            }

            const char* peptideID_cstr = peptideId.c_str();
            sml.update_VecNC(peptideIds_id, nPeptideIdChar, peptideID_cstr, peptideId.size());

            nPeptideIdChar += peptideId.size();

            const char* title_cstr = title.c_str();
            sml.update_VecNC(titles_id, nTitleChar, title_cstr, title.size());
            nTitleChar += title.size();

            sml.update_VecNC(precursorMz_id, m, &precursorMZ, 1);
            sml.update_VecNC(charge_id, m, &charge, 1);
            sml.update_VecNC(collisionEnergy_id, m, &collisionEnergy, 1);
            sml.update_VecNC(startTime_id, m, &startTime, 1);
            sml.update_VecNC(spectraGroup_j_id, nnz, js, spectrum_nnz, spectraGroup_id);
            sml.update_VecNC(spectraGroup_v_id, nnz, vs, spectrum_nnz, spectraGroup_id);

            nnz += spectrum_nnz;
            m++;

            sml.update_VecNC(peptideIdIndex_id, m, &nPeptideIdChar, 1);
            sml.update_VecNC(titleIndex_id, m, &nTitleChar, 1);
            sml.update_VecNC(spectraGroup_i_id, m, &nnz, 1, spectraGroup_id);

            //cout << m << "   " << peptideId << "   " << precursorMZ << ":" << collisionEnergy << endl;

            if (m % 10000 == 0) cout << m << " library spectra processed" << endl;
        }

        cout << m << " library spectra processed" << endl;

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

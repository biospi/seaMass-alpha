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
            "msp2sml [OPTIONS...] [MSP FILE]\n"
            "msp2sml <-m mz_scale> <file>"
        );

        general.add_options()
            ("help,h",
             "Produce this help message")
            ("file,f", po::value<string>(&filePathIn),
             "Input spectral library file in NIST msp or MassIVE sptxt format.")
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
        cout << "msp2sml : Copyright (C) 2018 - biospi Laboratory, University of Bristol, UK" << endl;
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
        int charge_id = sml.write_VecNC("charge", vector<char>(), NC_BYTE, 0, true);
        int peptideMz_id = sml.write_VecNC("peptideMz", vector<double>(), NC_DOUBLE, 0, true);
        int collisionEnergy_id = sml.write_VecNC("collisionEnergy", vector<float>(), NC_FLOAT, 0, true);

        ostringstream oss; oss << "mzScale=" << mzScale;
        int spectraGroup_id = sml.create_Group(oss.str());
        int spectraGroup_i_id = sml.write_VecNC("spectrumIndex", &zero, 1, NC_LONG, spectraGroup_id, true);
        int spectraGroup_j_id = sml.write_VecNC("binLocations", vector<long>(), NC_LONG, spectraGroup_id, true);
        int spectraGroup_v_id = sml.write_VecNC("binIntensities", vector<float>(), NC_FLOAT, spectraGroup_id, true);

        typedef boost::tokenizer< boost::escaped_list_separator<char> > so_tokenizer;

        ifstream msp(filePathIn, ios_base::in);
        ii m = 0;
        ii offset = ii(floor(log2(mzMin - PROTON_MASS) * (1L << mzScale))) - 1;
        ii n = (ii(ceil(log2(mzMax - PROTON_MASS) * (1L << mzScale))) + 1) - offset + 1;
        li nnz = 0;
        li nchar = 0;

        // output spectrum
        float* vs = new float[n];
        ii* js = new ii[n];

        Bspline bspline(3, 65536); // bspline basis function lookup table
        for(std::string line; std::getline(msp, line); )
        {
            boost::trim(line);
            if (line.size() == 0 || line[0] == '#')
                continue;

            if(line.compare(0, 5, "Name:") != 0)
                throw runtime_error("ERROR: MSP file malformed");

            string name = line.substr(5);
            boost::trim(name);

            // metadata
            string mods;
            float collisionEnergy = -1.0f;
            double parentMZ = 0.0;
            int nPeaks;

            // read header block
            for(std::string line; std::getline(msp, line); )
            {
                if(line.compare(0, 8, "Comment:") == 0)
                {
                    line = line.substr(8);
                    boost::trim(line);

                    so_tokenizer tok(line, boost::escaped_list_separator<char>("", " \t", "\"\'"));
                    for(so_tokenizer::iterator toki = tok.begin(); toki != tok.end(); ++toki)
                    {
                        if(toki->compare(0, 5, "Mods=") == 0)
                        {
                            mods = toki->substr(5);
                        }
                        else if(toki->compare(0, 7, "Parent=") == 0)
                        {
                            string s = toki->substr(7);
                            parentMZ = atof(s.c_str());
                        }
                        else if(energy < 0.0)
                        {
                            if(toki->compare(0, 4, "HCD=") == 0)
                            {
                                string s = toki->substr(4);
                                collisionEnergy = atof(s.c_str());
                            }
                        }
                        else
                        {
                            collisionEnergy = energy;
                        }
                    }
                }

                if(line.compare(0, 10, "Num peaks:") == 0)
                {
                    line = line.substr(10);
                    nPeaks = atoi(line.c_str());
                    break;
                }
            }

            if(energy < 0.0 && collisionEnergy < 0.0)
                throw runtime_error("ERROR: MSP does not contain collision energies. You must supply msp2sml with the '--energy' argument.");

            // read peaks
            for (ii j = 0; j < n; j++) vs[j] = 0.0f;
            for(std::string line; std::getline(msp, line); )
            {
                boost::trim(line);
                if (line.size() == 0)
                    break;

                so_tokenizer tok(line, boost::escaped_list_separator<char>("", " \t", "\"\'"));

                so_tokenizer::iterator toki = tok.begin();
                if (toki == tok.end())
                    break;

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

            // make sparse
            ii spectrum_nnz = 0;
            ii jOut = 0;
            for(ii jIn = 0; jIn < n; jIn++)
            {
                if(vs[jIn] > 0.0f)
                {
                    vs[jOut] = vs[jIn];
                    js[jOut] = jIn + offset;

                    spectrum_nnz++;
                    jOut++;
                }
            }

            size_t slashpos = name.find_first_of("/");
            string peptideSeq = name.substr(0, slashpos);
            string extra = name.substr(slashpos + 1, name.size() - slashpos - 1);
            size_t underscorepos = extra.find_first_of("_");
            if (underscorepos != string::npos)
                 extra = extra.substr(0, underscorepos);
            char charge = char(atoi(extra.c_str()));

            string peptideId = peptideSeq + ":" + mods;
            const char* peptideID_cstr = peptideId.c_str();
            sml.update_VecNC(peptideIds_id, nchar, peptideID_cstr, peptideId.size());

            nchar += peptideId.size();

            sml.update_VecNC(peptideMz_id, m, &parentMZ, 1);
            sml.update_VecNC(charge_id, m, &charge, 1);
            sml.update_VecNC(collisionEnergy_id, m, &collisionEnergy, 1);
            sml.update_VecNC(spectraGroup_j_id, nnz, js, spectrum_nnz, spectraGroup_id);
            sml.update_VecNC(spectraGroup_v_id, nnz, vs, spectrum_nnz, spectraGroup_id);

            nnz += spectrum_nnz;
            m++;

            sml.update_VecNC(peptideIdIndex_id, m, &nchar, 1);
            sml.update_VecNC(spectraGroup_i_id, m, &nnz, 1, spectraGroup_id);

            //cout << m << "   " << peptideID << "   " << peptideID.size() << "   " << parentMZ << ":" << collisionEnergy << endl;

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

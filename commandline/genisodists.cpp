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
#include "../io/FileNetcdf.hpp"
#include "../core/Bspline.hpp"
#include <limits>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>
using namespace std;
using namespace kernel;
namespace po = boost::program_options;


double PROTON_MASS = 1.007276466879;

double CARBON13_MASS = 13.0033548378;

void convolution(vector<double>& x, const vector<double>& a, const vector<double>& b)
{
    x.resize(a.size() + b.size() - 1, 0);
    for (ii i = 0; i < ii(a.size()); i++)
    {
        for (ii j = 0; j < ii(b.size()); j++)
        {
            x[i + j] += a[i] * b[j];
        }
    }
}


int main(int argc, const char * const * argv)
{
#ifdef NDEBUG
    try
#endif
    {
        string fileNameOut;
        double mz0;
        double mz1;
        int mzScale0;
        int mzScale1;
        int chargeStates;
        double carbonsPerDalton;
        fp threshold;
        int debugLevel;

        po::options_description general(
            "Usage\n"
            "-----\n"
            "Generates an isotope distribution database that seaMass needs to run.\n"
            "\n"
            "genisodists [OPTIONS...] <file>\n"
        );

        general.add_options()
            ("help,h", "Produce help message")
            ("file,f", po::value<string>(&fileNameOut),
             "Output file.")
            ("mz_min", po::value<double>(&mz0)->default_value(200.0),
             "Minimum m/z to consider [default=50.0]")
            ("mz_max", po::value<double>(&mz1)->default_value(2500.0),
             "Maximum m/z to consider [default=2500.0]")
            ("mz_scale_min", po::value<int>(&mzScale0)->default_value(10),
             "Minimum m/z scale to consider [default=10]")
            ("mz_scale_max", po::value<int>(&mzScale1)->default_value(16),
             "Maximum m/z scale to consider [default=20]")
            ("charge_states,z", po::value<int>(&chargeStates)->default_value(100),
             "Number of charge states to consider [default=100]")
            ("carbons,c", po::value<double>(&carbonsPerDalton)->default_value(0.03),
             "Number of carbon atoms per Dalton of mass. Set to 0 to generate no isotopes. [default=0.03]")
            ("threshold,t", po::value<fp>(&threshold)->default_value(0.00001),
             "Discard isotope if relative intensity < t. [default=0.00001]")
            ("debug,d", po::value<int>(&debugLevel)->default_value(0),
             "Debug level.")
        ;

        po::options_description desc;
        desc.add(general);

        po::positional_options_description pod;
        pod.add("file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
        po::notify(vm);

        cout << endl;
        cout << "genisodists : Copyright (C) 2017 - biospi Laboratory, University of Bristol, UK" << endl;
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

        FileNetcdf fileOut(fileNameOut, NC_NETCDF4);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // generate carbon isotope distributions
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        double massMax = pow(2.0, log2(mz1 - PROTON_MASS) + log2(chargeStates));

        vector< vector<double> > carbons;
        if (carbonsPerDalton == 0.0)
        {
            carbons.resize(1);
            carbons[0].resize(1);
            carbons[0][0] = 1.0;
        }
        else
        {
            auto maxCarbons = ii(ceil(massMax * carbonsPerDalton));
            cout << "maxCarbons=" << maxCarbons << endl;

            carbons.resize(maxCarbons);
            carbons[0].resize(2);
            carbons[0][0] = 0.9893;
            carbons[0][1] = 0.0107;

            for (ii i = 0; i < ii(carbons.size()) - 1; i++)
                convolution(carbons[i+1], carbons[i/2], carbons[i-i/2]);

            for (ii i = 0; i < ii(carbons.size()); i++)
            {
                bool trailing = false;
                for (ii j = 1; j < ii(carbons[i].size()); j++)
                {
                    if (carbons[i][j] < carbons[i][j-1])
                        trailing = true;

                    if (trailing && fp(carbons[i][j]) < threshold)
                    {
                        carbons[i].resize(j);
                        break;
                    }
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // For each charge state, Write At = m x n matrix where m is monoisotope centroid m/z and n is spectrum m/z.
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // our cubic b-spline kernel
        vector<double> hs(5);
        Bspline bspline(3, 65536);
        for (ii s = mzScale0; s <= mzScale1; s++)
        {
            cout << "mz_scale=" << s << endl;

            ostringstream oss;  oss << "s=" << setfill('0') << setw(2) << s;
            int groupId = fileOut.createGroup(oss.str());

            vector<ii> offset(2);
            // offset of monoisotope centroid m/z
            offset[0] = ii(floor(log2(mz0 - PROTON_MASS) * (1L << s)));
            // extent of monoisotope centroid m/z
            ii m = ii(floor(log2(mz1 - PROTON_MASS) * (1L << s))) - offset[0] + 1;
            // offset of spectrum m/z
            offset[1] = offset[0] - ii(hs.size() - 1) / 2;

            cout << "  offset=[" << offset[0] << "," << offset[1] << "]" << endl;
            cout << "  m=" << m << endl;

            fileOut.writeAttribute(offset, "offset", "", groupId);

            for (short z = 0; z < chargeStates; z++)
            {
                // offset to monoisotope centroid neutral mass
                double offsetMass = log2(double(z+1)) * (1L << s);

                // non-integer amount we need to shift this chargeState for chargeStates to line up
                double shift = offsetMass - ii(round(offsetMass));

                cout << " z=" << (z+1) << endl;
                cout << "  shift=" << fixed << shift << endl;

                vector<ii> is;
                vector<ii> js;
                vector<fp> vs;

                ii n = 0;
                for (ii i = offset[0]; i < offset[0] + m; i++)
                {
                    // monoisotopic neutral mass
                    double massMono = pow(2.0, (i + offsetMass) / (1L << s));

                    // number of carbons for this monoisotopic mass
                    ii nCarbons;
                    if (carbonsPerDalton == 0.0)
                        nCarbons = 1;
                    else
                        nCarbons = ii(ceil(massMono * carbonsPerDalton));

                    ii jMin = numeric_limits<ii>::max();
                    ii jMax = 0;
                    vector<double> feature(2 * m, 0); // just make it bigger than n probably is
                    for (ii p = 0; p < ii(carbons[nCarbons - 1].size()); p++)
                    {
                        double mzIsotope = (massMono + p * (CARBON13_MASS - 12.0)) / (z+1) + PROTON_MASS;

                        double offsetP = log2(mzIsotope - PROTON_MASS) * (1L << s) - shift;
                        auto offsetPi = ii(round(offsetP));
                        double offsetPf = offsetP - offsetPi;

                        hs[0] = bspline.ibasis(0.5 - offsetPf) - bspline.ibasis(0.0);
                        hs[1] = bspline.ibasis(1.5 - offsetPf) - bspline.ibasis(0.5 - offsetPf);
                        hs[2] = bspline.ibasis(2.5 - offsetPf) - bspline.ibasis(1.5 - offsetPf);
                        hs[3] = bspline.ibasis(3.5 - offsetPf) - bspline.ibasis(2.5 - offsetPf);
                        hs[4] = bspline.ibasis(4.0) - bspline.ibasis(3.5 - offsetPf);

                        for (ii k = 0; k < ii(hs.size()); k++)
                        {
                            ii j = offsetPi + k - ii(hs.size() - 1) / 2;

                            assert(j - offset[1] >= 0);
                            assert(j - offset[1] < ii(feature.size()));

                            if (hs[k] > 0.0)
                            {
                                jMin = jMin < j ? jMin : j;
                                jMax = jMax > j ? jMax : j;

                                feature[j - offset[1]] += hs[k] * carbons[nCarbons - 1][p];
                            }
                        }
                    }

                    for (ii j = jMin; j <= jMax; j++)
                    {
                        if (fp(feature[j - offset[1]]) >= threshold)
                        {
                            n = n > j ? n : j;
                            fp v = feature[j - offset[1]];

                            is.push_back(i - offset[0]);
                            js.push_back(j - offset[1]);
                            vs.push_back(v);
                        }
                    }
                }
                n++;

                cout << "  n=" << n << endl;

                MatrixSparse aTz;
                aTz.importFromCoo(m, n, ii(vs.size()), is.data(), js.data(), vs.data());

                ostringstream oss2;  oss2 << "z=" << setfill('0') << setw(4) << (z+1);
                fileOut.writeMatrixSparseCsr(aTz, oss2.str(), groupId);
            }
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

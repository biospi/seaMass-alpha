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


#include <limits>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>
#include "../core/DatasetSeamass.hpp"
#include "../peak/SMData.hpp"
#include "../peak/MathOperator.hpp"
#include "../peak/BsplineData.hpp"
#include "../peak/PeakOperator.hpp"
#include "../peak/PeakData.hpp"
#include "../peak/PeakManager.hpp"
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
        vector<char> scales(2);
        int shrinkageExponent;
        int toleranceExponent;
        int debugLevel;
		bool centroid;
		double threshold;

        // *******************************************************************

        po::options_description general(
            "Usage\n"
            "-----\n"
            "Peak detects 'seamass' output and writes the result to a new mzMLb or smb file.\n"
            "\n"
            "centroid [OPTIONS...] <file>\n"
        );

        general.add_options()
            ("help,h",
             "Produce this help message")
            ("file,f", po::value<string>(&filePathIn),
             "Input file in mzMLv or smv format produced by 'seamass'.")
			("threshold,t", po::value<double>(&threshold)->default_value(10.0),
             "Minimum ion counts in a peak. Default is 10.")
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
        cout << "seaMass-peak : Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
        cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
        cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
        cout << endl;
        initKernel(debugLevel);
        
		if(vm.count("help") || !vm.count("file"))
		{
			cout << desc << endl;
			return 0;
		}

        if(!vm.count("mz_scale"))
            scales[0] = numeric_limits<char>::max();

        if(!vm.count("st_scale"))
            scales[1] = numeric_limits<char>::max();

        string fileStemOut = boost::filesystem::path(filePathIn).stem().string();
        Dataset* dataset = FileFactory::createFileObj(filePathIn, fileStemOut, Dataset::WriteType::Input);
        if (!dataset)
            throw runtime_error("ERROR: Input file is missing or incorrect");

        Seamass::Input input;
        Seamass::Output output;
        string id;

        while (dataset->read(input, output, id))
        {
            if (getDebugLevel() % 10 == 0)
                cout << "Processing " << id << endl;

            // load back into Seamass
            Seamass seamassCore(input, output);

            // now centroid
            VecMat<double> mzPeak;
            VecMat<float> pkPeak;
            vector<size_t> mzpkVecSize;

            Seamass::ControlPoints contpts;
            seamassCore.getOutputControlPoints1d(contpts);
            uli dims[2];
            double mzRes;
            double rtRes;
            vector<ii> offset=contpts.offset;
            mzRes=double(contpts.scale[0]);
            if(contpts.scale.size() > 1)
                rtRes=double(contpts.scale[1]);
            else
                rtRes=0;

            VecMat<float> rawCoeff(contpts.extent[1],contpts.extent[0],contpts.coeffs);
            rawCoeff.getDims(dims);

            SMData2D<OpUnit> A(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
            SMData2D<OpNablaH> dhA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
            SMData2D<OpNabla2H> d2hA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);

            for(size_t i = 0; i < A.rt.size(); ++i)
            {
                A.rt[i] = input.startTimes[i];
                dhA.rt[i] = input.startTimes[i];
                d2hA.rt[i] = input.startTimes[i];
            }

            BsplineData<> bsData(A,dhA,d2hA);

            PeakManager<PeakData,BsplineData,Centroid2D> centriodPeak(bsData,threshold);
            centriodPeak.execute();
            centriodPeak.peak->getPeakMat(mzPeak, pkPeak, contpts.extent[1], mzpkVecSize);

            input.type = Seamass::Input::Type::Centroided;
            vector<double>().swap(input.locations);
            vector<fp>().swap(input.counts);
            vector<li>().swap(input.countsIndex);
            input.countsIndex.push_back(0);
            for (ii i = 0; i < mzpkVecSize.size(); i++)
            {
                input.locations.insert(input.locations.end(), mzPeak.v.begin(), mzPeak.v.end());
                input.counts.insert(input.counts.end(), pkPeak.v.begin(), pkPeak.v.end());
                input.countsIndex.push_back(input.counts.size());
            }
            dataset->write(input, id);

            if (getDebugLevel() % 10 == 0)
                cout << endl;
        }

        delete dataset;
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

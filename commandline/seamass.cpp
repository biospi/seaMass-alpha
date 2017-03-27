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

#include "../io/NetcdfWriter.hpp"
#include "../io/DatasetMzmlb.hpp"
#include "../io/MzMLb.hpp"


using namespace std;
namespace po = boost::program_options;


int main(int argc, const char * const* argv)
{
#ifdef NDEBUG
	try
#endif
	{
        cout << endl;
        SeamassCore::notice();

        string inFile;
        vector<short> scales(2);
        int shrinkageExponent;
        int toleranceExponent;
        int debugLevel;

        // *******************************************************************

        po::options_description general(
            "Usage\n"
            "-----\n"
            "seamass [OPTIONS...] [MZMLB]\n"
            "seamass <-m mz_scale> <-s st_scale> <-l lambda> <-t tol> <file>"
        );

        general.add_options()
            ("help,h",
             "Produce this help message")
            ("file,f", po::value<string>(&inFile),
             "Input file in mzMLb format. Use pwiz-mzmlb (https://github.com/biospi/mzmlb) to convert from mzML or vendor format.")
          /*("mz_range", po::value<string>(),
             "Input m/z range (in Thompsons) in the following format: [min,max]"
             "default: autodetect")
            ("st_range", po::value<string>(),
             "Input scan-time (in seconds) in the following format: [min,max]"
             "default: autodetect")*/
            ("mz_scale,m", po::value<short>(&scales[0]),
             "Output m/z resolution given as \"b-splines per Th = 2^mz_scale * 60 / 1.0033548378\". "
             "Use 0 or 1 for ToF (e.g. 1 is suitable for 30,000 resolution), 3 for Orbitrap. "
             "Default is to autodetect.")
            ("st_scale,s", po::value<short>(&scales[1]),
             "output scan-time resolution given as \"b-splines per second = 2^st_scale\". Use around 4. "
             "Default is to autodetect.")
            ("lambda,l", po::value<int>(&shrinkageExponent)->default_value(0),
             "Amount of denoising given as \"L1 lambda = 2^shrinkage\". Use around 0.")
            ("tol,t", po::value<int>(&toleranceExponent)->default_value(-10),
             "Convergence tolerance, given as \"gradient <= 2^tol\". Use around -10.")
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
        
		if(vm.count("help") || !vm.count("file"))
		{
			cout << desc << endl;
			return 0;
		}

        if(!vm.count("mz_scale"))
            scales[0] = numeric_limits<short>::max();

        if(!vm.count("st_scale"))
            scales[1] = numeric_limits<short>::max();

        cout << endl << getThreadInfo() << endl << endl;
        setDebugLevel(debugLevel);

        DatasetMzmlb msFile(inFile);
        OutmzMLb outmzMLb(inFile, msFile);
        SeamassCore::Input input;
        string id;
        double tolerance = pow(2.0, (double)toleranceExponent);
        double shrinkage = pow(2.0, (double)shrinkageExponent);
        while (msFile.next(input, id))
        {
            //if (id != "pos_411.25") continue;

            if (getDebugLevel() % 10 == 0) cout << "Processing " << id << endl;

            SeamassCore sm(input, scales, shrinkage, tolerance);

            do
            {
                if (getDebugLevel() >= 10)
                {
                     // for now, lets write out smo
                    ostringstream oss; oss << "." << setfill('0') << setw(4) << sm.getIteration();

                    // get 1D control points for centroiding and write to smo1D file
                    string smo1dFileName = boost::filesystem::change_extension(inFile, ".").string() + id + oss.str() + ".smo1D";
                    if (getDebugLevel() % 10 == 0) cout << "Writing " << smo1dFileName << endl;
                    NetcdfWriter netcdfWriter(smo1dFileName);
                    SeamassCore::ControlPoints controlPoints;
                    sm.getOutputControlPoints1d(controlPoints);
                    netcdfWriter.writeSmo(controlPoints);

                    // write smo 2D file if necessary
                    if (controlPoints.extent[1] > 1)
                    {
                        string smo2dFileName = boost::filesystem::change_extension(inFile, ".").string() + id + oss.str() + ".smo2D";
                        if (getDebugLevel() % 10 == 0) cout << "Writing " << smo2dFileName << endl;
                        NetcdfWriter netcdfWriter2d(smo2dFileName);
                        sm.getOutputControlPoints(controlPoints);
                        netcdfWriter2d.writeSmo(controlPoints);
                    }
                }
            }
            while (sm.step());

            // get 1D control points for centroiding and write to smo1D file
            string smo1dFileName = boost::filesystem::change_extension(inFile, ".").string() + id + ".smo1D";
            if (getDebugLevel() % 10 == 0) cout << "Writing " << smo1dFileName << endl;
            NetcdfWriter netcdfWriter(smo1dFileName);
            SeamassCore::ControlPoints controlPoints;
            sm.getOutputControlPoints1d(controlPoints);
            netcdfWriter.writeSmo(controlPoints);

            // write smo 2D file if necessary
            if (controlPoints.extent[1] > 1)
            {
                string smo2dFileName = boost::filesystem::change_extension(inFile, ".").string() + id + ".smo2D";
                if (getDebugLevel() % 10 == 0) cout << "Writing " << smo2dFileName << endl;
                NetcdfWriter netcdfWriter2d(smo2dFileName);
                sm.getOutputControlPoints(controlPoints);
                netcdfWriter2d.writeSmo(controlPoints);
            }

            // write seaMass outputBinCounts to new mzMLb file
            /*vector<fp> outputBinCounts(input.binCounts.size());
            cout << outputBinCounts.size() << endl;
            sm.getOutputBinCounts(outputBinCounts); // retrieve seaMass processed outputBinCounts
            // convert ion counts into ion density (counts per Th) and scale by exposures
            if (input.exposures.size() > 0)
            {
                if (input.spectrumIndex.size() > 0)
                {
                    // 2D data
                    for (li j = 0; j < (li)input.spectrumIndex.size() - 1; j++)
                    {
                        for (li i = input.spectrumIndex[j]; i < input.spectrumIndex[j + 1]; i++)
                        {
                            outputBinCounts[i] /= (fp) (input.binEdges[i + j + 1] - input.binEdges[i + j]) * input.exposures[j];
                        }
                    }
                }
                else
                {
                    // 1D data
                    for (li i = 0; i < (li)input.binCounts.size(); i++)
                    {
                        outputBinCounts[i] /= (fp) (input.binEdges[i + 1] - input.binEdges[i]) * input.exposures[0];
                    }
                }
            }
            outmzMLb.writeVecData(outputBinCounts); // write to mzMLb

            // write SMV file
            ostringstream oss;
            oss << boost::filesystem::change_extension(inFile, "").string() << "." << id << ".smv";
            NetcdfWriter smv(oss.str());
            vector<fp> originalBinCounts = input.binCounts; // save original input.binCounts
            sm.getOutputBinCounts(input.binCounts); // retrieve seaMass processed outputBinCounts
            for (size_t i = 0; i < input.binCounts.size(); i++) input.binCounts[i] = originalBinCounts[i] - input.binCounts[i]; // compute residuals
            smv.write_input(input); // write residuals to smv
            // write RTree
            SeaMass::Output output;
            sm.getOutput(output);
            smv.write_output(output, shrinkageExponent, toleranceExponent, 4096);*/

            // demonstration code to load from smv back into seaMass:Input and seaMass:Output structs
            // lets pretend Input struct and Output::baseline_size,baseline_scale,baseline_offset are already filled (as I'm not implementing a HDFReader class)
            // therefore we just need to load from the RTree ito Output::weights,scales,offsets
            /*SeaMass::Output loaded;
            loaded.baselineOffset = output.baselineOffset;
            loaded.baselineScale = output.baselineScale;
            loaded.baselineExtent = output.baselineExtent;
            RTreeReader rtree(oss.str());
            rtree.read(loaded);
            cout << "Number of saved bases: " << output.weights.size() << endl;
            cout << "Number of loaded bases: " << loaded.weights.size() << endl;*/

            if (getDebugLevel() % 10 == 0) cout << endl;
         }

        //vector<MetadataMzmlbSpectrum> *spcPtr = msFile.getSpectrumMetaData();
        //outmzMLb.writeXmlData(spcPtr);
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

//
// $Id$
//
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
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>

#include "../io/HDF5Writer.hpp"
#include "../io/RTreeReader.hpp"
#include "../io/MSFileData.hpp"
#include "../core/SeamassCore.hpp"
#include "../io/MzMLb.hpp"


using namespace std;
namespace po = boost::program_options;


int main(int argc, char **argv)
{
	SeamassCore::notice();

	string in_file;
	vector<int> scales(2);
	int shrinkageExponent;
	int toleranceExponent;
	int threads;
	int debugLevel;

	// *******************************************************************

	po::options_description general("Usage\n"
			"-----\n"
			"seamass [OPTIONS...] [MZMLB]\n"
			"seamass <-f in_file> <-m mz_scale> <-r st_scale> <-s shrinkage> <-l tol> <-t threads> <-o out_type>\n"
			"seamass -m 1 -r 4 -s -4 -l -9 -t 4 -o 0");

	general.add_options()
		("help,h", "Produce help message")
		("file,f", po::value<string>(&in_file),
			"Raw input file in seaMass Input format (mzMLb, csv etc.) "
			"guidelines: Use pwiz-seamass to convert from mzML or vendor format")
		("mz_scale,m", po::value<int>(&scales[0])->default_value(numeric_limits<short>::max()),
			"m/z resolution given as: \"b-splines per Th = 2^mz_scale * 60 / 1.0033548378\" "
			"guidelines: between 0 to 1 for ToF (e.g. 1 is suitable for 30,000 resolution), 3 for Orbitrap, "
			"default: auto")
		("st_scale,r", po::value<int>(&scales[1])->default_value(numeric_limits<short>::max()),
			"Scan time resolution given as: \"b-splines per second = 2^st_scale\" "
			"guidelines: around 4, "
			"default: auto")
		("shrinkage,s", po::value<int>(&shrinkageExponent)->default_value(0), ""
			"Amount of denoising given as: \"L1 lambda = 2^shrinkage\" "
			"guidelines: around 0, "
			"default: 0")
		("tolerance,t", po::value<int>(&toleranceExponent)->default_value(-10),
			"Convergence tolerance, given as: \"gradient <= 2^tol\" "
			"guidelines: around -10, "
			"default: -10")
		("debug_level,d", po::value<int>(&debugLevel)->default_value(0),
			"Debug level, "
			"guidelines: set to 1+ for convergence stats, 2+ for performance stats, 3+ to write intermediate iterations to disk, "
			"default: 0")
		("threads", po::value<int>(&threads)->default_value(0),
			"Number of OpenMP threads to use, "
			"guidelines: will automatically be set to amount of CPU cores, "
			"default: auto");

	po::options_description desc;
	desc.add(general);

	try
	{
		po::positional_options_description pod;
		pod.add("file", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
		po::notify(vm);
        
        setNumThreads(threads);

		if(vm.count("help"))
		{
			cout << desc << endl;
			return 0;
		}
		if(vm.count("file"))
		{
			cout << "Opening file: " << vm["file"].as<string>() << endl;
		}
		else
		{
			throw "Valid seamass input file was not given...";
		}
	}
	catch(exception& e)
	{
		cerr << "error: " << e.what() <<endl;
		cout << desc << endl;
		return 1;
	}
	catch(const char* msg)
	{
		cerr << "error: "<< msg <<endl;
		cout << desc << endl;
		return 1;
	}
	catch(...)
    {
		cerr<<"Exception of unknown type!\n";
	}
    
	mzMLbInputFile msFile(in_file);
    OutmzMLb outmzMLb(in_file,msFile);
	SeamassCore::Input input;
	string id;
	double tolerance = pow(2.0, toleranceExponent);
	double shrinkage = pow(2.0, shrinkageExponent);
	while (msFile.next(input, id))
	{
		cout << endl << "Processing " << id << ":" << endl;

		SeamassCore sm(input, scales, shrinkage, tolerance, debugLevel);

		do
		{
			if (debugLevel >= 3)
			{
				// create SMV file
				ostringstream oss;
				oss << boost::filesystem::change_extension(in_file, "").string() << "." << id << "." << setfill('0') << setw(4) << sm.getIteration() << ".smv";
				HDF5Writer smv(oss.str());

				// save back input but with bin_counts now containing the residuals
				vector<fp> originalBinCounts = input.binCounts;
				sm.getOutputBinCounts(input.binCounts);
				for (size_t i = 0; i < input.binCounts.size(); i++) input.binCounts[i] = originalBinCounts[i] - input.binCounts[i];
				smv.write_input(input);
				input.binCounts = originalBinCounts;

				// write RTree
				SeamassCore::Output output;
				sm.getOutput(output);
				smv.write_output(output, shrinkageExponent, toleranceExponent, 4096);

				// for now, lets also write out an smo
				ostringstream oss2;
				oss2 << boost::filesystem::change_extension(in_file, "").string() << "." << id << "." << setfill('0') << setw(4) << sm.getIteration() << ".smo";
				HDF5Writer smo(oss2.str());

				SeamassCore::ControlPoints controlPoints;
				sm.getOutputControlPoints(controlPoints);
				smo.write_output_control_points(controlPoints);
			}
		}
		while (sm.step());

		// write seaMass outputBinCounts to new mzMLb file 
		vector<fp> outputBinCounts; 
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
		oss << boost::filesystem::change_extension(in_file, "").string() << "." << id << ".smv";
		HDF5Writer smv(oss.str());
		vector<fp> originalBinCounts = input.binCounts; // save original input.binCounts
		sm.getOutputBinCounts(input.binCounts); // retrieve seaMass processed outputBinCounts 
		for (size_t i = 0; i < input.binCounts.size(); i++) input.binCounts[i] = originalBinCounts[i] - input.binCounts[i]; // compute residuals
		smv.write_input(input); // write residuals to smv
		// write RTree
		//SeaMass::Output output;
		//sm.getOutput(output);
		//smv.write_output(output, shrinkageExponent, toleranceExponent, 4096);

        // write SMO file
		ostringstream oss2;
		oss2 << boost::filesystem::change_extension(in_file, "").string() << "." << id << ".smo";
		HDF5Writer smo(oss2.str());
		SeamassCore::ControlPoints controlPoints;
		sm.getOutputControlPoints(controlPoints);
		smo.write_output_control_points(controlPoints);

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
	}

    vector<spectrumMetaData> *spcPtr = msFile.getSpectrumMetaData();
    outmzMLb.writeXmlData(spcPtr);

	return 0;
}

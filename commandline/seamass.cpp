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
#include <omp.h>

#include "../io/HDF5Writer.hpp"
#include "../io/RTreeReader.hpp"
#include "../io/MSFileData.hpp"
#include "../core/seaMass.hpp"


using namespace std;
namespace po = boost::program_options;


int main(int argc, char *argv[])
{
	SeaMass::notice();

	string in_file;
	vector<ii> scales(2);
	ii shrinkageExponent;
	ii toleranceExponent;
	ii threads;

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
		("mz_scale,m",po::value<ii>(&scales[0])->default_value(numeric_limits<short>::max()),
			"m/z resolution given as: \"b-splines per Th = 2^mz_scale * 60 / 1.0033548378\" "
			"guidelines: between 0 to 1 for ToF (e.g. 1 is suitable for 30,000 resolution), 3 for Orbitrap, "
			"default: auto")
		("st_scale,r", po::value<ii>(&scales[1])->default_value(numeric_limits<short>::max()),
			"Scan time resolution given as: \"b-splines per second = 2^st_scale\" "
			"guidelines: around 4, "
			"default: auto")
		("shrinkage,s",po::value<ii>(&shrinkageExponent)->default_value(-4),""
			"Amount of denoising given as: \"L1 shrinkage = 2^shrinkage\" "
			"guidelines: around -4, "
			"default: -4")
		("tolerance,l",po::value<ii>(&toleranceExponent)->default_value(-9),
			"Convergence tolerance, given as: \"gradient <= 2^tol\" "
			"guidelines: around -9, "
			"default: -9")
		("threads,t",po::value<ii>(&threads)->default_value(4),
			"Number of OpenMP threads to use, "
			"guidelines: set to amount of CPU cores or 4, whichever is smaller, "
			"default: 4");

	po::options_description desc;
	desc.add(general);

	try
	{
		po::positional_options_description pod;
		pod.add("file", 1);

		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
		po::notify(vm);

		if(vm.count("help"))
		{
			cout<<desc<<endl;
			return 0;
		}
			if(vm.count("threads"))
		{
			threads=vm["threads"].as<int>();
		}
		else
		{
			threads=omp_get_max_threads();
		}
		if(vm.count("file"))
		{
			cout<<"Opening file: "<<vm["file"].as<string>()<<endl;
		}
		else
		{
			throw "Valid seamass input file was not given...";
		}
	}
	catch(exception& e)
	{
		cerr<<"error: " << e.what() <<endl;
		cout<<desc<<endl;
		return 1;
	}
	catch(const char* msg)
	{
		cerr<<"error: "<<msg<<endl;
		cout<<desc<<endl;
		return 1;
	}
	catch(...)
		{
		cerr<<"Exception of unknown type!\n";
	}

	mzMLbInputFile msFile(in_file);
	SeaMass::Input input;
	string id;
	double shrinkage = pow(2.0, shrinkageExponent);
	double tolerance = pow(2.0, toleranceExponent);
	while(msFile.next(input, id))
	{
		cout << "Processing id=" << id << endl;
		SeaMass sm(input, scales);

		do
		{
			ostringstream oss;
			oss << boost::filesystem::change_extension(in_file, "").string() << "." << id << "." << setfill('0') << setw(4) << sm.getIteration() << ".smo";
			HDF5Writer smo(oss.str());

			SeaMass::ControlPoints controlPoints;
			sm.getOutputControlPoints(controlPoints);
			smo.write_output_control_points(controlPoints);

			vector<fp> binCounts;
			sm.getOutputBinCounts(binCounts);
			smo.write("outputBinCounts", binCounts);
		}
		while (sm.step(shrinkage, tolerance));

		{
			ostringstream oss;
			oss << boost::filesystem::change_extension(in_file, "").string() << "." << id << ".smo";
			HDF5Writer smo(oss.str());

			SeaMass::ControlPoints controlPoints;
			sm.getOutputControlPoints(controlPoints);
			smo.write_output_control_points(controlPoints);

			vector<fp> binCounts;
			sm.getOutputBinCounts(binCounts);
			smo.write("outputBinCounts", binCounts);
		}

		// save back input but with bin_counts now containing the residuals
		/*vector<fp> binCounts(input.binCounts.size(), 0.0);
		sm.getOutputBinCounts(binCounts);
		for (ii i = 0; i < input.binCounts.size(); i++) input.binCounts[i] -= binCounts[i];
		smo.write_input(input);*/

		// save seamass output as RTree
		/*SeaMass::Output output;
		sm.getOutput(output);
		smv.write_output(output, shrinkage, tolerance, 4096);

		// demonstration code to load from smv back into seaMass:Input and seaMass:Output structs
		// lets pretend Input struct and Output::baseline_size,baseline_scale,baseline_offset are already filled (as I'm not implementing a HDFReader class)
		// therefore we just need to load from the RTree ito Output::weights,scales,offsets
		SeaMass::Output loaded;
		loaded.baselineOffset = output.baselineOffset;
		loaded.baselineScale = output.baselineScale;
		loaded.baselineExtent = output.baselineExtent;
		RTreeReader rtree(oss.str());
		rtree.read(loaded);
		cout << "Number of saved bases: " << output.weights.size() << endl;
		cout << "Number of loaded bases: " << loaded.weights.size() << endl;*/
	}

	return 0;
}

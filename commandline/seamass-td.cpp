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

#include "../kernel/Matrix.hpp"
#include "../kernel/VecMat.hpp"
#include "../kernel/FileNetcdf.hpp"
#include "../topdown/SeamassTopdown.hpp"


using namespace std;
namespace po = boost::program_options;


int main(int argc, char **argv)
{
	//SeaMass::notice();

	string in_file;
	ii maxMass;
	ii binsPerDalton;
	ii shrinkageExponent;
	ii toleranceExponent;
	ii threads;
	ii debugLevel;

	// *******************************************************************

	po::options_description general("Usage\n"
		"-----\n"
		"seamass-td [OPTIONS...] [MZMLB]\n"
		"seamass-td <-f in_file> <-m mz_scale> <-r st_scale> <-s shrinkage> <-l tol> <-t threads> <-o out_type>\n"
		"seamass-td -m 1 -r 4 -s -4 -l -9 -t 4 -o 0");

	general.add_options()
		("help,h", "Produce help message")
		("file,f", po::value<string>(&in_file),
		"Raw input file in seaMass Input format (mzMLb, csv etc.) "
		"guidelines: Use pwiz-seamass to convert from mzML or vendor format")
		("max_mass,m", po::value<ii>(&maxMass)->default_value(100000),
		"Maximum output mass in Daltons, "
		"guidelines: depends on what you expect to see in the spectrum, "
		"default: 100000")
		("bpd,b", po::value<ii>(&binsPerDalton)->default_value(100),
		"guidelines: Number of output bins per Dalton,"
		"guidelines: this determines the precision of the result, "
		"default: 100")
		("shrinkage,s", po::value<ii>(&shrinkageExponent)->default_value(0), ""
		"Amount of denoising given as: \"L1 shrinkage = 2^shrinkage\" "
		"guidelines: around 0, "
		"default: 0")
		("tolerance,t", po::value<ii>(&toleranceExponent)->default_value(-10),
		"Convergence tolerance, given as: \"gradient <= 2^tol\" "
		"guidelines: around -10, "
		"default: -10")
		("debug_level,d", po::value<ii>(&debugLevel)->default_value(0),
		"Debug level, "
		"guidelines: set to 1 for debugging information, 2 to additionally write intermediate iterations to disk, "
		"default: 0")
		("threads", po::value<ii>(&threads)->default_value(4),
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

		if (vm.count("help"))
		{
			cout << desc << endl;
			return 0;
		}
		if (vm.count("threads"))
		{
			threads = vm["threads"].as<int>();
		}
		else
		{
		}
		if (vm.count("file"))
		{
			cout << "Opening file: " << vm["file"].as<string>() << endl;
		}
		else
		{
			throw "Valid seamass input file was not given...";
		}
	}
	catch (exception& e)
	{
		cerr << "error: " << e.what() << endl;
		cout << desc << endl;
		return 1;
	}
	catch (const char* msg)
	{
		cerr << "error: " << msg << endl;
		cout << desc << endl;
		return 1;
	}
	catch (...)
	{
		cerr << "Exception of unknown type!\n";
	}

	double tolerance = pow(2.0, toleranceExponent);
	double shrinkage = pow(2.0, shrinkageExponent);

	FileNetcdf inFile(in_file);
	VecMat<> in;
	inFile.read_MatNC("controlPoints", in);
	uli extent[2];
	in.getDims(extent);
	vector<int> offset, scale;
	int varid = inFile.read_VarIDNC("controlPoints");
	inFile.read_AttNC("scale", varid, scale);
	inFile.read_AttNC("offset", varid, offset);

	vector<fp> coeffs(extent[0], 0.0);
	for (li j = 0; j < extent[1]; j++)
	{
		for (li i = 0; i < extent[0]; i++)
		{
			coeffs[i] += (*in.m)[j + i * extent[1]];
		}
	}

	SeamassTopdown::Input input;
	input.binCounts.resize(extent[0] - 3, 0.0);
	for (li i = 0; i < input.binCounts.size(); i++)
	{
		input.binCounts[i]  = 0.04 * coeffs[i] + 0.46 * coeffs[i + 1] + 0.46 * coeffs[i + 2] + 0.04 * coeffs[i + 3];
	}
	input.scale = scale[0];
	input.offset = offset[0];

	FileNetcdf outFile("input.smr", NC_NETCDF4);
	outFile.write_VecNC("binCounts", input.binCounts, NC_FLOAT);

	SeamassTopdown sm(input, maxMass, binsPerDalton, shrinkage, tolerance, debugLevel);
	do
	{
	}
	while (sm.step());




	/*while (msFile.next(input, id))
	{
		cout << endl << "Processing " << id << ":" << endl;

		SeaMass sm(input, scales, shrinkage, tolerance, debugLevel);

		do
		{
			if (debugLevel > 1)
			{
				// create SMV file
				ostringstream oss;
				oss << boost::filesystem::change_extension(in_file, "").string() << "." << id << "." << setfill('0') << setw(4) << sm.getIteration() << ".smv";
				NetcdfWriter smv(oss.str());

				// save back input but with bin_counts now containing the residuals
				vector<fp> originalBinCounts = input.binCounts;
				sm.getOutputBinCounts(input.binCounts);
				for (ii i = 0; i < input.binCounts.size(); i++) input.binCounts[i] = originalBinCounts[i] - input.binCounts[i];
				smv.write_input(input);
				input.binCounts = originalBinCounts;

				// write RTree
				SeaMass::Output output;
				sm.getOutput(output);
				smv.write_output(output, shrinkageExponent, toleranceExponent, 4096);

				// for now, lets also write out an smo
				ostringstream oss2;
				oss2 << boost::filesystem::change_extension(in_file, "").string() << "." << id << "." << setfill('0') << setw(4) << sm.getIteration() << ".smo";
				NetcdfWriter smo(oss2.str());

				SeaMass::ControlPoints controlPoints;
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
		NetcdfWriter smv(oss.str());
		vector<fp> originalBinCounts = input.binCounts; // save original input.binCounts
		sm.getOutputBinCounts(input.binCounts); // retrieve seaMass processed outputBinCounts 
		for (ii i = 0; i < input.binCounts.size(); i++) input.binCounts[i] = originalBinCounts[i] - input.binCounts[i]; // compute residuals
		smv.write_input(input); // write residuals to smv
		// write RTree
		//SeaMass::Output output;
		//sm.getOutput(output);
		//smv.write_output(output, shrinkageExponent, toleranceExponent, 4096);

        // write SMO file
		ostringstream oss2;
		oss2 << boost::filesystem::change_extension(in_file, "").string() << "." << id << ".smo";
	    NetcdfWriter smo(oss2.str());
		SeaMass::ControlPoints controlPoints;
		sm.getOutputControlPoints(controlPoints);
		smo.write_output_control_points(controlPoints);
	}

    vector<MzmlbSpectrumMetadata> *spcPtr = msFile.getSpectrumMetaData();
    outmzMLb.writeXmlData(spcPtr);*/

	return 0;
}

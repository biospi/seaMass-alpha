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
#include "../io/SMVWriter.hpp"
#include "../io/MSFileData.hpp"
#include "../core/seaMass.hpp"


using namespace std;
namespace po = boost::program_options;


int main(int argc, char *argv[])
{
	seaMass::notice();

	string in_file;
	vector<ii> scales(2);
	ii _shrinkage;
	ii _tolerance;
	ii threads;

	// *******************************************************************

	po::options_description general("Usage\n"
			"-----\n"
			"seamass2controlpoints [OPTIONS...] [MZMLB]\n"
			"seamass2controlpoints <-f in_file>\n"
			"seamass2controlpoints");

	general.add_options()
		("help,h", "Produce help message")
		("file,f", po::value<string>(&in_file),
			"Raw input file in seaMass Input format (mzMLb, csv etc.) "
			"guidelines: Use pwiz-seamass to convert from mzML or vendor format");

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

	/*SMVFile smv(in_file);
	seaMass::Output output;
	smv.read(output);
	seaMass sm(output);

	ostringstream oss2;
	oss2 << boost::filesystem::change_extension(in_file, "").string() << ".smo";
	HDF5Writer smo(oss2.str());

	seaMass::ControlPoints control_points;
	sm.get_output_control_points(control_points);
	smo.write_output_control_points(control_points);*/

	return 0;
}

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
//#include <omp.h>

#include "../io/NetcdfWriter.hpp"
#include "../io/MSFileData.hpp"

using namespace std;
namespace po = boost::program_options;


int main(int argc, char *argv[])
{
	string in_file;

	// *******************************************************************

	po::options_description general("Usage\n"
			"-----\n"
			"mzmlb2smi [OPTIONS...] [MZMLB]\n"
			"mzmlb2smi <-f in_file>\n");

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

	mzMLbInputFile msFile(in_file);
	SeamassCore::Input input;
	string id;
	for (int i = 0; msFile.next(input, id); i++)
	{
		ostringstream oss;

		oss << boost::filesystem::change_extension(in_file, "").string() << "." << id << ".smi";
        NetcdfWriter netcdfWriter(oss.str());
		cout << "Writing file: " << oss.str() << endl;

		netcdfWriter.writeSmi(input);
	}

	return 0;
}

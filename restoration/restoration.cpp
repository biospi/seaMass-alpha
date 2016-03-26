//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> manchester.ac.uk>
//
// Copyright (C) 2013  CADET Laboratory for Medical Bioinformatics, University of Manchester, UK
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

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <H5Cpp.h>
#include <pugixml.hpp>
#include <map>

#include "seamass.hpp"
#include "NetCDFile.hpp"

using namespace std;
namespace xml = pugi;

struct spectrum
{
	size_t index;
	unsigned short preset_config;
	double precursor_mz;
	double scan_start_time;
	size_t scan_start_time_index;
	size_t count;
};


bool scan_start_time_order(const spectrum& lhs, const spectrum& rhs)
{
	return lhs.scan_start_time <= rhs.scan_start_time;
}


bool seamass_order(const spectrum& lhs, const spectrum& rhs)
{
	if (lhs.preset_config == rhs.preset_config)
	{
		if (lhs.precursor_mz == rhs.precursor_mz)
		{
			return lhs.scan_start_time <= rhs.scan_start_time;
		}
		else
		{
			return lhs.precursor_mz < rhs.precursor_mz;
		}
	}
	else
	{
		return lhs.preset_config < rhs.preset_config;
	}
}


void preSetScanConfig(vector<unsigned long> &scanConf)
{
	map<unsigned long, unsigned long> preSetScanConf;
	unsigned long idx=1;
	pair<map<unsigned long, unsigned long>::iterator, bool> ret;
	for(size_t i = 0; i < scanConf.size(); ++i)
	{
		ret = preSetScanConf.insert(pair<unsigned long, unsigned long>(scanConf[i],idx));
		if(ret.second == true)
		{
			++idx;
		}
	}
	if(preSetScanConf.size() == 1) preSetScanConf[scanConf[0]]=0;
	transform(scanConf.begin(),scanConf.end(),scanConf.begin(),
			[&preSetScanConf](unsigned long x){return preSetScanConf[x];} );
}

int main(int argc, char *argv[])
{
	H5::Exception::dontPrint();
	seamass::notice();
	if (argc != 8)
	{
		cout << "Usage" << endl;
		cout << "-----" << endl;
		cout << "restoration <in_file> <mz_res> <rt_res> <shrinkage> <tol> <threads> <out_type>" << endl;
		cout << endl;
		cout << "<in_file>:  Raw input file in seaMass Input format (smi)" << endl;
		cout << "            guidelines: Use pwiz-seamass to convert from mzML or vendor format" << endl;
		cout << "<mz_res>:   MS resolution given as: \"b-splines per Th = 2^mz_res * 60 / 1.0033548378\"" << endl;
		cout << "            guidelines: between 0 to 1 for ToF (e.g. 1 is suitable for 30,000 resolution), 3 for Orbitrap" << endl;
		cout << "<rt_res>:   LC resolution given as: \"b-splines per minute = 2^rt_res\"" << endl;
		cout << "            guidelines: around 4" << endl;
		cout << "<shrink>:   Amount of denoising given as: \"L1 shrinkage = 2^shrinkage\"" << endl;
		cout << "            guidelines: around -4" << endl;
		cout << "<tol>:      Convergence tolerance, given as: \"gradient <= 2^tol\"" << endl;
		cout << "            guidelines: around -9" << endl;
		cout << "<threads>:  Number of OpenMP threads to use" << endl;
		cout << "            guidelines: set to amount of CPU cores or 4, whichever is smaller" << endl;
		cout << "<out_type>: Type of output desired" << endl;
		cout << "            guidelines: 0 = just viz_client input; 1 = also smo; 2 = also debug" << endl;
		return 0;
	}
	string in_file = argv[1];
	int mz_res = atoi(argv[2]);
	int rt_res = atoi(argv[3]);
	int shrink = atoi(argv[4]);
	int tol = atoi(argv[5]);
	int threads = atoi(argv[6]);
	int out_type = atoi(argv[7]);

	int lastdot = in_file.find_last_of("."); 
	string id = (lastdot == string::npos) ? in_file : in_file.substr(0, lastdot); 

	// Open mzML file
	NetCDFile mzMLb3File(in_file);
	vector<char> mzMLBuff;
	vector<size_t> dataMatDim;
	vector<InfoGrpVar> dataSetList;

	// Load mzML metadata
	mzMLb3File.read_VecNC("mzML",mzMLBuff);
	size_t xmlSize=sizeof(char)*mzMLBuff.size();

	xml::xml_document docmzML;
	xml::xpath_node_set tools;
	xml::xml_parse_result result = docmzML.load_buffer_inplace(&mzMLBuff[0],xmlSize);

	hsize_t ns;
	istringstream(docmzML.child("mzML").child("run").child("spectrumList").attribute("count").value())>>ns;

	// query necessary metadata
    cout << "Querying metadata from " << ns << " spectra..." << endl;
	vector<double> start_times;
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/scanList/scan/cvParam[@accession='MS:1000016']");
	if(tools.empty())
	{
		cout<<"Cannot Find Start Time"<<endl;
		exit(1);
	}
	else
	{
		double rescale=1.0;
		if(string(tools.first().node().attribute("unitName").value()).compare("second") == 0)
			rescale = 1.0/60.0;

		for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			double scanRT;
			istringstream(itr->node().attribute("value").value()) >> scanRT;
			start_times.push_back(scanRT*rescale);
		}
	}

    vector<double> precursor_mzs(ns,0.0);
    tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/precursorList/precursor");
	if(!tools.empty())
	{
		for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			double preMZ;
			size_t idx;
			istringstream(itr->node().child("selectedIonList").child("selectedIon").
				find_child_by_attribute("accession","MS:1000744").attribute("value").value())>>preMZ;
			istringstream(itr->node().parent().parent().attribute("index").value())>>idx;
			precursor_mzs[idx]=preMZ;
		}
	}

    vector<unsigned long> config_indices(ns,0);
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/scanList/scan/cvParam[@accession='MS:1000616']");
	if(!tools.empty())
	{
		for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			unsigned long preConfig;
			size_t idx;
			istringstream(itr->node().attribute("value").value())>>preConfig;
			istringstream(itr->node().parent().parent().parent().attribute("index").value())>>idx;
			config_indices[idx]=preConfig;
		}
	}
	preSetScanConfig(config_indices);

	vector<size_t> specSize(ns,0);
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum");
	if(!tools.empty())
	{
		size_t idx=0;
		size_t scanLength=0;
		for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			istringstream(itr->node().attribute("defaultArrayLength").value())>>scanLength;
			istringstream(itr->node().attribute("index").value())>>idx;
			specSize[idx]=scanLength;
		}
	}

    unsigned long instrument_type = 1;
	if(string(docmzML.child("mzML").child("instrumentConfigurationList").child("instrumentConfiguration").child("componentList").
				child("analyzer").child("cvParam").attribute("name").value()).compare("orbitrap") == 0)
		instrument_type=2;

	mzMLb3File.search_Group("spectrum_MS_1000514");
	mzMLb3File.search_Group("spectrum_MS_1000515");
	dataSetList=mzMLb3File.get_Info();
	dataMatDim = mzMLb3File.read_DimNC(dataSetList[0].varName,dataSetList[0].grpid);

    vector<spectrum> spectra(ns);
    for (size_t i = 0; i < ns; i++)
    {
		spectra[i].index = i;
		spectra[i].preset_config = config_indices[i];
		spectra[i].precursor_mz = precursor_mzs[i];
		spectra[i].scan_start_time = start_times[i];
		spectra[i].count = specSize[i];
    }
	// determine start_scan_time order of spectra
	sort(spectra.begin(), spectra.end(), scan_start_time_order);
	for (size_t i = 0; i < spectra.size(); i++)
	{
		spectra[i].scan_start_time_index = i;
	}
	vector<double> scan_start_times(ns);
	for (size_t i = 0; i < spectra.size(); i++)
	{
		scan_start_times[i] = spectra[i].scan_start_time;
	}
	sort(spectra.begin(), spectra.end(), seamass_order);
	vector< std::vector<double> > mzs(ns-1);
	vector< std::vector<double> > intensities(ns-1);
	vector<size_t> hypIdx(2);
	vector<size_t> rdLen(2);
	hypIdx[1]=0; // Read from first Column.
	rdLen[0]=1; // Always 1 Row to read.
	int loaded = 0;
    bool precursor_mz_is_constant = true;
	for (size_t i = 0; i < spectra.size(); i++)
	{
        if (loaded > 1 && spectra[i].precursor_mz != spectra[i-1].precursor_mz)
            precursor_mz_is_constant = false;

		if(spectra[i].count > 0 && spectra[i].index != spectra.size()-1)
		{
			hypIdx[0]=spectra[i].index;
			rdLen[1]=spectra[i].count;
			mzMLb3File.read_HypVecNC(dataSetList[0].varName,mzs[spectra[i].scan_start_time_index],&hypIdx[0],&rdLen[0],dataSetList[0].grpid);
			mzMLb3File.read_HypVecNC(dataSetList[1].varName,intensities[spectra[i].scan_start_time_index],&hypIdx[0],&rdLen[0],dataSetList[1].grpid);
		}

        loaded++;
		if ((i == spectra.size()-1 ||
			spectra[i].preset_config != spectra[i+1].preset_config ||
            (spectra[i].precursor_mz != spectra[i+1].precursor_mz && spectra[i].precursor_mz == 0.0)))
		{
            if (loaded > 1 && precursor_mz_is_constant)
			{
                ostringstream oss; oss << spectra[i].preset_config << "_" << spectra[i].precursor_mz;

				// run seamass 2D
				seamass::process(id,
					oss.str().c_str(),
					instrument_type,
					scan_start_times, mzs, intensities,
					mz_res, mz_res,
					rt_res, rt_res,
					shrink, shrink,
					tol, tol,
					threads, out_type);
			}
			else
			{
				// run seamass 1D (todo)
			}
            loaded = 0;
            precursor_mz_is_constant = true;
        }
	}

	return 0;
}

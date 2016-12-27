//
// $Id$
//
//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Bristol, UK
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

#include <algorithm>
#include <map>
#include <pugixml.hpp>
#include "MSFileData.hpp"
#include <sstream>

namespace xml = pugi;

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


MassSpecFile* FileFactory::createFileObj(string fName)
{
	size_t lastdot = fName.find_last_of(".")+1;
	string ext = (lastdot == string::npos) ? fName : fName.substr(lastdot);

	if(ext == "mzMLb3")
		return new MSmzMLb3(fName);
	if(ext == "mzMLb")
		return new MSmzMLb(fName);
    /*
    if(ext == "csv")
		retun new MZcsv(fName);
    */
	return NULL;
}


MSmzMLb3::MSmzMLb3(string fileName)
{
	mzMLb3File.open(fileName);
	hypIdx.resize(2,0);
	rdLen.resize(2,0);
	hypIdx[1]=0; // Read from first Column.
	rdLen[0]=1; // Always 1 Row to read.
    instrument_type = 1;
	this->extractData();
}


void MSmzMLb3::extractData(void)
{
	// Load mzML metadata
	mzMLb3File.read_VecNC("mzML",mzMLBuff);
	size_t xmlSize=sizeof(char)*mzMLBuff.size();

	xml::xml_document docmzML;
	xml::xpath_node_set tools;
	xml::xml_parse_result result = docmzML.load_buffer_inplace(&mzMLBuff[0],xmlSize);

	uli ns;
	istringstream(docmzML.child("mzML").child("run").child("spectrumList").attribute("count").value())>>ns;

	// query necessary metadata
    cout << "Querying metadata from " << ns << " spectra..." << endl;

	// capture scan start times 
	vector<double> start_times;
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/scanList/scan/cvParam[@accession='MS:1000016']");
	if(tools.empty())
	{
		cout<<"Cannot Find Start Time"<<endl;
		exit(1);
	}
	else
	{
		double rescale = 1.0;
		if(string(tools.first().node().attribute("unitAccession").value()).compare("UO:0000031") == 0)
			rescale = 60.0;

		for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			double scanRT;
			istringstream(itr->node().attribute("value").value()) >> scanRT;
			start_times.push_back(scanRT * rescale);
		}
	}

	// capture polarity
	vector<bool> positive_polarity(ns, true);
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/cvParam[@accession='MS:1000129']");
	if (!tools.empty())
	{
		for (xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			size_t idx;
			istringstream(itr->node().parent().attribute("index").value()) >> idx;
			positive_polarity[idx] = false;
		}
	}

	// capture precursor m/z (if MS1 then set as 0.0)
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

	// capture preset scan configurations (note cannot be compared across files)
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

	if(string(docmzML.child("mzML").child("instrumentConfigurationList").child("instrumentConfiguration").child("componentList").
				child("analyzer").child("cvParam").attribute("name").value()).compare("orbitrap") == 0)
		instrument_type=2;

	mzMLb3File.search_Group("spectrum_MS_1000514");
	mzMLb3File.search_Group("spectrum_MS_1000515");
	dataSetList=mzMLb3File.get_Info();

	spectraMetaData.resize(ns);
    for (size_t i = 0; i < ns; i++)
    {
		spectraMetaData[i].index = i;
		spectraMetaData[i].positive_polarity = positive_polarity[i];
		spectraMetaData[i].preset_config = config_indices[i];
		spectraMetaData[i].precursor_mz = precursor_mzs[i];
		spectraMetaData[i].start_time = start_times[i];
		spectraMetaData[i].count = specSize[i];
    }
}

vector<spectrumMetaData>& MSmzMLb3::getSpectraMetaData()
{
	return spectraMetaData;
}

void MSmzMLb3::getScanMZs(vector<double> &mz, size_t index, size_t count)
{
	hypIdx[0]=index;
	rdLen[1]=count;
	mzMLb3File.read_HypVecNC(dataSetList[0].varName,mz,&hypIdx[0],&rdLen[0],dataSetList[0].grpid);
}

void MSmzMLb3::getScanIntensities(vector<double> &intensities, size_t index, size_t count)
{
	hypIdx[0]=index;
	rdLen[1]=count;
	mzMLb3File.read_HypVecNC(dataSetList[1].varName,intensities,&hypIdx[0],&rdLen[0],dataSetList[1].grpid);
}

unsigned long MSmzMLb3::getInstrument(void)
{
	return instrument_type;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



MSmzMLb::MSmzMLb(string fileName)
{
	mzMLbFile.open(fileName);
	hypIdx.resize(1,0);
	rdLen.resize(1,0);
	hypIdx[0]=0; // Index to read from.
	rdLen[0]=0; // Length of Hyperslab.
	instrument_type = 1;
	this->extractData();
}


void MSmzMLb::extractData(void)
{
	// Load mzML metadata
	mzMLbFile.read_VecNC("mzML",mzMLBuff);
	size_t xmlSize=sizeof(char)*mzMLBuff.size();

	xml::xml_document docmzML;
	xml::xpath_node_set tools;
	xml::xml_parse_result result = docmzML.load_buffer_inplace(&mzMLBuff[0],xmlSize);

	uli ns;
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
		if(string(tools.first().node().attribute("unitAccession").value()).compare("UO:0000031") == 0)
			rescale = 60.0;


		for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			double scanRT;
			istringstream(itr->node().attribute("value").value()) >> scanRT;
			start_times.push_back(scanRT*rescale);
		}
	}

	// capture if profile or centroid
	vector<bool> profile_mode(ns, false);
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/cvParam[@accession='MS:1000128']");
	if (!tools.empty())
	{
		for (xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			size_t idx;
			istringstream(itr->node().parent().attribute("index").value()) >> idx;
			profile_mode[idx] = true;
		}
	}

	// capture polarity
	vector<bool> positive_polarity(ns, true);
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/cvParam[@accession='MS:1000129']");
	if (!tools.empty())
	{
		for (xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			size_t idx;
			istringstream(itr->node().parent().attribute("index").value()) >> idx;
			positive_polarity[idx] = false;
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

	string spectrum_mz;
	string spectrum_intensities;
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/binaryDataArrayList/binaryDataArray");
	mzIdx.resize(ns);
	intensitiesIdx.resize(ns);
	if(!tools.empty())
	{
		size_t idx=0;
		size_t jdx=0;
		size_t scanMzIdx=0;
		size_t scanIntIdx=0;

		string first = tools[0].node().child("binary").attribute("externalDataset").value();
		string second = tools[1].node().child("binary").attribute("externalDataset").value();

		if(first.find("spectrum_MS_1000514") != string::npos) spectrum_mz = first;
		else if(first.find("spectrum_MS_1000515") != string::npos) spectrum_intensities = first;

		if(second.find("spectrum_MS_1000514") != string::npos) spectrum_mz = second;
		else if(second.find("spectrum_MS_1000515") != string::npos) spectrum_intensities = second;

		for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
		{
			if(itr->node().child("binary").attribute("externalDataset").value() == spectrum_mz)
			{
				istringstream(itr->node().child("binary").attribute("offset").value())>>scanMzIdx;
				mzIdx[idx]=scanMzIdx;
			}
			if(itr->node().child("binary").attribute("externalDataset").value() == spectrum_intensities)
			{
				istringstream(itr->node().child("binary").attribute("offset").value())>>scanIntIdx;
				intensitiesIdx[idx]=scanIntIdx;
			}
			if(jdx%2 == 1) ++idx;
			++jdx;
		}
	}

	if(string(docmzML.child("mzML").child("instrumentConfigurationList").child("instrumentConfiguration").child("componentList").
				child("analyzer").child("cvParam").attribute("name").value()).compare("orbitrap") == 0)
		instrument_type=2;

	mzMLbFile.search_Group(spectrum_mz);
	mzMLbFile.search_Group(spectrum_intensities);
	dataSetList=mzMLbFile.get_Info();

	spectraMetaData.resize(ns);
    for (size_t i = 0; i < ns; i++)
    {
		spectraMetaData[i].index = i;
		spectraMetaData[i].positive_polarity = positive_polarity[i];
		spectraMetaData[i].profile_mode = profile_mode[i];
		spectraMetaData[i].preset_config = config_indices[i];
		spectraMetaData[i].precursor_mz = precursor_mzs[i];
		spectraMetaData[i].start_time = start_times[i];
		spectraMetaData[i].count = specSize[i];
    }
}

vector<spectrumMetaData>& MSmzMLb::getSpectraMetaData()
{
	return spectraMetaData;
}

void MSmzMLb::getScanMZs(vector<double> &mz, size_t index, size_t count)
{
	hypIdx[0]=mzIdx[index];
	rdLen[0]=count;
	mzMLbFile.read_HypVecNC(dataSetList[0].varName,mz,&hypIdx[0],&rdLen[0],dataSetList[0].grpid);
}

void MSmzMLb::getScanIntensities(vector<double> &intensities, size_t index, size_t count)
{
	hypIdx[0]=intensitiesIdx[index];
	rdLen[0]=count;
	mzMLbFile.read_HypVecNC(dataSetList[1].varName,intensities,&hypIdx[0],&rdLen[0],dataSetList[1].grpid);
}

unsigned long MSmzMLb::getInstrument(void)
{
	return instrument_type;
}

/*
vector<double> MSmzMLb::getScanStartTimes(void)
{
	return scan_start_times;

}
*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool mzMLbInputFile::scan_start_time_order(const spectrumMetaData& lhs, const spectrumMetaData& rhs)
{
	return lhs.start_time < rhs.start_time;
}


bool mzMLbInputFile::seamass_order(const spectrumMetaData& lhs, const spectrumMetaData& rhs)
{
	if (lhs.preset_config == rhs.preset_config)
	{
		if (lhs.precursor_mz == rhs.precursor_mz)
		{
			if (lhs.positive_polarity == rhs.positive_polarity)
			{
				return lhs.start_time < rhs.start_time;
			}
			else
			{
				return lhs.positive_polarity < rhs.positive_polarity;
			}
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


mzMLbInputFile::mzMLbInputFile(string in_file) :
	i(0)
{
	msFile = FileFactory::createFileObj(in_file);
	spectraMetaData = msFile->getSpectraMetaData();
	instrument_type = msFile->getInstrument();

	int lastdot = in_file.find_last_of(".");
	string id = (lastdot == string::npos) ? in_file : in_file.substr(0, lastdot);

	if (spectraMetaData.size() > 1)
	{
		// determine start_scan_time order of spectra
		// set finish time as start time of next spectrum, must discard last spectrum as no finish time
		sort(spectraMetaData.begin(), spectraMetaData.end(), &mzMLbInputFile::scan_start_time_order);
		for (size_t i = 0; i < spectraMetaData.size() - 1; i++)
		{
			spectraMetaData[i].finish_time = spectraMetaData[i + 1].start_time;
		}
		spectraMetaData.resize(spectraMetaData.size() - 1);

		// remove centroided spectra as we do not process these
		vector<spectrumMetaData>::iterator iter = spectraMetaData.begin();
		while (iter != spectraMetaData.end())
		{
			if (iter->profile_mode)
			{
				++iter;
			}
			else
			{
				// erase returns the new iterator
				iter = spectraMetaData.erase(iter);
			}
		}
        
		// finally sort into contiguous blocks for seaMass
		sort(spectraMetaData.begin(), spectraMetaData.end(), &mzMLbInputFile::seamass_order);
	}
}


mzMLbInputFile::~mzMLbInputFile()
{
	delete msFile;
}

bool mzMLbInputFile::next(SeamassCore::Input& out, std::string& id)
{
	if (i >= spectraMetaData.size()) return false;

	ostringstream oss;
	oss << (spectraMetaData[i].positive_polarity ? "pos_" : "neg_") << spectraMetaData[i].precursor_mz;
	id = oss.str();

	vector< vector<double> > mzs;
	vector< vector<double> > intensities;
	vector<double> start_times;
	vector<double> finish_times;

	size_t loaded = 0;
	bool done = false;
	for (; !done; i++)
	{
		mzs.resize(loaded + 1);
		intensities.resize(loaded + 1);
		start_times.resize(loaded + 1);
		finish_times.resize(loaded + 1);

		msFile->getScanMZs(mzs[loaded], spectraMetaData[i].index, spectraMetaData[i].count);
		msFile->getScanIntensities(intensities[loaded], spectraMetaData[i].index, spectraMetaData[i].count);
		start_times[loaded] = spectraMetaData[i].start_time;
		finish_times[loaded] = spectraMetaData[i].finish_time;

		loaded++;

		if (i == spectraMetaData.size() - 1 ||
		    spectraMetaData[i].preset_config != spectraMetaData[i + 1].preset_config ||
			spectraMetaData[i].precursor_mz != spectraMetaData[i + 1].precursor_mz ||
			spectraMetaData[i].positive_polarity != spectraMetaData[i + 1].positive_polarity)
		{
			done = true;
		}
	}

	out.binCounts.resize(0);
	out.binEdges.resize(0);
	out.spectrumIndex.resize(intensities.size() > 1 ? intensities.size() + 1 : 0);
	out.startTimes.resize(intensities.size() > 1 ? intensities.size() : 0);
	out.finishTimes.resize(intensities.size() > 1 ? intensities.size() : 0);
	out.exposures.resize(intensities.size(), 1.0); // initialise exposures to default of 1 (unit exposure)
	size_t bck = 0;
	size_t blk = 0;
	for (size_t j = 0; j < intensities.size(); j++)
	{
		if (intensities.size() > 1)
		{
			out.spectrumIndex[j] = out.binCounts.size();
			out.startTimes[j] = start_times[j];
			out.finishTimes[j] = finish_times[j];
		}

		// This modifies the raw data for some limitations of the mzML spec and makes
		// sure the intensities are treated as binned between m/z datapoints.
		//
		// For ToF data, it interpolates the mzs to represent the edges of the bins rather
		//   than the centres
		// For FT data, it converts the sampled data to binned counts by treating the
		//   mz values as the bin edges, and using trapezoid rule to integrate intensity values
		//
		// This is all a bit rough at the moment, should be fitting splines to the data

		if (instrument_type == 1) // ToF
		{
			if (mzs[j].size() >= 2)
			{
				// dividing by minimum to get back to ion counts for SWATH data which appears to be automatic gain controlled to correct for dynamic range restrictions (hack!)
				double minimum = std::numeric_limits<double>::max();

				// we drop the first and last m/z datapoint as we don't know both their bin edges
				out.binCounts.resize(out.binCounts.size() + intensities[j].size() - 2);
				out.binEdges.resize(out.binEdges.size() + mzs[j].size() - 1);
				for (size_t k = 1; k < mzs[j].size(); k++)
				{
					if (intensities[j][k] > 0) minimum = minimum < intensities[j][k] ? minimum : intensities[j][k];

					// linear interpolation of mz extent (probably should be cubic)
					if (k < mzs[j].size() - 1) out.binCounts[bck++] = intensities[j][k];
					out.binEdges[blk++] = 0.5 * (mzs[j][k - 1] + mzs[j][k]);
				}

				// check to see if we can estimate the exposure i.e. is the minimum reasonable?
				if (minimum >= 1.0 && minimum <= 1000.0)
				{
					// correct bin_counts with derived exposures
					for (size_t k = bck - (intensities[j].size() - 2); k < bck; k++)
					{
						out.binCounts[k] /= minimum;
					}
					out.exposures[j] = 1.0 / minimum;
				}
			}
		}
		else // Orbitrap & FT-ICR
		{
			if (mzs[j].size() >= 2)
			{
				out.binCounts.resize(out.binCounts.size() + intensities[j].size() - 1);
				out.binEdges.resize(out.binEdges.size() + mzs[j].size());
				for (size_t k = 0; k < mzs[j].size(); k++)
				{
					if (k > 0) out.binCounts[bck++] = (mzs[j][k] - mzs[j][k - 1]) * 0.5 * (intensities[j][k] + intensities[j][k - 1]);
					out.binEdges[blk++] = mzs[j][k];
				}
			}
		}
	}
	if (intensities.size() > 1) out.spectrumIndex.back() = out.binCounts.size();

	return true;
}

vector<spectrumMetaData>* mzMLbInputFile::getSpectrumMetaData()
{
    return &this->spectraMetaData;
}

MassSpecFile* mzMLbInputFile::getGeometry()
{
    return msFile;
}

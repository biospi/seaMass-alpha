//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
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

#include <algorithm>
#include <map>
#include <pugixml.hpp>
#include "DatasetMzmlb.hpp"


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


// todo: use spectrum index
DatasetMzmlb::DatasetMzmlb(string& fileName) : spectrumIndex_(0)
{
	// Load mzML metadata
    vector<char> mzMLBuff;
    file_.open(fileName);
	file_.read_VecNC("mzML", mzMLBuff);
	size_t xmlSize = sizeof(char) * mzMLBuff.size();

	xml::xml_document docmzML;
	xml::xpath_node_set tools;
	xml::xml_parse_result result = docmzML.load_buffer_inplace(&mzMLBuff[0], xmlSize);

    // query necessary metadata
	uli ns;
	istringstream(docmzML.child("mzML").child("run").child("spectrumList").attribute("count").value()) >> ns;

	cout << "Querying metadata from " << ns << " spectra..." << endl;

    metadata_.resize(ns);
    for (size_t i = 0; i < ns; i++)
    {
        metadata_[i].index = i;
        metadata_[i].isProfileMode = false;
        metadata_[i].isPositivePolarity = true;
    }

	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/scanList/scan/cvParam[@accession='MS:1000016']");
	if(!tools.empty())
	{
		double rescale=1.0;
		if(string(tools.first().node().attribute("unitAccession").value()).compare("UO:0000031") == 0)
        {
            rescale = 60.0;
        }

        for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
        {
            size_t idx;
            istringstream(itr->node().parent().parent().parent().attribute("index").value()) >> idx;
            double scanRT;
            istringstream(itr->node().attribute("value").value()) >> scanRT;
            metadata_[idx].startTime = scanRT * rescale;
        }
	}

	// capture if profile or centroid
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/cvParam[@accession='MS:1000128']");
	if (!tools.empty())
	{
        for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
        {
            size_t idx;
            istringstream(itr->node().parent().attribute("index").value()) >> idx;
            metadata_[idx].isProfileMode = true;
        }
	}

	// capture polarity
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/cvParam[@accession='MS:1000129']");
	if (!tools.empty())
	{
        for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
        {
            size_t idx;
            istringstream(itr->node().parent().attribute("index").value()) >> idx;
            metadata_[idx].isPositivePolarity = false;
        }
	}

	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/precursorList/precursor");
	if(!tools.empty())
	{
        for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
        {
            size_t idx;
            istringstream(itr->node().parent().parent().attribute("index").value()) >> idx;
            double preMZ;
            istringstream(itr->node().child("selectedIonList").child("selectedIon").
                    find_child_by_attribute("accession", "MS:1000744").attribute("value").value()) >> preMZ;
            metadata_[idx].precursorMz = preMZ;
        }
	}

	vector<unsigned long> presetConfigs(ns, 0);
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/scanList/scan/cvParam[@accession='MS:1000616']");
	if(!tools.empty())
	{
        for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
        {
            unsigned long preConfig;
            size_t idx;
            istringstream(itr->node().attribute("value").value()) >> preConfig;
            istringstream(itr->node().parent().parent().parent().attribute("index").value()) >> idx;
            presetConfigs[idx] = preConfig;
        }
	}
	preSetScanConfig(presetConfigs);
    for (size_t i = 0; i < ns; i++)
    {
        metadata_[i].presetConfig = presetConfigs[i];
    }

	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum");
	if(!tools.empty())
	{
		size_t idx=0;
		size_t arrayLength=0;
        for(xml::xpath_node_set::const_iterator itr = tools.begin(); itr != tools.end(); ++itr)
        {
            istringstream(itr->node().attribute("defaultArrayLength").value()) >> arrayLength;
            istringstream(itr->node().attribute("index").value()) >> idx;
            metadata_[idx].arrayLength = arrayLength;
        }
	}

	string spectrum_mz;
	string spectrum_intensities;
	tools = docmzML.select_nodes("mzML/run/spectrumList/spectrum/binaryDataArrayList/binaryDataArray");
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
                metadata_[idx].mzIdx = scanMzIdx;
			}
			if(itr->node().child("binary").attribute("externalDataset").value() == spectrum_intensities)
			{
				istringstream(itr->node().child("binary").attribute("offset").value())>>scanIntIdx;
                metadata_[idx].intensitiesIdx = scanIntIdx;
			}
			if(jdx%2 == 1) ++idx;
			++jdx;
		}
	}

	if(string(docmzML.child("mzML").child("instrumentConfigurationList").child("instrumentConfiguration").child("componentList").
			child("analyzer").child("cvParam").attribute("name").value()).compare("orbitrap") == 0)
    {
        instrumentType_ = 2;
    }
    else
    {
        instrumentType_ = 1;
    }

	file_.search_Group(spectrum_mz);
	file_.search_Group(spectrum_intensities);
	dataSetList_=file_.get_Info();

    // sort spectra into appropriate order for next()
    if (metadata_.size() > 1)
    {
        // determine start time order of spectra
        // set finish time as start time of next spectrum, must discard last spectrum as no finish time
        sort(metadata_.begin(), metadata_.end(), &DatasetMzmlb::startTimeOrder);
        for (size_t i = 0; i < metadata_.size() - 1; i++)
        {
            metadata_[i].finishTime = metadata_[i + 1].startTime;
        }
        metadata_.resize(metadata_.size() - 1);

        // remove centroided spectra as we do not process these
        vector<MetadataMzmlbSpectrum>::iterator iter = metadata_.begin();
        while (iter != metadata_.end())
        {
            if (iter->isProfileMode)
            {
                ++iter;
            }
            else
            {
                // erase returns the new iterator
                iter = metadata_.erase(iter);
            }
        }

        // finally sort into contiguous blocks for seaMass
        sort(metadata_.begin(), metadata_.end(), &DatasetMzmlb::seamassOrder);
    }
}


DatasetMzmlb::~DatasetMzmlb()
{
}


bool DatasetMzmlb::next(SeamassCore::Input& out, std::string& id)
{
	if (spectrumIndex_ >= metadata_.size()) return false;

	ostringstream oss;
	oss << (metadata_[spectrumIndex_].isPositivePolarity ? "pos_" : "neg_") << metadata_[spectrumIndex_].precursorMz;
	id = oss.str();

    vector<fp>().swap(out.binCounts);
    vector<li>().swap(out.spectrumIndex);
    vector<double>().swap(out.binEdges);
    vector<double>().swap(out.startTimes);
    vector<double>().swap(out.finishTimes);
    vector<fp>().swap(out.exposures);

    size_t bck = 0;
    size_t blk = 0;
    bool done = false;
	for (; !done; spectrumIndex_++)
	{
        out.spectrumIndex.push_back((li)out.binCounts.size());
        out.startTimes.push_back(metadata_[spectrumIndex_].startTime);
        out.finishTimes.push_back(metadata_[spectrumIndex_].finishTime);
        out.exposures.push_back(1.0);

        vector<double> mzs, intensities;
        getScanMZs(mzs, metadata_[spectrumIndex_].index, metadata_[spectrumIndex_].arrayLength);
        getScanIntensities(intensities, metadata_[spectrumIndex_].index, metadata_[spectrumIndex_].arrayLength);

        // This modifies the raw data for some limitations of the mzML spec and makes
        // sure the intensities are treated as binned between m/z datapoints.
        //
        // For ToF data, it interpolates the mzs to represent the edges of the bins rather
        //   than the centres
        // For FT data, it converts the sampled data to binned counts by treating the
        //   mz values as the bin edges, and using trapezoid rule to integrate intensity values
        //
        // This is all a bit rough at the moment, should be fitting splines to the data
        if (instrumentType_ == 1) // ToF
        {
            if (mzs.size() >= 2)
            {
                // dividing by minimum to get back to ion counts for SWATH data which appears to be automatic gain controlled to correct for dynamic range restrictions (hack!)
                double minimum = std::numeric_limits<double>::max();

                // we drop the first and last m/z datapoint as we don't know both their bin edges
                out.binCounts.resize(out.binCounts.size() + intensities.size() - 2);
                out.binEdges.resize(out.binEdges.size() + mzs.size() - 1);
                for (size_t k = 1; k < mzs.size(); k++)
                {
                    if (intensities[k] > 0) minimum = minimum < intensities[k] ? minimum : intensities[k];

                    // linear interpolation of mz extent (probably should be cubic)
                    if (k < mzs.size() - 1) out.binCounts[bck++] = (fp)intensities[k];
                    out.binEdges[blk++] = 0.5 * (mzs[k - 1] + mzs[k]);
                }

                // check to see if we can estimate the exposure i.e. is the minimum reasonable?
                if (minimum >= 1.0 && minimum <= 1000.0)
                {
                    // correct bin_counts with derived exposures
                    for (size_t k = bck - (intensities.size() - 2); k < bck; k++)
                    {
                        out.binCounts[k] /= minimum;
                    }
                    out.exposures.back() = (fp)(1.0 / minimum);
                }
            }
        }
        else // Orbitrap & FT-ICR
        {
            if (mzs.size() >= 2)
            {
                out.binCounts.resize(out.binCounts.size() + intensities.size() - 1);
                out.binEdges.resize(out.binEdges.size() + mzs.size());
                for (size_t k = 0; k < mzs.size(); k++)
                {
                    if (k > 0) out.binCounts[bck++] = (fp)((mzs[k] - mzs[k - 1]) * 0.5 * (intensities[k] + intensities[k - 1]));
                    out.binEdges[blk++] = mzs[k];
                }
            }
        }

        // if at the end or any metadata is different in the next spectrum, we have come to the end of this channel
		if (spectrumIndex_ == metadata_.size() - 1 ||
		    metadata_[spectrumIndex_].presetConfig != metadata_[spectrumIndex_ + 1].presetConfig ||
			metadata_[spectrumIndex_].precursorMz != metadata_[spectrumIndex_ + 1].precursorMz ||
			metadata_[spectrumIndex_].isPositivePolarity != metadata_[spectrumIndex_ + 1].isPositivePolarity)
		{
			done = true;
		}
	}

	return true;
}


void DatasetMzmlb::getScanMZs(vector<double> &mz, size_t index, size_t count)
{
 	if(count > 0)
    {
        vector<size_t> hypIdx(1); hypIdx[0] = index;
        vector<size_t> rdLen(1); rdLen[0] = count;

        file_.read_HypVecNC(dataSetList_[0].varName, mz, &hypIdx[0], &rdLen[0], dataSetList_[0].grpid);
    }
}


void DatasetMzmlb::getScanIntensities(vector<double> &intensities, size_t index, size_t count)
{
	if(count > 0)
    {
        vector<size_t> hypIdx(1); hypIdx[0] = index;
        vector<size_t> rdLen(1); rdLen[0] = count;

        file_.read_HypVecNC(dataSetList_[1].varName, intensities, &hypIdx[0], &rdLen[0], dataSetList_[1].grpid);
    }
}


bool DatasetMzmlb::startTimeOrder(const MetadataMzmlbSpectrum &lhs, const MetadataMzmlbSpectrum &rhs)
{
	return lhs.startTime < rhs.startTime;
}


bool DatasetMzmlb::seamassOrder(const MetadataMzmlbSpectrum &lhs, const MetadataMzmlbSpectrum &rhs)
{
	if (lhs.presetConfig == rhs.presetConfig)
	{
		if (lhs.precursorMz == rhs.precursorMz)
		{
			if (lhs.isPositivePolarity == rhs.isPositivePolarity)
			{
				return lhs.startTime < rhs.startTime;
			}
			else
			{
				return lhs.isPositivePolarity < rhs.isPositivePolarity;
			}
		}
		else
		{
			return lhs.precursorMz < rhs.precursorMz;
		}
	}
	else
	{
		return lhs.presetConfig < rhs.presetConfig;
	}
}
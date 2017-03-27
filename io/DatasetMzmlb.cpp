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

#include "DatasetMzmlb.hpp"
#include <algorithm>
#include <iomanip>
#include <map>
#include <cassert>
#include <pugixml.hpp>


namespace xml = pugi;


DatasetMzmlb::DatasetMzmlb(string& fileName) : spectrumIndex_(0)
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "Opening " << fileName << " ..." << endl;

    file_.open(fileName);

    // Load "mzML_spectrumIndex"
    vector<li> mzML_spectrumIndex;
    file_.read_VecNC("mzML_spectrumIndex", mzML_spectrumIndex);

    // Load "mzML_chromatogramIndex"
    vector<li> mzML_chromatogramIndex;
    file_.read_VecNC("mzML_chromatogramIndex", mzML_chromatogramIndex);

    // Load "mzML" but without spectra
    vector<char> mzML;
    {
        vector<char> buffer;
        size_t offset = 0;
        size_t extent = mzML_spectrumIndex.front();
        file_.read_HypVecNC("mzML", buffer, &offset, &extent);
        mzML.insert(mzML.end(), buffer.begin(), buffer.end());

        offset = mzML_spectrumIndex.back();
        extent = mzML_chromatogramIndex.front() - offset;
        file_.read_HypVecNC("mzML", buffer, &offset, &extent);
        mzML.insert(mzML.end(), buffer.begin(), buffer.end());

        offset = mzML_chromatogramIndex.back();
        vector<size_t> dims = file_.read_DimNC("mzML");
        extent = dims[0] - offset;
        file_.read_HypVecNC("mzML", buffer, &offset, &extent);
        mzML.insert(mzML.end(), buffer.begin(), buffer.end());
    }

    // parse
	xml::xml_document mzmlDoc;
	xml::xpath_node_set nodes;
    xml::xml_parse_result result = mzmlDoc.load_buffer_inplace(&mzML[0], sizeof(char) * mzML.size());
    if (!result) throw runtime_error("Error: In mzMLb input - " + string(result.description()));

    // query all the detectors
    map<string, MetadataMzmlbSpectrum::DataType> dataTypes;
    nodes = mzmlDoc.select_nodes("mzML/instrumentConfigurationList/instrumentConfiguration");
    for(xml::xpath_node_set::const_iterator itr = nodes.begin(); itr != nodes.end(); ++itr)
    {
        string id;
        istringstream(itr->node().attribute("id").value()) >> id;

        MetadataMzmlbSpectrum::DataType dataType = MetadataMzmlbSpectrum::Unknown;
        xml::xpath_node_set detectors = itr->node().select_nodes("componentList/detector");
        for(xml::xpath_node_set::const_iterator itrDetectors = detectors.begin(); itrDetectors != detectors.end(); ++itrDetectors)
        {
            // is it an "electron multiplier" detector?
            xml::xpath_node_set detectorTypes = itrDetectors->node().select_nodes("cvParam[@accession='MS:1000253']");
            if(!detectorTypes.empty())
            {
                if (getDebugLevel() % 10 >= 1)
                    cout << getTimeStamp() << " Detector: electron multiplier (ion count)" << endl;
                if (dataType == MetadataMzmlbSpectrum::Unknown || dataType == MetadataMzmlbSpectrum::IonCount)
                    dataType = MetadataMzmlbSpectrum::IonCount;
                else
                    throw runtime_error("Error: Inconsistent <detector> types found in mzMLb");
            }
            // is it an "photomultiplier" detector?
            detectorTypes = itrDetectors->node().select_nodes("cvParam[@accession='MS:1000116']");
            if(!detectorTypes.empty())
            {
                if (getDebugLevel() % 10 >= 1)
                    cout << getTimeStamp() << " Detector: photomultiplier (ion count)" << endl;
                if (dataType == MetadataMzmlbSpectrum::Unknown || dataType == MetadataMzmlbSpectrum::IonCount)
                    dataType = MetadataMzmlbSpectrum::IonCount;
                else
                    throw runtime_error("Error: Inconsistent <detector> types found in mzMLb");
            }
            // is it an "induction" detector (i.e. orbitrap or FT-ICR)?
            detectorTypes = itrDetectors->node().select_nodes("cvParam[@accession='MS:1000624']");
            if(!detectorTypes.empty())
            {
                if (getDebugLevel() % 10 >= 1)
                    cout << getTimeStamp() << " Detector: inductive (ion current)" << endl;
                if (dataType == MetadataMzmlbSpectrum::Unknown || dataType == MetadataMzmlbSpectrum::IonCurrent)
                    dataType = MetadataMzmlbSpectrum::IonCurrent;
                else
                    throw runtime_error("Error: Inconsistent <detector> types found in mzMLb");
            }
        }

        dataTypes[id] = dataType;
    }

    // query default detector
    MetadataMzmlbSpectrum::DataType defaultDataType = dataTypes[mzmlDoc.select_node("mzML/run").node().attribute("defaultInstrumentConfigurationRef").value()];

    // query number of spectra to create "metadata_" structure
	li ns;
	istringstream(mzmlDoc.child("mzML").child("run").child("spectrumList").attribute("count").value()) >> ns;
    metadata_.resize(ns);

    if (getDebugLevel() % 10 >= 1)
 	    cout << getTimeStamp() << " Querying spectrum metadata ..." << endl;
    for (li i = 0; i < ns; i++)
    {
       // load each spectrum's mzML to get rest of info
        size_t offset = mzML_spectrumIndex[i];
        size_t extent = mzML_spectrumIndex[i + 1] - offset;
        file_.read_HypVecNC("mzML", mzML, &offset, &extent);

        //string str(mzML.begin(), mzML.end());
        //cout << str << endl;

        xml::xml_parse_result result = mzmlDoc.load_buffer_inplace(&mzML[0], sizeof(char) * mzML.size());
        if (!result) throw runtime_error("Error: In mzMLb input file - " + string(result.description()));

        // capture index
        istringstream(mzmlDoc.select_node("spectrum").node().attribute("index").value()) >> metadata_[i].mzmlSpectrumIndex;

        // capture mz and intensity array length
        istringstream(mzmlDoc.select_node("spectrum").node().attribute("defaultArrayLength").value()) >> metadata_[i].defaultArrayLength;

        // capture if profile mode spectrum
        nodes = mzmlDoc.select_nodes("spectrum/cvParam[@accession='MS:1000128']");
        if (!nodes.empty())
            metadata_[i].isProfileMode = true;
        else
            metadata_[i].isProfileMode = false;

        // capture polarity
        nodes = mzmlDoc.select_nodes("spectrum/cvParam[@accession='MS:1000129']");
        if (!nodes.empty())
            metadata_[i].id = "neg";
        else
        {
            nodes = mzmlDoc.select_nodes("spectrum/cvParam[@accession='MS:1000130']");
            if (!nodes.empty())
                metadata_[i].id = "pos";
            else
                metadata_[i].id = "unk";
        }

        nodes = mzmlDoc.select_nodes("spectrum/scanList/scan");
        if(!nodes.empty())
        {
            if (nodes.size() > 1)
                throw runtime_error("Error: We don't know how to process multiple <scan> per <spectrum> in mzMLb input file");

            // capture spectrum detector
            if (!nodes.first().node().attribute("instrumentConfigurationRef").empty())
                metadata_[i].dataType = dataTypes[nodes.first().node().attribute("instrumentConfigurationRef").value()];
            else
                metadata_[i].dataType = defaultDataType;

             // capture start time
            xml::xpath_node_set scantimes = nodes.first().node().select_nodes("cvParam[@accession='MS:1000016']");
            if(!scantimes.empty())
            {
                istringstream(scantimes.first().node().attribute("value").value()) >> metadata_[i].startTime;

                if(string(scantimes.first().node().attribute("unitAccession").value()).compare("UO:0000031") == 0)
                    metadata_[i].startTime *= 60.0;
            }
            else
                metadata_[i].startTime = -1.0;

            // capture spectrum window lower limit
            /*xml::xpath_node_set mz0s = nodes.first().node().select_nodes("scanWindowList/scanWindow/cvParam[@accession='MS:1000501']");
            if(!nodes.empty())
            {
                if (mz0s.size() > 1)
                    throw runtime_error("Error: We don't know how to process multiple <scanWindow> per <scanWindowList> in mzMLb input file");

                istringstream(mz0s.first().node().attribute("value").value()) >> metadata_[i].mz0;
            }
            else
                metadata_[i].mz0 = -1.0;

            // capture spectrum window upper limit
            xml::xpath_node_set mz1s = nodes.first().node().select_nodes("scanWindowList/scanWindow/cvParam[@accession='MS:1000500']");
            if(!nodes.empty())
            {
                if (mz1s.size() > 1)
                    throw runtime_error("Error: We don't know how to process multiple <scanWindow> per <scanWindowList> in mzMLb input file");

                istringstream(mz1s.first().node().attribute("value").value()) >> metadata_[i].mz1;
            }
            else
                metadata_[i].mz1 = -1.0;*/
        }
        else
            throw runtime_error("Error: <spectrum> missing <scan> in mzML input");

        nodes = mzmlDoc.select_nodes("spectrum/precursorList/precursor/selectedIonList/selectedIon/cvParam[@accession='MS:1000744']");
        if(!nodes.empty())
        {
            // capture precursor mz
            metadata_[i].id += "_";
            metadata_[i].id += nodes.first().node().attribute("value").value();

            // use precursor intensity as heuristic as to whether this is DDA (precursor intensity present) or DIA (not present)
            nodes = nodes.first().parent().select_nodes("cvParam[@accession='MS:1000042']");
            if(!nodes.empty())
            {
                metadata_[i].id += "_";
                metadata_[i].id += metadata_[i].mzmlSpectrumIndex;
            }
        }

        // capture dataset and offset of mzs
        nodes = mzmlDoc.select_nodes("spectrum/binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary");
        if(!nodes.empty())
        {
            istringstream(nodes.first().node().attribute("externalDataset").value()) >> metadata_[i].mzsDataset;
            istringstream(nodes.first().node().attribute("offset").value()) >> metadata_[i].mzsOffset;
        }
        else
            throw runtime_error("Error: No <binary> m/z data in mzMLb input file");

        // capture dataset and offset of intensities
        nodes = mzmlDoc.select_nodes("spectrum/binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary");
        if(!nodes.empty())
        {
            istringstream(nodes.first().node().attribute("externalDataset").value()) >> metadata_[i].intensitiesDataset;
            istringstream(nodes.first().node().attribute("offset").value()) >> metadata_[i].intensitiesOffset;
        }
        else
            throw runtime_error("Error: No <binary> intensity data in mzMLb input file");

        if (getDebugLevel() % 10 >= 5)
        {
            cout << getTimeStamp() << "  " << metadata_[i].mzmlSpectrumIndex << " id=" << metadata_[i].id << endl;
            cout << getTimeStamp() << "    Intensities dataset=" << metadata_[i].intensitiesDataset;
            cout << " offset=" << metadata_[i].intensitiesOffset;
            cout << " extent=" << metadata_[i].defaultArrayLength;
            cout << " profile=" << (metadata_[i].isProfileMode ? "true" : "false") << endl;
            cout << getTimeStamp() << "    Mzs dataset=" << metadata_[i].mzsDataset;
            cout << " offset=" << metadata_[i].mzsOffset;
            cout << " extent=" << metadata_[i].defaultArrayLength << endl;
            //cout << " window=[" << metadata_[i].mz0 << "," << metadata_[i].mz1 << "]Th" << endl;
            cout << getTimeStamp() << "    start_time=" << metadata_[i].startTime << "s" << endl;
        }

        // display progress update
        if (getDebugLevel() % 10 >= 1)
        {
            if ((i + 1) % 10000 == 0 || (i + 1) == ns)
            {
                cout << getTimeStamp() << "   " << setw(1 + (int)(log10((float)ns))) << (i+1) << "/" << ns << endl;
            }
        }
   }

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
    else
    {
        metadata_[0].finishTime = -1.0; // startTime and finishTime not used anyway for 1D seaMass
    }
}


DatasetMzmlb::~DatasetMzmlb()
{
}


bool DatasetMzmlb::next(SeamassCore::Input& out, std::string& id)
{
	if (spectrumIndex_ >= metadata_.size()) return false;

    id = metadata_[spectrumIndex_].id;
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << " Reading spectra with id=" << id << " ..." << endl;

    // determine next set of spectra
    bool done = false;
    size_t offset = spectrumIndex_;
    for (; !done; spectrumIndex_++)
    {
        // if at the end or any metadata is different in the next spectrum, we have come to the end of this channel
        if (spectrumIndex_ == metadata_.size() - 1 || metadata_[spectrumIndex_].id != metadata_[spectrumIndex_ + 1].id)
        {
            out.binCountsIndex.push_back((li)out.binCounts.size());
            done = true;
        }
    }
    size_t extent = spectrumIndex_ - offset;

    // read spectra
    vector< vector<double> > mzs(extent);
    vector< vector<fp> > intensities(extent);
	for (ii i = 0; i < extent; i++)
	{
       if(metadata_[offset + i].defaultArrayLength > 0)
        {
            size_t rdLen = metadata_[offset + i].defaultArrayLength;
            size_t hypIdx = metadata_[offset + i].mzsOffset;
            file_.read_HypVecNC(metadata_[offset + i].mzsDataset, mzs[i], &hypIdx, &rdLen);
            hypIdx = metadata_[offset + i].intensitiesOffset;
            file_.read_HypVecNC(metadata_[offset + i].intensitiesDataset, intensities[i], &hypIdx, &rdLen);
        }

        // display progress update
        if (getDebugLevel() % 10 >= 1)
        {
            if ((i + 1) % 1000 == 0 || (i + 1) == extent)
            {
                cout << getTimeStamp() << "   " << setw(1 + (int)(log10((float)extent))) << (i+1) << "/" << extent << endl;
            }
        }
    }

    // estimate of mz_range [mz0,mz1]
    // mzML spec does not guarentee this will exist (and hence do not know if regions have zero counts or were not scanned)
    double mz0 = numeric_limits<double>::max();
    double mz1 = 0.0;
    for (ii i = 0; i < (ii)mzs.size(); i++)
    {
        if (mzs[i].size() >= 2)
        {
            mz0 = mz0 < mzs[i].front() ? mz0 : mzs[i].front();
            mz1 = mz1 > mzs[i].back() ? mz1 : mzs[i].back();
        }
    }
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "   autodetected_mz_window=[" << mz0 << "," << mz1 << "]Th" << endl;
    // if this happens regularly, will have to use reported scan window or command line arguments
    if (mz1 <= mz0 || mz0 == numeric_limits<double>::max() || mz0 == 0.0)
        throw runtime_error("Error: Not enough datapoints in input mzMLb");

    // reinit out
    vector<li>(extent + 1).swap(out.binCountsIndex);
    vector<double>(extent).swap(out.startTimes);
    vector<double>(extent).swap(out.finishTimes);
    vector<fp>(extent).swap(out.exposures);
    vector<fp>().swap(out.binCounts);
    vector<double>().swap(out.binEdges);

    // This modifies the raw data for some limitations of the mzML spec and makes
    // sure the intensities are treated as binned between m/z datapoints.
    //
    // For IonCount data, it interpolates the mzs to represent the edges of the bins rather
    //   than the centres
    // For IonCurrent data, it converts the sampled data to binned counts by treating the
    //   mz values as the bin edges, and using trapezoid rule to integrate intensity values
    //
    size_t bck = 0;
    size_t blk = 0;
    for (ii i = 0; i < extent; i++)
    {
        out.binCountsIndex[i] = (li)out.binCounts.size();
        out.startTimes[i] = metadata_[offset + i].startTime;
        out.finishTimes[i] = metadata_[offset + i].finishTime;

        switch (metadata_[offset + i].dataType)
        {
            case MetadataMzmlbSpectrum::Unknown:
            case MetadataMzmlbSpectrum::IonCount: // ToF, Quad, Ion trap etc
            {
                if (intensities[i].size() == 0) // if empty we MUST add a zero bin!
                {
                    out.binEdges.push_back(mz0);
                    out.binCounts.push_back(0.0);
                    out.binEdges.push_back(mz1);

                    out.exposures[i] = 1.0;
                }
                else
                {
                    // dividing by minimum to get back to ion counts for SWATH data which appears to be automatic gain
                    // controlled to correct for dynamic range restrictions (hack!)
                    double minimum = std::numeric_limits<double>::max();
                    for (size_t k = 0; k < mzs[i].size(); k++)
                        if (intensities[i][k] > 0.0) minimum = minimum < intensities[i][k] ? minimum : intensities[i][k];
                    // check to see if we can estimate the exposure i.e. is the minimum reasonable?
                    if (minimum < 1.0 || minimum > 1000.0) minimum = 1.0;
                    out.exposures[i] = (fp)(1.0 / minimum);

                    if (intensities[i].front() == 0.0) // only use the first m/z if the intensity is zero
                    {
                        double frontEdge = mz0 < mzs[i].front() ? mz0 : mzs[i].front();
                        out.binEdges.push_back(frontEdge);
                        out.binCounts.push_back(0.0);
                    }

                    // use all the intensities except first and last
                    for (ii k = 1; k < (ii)intensities[i].size() - 1; k++)
                    {
                        if (intensities[i][k] != 0.0 || intensities[i][k - 1] != 0.0) // merge zeros
                        {
                            out.binEdges.push_back(0.5 * (mzs[i][k - 1] + mzs[i][k]));
                            out.binCounts.push_back(((fp)intensities[i][k]) * out.exposures[i]);
                        }
                    }
                    out.binEdges.push_back(0.5 * (mzs[i][mzs[i].size() - 2] + mzs[i].back()));

                    if (intensities[i].back() == 0.0) // only use the last m/z if the intensity is zero
                    {
                        out.binCounts.push_back(0.0);
                        double backEdge = mz1 > mzs[i].back() ? mz1 : mzs[i].back();
                        out.binEdges.push_back(backEdge);
                    }
                }
            } break;
            case MetadataMzmlbSpectrum::IonCurrent: // Orbitrap, FT-ICR etc
            {
                throw runtime_error("Error: Orbitrap support not implemented yet");
                /*if (mzs[i].size() >= 2)
                 {
                     out.binCounts.resize(out.binCounts.size() + intensities[i].size() - 1);
                     out.binEdges.resize(out.binEdges.size() + mzs[i].size());
                     for (size_t k = 0; k < mzs[i].size(); k++)
                     {
                         if (k > 0) out.binCounts[bck++] = (fp)((mzs[i][k] - mzs[i][k - 1]) * 0.5 * (intensities[i][k] + intensities[i][k - 1]));
                         out.binEdges[blk++] = mzs[i][k];
                     }
                 }*/
            } break;
        }

        vector<double>().swap(mzs[i]);
        vector<fp>().swap(intensities[i]);
    }
    out.binCountsIndex.back() = (li)out.binCounts.size();

    return true;
}


bool DatasetMzmlb::startTimeOrder(const MetadataMzmlbSpectrum &lhs, const MetadataMzmlbSpectrum &rhs)
{
	return lhs.startTime < rhs.startTime;
}


bool DatasetMzmlb::seamassOrder(const MetadataMzmlbSpectrum &lhs, const MetadataMzmlbSpectrum &rhs)
{
    if (lhs.id == rhs.id)
    {
        return lhs.startTime < rhs.startTime;
    }
    else
    {
        return lhs.id < rhs.id;
    }
}
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
#include <kernel.hpp>
#include <algorithm>
#include <iomanip>
#include <map>
using namespace kernel;
namespace xml = pugi;


DatasetMzmlb::DatasetMzmlb(const std::string& filePathIn, const std::string& filePathStemOut, std::vector<short>& scale, Dataset::WriteType writeType) : scale_(scale),
                        fileOut_(0), spectrumIndex_(0), lastSpectrumIndex_(-1000), spectrumListIdx_(0)
{
    if (filePathIn.empty())
        throw runtime_error("BUG: mzMLb/mzMLv file cannot be written without an mzMLb/mzMLv to readMatrixSparse from.");

    ostringstream oss;
    if (getDebugLevel() % 10 >= 1)
        oss << getTimeStamp();
    oss << "Querying " << filePathIn << " ...";
    info(oss.str());

    fileSml_ = NULL;

    fileIn_.open(filePathIn);

    // Load "mzML_spectrumIndex"
    vector<li> mzML_spectrumIndex;
    fileIn_.readVector(mzML_spectrumIndex, "mzML_spectrumIndex", 0);

    // Load "mzML_chromatogramIndex"
    vector<li> mzML_chromatogramIndex;
    fileIn_.readVector(mzML_chromatogramIndex, "mzML_chromatogramIndex", 0);

    // Load "mzML" but without spectra
    vector<char> mzML;
    {
        vector<char> buffer;
        size_t offset = 0;
        size_t extent = (size_t) mzML_spectrumIndex.front();
        fileIn_.read_HypVecNC("mzML", buffer, &offset, &extent);
        mzML.insert(mzML.end(), buffer.begin(), buffer.end());

        offset = (size_t) mzML_spectrumIndex.back();
        extent = mzML_chromatogramIndex.front() - offset;
        fileIn_.read_HypVecNC("mzML", buffer, &offset, &extent);
        mzML.insert(mzML.end(), buffer.begin(), buffer.end());

        offset = (size_t) mzML_chromatogramIndex.back();
        vector<size_t> dims;
        fileIn_.readExtent(dims, "mzML");
        extent = dims[0] - offset;
        fileIn_.read_HypVecNC("mzML", buffer, &offset, &extent);
        mzML.insert(mzML.end(), buffer.begin(), buffer.end());
    }

    // parse
    xml::xml_document mzmlDoc;
    xml::xpath_node_set nodes;
    xml::xml_parse_result result = mzmlDoc.load_buffer_inplace(&mzML[0], sizeof(char) * mzML.size());
    if (!result) throw runtime_error("Error: In mzMLb input - " + string(result.description()));

    map<string, SpectrumMetadata::DataType> dataTypes;
    nodes = mzmlDoc.select_nodes("mzML/instrumentConfigurationList/instrumentConfiguration");
    for(xml::xpath_node_set::const_iterator itr = nodes.begin(); itr != nodes.end(); ++itr)
    {
        string id;
        istringstream(itr->node().attribute("id").value()) >> id;

        SpectrumMetadata::DataType dataType = SpectrumMetadata::DataType::Unknown;
        xml::xpath_node_set detectors = itr->node().select_nodes("componentList/detector");
        for(xml::xpath_node_set::const_iterator itrDetectors = detectors.begin(); itrDetectors != detectors.end(); ++itrDetectors)
        {
            // is it an "electron multiplier" detector?
            xml::xpath_node_set detectorTypes = itrDetectors->node().select_nodes("cvParam[@accession='MS:1000253']");
            if(!detectorTypes.empty())
            {
                if (dataType == SpectrumMetadata::DataType::Unknown || dataType == SpectrumMetadata::DataType::IonCount)
                    dataType = SpectrumMetadata::DataType::IonCount;
                else
                    throw runtime_error("Error: Inconsistent <detector> types found in mzMLb");
            }
            // is it an "photomultiplier" detector?
            detectorTypes = itrDetectors->node().select_nodes("cvParam[@accession='MS:1000116']");
            if(!detectorTypes.empty())
            {
                if (dataType == SpectrumMetadata::DataType::Unknown || dataType == SpectrumMetadata::DataType::IonCount)
                    dataType = SpectrumMetadata::DataType::IonCount;
                else
                    throw runtime_error("Error: Inconsistent <detector> types found in mzMLb");
            }
            // is it an "induction" detector (i.e. orbitrap or FT-ICR)?
            detectorTypes = itrDetectors->node().select_nodes("cvParam[@accession='MS:1000624']");
            if(!detectorTypes.empty())
            {
                if (dataType == SpectrumMetadata::DataType::Unknown || dataType == SpectrumMetadata::DataType::IonCurrent)
                    dataType = SpectrumMetadata::DataType::IonCurrent;
                else
                    throw runtime_error("Error: Inconsistent <detector> types found in mzMLb");
            }
        }

        dataTypes[id] = dataType;
    }

    // query default detector
    SpectrumMetadata::DataType defaultDataType = dataTypes[mzmlDoc.select_node("mzML/run").node().attribute("defaultInstrumentConfigurationRef").value()];

    // query number of spectra to create "metadata_" structure
    li ns;
    istringstream(mzmlDoc.child("mzML").child("run").child("spectrumList").attribute("count").value()) >> ns;
    metadata_.resize(ns);

    for (li i = 0; i < ns; i++)
    {
       // load each spectrum's mzML to get rest of info
        size_t offset = mzML_spectrumIndex[i];
        size_t extent = mzML_spectrumIndex[i + 1] - offset;
        fileIn_.read_HypVecNC("mzML", mzML, &offset, &extent);

        xml::xml_parse_result result = mzmlDoc.load_buffer_inplace(&mzML[0], sizeof(char) * mzML.size());
        if (!result) throw runtime_error("Error: In mzMLb input file - " + string(result.description()));

        // capture index
        istringstream(mzmlDoc.select_node("spectrum").node().attribute("index").value()) >> metadata_[i].mzmlSpectrumIndex;

        // capture mz and intensity array length
        istringstream(mzmlDoc.select_node("spectrum").node().attribute("defaultArrayLength").value()) >> metadata_[i].defaultArrayLength;

        // capture if MS1 spectrum
        nodes = mzmlDoc.select_nodes("spectrum/cvParam[@accession='MS:1000579']");
        if (!nodes.empty())
            metadata_[i].config = "MS1";

        // capture if MSn spectrum
        nodes = mzmlDoc.select_nodes("spectrum/cvParam[@accession='MS:1000580']");
        if (!nodes.empty())
            metadata_[i].config = "MSn";

        // capture polarity
        nodes = mzmlDoc.select_nodes("spectrum/cvParam[@accession='MS:1000129']");
        if (!nodes.empty())
            metadata_[i].id = "n";
        else
        {
            nodes = mzmlDoc.select_nodes("spectrum/cvParam[@accession='MS:1000130']");
            if (!nodes.empty())
                metadata_[i].id = "p";
            else
                metadata_[i].id = "u";
        }

        // capture scan info (we can only process files with one scan per spectra, so for us scan = spectrum)
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
            xml::xpath_node_set scanStartTimes = nodes.first().node().select_nodes("cvParam[@accession='MS:1000016']");
            if(!scanStartTimes.empty())
            {
                metadata_[i].startTimeString = scanStartTimes.first().node().attribute("value").value();
                istringstream(metadata_[i].startTimeString) >> metadata_[i].startTime;

                if(string(scanStartTimes.first().node().attribute("unitAccession").value()).compare("UO:0000031") == 0)
                    metadata_[i].startTime *= 60.0;
            }
            else
                metadata_[i].startTime = -1.0;

            // capture preset scan config (will tell us if sciex data is dda or dia
            xml::xpath_node_set presetScanConfigs = nodes.first().node().select_nodes("cvParam[@accession='MS:1000616']");
            if(!presetScanConfigs.empty())
            {
                metadata_[i].config += "-";
                metadata_[i].config += presetScanConfigs.first().node().attribute("value").value();
            }
        }
        else
            throw runtime_error("Error: <spectrum> missing <scan> in mzML input");

        nodes = mzmlDoc.select_nodes("spectrum/precursorList/precursor/selectedIonList/selectedIon/cvParam[@accession='MS:1000744']");
        if(!nodes.empty())
        {
            // capture precursor mz
            string precursor = nodes.first().node().attribute("value").value();
            replace(precursor.begin(), precursor.end(), '.', '-');

            metadata_[i].id += "-";
            metadata_[i].id += precursor;
        }

        // capture if centroided spectrum
        nodes = mzmlDoc.select_nodes("spectrum/cvParam[@accession='MS:1000127']");
        if (!nodes.empty())
        {
            metadata_[i].dataType = SpectrumMetadata::DataType::Centroided;
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

        if (getDebugLevel() % 10 >= 3)
        {
            {
                ostringstream oss;
                oss << getTimeStamp() << "  " << metadata_[i].mzmlSpectrumIndex << " id=" << metadata_[i].id;
                info(oss.str());
            }
            {
                ostringstream oss;
                oss << getTimeStamp() << "    Intensities dataset=" << metadata_[i].intensitiesDataset;
                oss << " offset=" << metadata_[i].intensitiesOffset;
                oss << " extent=" << metadata_[i].defaultArrayLength;
                info(oss.str());
            }
            {
                ostringstream oss;
                oss << getTimeStamp() << "    Mzs dataset=" << metadata_[i].mzsDataset;
                oss << " offset=" << metadata_[i].mzsOffset;
                oss << " extent=" << metadata_[i].defaultArrayLength;
                info(oss.str());
            }
            {
                ostringstream oss;
                oss << getTimeStamp() << "    start_time=" << metadata_[i].startTime << "s";
                info(oss.str());
            }
        }

        // display progress update
        if (getDebugLevel() % 10 >= 1)
        {
            if ((i + 1) % 10000 == 0 || (i + 1) == ns)
            {
                ostringstream oss;
                oss << getTimeStamp() << "  " << setw(1 + (int)(log10((float)ns))) << (i+1) << "/" << ns;
                info(oss.str());
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
        // NOTE: KEEPING LAST SPECTRUM BECAUSE POSSIBLY GETS OUT OF SYNC BETWEEN READ AND WRITE
        metadata_[metadata_.size() - 1].finishTime = metadata_[metadata_.size() - 1].startTime + (metadata_[metadata_.size() - 1].startTime - metadata_[metadata_.size() - 2].startTime); // horrible guess
        //metadata_.resize(metadata_.size() - 1);

        // sort into contiguous blocks for seaMass
        sort(metadata_.begin(), metadata_.end(), &DatasetMzmlb::seamassOrder);
    }
    else
    {
        metadata_[0].finishTime = -1.0; // startTime and finishTime not used anyway for 1D seaMass
    }

    // check if spectra within each 'config' have different 'id', if so then DDA else assume DIA
    //bool dia = true;
    //bool done = false;
    li offset = 0;
    for (ii i = 1; i <= (ii)metadata_.size(); i++)
    {
        if (i == metadata_.size() || metadata_[i].config != metadata_[i - 1].config)
        {
            li extent = i - offset;

            bool isDia = extent > 1;
            for (ii j = offset + 1; j < offset + extent; j++)
            {
                if (metadata_[j].id != metadata_[j - 1].id)
                {
                    isDia = false;
                    break;
                }
            }

            if (!isDia)
            {
                for (ii j = offset; j < offset + extent; j++)
                {
                    replace(metadata_[j].startTimeString.begin(), metadata_[j].startTimeString.end(), '.', '-');

                    metadata_[j].id += "-";
                    metadata_[j].id += metadata_[j].startTimeString;
                }
            }

            if (getDebugLevel() % 10 >= 1)
            {
                ostringstream oss;
                oss << getTimeStamp() << "  config=" << metadata_[i - 1].config << " offset=" << offset << " extent=" << extent << " type=" << (isDia ? "DIA" : "DDA");
                info(oss.str());
            }

            offset = i;
        }
    }

    if (!filePathStemOut.empty())
    {
        // Setup and start saving mzMLb output file...
        if(mzML.size() > 0) vector<char>().swap(mzML);
        idxDataArrayOffSet_=0;
        vector<InfoGrpVar> dataSet;
        fileIn_.readVector(specIdx_, "mzML_spectrumIndex", 0);
        fileIn_.searchGroup("mzML");
        fileIn_.searchGroup("chromatogram_MS_1000595_float");
        fileIn_.searchGroup("chromatogram_MS_1000595_double");
        dataSet = fileIn_.get_Info();


        size_t loc=0;
        size_t len=specIdx_[0];
        fileIn_.read_HypVecNC("mzML",mzML,&loc,&len);

        fileOut_ = new FileNetcdf(filePathStemOut + (writeType == Dataset::WriteType::InputOutput ? ".mzMLv" : ".mzMLb"), NC_NETCDF4);

        if (dataSet.back().varName == "chromatogram_MS_1000595_double")
        {
            vector<double> chroMz;
            vector<fp> chroBinCounts;
            fileIn_.readVector(chroMz, "chromatogram_MS_1000595_double", 0);
            fileIn_.readVector(chroBinCounts, "chromatogram_MS_1000515_float", 0);
            fileOut_->writeVector(chroMz, "chromatogram_MS_1000595_double");
            fileOut_->writeVector(chroBinCounts, "chromatogram_MS_1000515_float");
        }
        else if (dataSet.back().varName == "chromatogram_MS_1000595_float")
        {
            vector<fp> chroMz;
            vector<fp> chroBinCounts;
            fileIn_.readVector(chroMz, "chromatogram_MS_1000595_float", 0);
            fileIn_.readVector(chroBinCounts, "chromatogram_MS_1000515_float", 0);
            fileOut_->writeVector(chroMz, "chromatogram_MS_1000595_float");
            fileOut_->writeVector(chroBinCounts, "chromatogram_MS_1000515_float");
        }

        fileOut_->write_DefHypVecNC<char>("mzML",NC_UBYTE);
        fileOut_->write_DefHypVecNC<double>("spectrum_MS_1000514_double",NC_DOUBLE);
        fileOut_->write_DefHypVecNC<float>("spectrum_MS_1000515_float",NC_FLOAT);
        fileOut_->write_DefHypVecNC<li>("mzML_spectrumIndex",NC_INT64);

        fileOut_->write_CatHypVecNC("mzML",mzML);
        newMzmlIndex_ = mzML.size();
        fileOut_->write_CatHypVecNC("mzML_spectrumIndex", &newMzmlIndex_,1);

        string s = "mzMLb 0.5";
        fileOut_->writeAttribute(s, "version", "mzML");

        mzML.clear();
    }

    //mzMin = 50.0;
    //mzMax = 3000.0;

    mzBin.mzMin = 50.0;
    mzBin.mzMax = 3000.0;

    ostringstream oss2;
    if (getDebugLevel() % 10 >= 1)
        oss2 << getTimeStamp() << " ";
    oss2 << "Processing ...";
    info(oss2.str());
}


DatasetMzmlb::~DatasetMzmlb()
{
    if (fileOut_)
        delete fileOut_;
}


bool DatasetMzmlb::read(std::string& filePathSml, std::string &id)
{
    // display total progress update
    if(spectrumIndex_ / 1000 > lastSpectrumIndex_ / 1000 || spectrumIndex_ >= metadata_.size())
    {
        if (getDebugLevel() % 10 == 0)
        {
            for (ii i = 0; i < 256; i++)
                cout << "\b";
            cout << setw(1 + (int)(log10((float)metadata_.size()))) << spectrumIndex_ << "/" << metadata_.size() << " " << flush;
        }
        else
        {
            ostringstream oss;
            oss << getTimeStamp() << " ";
            oss << setw(1 + (int)(log10((float)metadata_.size()))) << spectrumIndex_ << "/" << metadata_.size() << " " << flush;
            info(oss.str());
        }

        if (getDebugLevel() % 10 >= 0 && spectrumIndex_ >= metadata_.size())
            cout << endl;

        lastSpectrumIndex_ = spectrumIndex_;
    }
    
    // return if we are finished
    if (spectrumIndex_ >= metadata_.size()) return false;

    // Ranjeet todo: write into SML file rather than into Seamass::Input
    //enum class Type { Binned, Sampled, Centroided } type;
    std::vector<fp> counts;
    std::vector<li> countsIndex;
    std::vector<double> locations;
    std::vector<double> startTimes;
    std::vector<double> finishTimes;
    std::vector<fp> exposures;

    // determine next set of spectra
    id = metadata_[spectrumIndex_].id;
    bool done = false;
    li offset = spectrumIndex_;
    for (; !done; spectrumIndex_++)
    {
        // if at the end or any metadata is different in the next spectrum, we have come to the end of this channel
        if (spectrumIndex_ == metadata_.size() - 1 || metadata_[spectrumIndex_].id != metadata_[spectrumIndex_ + 1].id)
        {
            countsIndex.push_back((li)counts.size());
            done = true;
        }
    }
    li extent_ = spectrumIndex_ - offset;

    // readMatrixSparse mzs
    if ((extent_ > 1 && getDebugLevel() % 10 >= 1) || getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "  Reading m/z values for id=" << id << " ...";
        info(oss.str());
    }

    startTimes.resize(extent_);
    finishTimes.resize(extent_);
    vector< vector<double> > mzs(extent_);
    for (li i = 0; i < extent_; i++)
    {
        startTimes[i] = metadata_[offset + i].startTime;
        finishTimes[i] = metadata_[offset + i].finishTime;

        if(metadata_[offset + i].defaultArrayLength > 0)
        {
            size_t rdLen = metadata_[offset + i].defaultArrayLength;
            size_t hypIdx = metadata_[offset + i].mzsOffset;
            fileIn_.read_HypVecNC(metadata_[offset + i].mzsDataset, mzs[i], &hypIdx, &rdLen);
        }

        // display progress update
        if (extent_ > 1 && ((i + 1) % 1000 == 0 || (i + 1) == extent_))
        {
            if (extent_ > 1 && getDebugLevel() % 10 == 0)
            {
                cout << "." << flush;
            }
            else if (getDebugLevel() % 10 >= 1)
            {
                ostringstream oss;
                oss << getTimeStamp() << "     " << setw(1 + (int) (log10((float) extent_))) << (i + 1) << "/"
                    << extent_;
                info(oss.str());
            }
        }
    }

    // estimate of mz_range [mz0,mz1] as
    // mzML spec does not guarentee this will exist (and hence do not know if regions have zero counts or were not scanned)
    double mz0 = numeric_limits<double>::max();
    double mz1 = 0.0;
    if (metadata_[offset].dataType == SpectrumMetadata::DataType::IonCount ||
            metadata_[offset].dataType == SpectrumMetadata::DataType::IonCurrent)
    {
        if (extent_ > 1 && getDebugLevel() % 10 == 0)
        {
            cout << "." << flush;
        }
        else if ((extent_ > 1 && getDebugLevel() % 10 >= 1) || getDebugLevel() % 10 >= 2)
        {
            ostringstream oss;
            oss << getTimeStamp() << "  Converting from " << (metadata_[offset].dataType == SpectrumMetadata::DataType::IonCount ? "IonCount" : "IonCurrent") << " to bins ...";
            info(oss.str());
        }

        for (ii i = 0; i < extent_; i++)
        {
            if (mzs[i].size() >= 2)
            {
                mz0 = mz0 < mzs[i].front() ? mz0 : mzs[i].front();
                mz1 = mz1 > mzs[i].back() ? mz1 : mzs[i].back();
            }
        }

        if (extent_ > 1 && getDebugLevel() % 10 >= 1)
        {
            ostringstream oss;
            oss << getTimeStamp() << "   autodetected_mz_window=[" << fixed << setprecision(3) << mz0 << "," << mz1
                << "]Th";
            info(oss.str());
        }

        // if this happens regularly, will have to use reported scan window or command line arguments
        if (mz1 <= mz0 || mz0 == numeric_limits<double>::max() || mz0 == 0.0)
            throw runtime_error("Error: Not enough datapoints in input mzMLb");
    }

    // This modifies the raw data for some limitations of the mzML spec and makes
    // sure the intensities are treated as binned between m/z datapoints.
    //
    // For IonCount data, it interpolates the mzs to represent the edges of the bins rather
    //   than the centres
    // For IonCurrent data, it converts the sampled data to binned counts by treating the
    //   mz values as the bin edges, and using trapezoid rule to integrate intensity values
    //
    countsIndex.resize(extent_ + 1);
    for (ii i = 0; i < extent_; i++)
    {
        countsIndex[i] = (li) counts.size();

        // load intensities
        vector<fp> intensities;
        if(metadata_[offset + i].defaultArrayLength > 0)
        {
            size_t rdLen = metadata_[offset + i].defaultArrayLength;
            size_t hypIdx = metadata_[offset + i].intensitiesOffset;
            fileIn_.read_HypVecNC(metadata_[offset + i].intensitiesDataset, intensities, &hypIdx, &rdLen);
        }

        switch (metadata_[offset + i].dataType)
        {
            case SpectrumMetadata::DataType::Unknown: // Just save as profile mode sampled data
            {
                //type = Type::Sampled;

                locations.insert(locations.end(), mzs[i].begin(), mzs[i].end());
                counts.insert(counts.end(), intensities.begin(), intensities.end());
            } break;
            case SpectrumMetadata::DataType::Centroided: // Just save as centroided data
            {
                //type = Type::Centroided;

                locations.insert(locations.end(), mzs[i].begin(), mzs[i].end());
                counts.insert(counts.end(), intensities.begin(), intensities.end());
            } break;
            case SpectrumMetadata::DataType::IonCount: // ToF, Quad, Ion trap etc, is already binned but have to create bin edges and exposures
            {
                //type = Type::Binned;

                if (intensities.size() == 0) // if empty we MUST add a zero bin!
                {
                    locations.push_back(mz0);
                    counts.push_back(0.0);
                    locations.push_back(mz1);

                    exposures.push_back(1.0);
                }
                else
                {
                    // dividing by minimum to get back to ion counts for SWATH data which appears to be automatic gain
                    // controlled to correct for dynamic range restrictions (hack!)
                    double minimum = std::numeric_limits<double>::max();
                    for (size_t k = 0; k < mzs[i].size(); k++)
                        if (intensities[k] > 0.0)
                            minimum = minimum < intensities[k] ? minimum : intensities[k];
                    // check to see if we can estimate the exposure i.e. is the minimum reasonable?
                    if (minimum < 1.0 || minimum > 1000.0) minimum = 1.0;
                    exposures.push_back((fp) (1.0 / minimum));

                    if (intensities.front() == 0.0) // only use the first m/z if the intensity is zero
                    {
                        double frontEdge = mz0 < mzs[i].front() ? mz0 : mzs[i].front();
                        locations.push_back(frontEdge);
                        counts.push_back(0.0);
                    }

                    // use all the intensities except first and last
                    for (ii k = 1; k < (ii) intensities.size() - 1; k++)
                    {
                        if (intensities[k] != 0.0 || intensities[k - 1] != 0.0) // merge zeros
                        {
                            locations.push_back(0.5 * (mzs[i][k - 1] + mzs[i][k]));
                            counts.push_back(((fp) intensities[k]) * exposures[i]);
                        }
                    }
                    locations.push_back(0.5 * (mzs[i][mzs[i].size() - 2] + mzs[i].back()));

                    if (intensities.back() == 0.0) // only use the last m/z if the intensity is zero
                    {
                        counts.push_back(0.0);
                        double backEdge = mz1 > mzs[i].back() ? mz1 : mzs[i].back();
                        locations.push_back(backEdge);
                    }
                }
            } break;
            case SpectrumMetadata::DataType::IonCurrent: // Orbitrap, FT-ICR etc, bin this data
            {
                //type = Type::Binned;

                if (intensities.front() == 0.0 && mz0 < mzs[i].front())
                {
                    locations.push_back(mz0);
                    counts.push_back(0.0);
                }

                for (ii k = 0; k < (ii) intensities.size() - 1; k++)
                {
                    locations.push_back(mzs[i][k]);
                    counts.push_back((fp) ((mzs[i][k + 1] - mzs[i][k]) * 0.5 * (intensities[k + 1] + intensities[k])));
                }
                locations.push_back(mzs[i][mzs[i].size() - 1]);

                if (intensities.back() == 0.0 && mz1 > mzs[i].back())
                {
                    counts.push_back(0.0);
                    locations.push_back(mz1);
                }
           } break;
        }

        vector<double>().swap(mzs[i]);

        // display progress update
        if (extent_ > 1 && ((i + 1) % 1000 == 0 || (i + 1) == extent_))
        {
            if (extent_ > 1 && getDebugLevel() % 10 == 0)
            {
                cout << "." << flush;
            }
            else if (getDebugLevel() % 10 >= 1)
            {
                ostringstream oss;
                oss << getTimeStamp() << "     " << setw(1 + (int) (log10((float) extent_))) << (i + 1) << "/"
                    << extent_;
                info(oss.str());
            }
        }
    }
    countsIndex.back() = (li)counts.size();

    if (fileSml_ == NULL)
    {
        fileSml_ = new FileNetcdf(filePathSml, NC_NETCDF4);
        li zero=0;
        if (scale_[0] != -1 && scale_[1] != -1)
        {
            smlDataId_.resize(scale_[1]-scale_[0]+1);
            for (short i = scale_[0]; i < scale_[1]; ++i)
            {
                ostringstream oss;
                oss << "/mzScale=" << i;
                string groupName = oss.str();
                smlDataId_[i].matrixId = fileSml_->createGroup(groupName.c_str());
                smlDataId_[i].ijsId = fileSml_->write_VecNC("ijs", &zero, 1, NC_LONG,smlDataId[i].matrixId, true);
                smlDataId_[i].jsId = fileSml_->write_VecNC("js", vector<long int>(), NC_LONG, smlDataId[i].matrixId, true);
                smlDataId_[i].vsId = fileSml_->write_VecNC("vs", vector<float>(), NC_FLOAT, smlDataId[i].matrixId, true)
            }
        }
        else
        {
            smlDataId_.resize(1);
            ostringstream oss;
            oss << "/mzScale=" << (scale_[0] == -1 ? scale_[1]: scale_[0]);
            string groupName = oss.str();
            smlDataId_[0].matrixId = fileSml_->createGroup(groupName.c_str());
            smlDataId_[0].ijsId = fileSml_->write_VecNC("ijs", &zero, 1, NC_LONG,smlDataId[0].matrixId, true);
            smlDataId_[0].jsId = fileSml_->write_VecNC("js", vector<long int>(), NC_LONG, smlDataId[0].matrixId, true);
            smlDataId_[0].vsId = fileSml_->write_VecNC("vs", vector<float>(), NC_FLOAT, smlDataId[0].matrixId, true)
        }
    }

    if (fileSml_->getFileStatus())
    {
        genAxis(mzBin_,scale_[0]);
        rebinMZ(locations , counts, mzBin_);
    }

    return true;
}


void DatasetMzmlb::write(const std::string& filePathSml, const std::string &id)
{
    // Ranjeet todo: read from SML file rather than from Seamass::Input
    /*
     
    ///////// NOTE: 'id' IS CURRENTLY IGNORED, ASSUMES YOU ARE WRITING WHAT YOU'VE JUST READ ////////////
    li n = input.countsIndex.size() == 0 ? 1 : input.countsIndex.size() - 1;

    if ((n > 1 && getDebugLevel() % 10 >= 1) || getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "  Writing=" << id << " ...";
        info(oss.str());
    }

    li offset = spectrumIndex_ - extent_;
    for (ii i = 0; i < n; ++i)
    {
        fp exposure = input.exposures.size() > 0 ? input.exposures[i] : 1.0;

        vector<double> mzs;
        vector<fp> intensities;
        bool isCentroided = false;
        switch (input.type)
        {
            case Seamass::Input::Type::Binned:
                for (ii ci = input.countsIndex[i]; ci < input.countsIndex[i + 1]; ci++)
                {
                    ii li = ci + i;
                    mzs.push_back(0.5 * (input.locations[li] + input.locations[li + 1]));
                    intensities.push_back(input.counts[ci] / (exposure * (input.locations[li + 1] - input.locations[li])));
                }
                break;
            case Seamass::Input::Type::Centroided:
                isCentroided = true;
            case Seamass::Input::Type::Sampled:
                for (ii ci = input.countsIndex[i]; ci < input.countsIndex[i + 1]; ci++)
                {
                    mzs.push_back(input.locations[ci]);
                    intensities.push_back(input.counts[ci] / exposure);
                }
                break;
            default:
                throw runtime_error("BUG: input has no type");
        }

        writeXmlSpectrum(i + offset, mzs, intensities, isCentroided);
    }*/
}


bool DatasetMzmlb::read(std::string& filePathSml, Seamass::Output &output, std::string &id)
{
    throw runtime_error("BUG: not yet implemented!");
}


void DatasetMzmlb::write(const std::string& filePathSml, const Seamass &seamass, const std::string &id)
{
    throw runtime_error("BUG: not yet implemented!");
}


bool DatasetMzmlb::startTimeOrder(const SpectrumMetadata &lhs, const SpectrumMetadata &rhs)
{
    return lhs.startTime < rhs.startTime;
}


bool DatasetMzmlb::seamassOrder(const SpectrumMetadata &lhs, const SpectrumMetadata &rhs)
{
    if (lhs.config == rhs.config)
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
    else
    {
        return lhs.config < rhs.config;
    }
}


void DatasetMzmlb::rebinMZ(const vector<double>& rawMz, const vector<float>& rawSc, mzAxis& mzBin)
{
    for(li i = 0; i < rawMz.size(); ++i)
    {
        double dmz=fabs(rawMz[i+1]-rawMz[i]);
        double xn=(rawMz[i]-mzBin.mz[0])/mzBin.dmz;
        double xp=(rawMz[i+1]-mzBin.mz[0])/mzBin.dmz;
        li n = floor(xn);
        li p = floor(xp);

        //if(p > mzBin.xypxl.first-1) p=mzBin.xypxl.first-1;
        //if(n >= mzBin.xypxl.first-1) break;
        if(p > mzBin.idxMax-1) p=mzBin.idxMin-1;
        if(n >= mzBin.idxMin-1) break;
        // Overlapping image 1st box on boundary
        if(n < 0 && p >= 0)
        {
            n=0;
            for(li col=n; col <= p; ++col)
            {
                double boxBegin=0.0;
                double boxEnd=0.0;
                double scaleMZ=0.0;
                // Cut top overlap else its contained
                (mzBin.mz[col] < rawMz[i]) ? boxBegin=rawMz[i] : boxBegin=mzBin.mz[col];
                // Cut bottom overlap else its contained
                (mzBin.mz[col+1] < rawMz[i+1]) ? boxEnd=mzBin.mz[col+1] : boxEnd=rawMz[i+1];
                // Calculate MZ scaling
                scaleMZ=(boxEnd-boxBegin)/dmz;
                // Update Spectrum Count
                mzBin.sc[col]=mzBin.sc[col]+rawSc[i]*float(scaleMZ);
            }
        }
            // Multiple smaller image boxes within single raw data scanline
        else if(n < p && n >=0)
        {
            for(li col=n; col <= p; ++col)
            {
                double boxBegin=0.0;
                double boxEnd=0.0;
                double scaleMZ=0.0;
                // Cut top overlap else its contained
                (mzBin.mz[col] < rawMz[i]) ? boxBegin=rawMz[i] : boxBegin=mzBin.mz[col];
                // Cut bottom overlap else its contained
                (mzBin.mz[col+1] < rawMz[i+1]) ? boxEnd=mzBin.mz[col+1] : boxEnd=rawMz[i+1];
                // Calculate MZ scaling
                scaleMZ=(boxEnd-boxBegin)/dmz;
                // Update Spectrum Count
                mzBin.sc[col]=mzBin.sc[col]+rawSc[i]*float(scaleMZ);
            }
        }
            // Raw data scanline contained within a Image box hence no scale needed.
        else if(n == p && n > 0)
        {
            li col = n;
            double scaleMZ=1.0;
            // Update Sectural Count
            mzBin.sc[col]=mzBin.sc[col]+rawSc[i]*float(scaleMZ);
        }
    }
}

void DatasetMzmlb::genAxis(mzAxis& mz_, short scale)
{

    //ii offset = ii(floor(log2(min - PROTON_MASS) * (1L << mzScale))) - 1;
    //ii N = (ii(ceil(log2(max - PROTON_MASS) * (1L << mzScale))) + 1) - offset + 1;

    double mzScale = double (1L << scale);

    // make sure transform back to mz spans xMin -> xMax
    mz_.idxMin = ii(floor(log2(mz_.mzMin - PROTON_MASS) * mzScale));
    mz_.idxMax = ii(ceil(log2(mz_.mzMax - PROTON_MASS) * mzScale))+1;

    ii N = mz_.idxMax - mz_.idxMin + 1;

    mz_.mz.resize(N);
    mz_.dmz.clear();
    mz_.sc.resize(N-1, 0.0);

    for(ii i = mz_.idxMin; i <= mz_.idxMax; ++i)
    {
        mz_.mz[i] = pow(2,(i/mzScale)) + PROTON_MASS;
        if(i > mz_.idxMin)
            mz_.dmz.push_back(mz_.mz[i-1]-mz_.mz[i]);
    }
}

/*
void DatasetMzmlb::writeVecData(vector<float>& data_)
{
    fileOut_->write_CatHypVecNC("spectrum_MS_1000515_float",data_);
}


void DatasetMzmlb::writeXmlData()
{
    vector<char> mzML;

    if(spectrumIndex_ < metadata_.size() || spectrumIndex_ == 1)
    {
        li offset=spectrumIndex_ - extent_;
        for(size_t i = size_t(offset); i < spectrumIndex_; ++i)
        {
            if(mzML.size() > 0) vector<char>().swap(mzML);

            xml::xml_document mzmlXmlScan;
            vector<double> mzScan;
            size_t idx = metadata_[i].mzmlSpectrumIndex;
            size_t loc = size_t(specIdx_[idx]);
            size_t len = size_t(specIdx_[idx + 1] - specIdx_[idx]);
            fileIn_.read_HypVecNC("mzML",mzML,&loc,&len);

            size_t xmlSize = sizeof(char) * mzML.size();
            xml::xml_parse_result result = mzmlXmlScan.load_buffer_inplace(&mzML[0],xmlSize);

            size_t index = getXmlValue<size_t>(mzmlXmlScan,"spectrum","index");
            size_t arrayLen = getXmlValue<size_t>(mzmlXmlScan,"spectrum","defaultArrayLength");
            size_t mzOffSet = getXmlValue<size_t>(mzmlXmlScan,
                                                  "spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000514_double']",
                                                  "offset");

            fileIn_.read_HypVecNC("spectrum_MS_1000514_double", mzScan, &mzOffSet,
                                &arrayLen);

            //mzScan.erase(mzScan.begin());   // REMOVING FIRST POINT? TODO: DO THIS ONLY FOR BINNED DATA
            //mzScan.erase(mzScan.end() - 1); // REMOVING LAST POINT? TODO: DO THIS ONLY FOR BINNED DATA

            arrayLen = mzScan.size();

            setXmlValue<size_t>(mzmlXmlScan, "spectrum", "index", i);
            setXmlValue<size_t>(mzmlXmlScan, "spectrum", "defaultArrayLength", arrayLen);

            setXmlValue<size_t>(mzmlXmlScan,
                                "spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000514_double']",
                                "offset", idxDataArrayOffSet_);
            setXmlValue<size_t>(mzmlXmlScan,
                                "spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000515_float']",
                                "offset", idxDataArrayOffSet_);
            idxDataArrayOffSet_ += arrayLen;

            stringstream newmzML;
            mzmlXmlScan.print(newmzML);
            string output = newmzML.str();
            vector<char>().swap(mzML);
            mzML.assign(output.begin(), output.end());

            newSpecIdx_.push_back(newSpecIdx_[i] + mzML.size());
            fileOut_->write_CatHypVecNC("mzML", mzML);
            fileOut_->write_CatHypVecNC("spectrum_MS_1000514_double", mzScan);
        }
    }
    if(spectrumIndex_ >= metadata_.size())
    {
        writeChromatogramXmlEnd();
    }
}


void DatasetMzmlb::writePeakData(VecMat<double>& mzPeak_, VecMat<float>& pkPeak_,
                                 vector<size_t>& mzpkVecSize_)
{
    for(size_t i = 0; i < mzpkVecSize_.size(); ++i)
    {
        fileOut_->write_CatHypVecNC("spectrum_MS_1000514_double",mzPeak_.m[i],mzpkVecSize_[i]);
        fileOut_->write_CatHypVecNC("spectrum_MS_1000515_float",pkPeak_.m[i],mzpkVecSize_[i]);
    }
}


void DatasetMzmlb::writePeakXmlData(vector<size_t>& mzpkVecSize_)
{
    vector<char> mzML;

    if(spectrumIndex_ < metadata_.size() || spectrumIndex_ == 1)
    {
        li offset=spectrumIndex_ - extent_;
        size_t j=0;
        for(size_t i = size_t(offset); i < spectrumIndex_; ++i, ++j)
        {
            if(mzML.size() > 0) vector<char>().swap(mzML);

            xml::xml_document mzmlXmlScan;
            size_t idx = metadata_[i].mzmlSpectrumIndex;
            size_t loc = size_t(specIdx_[idx]);
            size_t len = size_t(specIdx_[idx + 1] - specIdx_[idx]);
            fileIn_.read_HypVecNC("mzML",mzML,&loc,&len);

            size_t xmlSize = sizeof(char) * mzML.size();
            xml::xml_parse_result result = mzmlXmlScan.load_buffer_inplace(&mzML[0],xmlSize);

            size_t arrayLen = mzpkVecSize_[j];

            setXmlValue<size_t>(mzmlXmlScan, "spectrum", "index", i);
            setXmlValue<size_t>(mzmlXmlScan, "spectrum", "defaultArrayLength", arrayLen);

            setXmlValue<size_t>(mzmlXmlScan,
                                "spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000514_double']",
                                "offset", idxDataArrayOffSet_);
            setXmlValue<size_t>(mzmlXmlScan,
                                "spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000515_float']",
                                "offset", idxDataArrayOffSet_);
            idxDataArrayOffSet_ += arrayLen;

            stringstream newmzML;
            mzmlXmlScan.print(newmzML);
            string output = newmzML.str();
            vector<char>().swap(mzML);
            mzML.assign(output.begin(), output.end());

            newSpecIdx_.push_back(newSpecIdx_[i] + mzML.size());
            fileOut_->write_CatHypVecNC("mzML", mzML);
        }
    }
    if(spectrumIndex_ >= metadata_.size())
    {
        writeChromatogramXmlEnd();
    }
}
*/

void DatasetMzmlb::writeXmlSpectrum(li offset_, vector<double> &mzs_,
                                    vector<fp> &intensities_, bool isCentroided_)
{
    vector<char> mzML;
    //if(spectrumIndex_ < metadata_.size() || spectrumIndex_ == 1)
    if(offset_ < metadata_.size() || offset_ == 1)
    {
        if(mzML.size() > 0) vector<char>().swap(mzML);

        xml::xml_document mzmlXmlScan;
        size_t idx = metadata_[offset_].mzmlSpectrumIndex;
        size_t loc = size_t(specIdx_[idx]);
        size_t len = size_t(specIdx_[idx + 1] - specIdx_[idx]);
        fileIn_.read_HypVecNC("mzML",mzML,&loc,&len);

        size_t xmlSize = sizeof(char) * mzML.size();
        mzmlXmlScan.load_buffer_inplace(&mzML[0],xmlSize);

        size_t arrayLen = mzs_.size();
        if (isCentroided_ == true)
        {
            setXmlValue<string>(mzmlXmlScan,"spectrum/cvParam[@accession='MS:1000128']","name","centroid spectrum");
            setXmlValue<string>(mzmlXmlScan,"spectrum/cvParam[@accession='MS:1000128']","accession","MS:1000127");
        }

        setXmlValue<size_t>(mzmlXmlScan, "spectrum", "index", spectrumListIdx_);
        spectrumListIdx_++;

        setXmlValue<size_t>(mzmlXmlScan, "spectrum", "defaultArrayLength", arrayLen);

        setXmlValue<size_t>(mzmlXmlScan,
                            "spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000514_double']",
                            "offset", idxDataArrayOffSet_);
        setXmlValue<size_t>(mzmlXmlScan,
                            "spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000515_float']",
                            "offset", idxDataArrayOffSet_);
        idxDataArrayOffSet_ += arrayLen;

        stringstream newmzML;
        mzmlXmlScan.print(newmzML);
        string output = newmzML.str();
        vector<char>().swap(mzML);
        mzML.assign(output.begin(), output.end());

        newMzmlIndex_ += mzML.size();
        fileOut_->write_CatHypVecNC("mzML", mzML);
        fileOut_->write_CatHypVecNC("spectrum_MS_1000514_double",mzs_);
        fileOut_->write_CatHypVecNC("spectrum_MS_1000515_float",intensities_);
        fileOut_->write_CatHypVecNC("mzML_spectrumIndex", &newMzmlIndex_,1);
    }

    if(offset_ >= metadata_.size() - 1)
    {
        writeChromatogramXmlEnd();
    }
}

void DatasetMzmlb::writeChromatogramXmlEnd()
{
    vector<char> mzML;
    vector<li> chromatogramIdx;
    vector<li> newChromatogramIdx;

    // Read in xml that is in between the /spectrum and chromatogram tag.
    fileIn_.readVector(chromatogramIdx, "mzML_chromatogramIndex", 0);
    size_t loc = specIdx_.back();
    size_t len = chromatogramIdx.front() - loc;
    fileIn_.read_HypVecNC("mzML", mzML, &loc, &len);
    fileOut_->write_CatHypVecNC("mzML",mzML);

    // Update New Chromatogram Index.
    vector<size_t> lenTotal;
    fileOut_->readExtent(lenTotal, "mzML");
    newChromatogramIdx.push_back(lenTotal[0]);

    // Read in old Chromatogram xml block.
    mzML.clear();
    loc = size_t(chromatogramIdx.front());
    len = size_t(chromatogramIdx.back() - chromatogramIdx.front());
    fileIn_.read_HypVecNC("mzML", mzML, &loc, &len);
    fileOut_->write_CatHypVecNC("mzML",mzML);

    // Update end of New Chromatogram Index.
    lenTotal.clear();
    fileOut_->readExtent(lenTotal, "mzML");
    newChromatogramIdx.push_back(lenTotal[0]);

    // Write end of mxML xml file tags.
    lenTotal.clear();
    fileIn_.readExtent(lenTotal, "mzML");
    loc = size_t(chromatogramIdx.back());
    len = lenTotal.front() - loc;
    mzML.clear();
    fileIn_.read_HypVecNC("mzML", mzML, &loc, &len);
    fileOut_->write_CatHypVecNC("mzML",mzML);

    // Write new Chromatogram Index...
    fileOut_->writeVector(newChromatogramIdx, "mzML_chromatogramIndex");
}

//
// Created by ranjeet.
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
	return NULL;
    /*
    if(ext == "mzMLb")
		return new MSMLb(fName);
    if(ext == "csv")
		retun new MZcsv(fName);
    */
}


MSmzMLb3::MSmzMLb3(string fileName)
{
	mzMLb3File.open(fileName);
	this->extractData();
	hypIdx.resize(2,0);
	rdLen.resize(2,0);
	hypIdx[1]=0; // Read from first Column.
	rdLen[0]=1; // Always 1 Row to read.
    instrument_type = 1;
}


void MSmzMLb3::extractData(void)
{
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

	if(string(docmzML.child("mzML").child("instrumentConfigurationList").child("instrumentConfiguration").child("componentList").
				child("analyzer").child("cvParam").attribute("name").value()).compare("orbitrap") == 0)
		instrument_type=2;

	mzMLb3File.search_Group("spectrum_MS_1000514");
	mzMLb3File.search_Group("spectrum_MS_1000515");
	dataSetList=mzMLb3File.get_Info();

	spectra.resize(ns);
    for (size_t i = 0; i < ns; i++)
    {
		spectra[i].index = i;
		spectra[i].preset_config = config_indices[i];
		spectra[i].precursor_mz = precursor_mzs[i];
		spectra[i].scan_start_time = start_times[i];
		spectra[i].count = specSize[i];
    }
}

vector<spectrum> MSmzMLb3::getSpectrum()
{
	return spectra;
}

void MSmzMLb3::getScanMZ(vector<double> &mz, size_t index, size_t count)
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



bool MSmzMLb3File::scan_start_time_order(const spectrum& lhs, const spectrum& rhs)
{
	return lhs.scan_start_time <= rhs.scan_start_time;
}


bool MSmzMLb3File::seamass_order(const spectrum& lhs, const spectrum& rhs)
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


MSFile::MSFile(string in_file) :
	i(0)
{
	msFile = FileFactory::createFileObj(in_file);
	spectra = msFile->getSpectrum();
	instrument_type = msFile->getInstrument();

	int lastdot = in_file.find_last_of(".");
	string id = (lastdot == string::npos) ? in_file : in_file.substr(0, lastdot);

	// determine start_scan_time order of spectra
	sort(spectra.begin(), spectra.end(), &MSFile::scan_start_time_order);
	for (size_t i = 0; i < spectra.size(); i++)
	{
		spectra[i].scan_start_time_index = i;
	}
	scan_start_times.resize(spectra.size());
	for (size_t i = 0; i < spectra.size(); i++)
	{
		scan_start_times[i] = spectra[i].scan_start_time;
	}
	sort(spectra.begin(), spectra.end(), &MSFile::seamass_order);
}


MSFile::~MSFile()
{
	delete msFile;
}


bool MSFile::next(seaMass::Input& output)
{
	out.mzs.assign(spectra.size() - 1, vector<double>());
	out.intensities.assign(spectra.size() - 1, vector<double>());
	int loaded = 0;
	bool precursor_mz_is_constant = true;
	for (; i < spectra.size(); i++)
	{
		if (loaded > 1 && spectra[i].precursor_mz != spectra[i - 1].precursor_mz)
			precursor_mz_is_constant = false;

		if (spectra[i].count > 0 && spectra[i].index != spectra.size() - 1)
		{
			msFile->getScanMZ(out.mzs[spectra[i].scan_start_time_index], spectra[i].index, spectra[i].count);
			msFile->getScanIntensities(out.intensities[spectra[i].scan_start_time_index], spectra[i].index, spectra[i].count);
		}

		loaded++;
		if ((i == spectra.size() - 1 ||
			spectra[i].preset_config != spectra[i + 1].preset_config ||
			(spectra[i].precursor_mz != spectra[i + 1].precursor_mz && spectra[i].precursor_mz == 0.0)))
		{
			bin_mzs_intensities();
			i++;
			return true;
		}
	}
	return false;
}


void MSFile::bin_mzs_intensities()
{
	// This modifies the raw data for some limitations of the mzML spec and makes
	// sure the intensities are treated as binned between m/z datapoints.
	//
	// For ToF data, it interpolates the mzs to represent the edges of the bins rather
	//   than the centres
	// For FT data, it converts the sampled data to binned counts by treating the
	//   mz values as the bin edges, and using trapezoid rule to integrate intensity values
	//
	// This is all a bit rough at the moment, should be fitting splines to the data

	// initialise exposures to default of 1 (unit exposure)
	out.exposures.assign(out.intensities.size(), 1);

	// if more than one spectrum, ignore last as we do not know its scan end time

	if (instrument_type == 1) // ToF
	{
		#pragma omp parallel for
		for (int j = 0; j < out.mzs.size(); j++)
		{
			if (out.mzs[j].size() >= 2)
			{
				// dividing by minimum to get back to ion counts for SWATH data which appears to be automatic gain controlled to correct for dynamic range restrictions (hack!)
				double minimum = std::numeric_limits<double>::max();

				// we drop the first and last m/z datapoint as we don't know both their bin edges
				for (size_t i = 1; i < out.mzs[j].size(); i++)
				{
					if (out.intensities[j][i] > 0) minimum = minimum < out.intensities[j][i] ? minimum : out.intensities[j][i];

					// linear interpolation of mz extent (probably should be cubic)
					out.mzs[j][i - 1] = 0.5 * (out.mzs[j][i - 1] + out.mzs[j][i]);
					out.intensities[j][i - 1] = out.intensities[j][i];
				}
				out.mzs[j].resize(out.mzs[j].size() - 1);
				out.intensities[j].resize(out.intensities[j].size() - 2);
				out.exposures[j] = 1.0 / minimum;
				for (size_t i = 0; i < out.intensities[j].size(); i++)
				{
					out.intensities[j][i] *= out.exposures[j];
					//fp baseline = (mzs[j][i+1] - mzs[j][i]) * 3.0/8.0 * 272.0;//68.0;//136.0;
					//intensities[j][i] = intensities[j][i] >= baseline ? intensities[j][i] : baseline;

					// if (intensities[j][i] == 0.0) intensities[j][i] = 3.0/8.0;
				}
			}
			else
			{
				out.mzs[j].resize(0);
				out.intensities[j].resize(0);
			}
		}
	}
	/*else if (instrument_type == 2) // Orbitrap - but disabled for now as doesn't play nicely with merge_bins function
	{
	#pragma omp parallel for
	for (ii j = 0; j < mzs.size(); j++)
	if (mzs[j].size() >= 2)
	{
	// for Orbitrap only, mark zeros as missing data
	for (ii i = 1; i < mzs[j].size(); i++)
	if (intensities[j][i-1] <= 0.0 || intensities[j][i] <= 0.0)
	{
	intensities[j][i-1] = -1.0;
	}
	else
	{
	intensities[j][i-1] = (mzs[j][i] - mzs[j][i-1]) * 0.5 * (intensities[j][i] + intensities[j][i-1]);
	}
	intensities[j].resize(intensities[j].size() - 1);
	}
	else
	{
	mzs[j].resize(0);
	intensities[j].resize(0);
	}
	}*/
	else // FT-ICR (Orbitrap is a type of FT-ICR)
	{
		#pragma omp parallel for
		for (int j = 0; j < out.mzs.size(); j++)
		{
			if (out.mzs[j].size() >= 2)
			{
				for (size_t i = 1; i < out.mzs[j].size(); i++)
				{
					out.intensities[j][i - 1] = (out.mzs[j][i] - out.mzs[j][i - 1]) * 0.5 * (out.intensities[j][i] + out.intensities[j][i - 1]);
				}
				out.intensities[j].resize(out.intensities[j].size() - 1);
			}
			else
			{
				out.mzs[j].resize(0);
				out.intensities[j].resize(0);
			}
		}
	}
}
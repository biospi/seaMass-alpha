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

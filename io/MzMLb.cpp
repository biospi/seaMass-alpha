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

#include "MzMLb.hpp"
#include "mzMLxml.hpp"

OutmzMLb::OutmzMLb(string _filename, DatasetMzmlb& inputFile) : filename(_filename)
{
    idxOffSet=0;
    //msInputFile = &inputFile;
    //specFile = msInputFile;
    //FileNetcdf input(filename);
    input.open(filename);
    //input.read_VecNC("mzML",mzMLBuff_);
    input.read_VecNC("mzML_spectrumIndex",specIdx);
    //input.read_VecNC("mzML_chromatogramIndex",chroIdx);
    input.read_VecNC("chromatogram_MS_1000595_double",chroMz);
    input.read_VecNC("chromatogram_MS_1000515_float",chroBinCounts);

    //input.read_VecNC("spectrum_MS_1000514_double",mz);
    input.search_Group("mzML");

	vector<InfoGrpVar> dataSet;
    dataSet = input.get_Info();
    input.read_AttNC("version",dataSet[0].varid,versionID,dataSet[0].grpid);
    size_t loc[1]={0};
    size_t len[1]={specIdx[0]};
    input.read_HypVecNC("mzML",mzMLBuff,loc,len);
    //input.close();


    //size_t xmlSize=sizeof(char)*mzMLBuff_.size();
    //xml::xml_parse_result result = mzMLXML.load_buffer_inplace(&mzMLBuff_[0],xmlSize);

    size_t lastdot = filename.find_last_of(".");
    string outFileName=filename.substr(0,lastdot)+".out.mzMLb";
    mzMLbFileOut.open(outFileName,NC_NETCDF4);
    //mzMLbFileOut.write_VecNC("mzML_spectrumIndex",specIdx,NC_UINT64);
    //mzMLbFileOut.write_VecNC("spectrum_MS_1000514_double",mz,NC_DOUBLE);
    //mzMLbFileOut.write_VecNC("mzML",mzMLBuff_,NC_CHAR);
    //mzMLbFileOut.write_AttNC("mzML","version",versionID,NC_CHAR);
    //mzMLbFileOut.write_VecNC("mzML_chromatogramIndex",chroIdx,NC_UINT64);
    mzMLbFileOut.write_VecNC("chromatogram_MS_1000595_double",chroMz,NC_DOUBLE);
    mzMLbFileOut.write_VecNC("chromatogram_MS_1000515_float",chroBinCounts,NC_FLOAT);
    //mzMLbFileOut.write_DefHypVecNC<char>("mzML",NC_CHAR);
    mzMLbFileOut.write_DefHypVecNC<char>("mzML",NC_UBYTE);
    mzMLbFileOut.write_DefHypVecNC<double>("spectrum_MS_1000514_double",NC_DOUBLE);
    mzMLbFileOut.write_DefHypVecNC<float>("spectrum_MS_1000515_float",NC_FLOAT);

    mzMLbFileOut.write_CatHypVecNC("mzML",mzMLBuff);
    newSpecIdx.push_back(mzMLBuff.size());
    mzMLBuff.clear();
}

OutmzMLb::~OutmzMLb()
{
    input.close();
    mzMLbFileOut.close();
}

void OutmzMLb::writeVecData(vector<float>& _data)
{
    mzMLbFileOut.write_CatHypVecNC("spectrum_MS_1000515_float",_data);
}

void OutmzMLb::writeXmlData(vector<DatasetMzmlb::MzmlbSpectrumMetadata> *spec)
{
    //vector<size_t> mzMLSize;
    for(size_t i = 0; i < spec->size(); ++i)
    {
        if(mzMLBuff.size() > 0) vector<char>().swap(mzMLBuff);

        xml::xml_document mzMLScan;
        vector<double> mzScan;

        size_t idx=(*spec)[i].mzmlSpectrumIndex;
        size_t loc[1] = {specIdx[idx]};
        size_t len[1] = {specIdx[idx+1]-specIdx[idx]};
        input.read_HypVecNC("mzML", mzMLBuff, loc, len);

        size_t xmlSize=sizeof(char)*mzMLBuff.size();
        xml::xml_parse_result result = mzMLScan.load_buffer_inplace(&mzMLBuff[0],xmlSize);

        size_t index=getXmlValue<size_t>(mzMLScan,"spectrum","index");
        size_t arrayLen=getXmlValue<size_t>(mzMLScan,"spectrum","defaultArrayLength");

        size_t mzOffSet=getXmlValue<size_t>(mzMLScan,"spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000514_double']","offset");
        size_t intenOffSet=getXmlValue<size_t>(mzMLScan,"spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000515_float']","offset");

        input.read_HypVecNC("spectrum_MS_1000514_double",mzScan,&mzOffSet,&arrayLen);

        arrayLen=mzScan.size();

        mzScan.erase(mzScan.begin());
        mzScan.erase(mzScan.end()-1);

        arrayLen=mzScan.size();

        setXmlValue<size_t>(mzMLScan,"spectrum","index",i);
        setXmlValue<size_t>(mzMLScan,"spectrum","defaultArrayLength",arrayLen);

        setXmlValue<size_t>(mzMLScan,"spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000514_double']","offset",idxOffSet);
        setXmlValue<size_t>(mzMLScan,"spectrum/binaryDataArrayList/binaryDataArray/binary[@externalDataset='spectrum_MS_1000515_float']","offset",idxOffSet);
        idxOffSet+=arrayLen;

        stringstream newmzML;
        mzMLScan.print(newmzML);
        //mzMLScan.save(newmzML);
        //mzMLScan.save(newmzML,"",pugi::format_raw);
        string output = newmzML.str();
        vector<char>().swap(mzMLBuff);
        mzMLBuff.assign(output.begin(),output.end());

        //mzMLSize.clear();
        //mzMLSize = mzMLbFileOut.read_DimNC("mzML");
        //newSpecIdx.push_back(mzMLSize[0]);
        newSpecIdx.push_back(newSpecIdx[i]+mzMLBuff.size());
        mzMLbFileOut.write_CatHypVecNC("mzML",mzMLBuff);
        mzMLbFileOut.write_CatHypVecNC("spectrum_MS_1000514_double",mzScan);
    }
    string subxml("</spectrumList>\n");

    if(mzMLBuff.size() > 0) vector<char>().swap(mzMLBuff);
    mzMLBuff.assign(subxml.begin(),subxml.end());
    mzMLbFileOut.write_CatHypVecNC("mzML",mzMLBuff);

    vector<size_t> lenTotal=input.read_DimNC("mzML");
    size_t loc = specIdx.back();
    size_t len = lenTotal[0]-loc;
    input.read_HypVecNC("mzML", mzMLBuff, &loc, &len);

    subxml.clear();
    subxml.assign(mzMLBuff.begin(),mzMLBuff.end());
    subxml = subxml.substr(subxml.find("<chromatogramList"));

    mzMLBuff.clear();
    mzMLBuff.assign(subxml.begin(),subxml.end());

    vector<size_t > pos;
    findVecString(mzMLBuff,pos,"<chromatogram ","</chromatogram>");
    lenTotal.clear();
    lenTotal = mzMLbFileOut.read_DimNC("mzML");
    newChroIdx.push_back(lenTotal[0]+pos[0]);
    newChroIdx.push_back(lenTotal[0]+pos[1]);

    mzMLbFileOut.write_CatHypVecNC("mzML",mzMLBuff);
    mzMLbFileOut.write_AttNC("mzML","version",versionID,NC_CHAR);

    mzMLbFileOut.write_VecNC("mzML_spectrumIndex",newSpecIdx,NC_INT64);
    mzMLbFileOut.write_VecNC("mzML_chromatogramIndex",newChroIdx,NC_INT64);

}

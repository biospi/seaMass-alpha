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

#include "MzMLb.hpp"

OutmzMLb::OutmzMLb(string _filename) : filename(_filename)
{
    NetCDFile input(filename);
    vector<InfoGrpVar> dataSet;
    vector<char> versionID;
    input.read_VecNC("mzML",mzMLBuff);
    input.read_VecNC("mzML_spectrumIndex",specIdx);
    input.read_VecNC("mzML_chromatogramIndex",chroIdx);
    input.read_VecNC("chromatogram_MS_1000595_double",chroMz);
    input.read_VecNC("chromatogram_MS_1000515_float",chroBinCounts);

    input.read_VecNC("spectrum_MS_1000514_double",mz);
    input.search_Group("mzML");
    dataSet = input.get_Info();
    input.read_AttNC("version",dataSet[0].varid,versionID,dataSet[0].grpid);
    input.close();

    size_t xmlSize=sizeof(char)*mzMLBuff.size();
    xml::xml_parse_result result = mzMLXML.load_buffer_inplace(&mzMLBuff[0],xmlSize);

    size_t lastdot = filename.find_last_of(".");
    string outFileName=filename.substr(0,lastdot)+"_new.mzMLb";
    mzMLbFileOut.open(outFileName,NC_NETCDF4);
    mzMLbFileOut.write_VecNC("mzML_spectrumIndex",specIdx,NC_UINT64);
    mzMLbFileOut.write_VecNC("spectrum_MS_1000514_double",mz,NC_DOUBLE);
    mzMLbFileOut.write_VecNC("mzML",mzMLBuff,NC_CHAR);
    mzMLbFileOut.write_AttNC("mzML","version",versionID,NC_CHAR);
    mzMLbFileOut.write_VecNC("mzML_chromatogramIndex",chroIdx,NC_UINT64);
    mzMLbFileOut.write_VecNC("chromatogram_MS_1000595_double",chroMz,NC_DOUBLE);
    mzMLbFileOut.write_VecNC("chromatogram_MS_1000515_float",chroBinCounts,NC_FLOAT);
}

OutmzMLb::~OutmzMLb()
{
    mzMLbFileOut.close();
}

void OutmzMLb::writeData(vector<float>& _data)
{
    mzMLbFileOut.write_VecNC("spectrum_MS_1000515_float",_data,NC_FLOAT);
}

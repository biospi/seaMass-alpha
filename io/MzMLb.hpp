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

#ifndef SEAMASS_MZMLB_HPP
#define SEAMASS_MZMLB_HPP

#include <pugixml.hpp>
#include "NetCDFile.hpp"

namespace xml = pugi;

class OutmzMLb
{
public:
    OutmzMLb(string _filename);
    ~OutmzMLb();
    void writeData(vector<float>& _data);
private:
    string filename;
    vector<char> mzMLBuff;
    vector<double> mz;
    vector<unsigned long long int> specIdx;
    vector<unsigned long long int> chroIdx;
    vector<float> chroBinCounts;
    vector<double> chroMz;
    NetCDFile mzMLbFileOut;
    xml::xml_document mzMLXML;
};

#endif //SEAMASS_MZMLB_HPP

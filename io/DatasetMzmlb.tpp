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

#ifndef SEAMASS_DATASETMZMLB_TPP
#define SEAMASS_DATASETMZMLB_TPP

#include "DatasetMzmlb.hpp"

template<typename T>
T DatasetMzmlb::getXmlValue(xml::xml_document &scan, string xpath, string attrib)
{
	T value;
	xml::xpath_node_set tool;
	tool = scan.select_nodes(xpath.c_str());
	istringstream(tool.first().node().attribute(attrib.c_str()).value())>>value;
	return value;
}

template<typename T>
void DatasetMzmlb::setXmlValue(xml::xml_document &scan, string xpath, string attrib, T value)
{
	xml::xpath_node_set tool;
	tool = scan.select_nodes(xpath.c_str());

	ostringstream buff;
	buff << value;
	string newVal(buff.str());
	tool.first().node().attribute(attrib.c_str()).set_value(newVal.c_str());
}

#endif //SEAMASS_DATASETMZMLB_TPP

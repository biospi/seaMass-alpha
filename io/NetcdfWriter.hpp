//
// $Id$
//
//
// Original author: Ranjeet Bhamber <ranjeet.bhamber <a.t> bristol.ac.uk>
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

#ifndef SEAMASS_NETCDFWRITER_HPP
#define SEAMASS_NETCDFWRITER_HPP

#include <iostream>
#include "../kernel/NetcdfFile.hpp"
#include "../core/SeamassCore.hpp"


class NetcdfWriter
{
protected:
    std::string filename;
    NetCDFile fileout;
    template<typename T>
    void write(const std::string& objectname, std::vector<T>& cdata, nc_type xtype);

public:
    NetcdfWriter(const std::string& filename);
    ~NetcdfWriter();

    void write_input(SeamassCore::Input& input);
    void write_output(SeamassCore::Output& output, ii shrinkage, ii tolerance, ii page_size);
    void write_output_control_points(SeamassCore::ControlPoints& control_points);

    void write(const std::string& objectname, std::vector<unsigned char>& cdata);
    void write(const std::string& objectname, std::vector<float>& cdata);
    void write(const std::string& objectname, std::vector<double>& cdata);
    void write(const std::string& objectname, std::vector<long>& cdata);
    void write(const std::string& objectname, std::vector<long long>& cdata);
};

template<typename T>
void NetcdfWriter::write(const std::string& objectname, std::vector<T>& cdata, nc_type xtype)
{
    fileout.write_VecNC(objectname,cdata,xtype);
}


#endif //SEAMASS_NETCDFWRITER_HPP

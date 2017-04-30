//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> bristol.ac.uk>
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


#include "ObserverMatrix.hpp"
#include "FileNetcdf.hpp"
#include <iomanip>
using namespace std;


ObserverMatrix::ObserverMatrix()
{
}


ObserverMatrix::~ObserverMatrix()
{
}


void ObserverMatrix::notice(const std::string& message, const Matrix* a)
{
    if (a)
    {
        istringstream iss(message);
        char leftBracket;
        int id;
        iss >> leftBracket >> id;

        ostringstream oss;
        oss << setw(8) << setfill('0') << id << ".coo";
        FileNetcdf file(oss.str(), NC_NETCDF4);

        file.write(*a, "a");
    }
}


void ObserverMatrix::warning(const std::string& message, const Matrix* a)
{
    if (a)
    {
        istringstream iss(message);
        char leftBracket;
        int id;
        iss >> leftBracket >> id;

        ostringstream oss;
        oss << setw(8) << setfill('0') << id << ".coo";
        FileNetcdf file(oss.str(), NC_NETCDF4);

        file.write(*a, "a");
    }
}


void ObserverMatrix::error(const std::string& message, const Matrix* a)
{
    if (a)
    {
        istringstream iss(message);
        char leftBracket;
        int id;
        iss >> leftBracket >> id;

        ostringstream oss;
        oss << setw(8) << setfill('0') << id << ".coo";
        FileNetcdf file(oss.str(), NC_NETCDF4);

        file.write(*a, "a");
    }
}

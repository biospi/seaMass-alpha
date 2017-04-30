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


#include "Callback.hpp"
#include "kernel.hpp"
using namespace std;
using namespace kernel;


void Callback::addCallback(Callback* callback)
{
    callbacks_.push_back(callback);
}


Callback::Callback()
{
}


Callback::~Callback()
{
}


bool Callback::isCallback() const
{
    return callbacks_.size() > 0;
}


void Callback::notice(const string& message) const
{
    for (ii i = 0; i < ii(callbacks_.size()); i++)
        cout << getTimeStamp() << message << endl;
}


std::vector<Callback*> Callback::callbacks_;

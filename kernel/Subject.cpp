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


#include "Subject.hpp"
#include "iostream"
using namespace std;


Observer::Observer()
{
}


Observer::~Observer()
{
}


void Observer::notify(const std::string& message)
{
    cout << message << endl;
}



void Subject::registerObserver(Observer *observer)
{
    observers_.push_back(observer);
}


Subject::Subject()
{
}


Subject::~Subject()
{
}


bool Subject::areObservers() const
{
    return observers_.size() > 0;
}


void Subject::notifyObservers(const string& message) const
{
    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->notify(message);

    //for (ii i = 0; i < ii(observers_.size()); i++)
    //    cout << getTimeStamp() << message << endl;
}


std::vector<Observer*> Subject::observers_;

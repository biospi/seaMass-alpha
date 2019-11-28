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


void Subject::setDebugLevel(int debugLevel)
{
    debugLevel_ = debugLevel;
}


int Subject::getDebugLevel()
{
    return observers_.size() > 0 ? debugLevel_ : 0;
}


void Subject::info(const string &message) const
{
    int debugLevel = getDebugLevel();
    setDebugLevel(0);

    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->info(message);

    setDebugLevel(debugLevel);
}


void Subject::warning(const string &message) const
{
    int debugLevel = getDebugLevel();
    setDebugLevel(0);

    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->warning(message);

    setDebugLevel(debugLevel);
}


void Subject::error(const string &message) const
{
    int debugLevel = getDebugLevel();
    setDebugLevel(0);

    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->error(message);

    setDebugLevel(debugLevel);
}


int Subject::debugLevel_ = 0;


std::vector<Observer*> Subject::observers_;

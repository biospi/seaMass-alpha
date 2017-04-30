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


#include "SubjectMatrixSparse.hpp"
#include <MatrixSparse.hpp>
#include <iostream>
using namespace std;


void SubjectMatrixSparse::registerObserver(ObserverMatrixSparse *observer)
{
    observers_.push_back(observer);
}


SubjectMatrixSparse::SubjectMatrixSparse()
{
}


SubjectMatrixSparse::~SubjectMatrixSparse()
{
}


void SubjectMatrixSparse::info(const string &message, const MatrixSparse *a) const
{
    int debugLevel = getDebugLevel();
    setDebugLevel(0);

    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->notice(message, a);

    setDebugLevel(debugLevel);

    Subject::info(message);
}


void SubjectMatrixSparse::warning(const string &message, const MatrixSparse* a) const
{
    int debugLevel = getDebugLevel();
    setDebugLevel(0);

    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->warning(message, a);

    setDebugLevel(debugLevel);

    Subject::warning(message);
}


void SubjectMatrixSparse::error(const string &message, const MatrixSparse* a) const
{
    int debugLevel = getDebugLevel();
    setDebugLevel(0);

    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->error(message, a);

    setDebugLevel(debugLevel);

    Subject::error(message);
}


std::vector<ObserverMatrixSparse*> SubjectMatrixSparse::observers_;

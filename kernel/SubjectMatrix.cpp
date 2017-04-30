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


#include "SubjectMatrix.hpp"
#include <Matrix.hpp>
#include <iostream>
using namespace std;


void SubjectMatrix::registerObserver(ObserverMatrix *observer)
{
    observers_.push_back(observer);
}


SubjectMatrix::SubjectMatrix()
{
}


SubjectMatrix::~SubjectMatrix()
{
}


void SubjectMatrix::notice(const string &message, const Matrix* a) const
{
    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->notice(message, a);

    Subject::notice(message);
}


void SubjectMatrix::warning(const string &message, const Matrix* a) const
{
    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->warning(message, a);

    Subject::warning(message);
}


void SubjectMatrix::error(const string &message, const Matrix* a) const
{
    for (int i = 0; i < int(observers_.size()); i++)
        observers_[i]->error(message, a);

    Subject::error(message);
}


std::vector<ObserverMatrix*> SubjectMatrix::observers_;

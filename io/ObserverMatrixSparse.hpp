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


#ifndef SEAMASS_KERNEL_OBSERVERMATRIXSPARSE_HPP
#define SEAMASS_KERNEL_OBSERVERMATRIXSPARSE_HPP


#include "../kernel/Observer.hpp"
#include <string>
#include <vector>


class MatrixSparse;


class ObserverMatrixSparse : public Observer
{
public:
    ObserverMatrixSparse();
    virtual ~ObserverMatrixSparse();

    using Observer::info;
    using Observer::warning;
    using Observer::error;
    
    virtual void info(const std::string& message, const MatrixSparse* a);
    virtual void warning(const std::string& message, const MatrixSparse* a);
    virtual void error(const std::string& message, const MatrixSparse* a);
};


#endif

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

#ifndef SEAMASS_IO_VECMAT_HPP
#define SEAMASS_IO_VECMAT_HPP


#include <vector>
#include <cmath>
using namespace std;


typedef long long int lli;
typedef unsigned long long uli;


template<typename T = float>
struct VecMat
{
    VecMat(void);
    VecMat(uli _r, uli _c, const vector<T> &_vec);
    VecMat(uli _r, uli _c);
    vector<T> v; // Vector of Matrix data.
    T** m; // Data Matrix
    void set(uli _r, uli _c, const vector<T> &_vec);
    void set(uli _r, uli _c);
    void getDims(uli dims[]) const;
    void clear(void);
private:
    vector<T*> matIdx;
    uli row;
    uli col;
};


#include "VecMat.tpp"


#endif

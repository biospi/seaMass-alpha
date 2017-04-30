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

#ifndef SEAMASS_IO_VECMAT_TPP
#define SEAMASS_IO_VECMAT_TPP


#include "VecMat.hpp"


template<typename T>
VecMat<T>::VecMat(void):m(0),row(0),col(0){}

template<typename T>
VecMat<T>::VecMat(uli _r, uli _c, const vector<T> &_vec):v(_vec),row(_r), col(_c)
{
    matIdx.resize(row,0);
    for(uli i=0; i < row; ++i)
        matIdx[i]=&v[i*col];
    m=&matIdx[0];
}

template<typename T>
VecMat<T>::VecMat(uli _r, uli _c) : row(_r), col(_c)
{
    v.resize(row*col);
    matIdx.resize(row,0);
    for(uli i=0; i < row; ++i)
        matIdx[i]=&v[i*col];
    m=&matIdx[0];
}

template<typename T>
void VecMat<T>::set(uli _r, uli _c, const vector<T> &_vec)
{
    this->clear();
    row=_r;
    col=_c;
    v=_vec;
    matIdx.resize(row,0);
    for(uli i=0; i < row; ++i)
        matIdx[i]=&v[i*col];
    m=&matIdx[0];
}

template<typename T>
void VecMat<T>::set(uli _r, uli _c)
{
    this->clear();
    row=_r;
    col=_c;
    v.resize(row*col,0);
    matIdx.resize(row,0);
    for(uli i=0; i < row; ++i)
        matIdx[i]=&v[i*col];
    m=&matIdx[0];
}

template<typename T>
void VecMat<T>::getDims(uli dims[]) const
{
    dims[0]=row;
    dims[1]=col;
}

template<typename T>
void VecMat<T>::clear(void)
{
    row=0;
    col=0;
    vector<T>().swap(this->v);
    vector<T*>().swap(this->matIdx);
    m=0;
}


#endif

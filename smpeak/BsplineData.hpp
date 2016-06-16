//
// $Id$
//
//
// Author: Ranjeet Bhamber <ranjeet <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Liverpool, UK
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

#ifndef SMPEAK_BSPLINEDATA_HPP_
#define SMPEAK_BSPLINEDATA_HPP_

#include "peakcore.hpp"
#include "SMData.hpp"
#include "core.hpp"

template<typename T = float>
struct BasisPatch
{
public:
	BasisPatch();
	void set(vector<T> &x, vector<int> &t);
	VecMat<T> b;
};

template<typename R = double, typename T = float>
class BsplineData
{
public:
	BsplineData(DataAxis<T,R> &bs, DataAxis<T,R> &dbs, DataAxis<T,R> &d2bs);
	vector<DataAxis<T,R>* > get(void);
private:
	vector<DataAxis<T,R>* > bspObjP;
};


template<typename R = double, typename T = float>
class BsplineBasisData
{
public:
	BsplineBasisData(DataAxis<T,R> &bs,
			DataAxis<T,R> &dhbs, DataAxis<T,R> &d2hbs,
			DataAxis<T,R> &dvbs, DataAxis<T,R> &d2vbs,
			BasisPatch<T> &bp);
	BsplineBasisData(DataAxis<T,R> &bs,
			DataAxis<T,R> &dhbs, DataAxis<T,R> &d2hbs,
			DataAxis<T,R> &dvbs, BasisPatch<T> &bp);
	void get(vector<DataAxis<T,R>* > &bsDat, BasisPatch<T> *&bp);
	void dumpData(string filename, const H5::DataType &data_type_id = H5::PredType::NATIVE_FLOAT);
private:
	vector<DataAxis<T,R>* > bspObjP;
	BasisPatch<T> *bPat;
};

#include "BsplineData.tpp"

#endif /* SMPEAK_BSPLINEDATA_HPP_ */

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

#ifndef SMPEAK_PEAKMANAGER_HPP_
#define SMPEAK_PEAKMANAGER_HPP_

#include "../io/VecMat.hpp"


template
<
	template<class Peak> class PeakCon,
	template<class Bspline1,class Bspline2> class BsplineCon,
	template
	<
		template<class SubOperator1> class Operator1,
		template<class SubOperator1, class SubOperator2> class Operator2,
		class OpT,
		class OpR
	> class PeakOp,
	typename R = double,
	typename T = float
>
class PeakManager : public PeakOp<PeakCon,BsplineCon, T, R>
{
public:
	PeakManager(BsplineCon<R,T> &_data, T _threshold = 0) : threshold(_threshold)
	{
		peak=new PeakCon<T>();
		data=&_data;
	};
	void execute()
	{
		this->calculate(this->peak, this->data, threshold);
	};
	PeakCon<T> *peak;
	~PeakManager() {delete this->peak;};
private:
	BsplineCon<R,T> *data;
	T threshold;
};


#endif /* SMPEAK_PEAKMANAGER_HPP_ */

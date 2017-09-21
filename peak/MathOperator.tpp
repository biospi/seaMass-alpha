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

#ifndef SMPEAK_MATHOPERATOR_TPP_
#define SMPEAK_MATHOPERATOR_TPP_

#include <BasisBsplineMz.hpp>
#include "MathOperator.hpp"

template<class T>
void OpUnitS<T>::axisRT(li dims, int _offset, double rt_res, vector<double> &_rt)
{
	_rt.resize(dims);
	double offset = double(_offset);
	double ppbrt = 1.0 / (pow(2.0, rt_res));
	for (li i = 0; i < _rt.size(); ++i) {
		_rt[i] = (offset + i) * ppbrt;
	}
}

template<class T>
void OpUnitS<T>::axisMZ(li dims, int _offset, double mz_res, vector<double> &_mz)
{
	_mz.resize(dims);
	double offset = double(_offset);
	//double ppbmz = 1.0033548378 / (pow(2.0, mz_res) * 60.0);
	//for (li i = 0; i < _mz.size(); ++i) {
	//	_mz[i] = (offset + i) * ppbmz;
	//}
    for (li i = 0; i < _mz.size(); ++i)
    {
        _mz[i] = pow(2.0, (offset + i) / double(1L << ii(mz_res) )) + BasisBsplineMz::PROTON_MASS;
    }
}


template<class T>
void OpNablaHS<T>::apply(li row, li col, T** alpha)
{
	for(li i=0; i < row; ++i)
	{
		T Nm1=alpha[i][0];
		for(li j=1; j < col; ++j)
		{
			T N=alpha[i][j];
			alpha[i][j]=alpha[i][j]-Nm1;
			Nm1=N;
		}
	}
	for(li i=0; i < row; ++i)
		alpha[i][0]=0.0;
}

template<class T>
void OpNablaHS<T>::axisMZ(li dims, int _offset, double mz_res, vector<double> &_mz)
{
	_mz.resize(dims);
	double offset = double(_offset) - 0.5;
    //double offset = double(_offset) - 0.5;
	//double ppbmz = 1.0033548378 / (pow(2.0, mz_res) * 60.0);
	//for (li i = 0; i < _mz.size(); ++i) {
	//	_mz[i] = (offset + i) * ppbmz;
	//}
    for (li i = 0; i < _mz.size(); ++i)
    {
        _mz[i] = pow(2.0, (offset + i) / double(1L << ii(mz_res) )) + BasisBsplineMz::PROTON_MASS;
    }
}


template<class T>
void OpNabla2HS<T>::apply(li row, li col, T** alpha)
{
	for(li i=0; i < row; ++i)
	{
		T Nm1=alpha[i][1];
		T Nm2=alpha[i][0];
		for(li j=2; j < col; ++j)
		{
			T N=alpha[i][j];
			alpha[i][j]=alpha[i][j]-2.0*Nm1+Nm2;
			Nm2=Nm1;
			Nm1=N;
		}
	}
	for(li j=0; j < 2; ++j)
		for(li i=0; i < row; ++i)
			alpha[i][j]=0.0;
}

template<class T>
void OpNabla2HS<T>::axisMZ(li dims, int _offset, double mz_res, vector<double> &_mz)
{
	_mz.resize(dims);
	double offset = double(_offset) - 1;
    //double offset = double(_offset) - 1;
	//double ppbmz = 1.0033548378 / (pow(2.0, mz_res) * 60.0);
	//for (li i = 0; i < _mz.size(); ++i) {
	//	_mz[i] = (offset + i) * ppbmz;
	//}
    for (li i = 0; i < _mz.size(); ++i)
    {
        _mz[i] = pow(2.0, (offset + i) / double(1L << ii(mz_res) )) + BasisBsplineMz::PROTON_MASS;
    }
}


template<class T>
void OpUnit<T>::axisRT(li dims, int _offset, double rt_res, vector<double> &_rt)
{
	_rt.resize(dims);
	double offset = double(_offset);
	double ppbrt = 1.0 / (pow(2.0, rt_res));
	//#pragma omp parallel for
	for (li i = 0; i < _rt.size(); ++i) {
		_rt[i] = (offset + i) * ppbrt;
	}
}

template<class T>
void OpUnit<T>::axisMZ(li dims, int _offset, double mz_res, vector<double> &_mz)
{
	_mz.resize(dims);
	double offset = double(_offset);
    //double offset = double(_offset);
	//double ppbmz = 1.0033548378 / (pow(2.0, mz_res) * 60.0);
	//#pragma omp parallel for
	//for (li i = 0; i < _mz.size(); ++i) {
	//	_mz[i] = (offset + i) * ppbmz;
	//}
    for (li i = 0; i < _mz.size(); ++i)
    {
        _mz[i] = pow(2.0, (offset + i) / double(1L << ii(mz_res) )) + BasisBsplineMz::PROTON_MASS;
    }
}


template<class T>
void OpNablaH<T>::apply(li row, li col, T** alpha)
{
	//#pragma omp parallel for
	for(li i=0; i < row; ++i)
	{
		T Nm1=alpha[i][0];
		for(li j=1; j < col; ++j)
		{
			T N=alpha[i][j];
			alpha[i][j]=alpha[i][j]-Nm1;
			Nm1=N;
		}
	}
	//#pragma omp parallel for
	for(li i=0; i < row; ++i)
		alpha[i][0]=0.0;
}

template<class T>
void OpNablaH<T>::axisMZ(li dims, int _offset, double mz_res, vector<double> &_mz)
{
	_mz.resize(dims);
	double offset = double(_offset) - 0.5;
    //double offset = double(_offset) - 0.5;
	//double ppbmz = 1.0033548378 / (pow(2.0, mz_res) * 60.0);
	//#pragma omp parallel for
	//for (li i = 0; i < _mz.size(); ++i) {
	//	_mz[i] = (offset + i) * ppbmz;
	//}
    for (li i = 0; i < _mz.size(); ++i)
    {
        _mz[i] = pow(2.0, (offset + i) / double(1L << ii(mz_res) )) + BasisBsplineMz::PROTON_MASS;
    }
}


template<class T>
void OpNabla2H<T>::apply(li row, li col, T** alpha)
{
	//#pragma omp parallel for
	for(li i=0; i < row; ++i)
	{
		T Nm1=alpha[i][1];
		T Nm2=alpha[i][0];
		for(li j=2; j < col; ++j)
		{
			T N=alpha[i][j];
			alpha[i][j]=alpha[i][j]-2.0*Nm1+Nm2;
			Nm2=Nm1;
			Nm1=N;
		}
	}
	//#pragma omp parallel for
	for(li j=0; j < 2; ++j)
		for(li i=0; i < row; ++i)
			alpha[i][j]=0.0;
}

template<class T>
void OpNabla2H<T>::axisMZ(li dims, int _offset, double mz_res, vector<double> &_mz)
{
	_mz.resize(dims);
	double offset = double(_offset) - 1;
    //double offset = double(_offset) - 1;
	//double ppbmz = 1.0033548378 / (pow(2.0, mz_res) * 60.0);
	//#pragma omp parallel for
	//for (li i = 0; i < _mz.size(); ++i) {
	//	_mz[i] = (offset + i) * ppbmz;
	//}
    for (li i = 0; i < _mz.size(); ++i)
    {
        _mz[i] = pow(2.0, (offset + i) / double(1L << ii(mz_res) )) + BasisBsplineMz::PROTON_MASS;
    }
}


template<class T>
void OpNablaV<T>::apply(li row, li col, T** alpha)
{
	//#pragma omp parallel for
	for(li i=0; i < col; ++i)
	{
		T Nm1=alpha[0][i];
		for(li j=1; j < row; ++j)
		{
			T N=alpha[j][i];
			alpha[j][i]=alpha[j][i]-Nm1;
			Nm1=N;
		}
	}
	//#pragma omp parallel for
	for(li i=0; i < col; ++i)
		alpha[0][i]=0.0;
}

template<class T>
void OpNablaV<T>::axisRT(li dims, int _offset, double rt_res, vector<double> &_rt)
{
	_rt.resize(dims);
	double offset = double(_offset) - 0.5;
	double ppbrt = 1.0 / (pow(2.0, rt_res));
	//#pragma omp parallel for
	for (li i = 0; i < _rt.size(); ++i) {
		_rt[i] = (offset + i) * ppbrt;
	}
}


template<class T>
void OpNabla2V<T>::apply(li row, li col, T** alpha)
{
	//#pragma omp parallel for
	for(li i=0; i < col; ++i)
	{
		T Nm1=alpha[1][i];
		T Nm2=alpha[0][i];
		for(li j=2; j < row; ++j)
		{
			T N=alpha[j][i];
			alpha[j][i]=alpha[j][i]-2.0*Nm1+Nm2;
			Nm2=Nm1;
			Nm1=N;
		}
	}
	//#pragma omp parallel for
	for(li i=0; i < 2; ++i)
		for(li j=0; j < col; ++j)
			alpha[i][j]=0.0;
}

template<class T>
void OpNabla2V<T>::axisRT(li dims, int _offset, double rt_res, vector<double> &_rt)
{
	_rt.resize(dims);
	double offset = double(_offset) - 1;
	double ppbrt = 1.0 / (pow(2.0, rt_res));
	//#pragma omp parallel for
	for (li i = 0; i < _rt.size(); ++i) {
		_rt[i] = (offset + i) * ppbrt;
	}
}

#endif /* SMPEAK_MATHOPERATOR_TPP_ */

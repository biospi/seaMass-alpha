#ifndef SMPEAK_MATHOPERATOR_HPP_
#define SMPEAK_MATHOPERATOR_HPP_

#include "peakcore.hpp"

template<class T>
class OpUnit
{
protected:
	void apply(lli row, lli col, T** alpha){};
	void axisRT(hsize_t dims, int _offset, double rt_res, vector<double> &_rt);
	void axisMZ(hsize_t dims, int _offset, double mz_res, vector<double> &_mz);
	~OpUnit(){};
};

template<class T>
class OpNablaH : public OpUnit<T>
{
protected:
	void apply(lli row, lli col, T** alpha);
	void axisMZ(hsize_t dims, int _offset, double mz_res, vector<double> &_mz);
	~OpNablaH(){};
};

template<class T>
class OpNabla2H : public OpUnit<T>
{
protected:
	void apply(lli row, lli col, T** alpha);
	void axisMZ(hsize_t dims, int _offset, double mz_res, vector<double> &_mz);
	~OpNabla2H(){};
};

template<class T>
class OpNablaV : public OpUnit<T>
{
protected:
	void apply(lli row, lli col, T** alpha);
	void axisRT(hsize_t dims, int _offset, double rt_res, vector<double> &_rt);
	~OpNablaV(){};
};


template<class T>
class OpNabla2V : public OpUnit<T>
{
protected:
	void apply(lli row, lli col, T** alpha);
	void axisRT(hsize_t dims, int _offset, double rt_res, vector<double> &_rt);
	~OpNabla2V(){};
};

#include"MathOperator.tpp"

#endif /* SMPEAK_MATHOPERATOR_HPP_ */

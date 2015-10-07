#ifndef SMPEAK_PEAKOPERATOR_HPP_
#define SMPEAK_PEAKOPERATOR_HPP_

#include "peakcore.hpp"


template<typename T>
class MathOp
{
protected:
	void calMidPoint(lli rtIdx, lli mzIdx, T** alpha,
					const vector<double> &mza, double &mz1, double &a1);
	double calT(const double y0, const double y1, const double y2);
	double calX(double t, double x0, double x1, double x2);
	vector<T> cal3rdMidPoint(lli rtIdx, lli mzIdx, T **P);
	T calPeakCount(vector<T> &ry, double t);
	void calPeakWidth(lli rtIdx,lli mzIdx, T** alpha, const vector<double> d2mz,
					double &mzlhs, double &mzrhs);
	void calPeakMZ(DataAxis<T> const *bs,DataAxis<T> const *dbs,DataAxis<T> const *d2bs,
					lli i, lli j, double &mzPeak, T &countMax,
					double &mzlhs, double &mzrhs, double &t0, lli &falsePeak);
};

template
<
	template <class Peak> class pPeak,
	template <class Data> class pData,
	typename T = float
>
class Centroid : public MathOp<T>
{
public:
	void calculate(pPeak<T> *peak, pData<T> *data);
protected:
	~Centroid(){};
};


template
<
	template <class Peak> class pPeak,
	template <class Data> class pData,
	typename T = float
>
class ExtractPeak : public MathOp<T>
{
public:
	void calculate(pPeak<T> *peak, pData<T> *data);
protected:
	~ExtractPeak(){};
private:
	void calPeakMZ(DataAxis<T> const *bs,DataAxis<T> const *dbs,DataAxis<T> const *d2bs,
					lli i, lli j, double &mzPeak, T &countMax,
					double &mzlhs, double &mzrhs, double &t0, lli &falsePeak, bool &peakRT);
};

#include "PeakOperator.tpp"

#endif /* SMPEAK_PEAKOPERATOR_HPP_ */

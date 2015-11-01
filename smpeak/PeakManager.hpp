#ifndef SMPEAK_PEAKMANAGER_HPP_
#define SMPEAK_PEAKMANAGER_HPP_

#include "peakcore.hpp"


template
<
	template<class Peak> class PeakCon,
	template<class Bspline1,class Bspline2> class BsplineCon,
	template
	<
		template<class SubOperator1> class Operator1,
		template<class SubOperator1, class SubOperator2> class Operator2,
		typename OpT,
		typename OpR
	> class PeakOp,
	typename T = float,
	typename R = double
>
class PeakManager : public PeakOp<PeakCon<T>,BsplineCon<T,R>, T, R>
{
public:
	PeakManager(BsplineCon<T,R> &_data)
	{
		peak=new PeakCon<T>();
		data=&_data;
	};
	void execute()
	{
		calculate(this->peak, this->data);
	};
	PeakCon<T> *peak;
	~PeakManager() {delete this->peak;};
	private:
	BsplineCon<T,R> *data;

};



#endif /* SMPEAK_PEAKMANAGER_HPP_ */

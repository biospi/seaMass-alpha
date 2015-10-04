#ifndef SMPEAK_PEAKCMD_HPP_
#define SMPEAK_PEAKCMD_HPP_

#include"peakcore.hpp"
//#include"PeakOperator.hpp"

template
<
	template<class Peak> class PeakCon,
	template<class Bspline> class BsplineCon,
	template<
		template<class SubOperator1> class Operator1,
		template<class SubOperator2> class Operator2,
		class Op> class PeakOp,
	typename T = float
>
class PeakCmd : public PeakOp<PeakCon<T>, BsplineCon<T>, T>
{
public:
	PeakCmd(BsplineCon<T> &_data)
	{
		peak=new PeakCon<T>();
		data=&_data;
	};
	void execute()
	{
		calculate(this->peak, this->data);
	};
	PeakCon<T> *peak;
	~PeakCmd() {delete this->peak;};
	private:
	BsplineCon<T> *data;

};



#endif /* SMPEAK_PEAKCMD_HPP_ */

#ifndef SMPEAK_PEAKOPERATOR_TPP_
#define SMPEAK_PEAKOPERATOR_TPP_

template<template<class> class pPeak, template<class> class pData, typename T>
void Centroid<pPeak,pData,T>::calculate(pPeak<T> *peak, pData<T> *data)
{
	vector<DataAxis<T>* > bsData =  data->get();

	DataAxis<T>* bs=bsData[0];
	DataAxis<T>* dbs=bsData[1];
	DataAxis<T>* d2bs=bsData[2];

	hsize_t dims[2];
	bs->alpha->getDims(dims);
	lli row = dims[0];
	lli col = dims[1];

	lli falsePeak=0;
	lli falseWidth=0;
	// Find Peaks and exact MZ values.
	cout<<"Preform Centroid of Mass Spec DATA along MZ!"<<endl;
	cout<<"Find Peaks along MZ axis"<<endl;

	#pragma omp parallel for reduction(+:falsePeak,falseWidth)
	for(lli i = 0; i < row; ++i)
	{
		for(lli j = 2; j < col-2; ++j)
		{
			if((dbs->alpha->m[i][j] > 0) && (dbs->alpha->m[i][j+1] < 0) )
			{
				double pa1=0.0;
				double pmz1=0.0;
				calMidPoint(i,j,dbs->alpha->m,dbs->mz,pmz1,pa1);
				if(pa1 < 0){
					double pa0=0.0;
					double pmz0=0.0;
					vector<T> ry;
					calMidPoint(i,j-1,dbs->alpha->m,dbs->mz,pmz0,pa0);
					double t0 = calT(pa0,double(dbs->alpha->m[i][j]),pa1);
					if(t0>=0)
					{
						double mzPeak=calX(t0,pmz0,dbs->mz[j],pmz1);
						double mzlhs=0.0;
						double mzrhs=0.0;
						ry = cal3rdMidPoint(i,j,bs->alpha->m);
						T countMax = calPeakCount(ry,t0);
						calPeakWidth(i,j,d2bs->alpha->m,d2bs->mz,mzlhs,mzrhs);
						if(mzlhs >= 0 && mzrhs >= 0)
						{
							#pragma omp critical(peak)
							{
								peak->addPeak(mzPeak,bs->rt[i],countMax,make_pair(mzlhs,mzrhs),
									make_pair(bs->rt[i],bs->rt[i]),t0,i,j);
							}
						}
						else
						{
							++falseWidth;
						}
					}
					else
					{
						++falsePeak;
					}
				}
				else
				{
					double pa2=0.0;
					double pmz2=0.0;
					vector<T> ry;
					calMidPoint(i,j+1,dbs->alpha->m,dbs->mz,pmz2,pa2);
					double t0 = calT(pa1,double(dbs->alpha->m[i][j+1]),pa2);
					if (t0>=0)
					{
						double mzPeak=calX(t0,pmz1,dbs->mz[j+1],pmz2);
						double mzlhs=0.0;
						double mzrhs=0.0;
						ry = cal3rdMidPoint(i,j,bs->alpha->m);
						T countMax = calPeakCount(ry,t0);
						calPeakWidth(i,j,d2bs->alpha->m,d2bs->mz,mzlhs,mzrhs);
						if(mzlhs >= 0 && mzrhs >= 0)
						{
							#pragma omp critical(peak)
							{
								peak->addPeak(mzPeak,bs->rt[i],countMax,make_pair(mzlhs,mzrhs),
									make_pair(bs->rt[i],bs->rt[i]),t0,i,j);
							}
						}
						else
						{
							++falseWidth;
						}
					}
					else
					{
						++falsePeak;
					}
				}
			}
		}
	}

	cout<<"Found ["<<peak->numOfPeaks()<<"] Peaks."<<endl;
	if(falsePeak > 0)
		cout<<"Warning !!! Found ["<<falsePeak<<"] Insignificant False Peaks Detected - Peaks Ignored"<<endl;
	if(falseWidth > 0)
		cout<<"Warning !!! Found ["<<falseWidth<<"] False Peaks Detected with Incorrect Peak Widths - Peaks Ignored"<<endl;
}

template<template<class> class pPeak, template<class> class pData, typename T>
void Centroid<pPeak,pData,T>::calMidPoint(lli rtIdx, lli mzIdx, T** alpha,
		vector<double> &mza, double &mz1, double &a1)
{
	mz1=0.5*(mza[mzIdx]+mza[mzIdx+1]);
	a1=double(0.5*(alpha[rtIdx][mzIdx]+alpha[rtIdx][mzIdx+1]));
}

template<template<class> class pPeak, template<class> class pData, typename T>
double Centroid<pPeak,pData,T>::calT(const double y0, const double y1, const double y2)
{
	double a=y1-y0;
	double b=y0-2.0*y1+y2;
	double g=a/b;
	double r=sqrt(g*g-y0/b);
	double tp=-g+r;
	double tn=-g-r;

	if(tp >= 0 && tp <= 1) return tp;
	else if(tn >= 0 && tn <= 1)	return tn;
	else
	{
		//cout<<"Insignificant False Peak Detected - Peak Ignored"<<endl;
		return -1.0;
	}
}

template<template<class> class pPeak, template<class> class pData, typename T>
double Centroid<pPeak,pData,T>::calX(double t, double x0, double x1, double x2)
{
	double c=(1.0-t);
	return c*c*x0+2.0*c*t*x1+t*t*x2;
}

template<template<class> class pPeak, template<class> class pData, typename T>
vector<T> Centroid<pPeak,pData,T>::cal3rdMidPoint(lli row, lli col, T **P)
{
	vector<T> ry(4,0.0);
	T p0=P[row][col-1];
	T p1=P[row][col];
	T p2=P[row][col+1];
	T p3=P[row][col+2];
	T u=1.0/3.0;
	T v=2.0/3.0;
	T w=1.0/2.0;
	T up=(1-u);
	T vp=(1-v);
	T wp=(1-w);

	ry[0]=wp*(vp*p0+v*p1)+w*(up*p1+u*p2);
	ry[1]=up*p1+u*p2;
	ry[2]=vp*p1+v*p2;
	ry[3]=wp*(vp*p1+v*p2)+w*(up*p2+u*p3);

	return ry;
}


template<template<class> class pPeak, template<class> class pData, typename T>
T Centroid<pPeak,pData,T>::calPeakCount(vector<T> &ry, double t)
{
	double tp=1-t;
	T p0=tp*ry[0]+t*ry[1];
	T p1=tp*ry[1]+t*ry[2];
	T p2=tp*ry[2]+t*ry[3];
	T q0=tp*p0+t*p1;
	T q1=tp*p1+t*p2;
	return tp*q0+t*q1;
}

template<template<class> class pPeak, template<class> class pData, typename T>
void Centroid<pPeak,pData,T>::calPeakWidth(lli rtIdx,lli mzIdx, T** alpha, vector<double> d2mz,
		double &mzlhs, double &mzrhs)
{
	lli n=d2mz.size();
	lli lhsIdx=mzIdx+1;
	lli rhsIdx=mzIdx+1;

	do
	{
		--lhsIdx;
		if(lhsIdx < 0) break;
	} while( alpha[rtIdx][lhsIdx] < 0);
	do
	{
		++rhsIdx;
		if(rhsIdx > n) break;
	} while( alpha[rtIdx][rhsIdx] < 0);

	if(lhsIdx > 0 && rhsIdx < n ){
		double x1,x2,y1,y2,m;

		x1=d2mz[lhsIdx];
		y1=alpha[rtIdx][lhsIdx];
		x2=d2mz[lhsIdx+1];
		y2=alpha[rtIdx][lhsIdx+1];
		m=(y2-y1)/(x2-x1);
		mzlhs=(-y1/m)+x1;

		x1=d2mz[rhsIdx];
		y1=alpha[rtIdx][rhsIdx];
		x2=d2mz[rhsIdx-1];
		y2=alpha[rtIdx][rhsIdx-1];
		m=(y2-y1)/(x2-x1);
		mzrhs=(-y1/m)+x1;

		/*
		if (mzlhs < d2mz[0])
		{
			cout<<"WTF LHS is out of range!!"<<endl;
		}
		if (mzrhs > d2mz[n-1])
		{
			cout<<"WTF RHS is out of range!!"<<endl;
		}
		*/
	}
	else
	{
		mzlhs=-1;
		mzrhs=-1;
	}

}

#endif /* SMPEAK_PEAKOPERATOR_TPP_ */
